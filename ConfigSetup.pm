package ConfigSetup;
use strict; 
use File::Path qw/make_path remove_tree/;
use Common; 
use Data::Dumper;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/config_setup/; 

sub config_setup {
    my ($f_cfg, $org, $f_fas, $f_gff, $cutoff_e) = @_;
   
    print "setting up environment variables\n";
    open(FH, "<$f_cfg") || die "config file $f_cfg is not there\n";
    while(<FH>) {
        chomp;
        next unless $_;
        next if /^\#/;
        $_ =~ s/\s//g;
        my ($k, $v) = split "=";
        while( $v =~ /\$\{(\w+)\}/g ) {
            die "no env variable named $1\n" unless exists $ENV{$1};
            my $rep = $ENV{$1};
            $v =~ s/\$\{$1\}/$rep/;
        }
        if($k eq "SPADA_ORG" && $org) {
            $ENV{$k} = $org;
        } elsif($k eq "SPADA_FAS" && $f_fas) {
            $ENV{$k} = $f_fas;
        } elsif($k eq "SPADA_GFF" && $f_gff) {
            $ENV{$k} = $f_gff;
        } elsif($k eq "evalue" && $cutoff_e) {
            $ENV{$k} = $cutoff_e;
        } elsif($k eq "method") {
            $ENV{$k} = { map {$_=>0} split(";", $v) };
        } else {
            $ENV{$k} = $v;
        }
    }

    my @keys = qw/SPADA_SOURCE SPADA_PROFILE SPADA_DATA SPADA_ORG SPADA_FAS
        ClustalO GeneWise SplicePredictor SignalP HMMER Augustus/;
    for my $key (@keys) {
        exists $ENV{$key} || die "$key not defined\n";
    }
    
    # check availability of called programs
    for my $soft (keys %{$ENV{"method"}}) {
        my @f_bins;
        if($soft eq "Augustus_evidence") {
            push @f_bins, $ENV{"Augustus"}."/bin/augustus";
        } elsif($soft eq "Augustus_de_novo") {
            push @f_bins, $ENV{"Augustus"}."/bin/augustus";
        } elsif($soft eq "SPADA") {
            push @f_bins, $ENV{"GeneWise"}."/bin/genewise";
            push @f_bins, $ENV{"SplicePredictor"}."/bin/SplicePredictor";
        } elsif($soft eq "GeneMark") {
            push @f_bins, $ENV{"GeneMark"}."/gmhmme3";
        } elsif($soft eq "GlimmerHMM") {
            push @f_bins, $ENV{"GlimmerHMM"}."/bin/glimmerhmm";
        } elsif($soft eq "GeneID") {
            push @f_bins, $ENV{"GeneID"}."/bin/geneid";
        }

        my $tag = 1;
        for my $f_bin (@f_bins) { $tag = 0 unless -s $f_bin; }
        if($tag == 1) {
            printf "\twill run %s\n", $soft;
            $ENV{"method"}->{$soft} = 1;
        } else {
            printf "\twill NOT run %s\n", $soft;
        }
    }

    $ENV{"TMP_DIR"} = $ENV{"SPADA_DATA"}."/".$ENV{"SPADA_ORG"};
    make_path($ENV{"TMP_DIR"});
    
    push @INC, $ENV{"SPADA_SOURCE"};
    $ENV{'PATH'} = join(":", $ENV{"SPADA_SOURCE"}, $ENV{'PATH'});
}

1;
__END__
