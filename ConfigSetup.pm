package ConfigSetup;
use strict; 
use File::Path qw/make_path remove_tree/;
use Common; 
use Data::Dumper;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/config_setup config_setup_simple/; 

sub read_cfg_hash {
    my ($f_cfg) = @_;
    my $h;
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
        $h->{$k} = $v;
    }
    return $h;
}
sub config_setup_simple {
    my ($f_cfg, $dir_hmm) = @_;
   
    print "=====  setting up environment variables  =====\n";
    my $h = read_cfg_hash($f_cfg);
    for my $k (keys %$h) { $ENV{$k} = $h->{$k}; }

    $ENV{"SPADA_HMM_DIR"} = $dir_hmm if defined $dir_hmm;
    make_path($ENV{"SPADA_HMM_DIR"}) if ! -d $ENV{"SPADA_HMM_DIR"};
    $ENV{"TMP_DIR"} = $ENV{"SPADA_HMM_DIR"};

    my @keys = qw/ClustalO trimAl HMMER/;
    for my $key (@keys) {
        exists $ENV{$key} || die "$key not defined\n";
    }
} 

sub config_setup {
    my ($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, $org, $sp, $e, $methods) = @_;
   
    print "=====  setting up environment variables  =====\n";
    my $h = read_cfg_hash($f_cfg);
    for my $k (keys %$h) { $ENV{$k} = $h->{$k}; }

    $ENV{"SPADA_OUT_DIR"} = $dir if defined $dir;
    $ENV{"SPADA_HMM_DIR"} = $dir_hmm if defined $dir_hmm;
    $ENV{"SPADA_FAS"} = $f_fas if defined $f_fas;
    $ENV{"SPADA_GFF"} = "" if defined $f_fas;
    $ENV{"SPADA_GFF"} = $f_gff if defined $f_gff;
    $ENV{"SPADA_ORG"} = $org if defined $org;
    $ENV{"evalue"} = $e if defined $e;
    $ENV{"eval_sp"} = $sp if defined $sp;
    $ENV{"methods"} = $methods if defined $methods;

    my @keys = qw/SPADA_SRC_DIR SPADA_OUT_DIR SPADA_HMM_DIR SPADA_ORG SPADA_FAS
        ClustalO SignalP HMMER/;
    for my $key (@keys) {
        exists $ENV{$key} || die "$key not defined\n";
    }

    $dir = $ENV{"SPADA_OUT_DIR"};
    $dir_hmm = $ENV{"SPADA_HMM_DIR"};
    $f_fas = $ENV{"SPADA_FAS"};
    $f_gff = $ENV{"SPADA_GFF"};
    
    printf "  using %s matrix\n", $ENV{"SPADA_ORG"};

    # check HMM / profile directory
    -s "$dir_hmm/21_all.hmm" || die "$dir_hmm/21_all.hmm is not there\n";

    # check availability of called programs
    $ENV{"method"} = { map {$_=>{}} split(";", $ENV{"method"}) };
    for my $soft (keys %{$ENV{"method"}}) {
        my $hb = $ENV{"method"}->{$soft};
        if($soft eq "Augustus_evidence") {
            $hb->{"Augustus"} = "bin/augustus";
        } elsif($soft eq "Augustus_de_novo") {
            $hb->{"Augustus"} = "bin/augustus";
        } elsif($soft eq "GeneWise_SplicePredictor") {
            $hb->{"GeneWise"} = "bin/genewise";
            $hb->{"SplicePredictor"} = "bin/SplicePredictor";
        } elsif($soft eq "GeneMark") {
            $hb->{"GeneMark"} = "gmhmme3";
        } elsif($soft eq "GlimmerHMM") {
            $hb->{"GlimmerHMM"} = "bin/glimmerhmm";
        } elsif($soft eq "GeneID") {
            $hb->{"GeneID"} = "bin/geneid";
        }
        
        for my $key (keys %$hb) {
            exists $ENV{$key} || die "$key not defined\n";
            my $fb = $ENV{$key}."/".$hb->{$key};
            -s $fb || die "$key: $fb is not there\n";
        }
        printf "\twill run %s\n", $soft;
    }

    
    make_path($ENV{"SPADA_OUT_DIR"}) if ! -d $ENV{"SPADA_OUT_DIR"};
    $ENV{"TMP_DIR"} = $ENV{"SPADA_OUT_DIR"};
    
    $ENV{"BLOSUM62"} = $ENV{"SPADA_SRC_DIR"}."/BLOSUM62";
    push @INC, $ENV{"SPADA_SRC_DIR"};
    $ENV{'PATH'} = join(":", $ENV{"SPADA_SRC_DIR"}, $ENV{'PATH'});
}

1;
__END__
