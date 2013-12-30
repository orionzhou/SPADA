#!/usr/bin/perl -w
=pod BEGIN
  
=head1 NAME
  
  gfffix.pl: format a GFF3 file into standard GFF3 format.

=head1 SYNOPSIS
  
 gfffix.pl [-help] -p <option> -i <input-file> -o <output-file>

 Options:
    -help    brief help message
    -in      input file
    -out     output file
    -opt     input GFF file format option 

=cut

use strict;
use FindBin;
use lib "$FindBin::Bin";

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path/;
use Common;
use Gff;

my $help_flag;
my ($opt, $fi, $fo);
my %options = (
    "help|h" => \$help_flag,
    "in|i=s" => \$fi,
    "out|o=s" => \$fo
    "opt|p=s" => \$opt,
);

sub format_gff_tair {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") || die "cannot read $fi\n";
    open(FHO, ">$fo") || die "cannot write to $fo\n";
    print FHO "##gff-version 3\n";
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        my $type = $ps[2];
        if($type =~ /^(protein)|(chromosome)|(transposon_fragment)|(transposable_element)$/) {
            next;
        } elsif($type eq "pseudogenic_transcript") {
            $ps[2] = "mRNA";
        } elsif($type eq "pseudogenic_exon") {
            $ps[2] = "exon";
        } elsif($type eq "CDS") {
            $ps[8] =~ /Parent=([\w\.]+)/;
            $ps[8] = "Parent=$1";
        }
        print FHO join("\t", @ps)."\n";
    }
    close FHI;
    close FHO;
}
sub format_gff_phytozome {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") or die "cannot open $fi for reading\n";
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO "##gff-version 3\n";
    my $h;
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        my $ht = parse_gff_tags($ps[8]);
        if($ps[2] eq "mRNA") {
            my ($id, $name) = ($ht->{"ID"}, $ht->{"Name"});
            $h->{$id} = $ht->{"Name"};
            $ht->{"ID"} = $name;
            delete $ht->{"Name"};
        } elsif($ps[2] ne "gene") {
            my $pa = $ht->{"Parent"};
            die "cannot find Name for ID[$pa]\n" unless exists $h->{$pa};
            $ht->{"Parent"} = $h->{$pa};
        }
        my $tagStr = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
        $ps[8] = $tagStr;
        print FHO join("\t", @ps)."\n";
    }
    close FHI;
    close FHO;
}
sub format_gff_ensembl {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") || die "cannot read $fi\n";
    open(FHO, ">$fo") || die "cannot write to $fo\n";
    print FHO "##gff-version 3\n";
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        my ($seqid, $type) = @ps[0,2];
        next if $type eq "chromosome";
        if($seqid =~ /^\d+$/) {
            $ps[0] = "chr$seqid";
        } elsif($seqid eq "UNKNOWN") {
            $ps[0] = "chrU";
        }
        print FHO join("\t", @ps)."\n";
    }
    close FHI;
    close FHO;
}
sub format_gff_cufflinks {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") || die "cannot read $fi\n";
    open(FHO, ">$fo") || die "cannot write to $fo\n";
    print FHO "##gff-version 3\n";
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        if($ps[8] =~ /\"(.+)\"/) {
            my $str = $1;
            my $offset = $-[1];
            $str =~ s/\=/\:/g;
            $str =~ s/\;/\|/g;
            substr($ps[8], $offset, length($str), $str);
        }
        my $ht = parse_gff_tags($ps[8]);
        if(exists $ht->{"Parent"}) {
            my @parents = split(",", $ht->{"Parent"});
            if(@parents > 1) {
                for my $i (0..$#parents) {
                    my $htn = { map {$_=>$ht->{$_}} keys(%$ht) };
                    $htn->{"ID"} = $ht->{"ID"}."_".($i+1) if exists $htn->{"ID"};
                    $htn->{"Parent"} = $parents[$i];
                    my $note = join(";", map {$_."=".$htn->{$_}} keys(%$htn))."\n";
                    print FHO join("\t", @ps[0..7], $note)."\n";
                }
            } else {
                print FHO join("\t", @ps)."\n";
            }
        } else {
            print FHO join("\t", @ps)."\n";
        }
    }
    close FHI;
    close FHO;
}
sub format_gff_jcvi {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") || die "cannot read $fi\n";
    open(FHO, ">$fo") || die "cannot write to $fo\n";
    print FHO "##gff-version 3\n";
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        $ps[2] = "transposable_element_gene" if $ps[2] eq "transposable_element";
        my $ht = parse_gff_tags($ps[8]);
        if(exists($ht->{"conf_class"})) {
            $ht->{"Note"} = sprintf "[%s]%s", $ht->{"conf_class"}, $ht->{"Note"};
            $ps[8] = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
        }
        print FHO join("\t", @ps)."\n";
    }
    close FHI;
    close FHO;
}
sub format_gff_lj {
    my ($fi, $fo) = @_;
    open(FHI, "<$fi") || die "cannot read $fi\n";
    open(FHO, ">$fo") || die "cannot write to $fo\n";
    print FHO "##gff-version 3\n";
    while(<FHI>) {
        chomp;
        next if !$_ || /^\#/;
        my @ps = split "\t";
        my $ht = parse_gff_tags($ps[8]);
        for my $k (keys(%$ht)) {
            my $v = $ht->{$k};
            my $l = length($v);
            if(substr($v, 0, 1) eq "\"" && substr($v, $l-1, 1) eq "\"") {
                $ht->{$k} = substr($v, 1, $l-2) ;
                $ps[8] = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
            }
        }
        next if $ps[0] !~ /^chr/i || $ps[2] eq "tRNA" || (exists $ht->{"Note"} && $ht->{"Note"} =~ /tRNA/);
        
        if($ps[2] eq "gene") {
            print FHO join("\t", @ps)."\n";
            
            $ps[2] = "mRNA";
            $ht->{"Parent"} = $ht->{"ID"};
            $ht->{"ID"} =~ s/gene/mRNA/g;
            $ps[8] = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
            print FHO join("\t", @ps)."\n";
        } else {
            $ht->{"Parent"} =~ s/gene/mRNA/g;
            $ps[8] = join(";", map {$_."=".$ht->{$_}} keys(%$ht));
            print FHO join("\t", @ps)."\n";
        }
    }
    close FHI;
    close FHO;
}

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag || !$opt || !$fi || !$fo;
pod2usage({-message => "$fi is not there\n", -exitval => 2}) if ! -s $fi;

if($opt eq "tair") {
    format_gff_tair($fi, $fo);
} elsif($opt eq "phytozome") {
    format_gff_phytozome($fi, $fo);
} elsif($opt eq "ensembl") {
    format_gff_ensembl($fi, $fo);
} elsif($opt eq "cufflinks") {
    format_gff_cufflinks($fi, $fo);
} elsif($opt eq "jcvi") {
    format_gff_jcvi($fi, $fo);
} elsif($opt eq "lj") {
    format_gff_lj($fi, $fo);
} else {
    print "option [$opt] not supported\n";
    pod2usage(1, -exitval=>2);
}




