#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtbcheckphase.pl - check and fix a Gtb file

=head1 SYNOPSIS
  
  gtbcheckphase.pl [-help] [-in input-Gtb] [-seq refseq-fasta] [-out output-Gtb]

  Options:
      -help   brief help message
      -in     input Gtb file
      -out    output Gtb file
      -seq    refseq fasta

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fs;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $t = readTable(-inh=>$fhi, -header=>1);
close $fhi;

my $cntF = 0;
print $fho join("\t", @HEAD_GTB)."\n";
for my $i (0..$t->lastRow) {
    my $ts = $t->subTable([$i]);
    my $gene = Gene->new( -gtb=>$ts );
    my ($rna) = $gene->get_rna(); 
    $cntF += $rna->check_phase($fs);
    print $fho $rna->to_gtb()."\n";
    printf "%5d: %5d non-0 phase fixed\n", $i+1, $cntF if ($i+1) % 1000 == 0;
}
close $fho;

__END__