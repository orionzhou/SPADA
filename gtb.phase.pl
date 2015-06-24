#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.phase.pl - check and fix phases for a Gtb file

=head1 SYNOPSIS
  
  gtb.phase.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (Gtb) file
    -o (--out)    output (Gtb) file
    -s (--seq)    refseq fasta

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gtb;
use Rna;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
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
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-inh => $fhi, -header => 1);
close $fhi;

print $fho join("\t", @HEAD_GTB)."\n";
my $cntF = 0;
for my $i (0..$t->lastRow) {
  my $rna = Rna->new(-gtb => $t->rowRef($i));
  $cntF += $rna->check_phase($fs);
  print $fho $rna->to_gtb()."\n";
#  printf "%5d: %5d non-0 phase fixed\n", $i+1, $cntF if ($i+1) % 1000 == 0;
}
close $fho;

__END__
