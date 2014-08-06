#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gff2gtb.pl - convert a Gff file to Gtb file and validates

=head1 SYNOPSIS
  
  gff2gtb.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (Gff)
    -o (--out)    output file (Gtb)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gff;
use Gtb;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my ($fhi, $fho);
open ($fhi, "<$fi") || die "cannot read $fi\n";

if ($fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

print $fho join("\t", @HEAD_GTB)."\n";

my $it = parse_gff($fhi);
my ($cntR, $cntG) = (1, 1);
while(my $gene = $it->()) {
  for my $rna ($gene->get_rna) {
    print $fho $rna->to_gtb()."\n";
    printf "%5d RNA | %5d gene...\n", $cntR, $cntG if $cntR % 1000 == 0;
    $cntR ++;
  }
  $cntG ++;
}
close $fhi;
close $fho;

__END__
