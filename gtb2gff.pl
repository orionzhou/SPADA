#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2gff.pl - convert a Gtb file to GFF3 format

=head1 SYNOPSIS
  
  gtb2gff.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file (Gtb)
      -out    output file (Gff)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Rna;
use List::Util qw/min max sum/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi;

my ($fhi, $fho);
if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-in => $fi, -header => 1);
$t->sort("par", 1, 0, "id", 1, 0);
my $ref = group($t->colRef("par"));
my @gids = sort {$ref->{$a}->[0] <=> $ref->{$b}->[0]} keys %$ref;
printf "%d genes : %d models\n", scalar(@gids), $t->nofRow;

print $fho "##gff-version 3\n";
my ($cntG, $cntR) = (1, 1);
for my $gid (@gids) {
  my ($idxB, $cnt) = @{$ref->{$gid}};

  my @rnas;
  for my $i ($idxB..$idxB+$cnt-1) {
    push @rnas, Rna->new(-gtb => $t->rowRef($i));
  }

  my $seqid = $rnas[0]->seqid;
  my $beg = min(map {$_->beg} @rnas);
  my $end = max(map {$_->end} @rnas);
  my $srd = $rnas[0]->strand;
  my $tag = "ID=".$rnas[0]->parent;
  print $fho join("\t", $seqid, '.', 'gene', $beg, $end, '.', $srd, '.', $tag)."\n";;

  for my $rna (@rnas) {
    print $fho $rna->to_gff()."\n";
  }
}
close $fho;



__END__
