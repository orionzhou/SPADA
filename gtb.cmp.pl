#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb.cmp.pl - compares a Qry-Gtb with a Tgt-Gtb 

=head1 SYNOPSIS
  
  gtb.cmp.pl [-help] [-qry query-Gtb] [-tgt target-Gtb] [-out cmp-result]

  Options:
    -h (--help)  brief help message
    -q (--qry)   query Gtb file
    -t (--tgt)   target Gtb file
    -o (--out)   cmp result output 

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Common;
use Gtb;
use Location;

my ($fq, $ft, $fo) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "qry|q=s" => \$fq,
  "tgt|t=s" => \$ft,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fq || !$ft || !$fo;

my $tq = readTable(-in => $fq, -header => 1);
my $tt = readTable(-in => $ft, -header => 1);

my $h;
for my $i (0..$tt->lastRow) {
  my ($id, $par, $chr, $beg, $end, $srd, $eloc, $iloc, $cloc) = $tt->row($i);
  $h->{$chr}->{$srd} ||= [];
  push @{$h->{$chr}->{$srd}}, [$beg, $end, $i];
}
  
open(my $fho, ">$fo") || die "cannot write $fo\n";
print $fho join("\t", qw/qid tid tag olen qlen tlen oexon qexon texon/)."\n";
for my $i (0..$tq->nofRow-1) {
  my ($qid, $chr, $qb, $qe, $srd, $qlocS, $qphase) = 
    map {$tq->elm($i, $_)} qw/id chr beg end srd cloc phase/;
  $qlocS ne "" || die "no CDS for $qid\n";
  my $rqloc = locStr2Ary($qlocS);
  my $qloc = $srd eq "+" ? [ map {[$qb+$_->[0]-1, $qb+$_->[1]-1]} @$rqloc ]
    : [ map {[$qe-$_->[1]+1, $qe-$_->[0]+1]} @$rqloc ];
  
  my $idxs = get_ovlp_idxs($chr, $srd, $qloc, $h);
  my $t2 = $tt->subTable($idxs);
#  my $t2 = get_ovlp_gtb(-loc=>$qloc, -chr=>$chr, -srd=>$srd, -tgt=>$tt, -opt=>1);
  my @stats;
  if($t2->nofRow == 0) {
    push @stats, ["", 9, 0, locAryLen($qloc), 0, 0, scalar(@$qloc), 0];
  } else {
    for my $j (0..$t2->nofRow-1) {
      my ($tid, $tb, $te, $tlocS, $tphase) = 
        map {$t2->elm($j, $_)} qw/id beg end cloc phase/;
      my $rtloc = locStr2Ary($tlocS);
      my $tloc = $srd eq "+" ? [ map {[$tb+$_->[0]-1, $tb+$_->[1]-1]} @$rtloc ]
        : [ map {[$te-$_->[1]+1, $te-$_->[0]+1]} @$rtloc ];
      my ($tag, $lenO, $len1, $len2, $exonO, $exon1, $exon2) =
        cmp_cds($qloc, $qphase, $tloc, $tphase, $srd);
      push @stats, [$tid, $tag, $lenO, $len1, $len2, $exonO, $exon1, $exon2];
    }
  }
  @stats = sort {$a->[1]<=>$b->[1] || $b->[2]<=>$a->[2] || $a->[3]<=>$b->[3] || $a->[4]<=>$b->[4]} @stats;
  print $fho join("\t", $qid, @{$stats[0]})."\n";
#    printf "  comparing gene models [%5d / %5d]\n", $i+1, $tq->nofRow;
}
close $fho;

sub get_ovlp_idxs {
  my ($chr, $srd, $qloc, $h) = @_;
  my $loc1 = [ map {[@{$qloc->[$_]}, "q".$_]} (0..@$qloc-1) ];
  my $loc2 = exists $h->{$chr}->{$srd} ? $h->{$chr}->{$srd} : [];
  my $ref = posSplit([@$loc1, @$loc2]);
  my @idxs;
  for (@$ref) {
    my ($beg, $end, $idxs) = @$_;
    my @idxs1 = grep /^\d+$/, @$idxs;
    my @idxs2 = grep /^q\d+$/, @$idxs;
    push @idxs, @idxs1 if @idxs1 > 0 && @idxs2 > 0;
  }
  return \@idxs;
}

#model_eval($f_tmp, $ft, $fc);
#my ($sn_nt, $sp_nt, $sn_ex, $sp_ex) = get_sn_sp($fc);
#printf "sn_nt=%.03f sp_nt=%.03f sn_ex=%.03f sp_ex=%.03f\n", $sn_nt, $sp_nt, $sn_ex, $sp_ex;
#system("rm $f_tmp");

__END__
