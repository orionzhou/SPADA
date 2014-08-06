package CompareModel;
use strict; 
use Location; 
use Common; 
use Data::Dumper;
use Gtb;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/model_eval get_sn_sp
    extract_compatible_models/;

sub model_eval {
  my ($f_cmp, $f_tgt, $fo) = @_;
  my $tc = read_gtb($f_cmp);
  my $tt = read_gtb($f_tgt);

  my $ht;
  for my $i (0..$tt->nofRow-1) {
    my ($idT, $locS) = map {$tt->elm($i, $_)} qw/id locC/;
    my $loc = locStr2Ary($locS);
    my $len = locAryLen($loc);
    my $n_cds = @$loc;
    $ht->{$idT} = [0, $len, $n_cds];
  }

  open(FH, ">$fo") || die "cannot open $fo for writing\n";
  print FH join("\t", qw/idQ tag idT lenTP lenFP lenFN exonTP exonFP exonFN/)."\n";
  for my $i (0..$tc->nofRow-1) {
    my ($idQ, $idT, $tag, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN) = $tc->row($i);
    $tag = 5 if $tag == 7 || $tag == 8 || ($tag==2 && $lenFP+$lenFN>30);

    print FH join("\t", $idQ, $tag, $idT, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN)."\n";
    $ht->{$idT}->[0] ++ if $idT ne "";
  }
  for my $idT (keys(%$ht)) {
    my ($cnt, $len, $n_cds) = @{$ht->{$idT}};
    if($cnt == 0) {
      print FH join("\t", '', 10, $idT, 0, 0, $len, 0, 0, $n_cds)."\n";
    } elsif($cnt > 1) {
      print "  $idT hit $cnt times\n";
    }
  }
  close FH;
}

sub get_sn_sp {
  my ($fi) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  my $sn_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFN")) );
  my $sp_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFP")) );
  my $sn_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFN")) );
  my $sp_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFP")) );
  return ($sn_nt, $sp_nt, $sn_ex, $sp_ex);
}

sub compare_models_2 {
  my ($fi, $f_gtb, $fo) = rearrange(['in', 'gtb', 'out'], @_);
  my $ti = readTable(-in=>$fi, -header=>1);
  my $tg = read_gtb($f_gtb);

  open(FH, ">$fo") || die "cannot write $fo\n";
  print FH join("\t", qw/id gene tag lenC lenI len5 len3 lenO/)."\n";
  for my $i (0..$ti->nofRow-1) {
    my ($id, $chr, $srd, $locS) = map {$ti->elm($i, $_)} qw/id chr strand loc/;
    my $locQ = locStr2Ary($locS);
    my $phaseQ = join(",", getPhase($locQ, $srd));

    my $t2 = get_ovlp_models(-loc=>$locQ, -chr=>$chr, -strand=>$srd, -tgt=>$tg, -opt=>1);
    my @stats;
    if($t2->nofRow == 0) {
        push @stats, ["", 9, ("") x 4, locAryLen($locQ), 0];
    } else {
        for my $j (0..$t2->nofRow-1) {
            my ($gene, $srdT, $locTStr, $phaseT) = map {$t2->elm($j, $_)} qw/id strand locC phase/;
            my $locT = locStr2Ary($locTStr);
            push @stats, [$gene, qry_one_model($chr, $locQ, $srd, $phaseQ, $t2->rowRef($j))];
        }
    }
    @stats = sort {$a->[1]<=>$b->[1] || $b->[2]<=>$a->[2] || $a->[6]<=>$b->[6] || $a->[7]<=>$b->[7]} @stats;
    my ($gene, $tag, $lenC, $lenI, $len5, $len3, $lenO, $lenM) = @{$stats[0]};
    print FH join("\t", $id, $gene, $tag, $lenC, $lenI, $len5, $len3, $lenO)."\n";
    printf "  comparing gene models [%5d / %5d]\n", $i+1, $ti->nofRow if ($i+1) % 1000 == 0;
  }
  close FH;
}
sub qry_one_model {
  my ($chr, $locQ, $strand, $phaseStrQ, $row) = @_;
  $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
  my ($beg, $end) = ($locQ->[0]->[0], $locQ->[-1]->[1]);
  my $lenQ = locAryLen($locQ);
  
  my ($begG, $endG, $locSC, $locSI, $locS5, $locS3, $phaseStrT) = @$row[3,4,8,7,9,10,11];
  my ($locC, $locI, $loc5, $loc3) = map {locStr2Ary($_)} ($locSC, $locSI, $locS5, $locS3);
  
  my ($oPosAry, $sPosAry1, $sPosAry2, $oLen, $sLen1, $sLen2) = posCmp($locQ, $locC);
  my ($refC, $lenC, $lenM) = ($oPosAry, $oLen, $sLen2);
  my $tag_exact = ($sLen1 == 0 && $sLen2 == 0) ? 1 : 0;
  
  my ($refI, $lenI) = posOvlp($locQ, $locI);
  my ($ref5, $len5) = posOvlp($locQ, $loc5);
  my ($ref3, $len3) = posOvlp($locQ, $loc3);
  
  my $locM = [[$begG, $endG]];
  ($oPosAry, $sPosAry1, $sPosAry2, $oLen, $sLen1, $sLen2) = posCmp($locQ, $locM);
  my $lenO = $sLen1;
  die "length do not add up to $lenQ\n".join("\t", @$row, "\n", $lenC, $lenI, $len5, $len3, $lenO)."\n" unless $lenC+$lenI+$len5+$len3+$lenO == $lenQ;

  my $tag = "";
  if($lenC == 0) {
    $tag = 8;
  } else {
    $refC = [ sort {$a->[0] <=> $b->[0]} @$refC ];
    my $origin = $strand eq "-" ? $refC->[-1]->[-1] : $refC->[0]->[0];

    my @locTs = sort {$a->[0] <=> $b->[0]} @$locC;
    @locTs = reverse @locTs if is_revsrd($strand, "+");
    my $idx = first_index {$origin >= $_->[0] && $origin <= $_->[1]} @locTs;
    die "origin[$origin] not within CDS\n" if $idx == -1;
    my ($begT, $endT) = @{$locTs[$idx]};
    my @phasesT = split(",", $phaseStrT);
    my $phaseT = $phasesT[$idx];
    printf "no phase for %dth CDS: %s\n".join("\t", @$row)."\n", $idx+1, $phaseStrT if !defined($phaseT);
    my $disT = ($strand eq "-") ? $endT - $phaseT - $origin + 3 : $origin - $begT - $phaseT + 3;

    my @locQs = sort {$a->[0] <=> $b->[0]} @$locQ;
    @locQs = reverse @locQs if $strand eq "-";
    $idx = first_index {$origin >= $_->[0] && $origin <= $_->[1]} @locQs;
    die "origin[$origin] not within CDS\n" if $idx == -1;
    my ($begQ, $endQ) = @{$locQs[$idx]};
    my @phasesQ = split(",", $phaseStrQ);
    my $phaseQ = $phasesQ[$idx];
    printf "no phase for %dth CDS: %s\n".Dumper($locQ), $idx+1, $phaseStrQ if !defined($phaseQ);
    my $disQ = ($strand eq "-") ? $endQ - $phaseQ - $origin + 3 : $origin - $begQ - $phaseQ + 3;
    
    if($disQ % 3 == $disT % 3) {
      if($tag_exact == 1) {
        $tag = 0.5;
      } elsif($lenC + 9 >= $lenQ) {
        $tag = 1;
      } else {
        $tag = 2;  
      }
    } else {
      $tag = 7; 
    }
  }
  return ($tag, $lenC, $lenI, $len5, $len3, $lenO, $lenM);
}

sub extract_compatible_models {
  my ($fi, $f_gtb, $fo) = @_;
  my $ti = readTable(-in=>$fi, -header=>1);
  my $h = read_gtb_hash($f_gtb);
  
  open(FH, ">$fo") || die "cannot write $fo\n";
  my $cnt = 0;
  for my $i (0..$ti->nofRow-1) {
    my ($id, $gene, $tag) = map {$ti->elm($i, $_)} qw/id gene tag/;
    next unless $tag =~ /^[12]$/;
    die "cannot find gene $gene\n" unless exists $h->{$gene};
    my @row = @{$h->{$gene}};
    @row[0,1,-1] = ("$id.1", $id, $gene);

    my ($locCStr, $strand) = @row[8,5];
    my $locCAry = locObj2Ary(locStr2Obj($locCStr));
    $locCAry = [ sort {$a->[0] <=> $b->[0]} @$locCAry ];
    
    my @phases = split(",", $row[11]);
    if($phases[0] > 0) {
      if($strand == -1) {
        $locCAry->[-1]->[-1] -= $phases[0];
      } else {
        $locCAry->[0]->[0] += $phases[0];
      }
      $phases[0] = 0;
      $row[11] = join(",", @phases);
      $locCStr = locAry2Str($locCAry, $strand);
    }

    my ($begM, $endM) = ($locCAry->[0]->[0], $locCAry->[-1]->[-1]);
    my ($locIAry) = posDiff([[$begM, $endM]], $locCAry);
    my $locIStr = locAry2Str($locIAry, $strand);

    @row[3,4,6,7,8,9,10] = ($begM, $endM, $locCStr, $locIStr, $locCStr, '', '');

    print FH join("\t", @row)."\n";
    $cnt ++;
  }
  printf "%d compatible models extracted (out of %d)\n", $cnt, $ti->nofRow;
  close FH;
}
  


__END__
