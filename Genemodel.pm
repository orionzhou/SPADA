package Genemodel;
use strict; 
use Common; 
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/compare_models compare_2_model
    compare_models_2
    extract_compatible_models/;

sub get_ovlp_models { # opt_strand = 1 (strand sensitive); 2 (strand Insensitive)
    my ($locQ, $chr, $strand, $t, $opt_strand) = rearrange(['loc', 'chr', 'strand', 'tgt', 'opt'], @_);
    $opt_strand ||= 1;
    die "no strand provided\n" if $opt_strand == 1 && !defined($strand);
    $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
    my ($beg, $end) = ($locQ->[0]->[0], $locQ->[-1]->[1]);

    my $idx_chr = first_index {$_ eq 'chr'} $t->header;
    my $idx_beg = first_index {$_ eq 'beg'} $t->header;
    my $idx_end = first_index {$_ eq 'end'} $t->header;
    my $t2 = $t->match_pattern("\$_->[$idx_chr] eq '$chr' && \$_->[$idx_beg] <= $end && \$_->[$idx_end] >= $beg");
    if($opt_strand == 1) {
        my $idx_str = first_index {$_ eq 'strand'} $t->header;
        $t2 = $t2->match_pattern("\$_->[$idx_str] eq '$strand'");
    }
    return $t2;
}


sub exonCmp {
    my ($loc1, $loc2) = @_;
    my ($h1, $h2);
    for my $i (0..@$loc1-1) {
        my ($b1, $e1) = @{$loc1->[$i]};
        my $ho;
        for my $j (0..@$loc2-1) {
            my ($b2, $e2) = @{$loc2->[$j]};
            my $lenO = max(0, min($e1, $e2) - max($b1, $b2) + 1);
            $ho->{$j} = $lenO;
        }
        my @js = sort {$ho->{$b} <=> $ho->{$a}} keys(%$ho);
        my $j = $js[0];
        if($ho->{$j} == 0) {
            $h1->{$i} = 2;
        } elsif($loc2->[$j]->[0] == $b1 && $loc2->[$j]->[1] == $e1) {
            $h1->{$i} = 1;
            $h2->{$j} = 1;
        } else {
            $h1->{$i} = 4;
            $h2->{$j} = 4;
        }
    }
    for my $j (0..@$loc2-1) {
        $h2->{$j} ||= 3;
    }
    my $exonO = scalar(grep {$h1->{$_} == 1} keys(%$h1));
    my $exon1 = scalar(grep {$h1->{$_} == 2} keys(%$h1));
    my $exon2 = scalar(grep {$h2->{$_} == 3} keys(%$h2));
    return ($exonO, $exon1, $exon2);
}
sub compare_2_model {
    my ($locQ, $phaseStrQ, $locT, $phaseStrT, $srd) = @_;
    my ($lenQ, $lenT) = map {locAryLen($_)} ($locQ, $locT);
    my ($refO, $sPosAry1, $sPosAry2, $lenO, $len1, $len2) = posCmp($locQ, $locT);
    my ($exonO, $exon1, $exon2) = exonCmp($locQ, $locT);
    
    my $tag;
    if($len1 == 0 && $len2 == 0) {
        $tag = 1;
    } elsif($lenO == 0) {
        $tag = 8;
    } else {
        my @locTs = sort {$a->[0] <=> $b->[0]} @$locT;
        @locTs = reverse @locTs if $srd =~ /^\-1?$/;
        my @locQs = sort {$a->[0] <=> $b->[0]} @$locQ;
        @locQs = reverse @locQs if $srd =~ /^\-1?$/;
        
        $tag = 2;
        for (sort {$a->[0] <=> $b->[0]} @$refO) {
            my ($begO, $endO) = @$_;
      
            my $idxQ = first_index {$_->[0] <= $begO && $_->[1] >= $endO} @locQs;
            die "Overlap[$begO-$endO] not within locQ\n".Dumper($locQ) if $idxQ == -1;
            my ($begQ, $endQ) = @{$locQs[$idxQ]};
            my $phaseQ = [split(",", $phaseStrQ)]->[$idxQ];
            die "no phase for locQ: $phaseStrQ\n".Dumper($locQ) if !defined($phaseQ);
        
            my $idxT = first_index {$_->[0] <= $begO && $_->[1] >= $endO} @locTs;
            die "Overlap[$begO-$endO] not within locT\n".Dumper($locT) if $idxT == -1;
            my ($begT, $endT) = @{$locTs[$idxT]};
            my $phaseT = [split(",", $phaseStrT)]->[$idxT];
            die "no phase for locT: $phaseStrT\n".Dumper($locT) if !defined($phaseT);
            
            my ($pQ, $pT);
            if($srd =~ /^[\+1]$/) {
                $pQ = ($phaseQ - ($begO - $begQ)) % 3;
                $pT = ($phaseT - ($begO - $begT)) % 3;
            } else {
                die "unknown strand: $srd\n" unless $srd =~ /^\-1?$/;
                $pQ = ($phaseQ - ($endQ - $endO)) % 3;
                $pT = ($phaseT - ($endT - $endO)) % 3;
            }
            $tag = 7 if($pQ != $pT);
        }
    }
    return ($tag, $lenO, $len1, $len2, $exonO, $exon1, $exon2);
}
sub compare_models {
    my ($f_qry, $f_tgt, $fo) = @_;
    my $tq = readTable(-in=>$f_qry, -header=>1);
    my $tt = readTable(-in=>$f_tgt, -header=>1);
  
    open(FH, ">$fo");
    print FH join("\t", qw/id gene tag lenO len1 len2 exonO exon1 exon2/)."\n";
    for my $i (0..$tq->nofRow-1) {
        my ($id, $chr, $strand, $locQStr, $phaseQ) = map {$tq->elm($i, $_)} qw/id chr strand locC phase/;
        my $locQ = locStr2Ary($locQStr);
        if(!$locQStr) {
            $locQ = locStr2Ary($tq->elm($i, "locE"));
            $phaseQ = join(",", getPhase($locQ, $strand));
        }

        my $t2 = get_ovlp_models(-loc=>$locQ, -chr=>$chr, -strand=>$strand, -tgt=>$tt, -opt=>1);
        my @stats;
        if($t2->nofRow == 0) {
            push @stats, ["", 9, 0, locAryLen($locQ), 0, 0, scalar(@$locQ), 0];
        } else {
            for my $j (0..$t2->nofRow-1) {
                my ($gene, $locTStr, $phaseT) = map {$t2->elm($j, $_)} qw/id locC phase/;
                my $locT = locStr2Ary($locTStr);
                if(!$locTStr) {
                    $locT = locStr2Ary($t2->elm($j, "locE"));
                    $phaseT = join(",", getPhase($locT, $strand));
                }
                my ($tag, $lenO, $len1, $len2, $exonO, $exon1, $exon2) =
                    compare_2_model($locQ, $phaseQ, $locT, $phaseT, $strand);
                push @stats, [$gene, $tag, $lenO, $len1, $len2, $exonO, $exon1, $exon2];
            }
        }
        @stats = sort {$a->[1]<=>$b->[1] || $b->[2]<=>$a->[2] || $a->[3]<=>$b->[3] || $a->[4]<=>$b->[4]} @stats;
        print FH join("\t", $id, @{$stats[0]})."\n";
        printf "  comparing gene models... %5d out of %d done\r", $i+1, $tq->nofRow;
    }
    print "\n";
    close FH;
}

sub compare_models_2 {
    my ($fi, $f_gtb, $fo) = rearrange(['in', 'gtb', 'out'], @_);
    my $ti = readTable(-in=>$fi, -header=>1);
    my $tg = readTable(-in=>$f_gtb, -header=>1);
  
    open(FH, ">$fo");
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
        printf "  comparing gene models... %5d out of %d done\r", $i+1, $ti->nofRow;
    }
    print "\n";
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
        @locTs = reverse @locTs if is_opposite_strands($strand, "+");
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
    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $h = { map {$tg->elm($_, "id") => $tg->rowRef($_)} (0..$tg->nofRow-1) };
    
    open(FH, ">$fo");
    print FH join("\t", $tg->header)."\n";
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
