package Location;
use strict; 
use Data::Dumper;
use Clone qw/clone/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use Common;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/
    locStr2Ary locAry2Str locAryLen trimLoc cropLoc cropLoc_cds
    posOvlp posCmp posMerge posSplit posDiff posMergeDeep
    tiling
    coordTransform coordTransform_rough coordTransform_itv/;
@EXPORT_OK = qw//;

sub locStr2Ary {
    my ($locS) = @_;
    my @locA = ();
    return \@locA unless defined $locS;
    my @locs = split(",", $locS);
    for (@locs) {
        if(/^(\d+)\-(\d+)$/) {
            push @locA, [$1, $2];
        } else {
            die "malformat location string: $locS\n";
        }
    }
    return \@locA;
}
sub locAry2Str {
    my ($locA) = @_;
    $locA ||= [];
    my $locS = join(",", map {$_->[0]."-".$_->[1]} @$locA);
    return $locS;
}
sub locAryLen {
    my ($loc) = @_;
    my $len = 0;
    for (@$loc) {
        my ($b, $e) = @$_;
        ($b, $e) = ($e, $b) if $b > $e;
        $len += $e - $b + 1;
    }
    return $len;
}
sub trimLoc {
    my ($loc, $beg, $end) = @_;
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    $beg = $loc->[0]->[0] unless defined $beg;
    $end = $loc->[-1]->[1] unless defined $end;
    ($beg, $end) = ($end, $beg) if $beg > $end;
    my $locO = [];
    for (@$loc) {
        my ($b, $e) = @$_;
        ($b, $e) = ($e, $b) if $b > $e;
        my $begO = max($beg, $b);
        my $endO = min($e, $end);
        push @$locO, [$begO, $endO] if $begO <= $endO;
    }
    return $locO;
}
sub cropLoc {
    my ($loc, $len, $opt) = @_;
    die "unknown opt: $opt\n" if $opt !~ /^[12]$/;
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    $loc = [ reverse @$loc ] if $opt == 2;
    my ($locO, $locC) = ([], []);
    my $lenC = 0;
    for (@$loc) {
        my ($beg, $end) = @$_;
        my $lenT = $end - $beg + 1;
        if($lenC >= $len) {
            push @$locO, [$beg, $end];
        } elsif($lenC + $lenT > $len) {
            my ($begO, $endO, $begC, $endC);
            if($opt == 1) {
                $begO = $beg + $len - $lenC;
                $endO = $end;
                $begC = $beg;
                $endC = $begO - 1;
            } else {
                $begO = $beg;
                $endO = $end - ($len - $lenC);
                $begC = $endO + 1;
                $endC = $end;
            }
            $lenC = $len;
            push @$locO, [$begO, $endO];
            push @$locC, [$begC, $endC];
        } else {
            $lenC += $lenT;
            push @$locC, [$beg, $end];
        }
    }
    return ($locO, $locC);
}
sub cropLoc_cds {
    my ($loc, $srd, $phase) = @_;
    die "unknown strand: $srd\n" if $srd !~ /^[\+\-]$/;
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    $loc = [ reverse @$loc ] if $srd eq "-";
    my @phases = $phase ? split(",", $phase) : getPhase($loc, $srd);
    my $locO = [];
    for my $i (0..@$loc-1) {
        my ($beg, $end) = @{$loc->[$i]};
        my $phase = $phases[$i];
        if($end - $beg + 1 > $phase) {
            if($srd eq "+") {
                $beg += $phase;
            } else {
                $end -= $phase;
            }
            my $res = ($end - $beg + 1) % 3;
            if($end - $beg + 1 > $res) {
                if($srd eq "+") {
                    $end -= $res;
                } else {
                    $beg += $res;
                }
                push @$locO, [$beg, $end];
            }
        }
    }
    return $locO;
}

sub posMerge {
    my ($locI) = @_;
    return undef if !$locI;
    my @locs = @$locI;
    @locs = grep {@$_ >= 2} @locs;
    return undef if @locs == 0;
    @locs = map {[$locs[$_]->[0], $locs[$_]->[1], $_]} (0..$#locs) if @{$locs[0]} == 2;
    @locs = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @locs;
    my @refs;
    my ($prevS, $prevE, $hi);
    for my $i (0..$#locs) {
        my ($s, $e, $idx) = @{$locs[$i]};
        if($i == 0) {
            ($prevS, $prevE) = ($s, $e);
            $hi->{$idx} = 1;
        }
        if($prevE >= $s-1) {
            $prevE = $prevE<$e ? $e : $prevE;
            $hi->{$idx} ||= 0;
            $hi->{$idx} ++;
        } else {
            push @refs, [$prevS, $prevE, $hi];
            ($prevS, $prevE) = ($s, $e);
            $hi = {$idx => 1};
        }
    }
    push @refs, [$prevS, $prevE, $hi] if $prevE;
    @refs = map {[ $_->[0], $_->[1], [keys(%{$_->[2]})] ]} @refs;
    return \@refs;
}
sub posSplit {
    my ($locI) = @_;
    $locI = @{$locI->[0]} == 3 ? $locI : [map {[@{$locI->[$_]}, $_]} (0..@$locI-1)];
    $locI = [ sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @$locI ];
    my @ss = sort {$a<=>$b} uniq( map {@$_[0..1]} @$locI );
    my (@posSs, @posEs, @idxs);
    for my $i (0..$#ss) {
        my $s = $ss[$i];
        push @posSs, $s;
        push @posEs, $s;
        push @idxs, [];
        if($i != $#ss) {
            my $posNext = $ss[$i+1];
            if($posNext > $s+1) {
                push @posSs, $s + 1;
                push @posEs, $posNext - 1;
                push @idxs, [];
            }
        }
    }
    for my $i (0..@$locI-1) {
        my ($s1, $e1, $idx) = @{$locI->[$i]};
        my $eP = 0;
        my $tmpIdx = bsearch(\@posSs, $s1);
        for my $i ($tmpIdx..$#posSs) {
            my ($s2, $e2, $idxAry) = ($posSs[$i], $posEs[$i], $idxs[$i]);
            if($s1<=$s2 && $e1>=$e2) {
                push @$idxAry, $idx;
            }
            $eP = $e2;
            last if $eP > $e1;
        }
    }
    my $ref = {};
    my $eP = 0;
    for my $i (0..$#posSs) {
        my ($s1, $e1, $idxAry1) = ($posSs[$i], $posEs[$i], $idxs[$i]);
        next if $s1 <= $eP;
        $eP = $e1;
        for my $j ($i+1..$#posSs) {
            my ($s2, $e2, $idxAry2) = ($posSs[$j], $posEs[$j], $idxs[$j]);
            my ($aryS, $ary1, $ary2) = aryCmp($idxAry1, $idxAry2);
            if(@$ary1 == 0 && @$ary2 == 0) {
                $eP = $e2;
            } else{
                last;
            }
        }
        $ref->{$s1} = [$eP, [uniq(@$idxAry1)]] if @$idxAry1 > 0;
    }
    my $res = [ map {[$_, @{$ref->{$_}}]} sort {$a<=>$b} keys %$ref ];
    return $res;
}
sub posOvlp {
    my ($posAry1, $posAry2) = @_;
    my $rst = [];
    my $oLen = 0;
    for my $posR1 (@$posAry1) {
        my ($beg1, $end1) = @{$posR1}[0..1];
        ($beg1, $end1) = ($end1, $beg1) if $beg1 > $end1;
        for my $posR2 (@$posAry2) {
            my ($beg2, $end2) = @{$posR2}[0..1];
            ($beg2, $end2) = ($end2, $beg2) if $beg2 > $end2;
            next if $beg1 > $end2 || $beg2 > $end1;
            my $oBeg = $beg1 < $beg2 ? $beg2 : $beg1;
            my $oEnd   = $end1   > $end2   ? $end2   : $end1;
            push (@$rst, [$oBeg, $oEnd]);
            $oLen += $oEnd - $oBeg + 1;
        }
    }
    return ($rst, $oLen);
}
sub posDiff {
    my ($posGAry, $posLAry) = @_;
    my $posDAry = clone($posGAry);
    my $dLen = 0;
    for my $posL (@$posLAry) {
        my ($begL, $endL) = ($posL->[0], $posL->[1]);
        die "error in posDiff[$begL > $endL]\n".Dumper($posGAry).Dumper($posLAry) unless $begL <= $endL;
        my $aryIndex = -1;
        for my $posG (@$posDAry) {
            $aryIndex ++;
            next if !defined($posG->[0]);
            my ($begG, $endG) = ($posG->[0], $posG->[1]);
            die("Error: begG[$begG] > endG[$endG]\n") unless $begG <= $endG;
            if ($begG<=$begL && $endG>=$endL) {
                push (@$posDAry, [$begG, $begL-1]) if ($begG<$begL);
                push (@$posDAry, [$endL+1, $endG]) if ($endG>$endL);
                delete $posDAry->[$aryIndex];
            } elsif($endG<$begL || $begG>$endL) {
                next;
            } else {
                print Dumper($posGAry, $posLAry);
                die "[$begL - $endL] not in [$begG - $endG]\n";
            }
        }
    }
    my $ref = [];
    for my $posPair (@$posDAry) {
        next unless $posPair;
        next unless @$posPair;
        push @$ref, [ (@$posPair) ];
        $dLen += $posPair->[1] - $posPair->[0] + 1;
    }
    $ref = [ sort {$a->[0] <=> $b->[0]} @$ref ];
    return ($ref, $dLen);
}
sub tiling { #opt = 1(minimal) / 2(maximal)
    my ($locAry, $statAry, $opt) = @_;
    $opt ||= 1;
    my $r = posSplit($locAry);
    my (@refs, @locs);
    for (@$r) {
        my ($s, $e, $idxAry) = @$_;
        my $idxExt = 0;
        my $tmp = "";
        for my $i (0..@$idxAry-1) {
            my $stat = $statAry->[$idxAry->[$i]];
            $tmp = $stat if $tmp eq "";
            if(($opt == 1 && $stat < $tmp) || ($opt == 2 && $stat > $tmp)) {
                $idxExt = $i;
                $tmp = $stat;
            }
        }
        push @locs, [$s, $e, $idxAry->[$idxExt], $idxAry];
    }
    @locs = sort {$a->[0] <=> $b->[0]} @locs;
    
    my $prevE = 0;
    my $idxAry = [];
    for my $i (0..$#locs) {
        my ($s1, $e1, $idx1, $idxAry1) = @{$locs[$i]};
        next if $s1 <= $prevE;
        ($prevE, $idxAry) = ($e1, $idxAry1);
        for my $j ($i+1..$#locs) {
            my ($s2, $e2, $idx2, $idxAry2) = @{$locs[$j]};
            if($prevE == $s2-1 && $idx1 == $idx2) {
                $prevE = $e2;
                push @$idxAry, @$idxAry2;
            } else{
                last;
            }
        }
        push @refs, [$s1, $prevE, $idx1, [uniq(@$idxAry)]];
        #print "\t".join("\t", $start1, @{$loc2->{$start1}})."\n";
    }
    return \@refs;
}
sub posMergeDeep {
    my ($locAs) = @_;
    my $ref1 = [];
    for my $i (0..@$locAs-1) {
        my $locA = $locAs->[$i];
        push @$ref1, map {[@$_, $i]} @$locA;
    }
    my $ref2 = posMerge($ref1);
    my @idxss = map {$_->[2]} @$ref2;
    my @idxMap;
    for my $i (0..$#idxss) {
        my $idxs = $idxss[$i];
        push @idxMap, map {[$i, $_]} @$idxs;
    }
    @idxMap = sort {$a->[1] <=> $b->[1]} @idxMap;
    my $ref3 = group([ map {$_->[1]} @idxMap ]);
    for my $idxP (sort {$a<=>$b} keys %$ref3) {
        my ($idx, $cnt) = @{$ref3->{$idxP}};
        $ref3->{$idxP} = [ sort {$a<=>$b} map {$_->[0]} @idxMap[$idx..$idx+$cnt-1] ];
    }
    my $g = Graph::Undirected->new();
    for my $i (keys %$ref3) {
        my $js = $ref3->{$i};
        if(@$js > 1) {
            for my $j (@$js[1..@$js-1]) {
                $g->add_edge($js->[0], $j);
            }
        } else {
            $g->add_vertex($js->[0]);
        }
    }
    my @ccs = $g->connected_components();
    my ($ref4, $ref5) = ({}, {});
    for (0..$#ccs) {
        @$ref4{@{$ccs[$_]}} = ($_) x @{$ccs[$_]};
    }
    for my $i (keys %$ref3) {
        my @js = @{$ref3->{$i}};
        my @ks = uniq( map {$ref4->{$_}} @js );
        die "internal error at posMergeDeep\n" unless @ks == 1;
        $ref5->{$i} = $ks[0];
    }
    my $ref = [];
    for my $k (uniq(values %$ref4)) {
        my @is1 = sort {$ref2->[$a]->[0] <=> $ref2->[$b]->[0]} grep {$ref4->{$_} == $k} keys %$ref4;
        my @is2 = sort {$a<=>$b} grep {$ref5->{$_} == $k} keys %$ref5;
        push @$ref, [ [map {[@$_[0..1]]} @{$ref2}[@is1] ], \@is2 ];
    }
    return $ref;
}

sub posCmp {
    my ($posAryO1, $posAryO2) = @_;
    my ($posAry1, $posAry2) = (clone($posAryO1), clone($posAryO2));
    ($posAry1, $posAry2) = (posMerge($posAry1), posMerge($posAry2));
    my ($oPosAry, $oLen) = posOvlp($posAry1, $posAry2);
    my ($sPosAry1, $sLen1) = posDiff($posAry1, $oPosAry);
    my ($sPosAry2, $sLen2) = posDiff($posAry2, $oPosAry);
    return ($oPosAry, $sPosAry1, $sPosAry2, $oLen, $sLen1, $sLen2);
}

sub find_interval {
    my ($pos, $loc) = @_;
    my @dists;
    for my $i (0..@$loc-1) {
        my ($beg, $end) = @{$loc->[$i]};
        my ($dist1, $dist2) = map {$_ - $pos} ($beg, $end);
        my $dist = ($dist1 <= 0 && $dist2 >= 0) ? 0 : min(abs($dist1), abs($dist2));
        push @dists, [$i, $dist1, $dist2, $dist];
    }
    @dists = sort {$a->[3] <=> $b->[3]} @dists;

    my $idx = $dists[0]->[0];
    if($dists[0]->[3] != 0) {
        my ($beg, $end) = @{$loc->[$idx]};
        $pos = $pos < $beg ? $beg : $end;
    }
    return ($idx, $pos);
}
sub coordTransform {
    my ($pos, $locI, $srdI, $locO, $srdO) = @_;
    $locI = [ sort {$a->[0] <=> $b->[0]} @$locI ];
    $locI = [ reverse @$locI ] if $srdI =~ /^\-1?$/;
    $locO = [ sort {$a->[0] <=> $b->[0]} @$locO ];
    $locO = [ reverse @$locO ] if $srdO =~ /^\-1?$/;
    
    my $idx;
    ($idx, $pos) = find_interval($pos, $locI);
      
    my ($rangeI, $rangeO) = map {$_->[$idx]} ($locI, $locO);
    my $pos2;
    if(is_opposite_strands($srdI, $srdO)) {
        $pos2 = $rangeO->[1] - ($pos - $rangeI->[0]);
    } else {
        $pos2 = ($pos - $rangeI->[0]) + $rangeO->[0];
    }
    return $pos2;
}
sub coordTransform_itv {
    my ($pos1, $pos2, $locI, $srdI, $locO, $srdO) = @_;
    $locI = [ sort {$a->[0] <=> $b->[0]} @$locI ];
    $locO = [ sort {$a->[0] <=> $b->[0]} @$locO ];
    $locO = [ reverse @$locO ] if is_opposite_strands($srdI, $srdO);
    ($pos1, $pos2) = ($pos2, $pos1) if $pos1 > $pos2;
    
    my ($idx1, $idx2);
    ($idx1, $pos1) = find_interval($pos1, $locI);
    ($idx2, $pos2) = find_interval($pos2, $locI);
    die "idx1[$idx1] > idx2[$idx2]\n" if $idx1 > $idx2;

    my $loc = [];
    for my $i ($idx1..$idx2) {
        my ($begI, $endI) = @{$locI->[$i]};
        my ($begO, $endO) = @{$locO->[$i]};
        die "locI[$begI-$endI] <> locO[$begO-$endO]\n".Dumper($locI).Dumper($locO) if $endI-$begI != $endO-$begO;
        my ($bl, $el);
        if($i == $idx1 && $i == $idx2) {
            ($bl, $el) = ($pos1, $pos2);
        } elsif($i == $idx1 && $i < $idx2) {
            ($bl, $el) = ($pos1, $endI);
        } elsif($i == $idx2 && $i > $idx1) {
            ($bl, $el) = ($begI, $pos2);
        } else {
            ($bl, $el) = ($begI, $endI);
        }
        my $bg = is_opposite_strands($srdI, $srdO) ? $endO - ($el - $begI) : $begO + ($bl - $begI);
        my $eg = is_opposite_strands($srdI, $srdO) ? $endO - ($bl - $begI) : $begO + ($el - $begI);
        die "bg[$bg] > eg[$eg]\n" if $bg > $eg;
        push @$loc, [$bg, $eg];
    }
    return $loc;
}
sub coordTransform_rough {
    my ($posI, $locI, $srdI, $locO, $srdO) = @_;
    my $itvI = 0;
    if($srdI =~ /^\-1?/) {
        for (@$locI) {
            my ($begI, $endI) = @$_;
            if($posI >= $begI && $posI <= $endI) {
                $itvI += $endI - $posI + 1;
            } elsif($posI < $begI) {
                $itvI += $endI - $begI + 1;
            }
        }
    } else {
        for (@$locI) {
            my ($begI, $endI) = @$_;
            if($posI >= $begI && $posI <= $endI) {
                $itvI += $posI - $begI + 1;
            } elsif($posI > $endI) {
                $itvI += $endI - $begI + 1;
            }
        }
    }

    my $itvO = sprintf "%.0f", ($itvI - 1) * (locAryLen($locO)-1) / (locAryLen($locI)-1) + 1;
    $locO = [ sort {$a->[0] <=> $b->[0]} @$locO ];
    my $posO;
    if($srdO =~ /^\-1?/) {
        for (reverse @$locO) {
            my ($begO, $endO) = @$_;
            my $lenO = $endO - $begO + 1;
            if($lenO < $itvO) {
                $itvO -= $lenO;
            } else {
                $posO = $endO - $itvO + 1;
                last;
            }
        }
    } else {
        for (@$locO) {
            my ($begO, $endO) = @$_;
            my $lenO = $endO - $begO + 1;
            if($lenO < $itvO) {
                $itvO -= $lenO;
            } else {
                $posO = $begO + $itvO - 1;
                last;
            }
        }
    }
    return $posO;
}

1;
__END__
