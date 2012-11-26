package Hits;
use strict; 
use Common; 
use Seq;
use Bio::Seq;
use Align;
use Log::Log4perl;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/pipe_hit hit2Gff/;

sub prepare_for_tiling {
    my ($fis, $fo1, $fo2) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("preparing for hit tiling");

    my $t = readTable(-in=>$fis->[0], -header=>1);
    my $id = 1;
    for my $i (1..@$fis-1) {
        next unless -s $fis->[$i];
        my $ti = readTable(-in=>$fis->[$i], -header=>1);
        for my $j (0..$ti->nofRow-1) {
            $ti->setElm($j, "id", $id++);
            $t->addRow( $ti->rowRef($j) );
        }
    }

    open(FH, ">", $fo1) || die "Can't open file $fo1: $!\n";
    print FH $t->tsv(1);
    close FH;
    
    $t->delCols([qw/score locH locQ locA alnQ alnH alnP/]);
    $t->addCol([('')x$t->nofRow], "note");
    open(FH, ">", $fo2) || die "Can't open file $fo2: $!\n";
    print FH $t->tsv(1);
    close FH;
}
sub hit_tiling {
    my ($fi, $fo, $f_hit) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("tiling HMM hits");

    my $hh;
    my $th = readTable(-in=>$f_hit, -header=>1);
    for my $i (0..$th->nofRow-1) {
        my ($id, $idQ, $srdQ, $locQS, $idH, $srdH, $locHS, $src) =
            map {$th->elm($i, $_)} qw/id idQ srdQ locQ idH srdH locH source/;
        my $locQ = locStr2Ary($locQS);
        my $locH = locStr2Ary($locHS);
        $hh->{$id} = [$idQ, $locQ, $srdQ, $idH, $locH, $srdH, $src];
    }

    my $t = readTable(-in=>$fi, -header=>1);
    my $idx = first_index {$_ eq "domains"} $t->header;
    my $colD = $idx == -1 ? "hits" : "domains";
    my $refH = {};
    for my $i (0..$t->nofRow-1) {
        my ($id, $idQ, $begQ, $endQ, $srdQ, $idH, $begH, $endH, $srdH, $e, $score) = $t->row($i);
        $log->error_die("qry not + strand: ".join("\t", $t->row($i))."\n") unless $srdQ eq "+";
        $refH->{$idH} = [] unless exists $refH->{$idH};
        push @{$refH->{$idH}}, [$begH, $endH, $e, $id];
    }
    
    open(FH, ">", $fo) || die "Can't open file $fo: $!\n";
    print FH join("\t", qw/id idH begH endH srdH idQ begQ endQ srdQ e source note/)."\n";
    for my $idH (sort(keys(%$refH))) {
        my $ref = $refH->{$idH};
        $ref = [ sort {$a->[0] <=> $b->[0]} @$ref ];
        my $locs1 = [ map {[$_->[0], $_->[1]]} @$ref ];
        my $stats = [ map {$_->[2]} @$ref ];
        my $locs2 = tiling($locs1, $stats, 1);
        for (@$locs2) {
            my ($begH, $endH, $idx, $idxs) = @$_;
            my $hs = {};
            for (@$idxs) {
                my ($e, $id) = @{$ref->[$_]}[2,3];
                die "no stat for $id\n" unless exists $hh->{$id};
                my ($idQ, $locQ, $srdQ, $idH2, $locH, $srdH, $src) = @{$hh->{$id}};
                my ($begQ, $endQ) = map {coordTransform($_, $locH, $srdH, $locQ, $srdQ)} ($begH, $endH);
                ($begQ, $endQ) = ($endQ, $begQ) if $begQ > $endQ;
                $hs->{$id} = [$idQ, $begQ, $endQ, $srdQ, $e, $srdH, $src];
            }
            my $id = $ref->[$idx]->[3];
            my @notes;
            for (values(%$hs)) {
                my ($idQ, $begQ, $endQ, $srdQ, $e, $srdH, $src) = @$_;
                my $str = join("_", $idQ, "$begQ-$endQ", $srdQ, $e, $src);
                push @notes, $str;
            }
            my $note = join(" ", @notes);
            my ($idQ, $begQ, $endQ, $srdQ, $e, $srdH, $src) = @{$hs->{$id}};
            print FH join("\t", $id, $idH, $begH, $endH, $srdH, $idQ, $begQ, $endQ, $srdQ, $e, $src, $note)."\n";
        }
    }
    close FH;
}
sub hit_noise_reduction {
    my ($fi, $fo, $p) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("removing noise hits");
    
    my ($min_len, $min_e) = map {$p->{$_}} qw/min_len min_e/;
    $min_len ||= 10;
    $min_e ||= 10;
    
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("idH", 1, 0, "begH", 0, 0);
    my @idxs_rm;

    my $ref = group([$t->col("idH")]);
    for my $idH (sort(uniq(keys %$ref))) {
        my ($idxB, $cnt) = @{$ref->{$idH}};
        my ($sid, $srdH, $rf) = reverse map {scalar reverse} split("_", reverse($idH), 3);
        die join("\t", $idH)."\n" if $srdH eq "1";
        my @locs1 = map {[$t->elm($_, "begH"), $t->elm($_, "endH")]} ($idxB..$idxB+$cnt-1);
        my $locs2 = posMerge(\@locs1);
        for (@$locs2) {
            my ($begHR, $endHR, $idxsR) = @$_;
            my @idxs = map {$idxB + $_} @$idxsR;

            my @stats = map { [ $_, $t->elm($_, "e"), $t->elm($_, "endH") - $t->elm($_, "begH") + 1 ] } @idxs;
            @stats = sort {$a->[1] <=> $b->[1] || $b->[2] <=> $a->[2]} @stats;
            my $idxM = $stats[0]->[0];

            if(@stats > 1) {
                my @idxs_nonM = map {$_->[0]} @stats[1..$#stats];
                push @idxs_rm, @idxs_nonM;
            }

            my ($begH, $endH, $e) = map {$t->elm($idxM, $_)} qw/begH endH e/;
            push @idxs_rm, $idxM if $endH-$begH+1 < $min_len || $e > $min_e;
        }
    }

    $t->delRows(\@idxs_rm);
    $log->info(sprintf("  %d removed / %d passed", scalar(@idxs_rm), $t->nofRow));
    open(FH, ">", $fo) || die "Can't open file $fo: $!\n";
    print FH $t->tsv(1);
    close FH;
}
sub hit_resolve_rf {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("removing hits on wrong reading frames");
    
    my $t = readTable(-in=>$fi, -header=>1);

    my @chrs = map { [reverse map {scalar reverse} split("_", reverse($_), 3)]->[-3] } $t->col("idH");
    $t->addCol(\@chrs, "chr");
    $t->sort("chr", 1, 0, "begH", 0, 0);
    my $ref = group($t->colRef("chr"));

    my @idxs_rm;
    for my $chr (sort(uniq(keys %$ref))) {
        my ($idxB, $cnt) = @{$ref->{$chr}};
        my @locs1 = map {[$t->elm($_, "begH"), $t->elm($_, "endH")]} ($idxB..$idxB+$cnt-1);
        my $locs2 = posMerge(\@locs1);
        for (@$locs2) {
            my ($begHR, $endHR, $idxsR) = @$_;
            my @idxs = map {$idxB + $_} @$idxsR;

            my @stats = map { [ $_, $t->elm($_, "e") ] } @idxs;
            @stats = sort {$a->[1] <=> $b->[1]} @stats;
            my $idxM = $stats[0]->[0];

            if(@stats > 1) {
                my @idxs_nonM = map {$_->[0]} @stats[1..$#stats];
#        print join("\t", map {$t->elm($_, "id")} @idxs)."\n";
                push @idxs_rm, @idxs_nonM;
            }
        }
    }

    $t->delRows(\@idxs_rm);
    my $chrs = $t->delCol("chr");
    $t->delCol("idH");
    $t->addCol($chrs, "idH", 1);
    $log->info(sprintf("  %d removed / %d passed", scalar(@idxs_rm), $t->nofRow));
    open(FH, ">", $fo) || die "Can't open file $fo: $!\n";
    print FH $t->tsv(1);
    close FH;
}

sub unifyHits {
    my ($locAs, $srds, $group) = @_;
    my $n = @$locAs;
    my $r = [];
    my ($h, @loc);
    for my $i (0..$n-1) {
        my ($locA, $srd) = ($locAs->[$i], $srds->[$i]);
        for my $j (0..@$locA-1) {
            push @loc, [$locA->[$j]->[0], $locA->[$j]->[1], $i];
        }
        $h->{$i} = ['+0', $locA, $srd];
    }
    my $ref = posSplit(\@loc);
    my $ref2 = {};
    for my $j (0..@$ref-1) {
        my ($s, $e, $idxs) = @{$ref->[$j]};
        for my $i (@$idxs) {
            push @{$ref2->{$i}}, $j;
        }
    }
    for my $i (0..$n-1) {
        my $idxPs = $ref2->{$i};
        my $srd = $h->{$i}->[2];
        $idxPs = [sort {$ref->[$a]->[0] <=> $ref->[$b]->[0]} @$idxPs];
        $idxPs = [reverse @$idxPs] if $srd =~ /^\-1?$/;
        my $locAry = [ map { [$ref->[$_]->[0], $ref->[$_]->[1]] } @$idxPs ];
        my @phases = getPhase($locAry, $srd);
        $ref2->{$i} = { map {$idxPs->[$_] => $phases[$_]} (0..$#phases) };
    }
    for my $j (0..@$ref-1) {
        my ($s, $e, $idxs) = @{$ref->[$j]};
        my $pH;
        for my $i (@$idxs) {
            my $phase = $ref2->{$i}->{$j};
            $pH->{$phase} = [] unless exists $pH->{$phase};
            push @{$pH->{$phase}}, $i;
        }
        push @$r, [$s, $e, $pH];
    }
    return $r;
}
sub group_pick_hits {
    my ($fi, $fo1, $fo2) = rearrange(['in', 'out1', 'out2'], @_);
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("sorting hits into groups");
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("chr", 1, 0);
    my $ref = group($t->colRef("chr"));

    open(FH1, ">$fo1");
    print FH1 join("\t", qw/id1 id2/, $t->header)."\n";
    open(FH2, ">$fo2");
    print FH2 join("\t", $t->header)."\n";
    
    my $id1 = 0;
    for my $chr (sort keys %$ref) {
        my ($idxS, $n) = @{$ref->{$chr}};
        my $t1 = $t->subTable([$idxS..$idxS+$n-1]);

        my @locAs = map {[split("-", $_)]} $t1->col("loc");
        my $ref = posMerge(\@locAs);
        for (@$ref) {
            my ($beg, $end, $idxs) = @$_;
            $id1 ++;

            for my $i (0..@$idxs-1) {
                my $idx = $idxs->[$i];
                print FH1 join("\t", $id1, $i+1, $t1->row($idx))."\n";
            }
            
            my @locAs = map {locStr2Ary($t1->elm($_, "loc"))} @$idxs;
            my @srds = map {$t1->elm($_, "srd")} @$idxs;
            my ($ref2) = unifyHits(\@locAs, \@srds, $id1);
            my @tmp = grep {scalar keys %{$_->[2]} > 1} @$ref2;
            my $tag = @tmp > 0 ? 1 : 0;
#      print "Possible pseudogene: $id1\n" if $tag == 1;

            my @es = map {$t1->elm($_, "e")} @$idxs;
            my $idx_min = first_index {$_ == min(@es)} @es;
            my $id2 = $idx_min + 1;
            print FH2 join("\t", $t1->row($idxs->[$idx_min]))."\n";
        }
    }
    $log->info(sprintf "  %d groups in total", $id1);
    close FH1;
    close FH2;
}

sub note2hash {
    my ($str) = @_;
    my $h;
    my @ps = split " ", $str;
    for (@ps) {
        my ($fam, $locQS, $srdQ, $e, $src) = split("_", $_);
        my $locQ = locStr2Ary($locQS);
        $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
        my ($begQ, $endQ) = ($locQ->[0]->[0], $locQ->[-1]->[1]);
        $h->{$fam} = [$begQ, $endQ, $srdQ, $e, $src, $locQ] unless exists $h->{$fam};
        $h->{$fam} = [$begQ, $endQ, $srdQ, $e, $src, $locQ] if $h->{$fam}->[3] > $e;
    }
    return $h;
}
sub hash2note {
    my ($h) = @_;
    my @notes;
    for my $key (sort(keys(%$h))) {
        my ($beg, $end, $srd, $e, $src, $loc) = @{$h->{$key}};
        push @notes, join("_", $key, locAry2Str($loc), $srd, $e, $src);
    }
    return join(" ", @notes);
}
sub reformat_hit_info {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("re-formatting hit info");
    my $t = readTable(-in=>$fi, -header=>1);
    
    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end srd loc family begQ endQ srdQ locQ e source note/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $srd, $fam, $begQ, $endQ, $srdQ, $e, $source, $note) = $t->row($i);
        die "$id idQ[$fam] not + strand\n".join("\t", $t->row($i))."\n" if $srdQ ne "+";
        my $h = note2hash($note);
        my @fams = sort {$h->{$a}->[3] <=> $h->{$b}->[3]} keys(%$h);
        @fams = @fams[0..2] if @fams > 3;
        my $hN = { map {$_ => $h->{$_}} @fams };
        my $noteN = hash2note($hN);
        print FH join("\t", $chr, $beg, $end, $srd, "$beg-$end", $fam, $begQ, $endQ, $srdQ, "$begQ-$endQ", $e, $source, $noteN)."\n";
    }
    close FH;
}

sub merge_fam_stat {
    my ($h1, $h2) = @_;
    my $h;
    for my $fam (uniq(keys(%$h1), keys(%$h2))) {
        if(exists $h1->{$fam} && exists $h2->{$fam}) {
            my ($beg1, $end1, $srd1, $e1, $src1, $loc1) = @{$h1->{$fam}}; 
            my ($beg2, $end2, $srd2, $e2, $src2, $loc2) = @{$h2->{$fam}};
            my $loc = [ sort {$a->[0] <=> $b->[0]} (@$loc1, @$loc2) ];
            my ($beg, $end, $srd) = ($loc->[0]->[0], $loc->[-1]->[1], $srd1);
            my ($e, $src) = ($e1<$e2) ? ($e1, $src1) : ($e2, $src2);
            $h->{$fam} = [$beg, $end, $srd, $e, $src, $loc];
#    } elsif(exists $h1->{$fam}) {
#      $h->{$fam} = [@{$h1->{$fam}}];
#    } else {
#      $h->{$fam} = [@{$h2->{$fam}}];
        }
    }
    return $h;
}

sub remove_noisy_hits {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("removing noise hits between multi-part hits");
    my ($lenO1, $lenO2, $dist_max) = (45, 15, 5000);

    my @idxs_rm;
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (2..$t->nofRow-1) {
        my ($chr, $beg, $end, $srd) = map {$t->elm($i, $_)} qw/chr beg end srd/;
        my ($loc, $h) = (locStr2Ary($t->elm($i, "loc")), note2hash($t->elm($i, "note")));

        my $idxP = $i - 1;
        my ($chrP, $begP, $endP, $srdP) = map {$t->elm($idxP, $_)} qw/chr beg end srd/;
        my ($locP, $hP) = (locStr2Ary($t->elm($idxP, "loc")), note2hash($t->elm($idxP, "note")));

        my @fams_all = uniq(keys(%$hP), keys(%$h));
        my @fams = grep { exists($h->{$_}) && exists($hP->{$_}) && 
            ( ($srd eq "+" && $h->{$_}->[0] >= $hP->{$_}->[1] - $lenO1 && $h->{$_}->[0] > $hP->{$_}->[0] + $lenO2) ||
                ($srd eq "-" && $hP->{$_}->[0] >= $h->{$_}->[1] - $lenO1 && $hP->{$_}->[0] > $h->{$_}->[0] + $lenO2)
            ) } @fams_all;

        if($chr eq $chrP && $srd eq $srdP && ($beg-$endP) < $dist_max && @fams > 0) {
            next;
        }

        $idxP = $i - 2;
        ($chrP, $begP, $endP, $srdP) = map {$t->elm($idxP, $_)} qw/chr beg end srd/;
        ($locP, $hP) = (locStr2Ary($t->elm($idxP, "loc")), note2hash($t->elm($idxP, "note")));
        my $eN = $t->elm($i-1, "e");

        @fams_all = uniq(keys(%$hP), keys(%$h));
        @fams = grep { exists($h->{$_}) && exists($hP->{$_}) && 
            ( ($srd eq "+" && $h->{$_}->[0] >= $hP->{$_}->[1] - $lenO1 && $h->{$_}->[0] > $hP->{$_}->[0] + $lenO2) ||
                ($srd eq "-" && $hP->{$_}->[0] >= $h->{$_}->[1] - $lenO1 && $hP->{$_}->[0] > $h->{$_}->[0] + $lenO2)
            ) } @fams_all;
        if($chr eq $chrP && $srd eq $srdP && ($beg-$endP) < $dist_max && @fams > 0 && $eN > 1E-4) {
            push @idxs_rm, $i-1;
        }
    }
    $t->delRows(\@idxs_rm);
    $log->info(sprintf "  %d removed / %d passed", scalar(@idxs_rm), $t->nofRow);
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub merge_hits_multi {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("merging multi-part hits");
    my ($lenO1, $lenO2, $dist_max) = (45, 15, 5000);

    my @idxs_rm;
    my $idxP = 0;
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (1..$t->nofRow-1) {
        my ($chr, $beg, $end, $srd) = map {$t->elm($i, $_)} qw/chr beg end srd/;
        my ($loc, $h) = (locStr2Ary($t->elm($i, "loc")), note2hash($t->elm($i, "note")));

        my ($chrP, $begP, $endP, $srdP) = map {$t->elm($idxP, $_)} qw/chr beg end srd/;
        my ($locP, $hP) = (locStr2Ary($t->elm($idxP, "loc")), note2hash($t->elm($idxP, "note")));
      
        my @fams_all = uniq(keys(%$hP), keys(%$h));
        my @fams = grep { exists($h->{$_}) && exists($hP->{$_}) && 
            ( ($srd eq "+" && $h->{$_}->[0] >= $hP->{$_}->[1] - $lenO1 && $h->{$_}->[0] > $hP->{$_}->[0] + $lenO2) ||
                ($srd eq "-" && $hP->{$_}->[0] >= $h->{$_}->[1] - $lenO1 && $hP->{$_}->[0] > $h->{$_}->[0] + $lenO2)
            ) } @fams_all;
        if($chr eq $chrP && $srd eq $srdP && ($beg-$endP) < $dist_max && @fams > 0) {
            my $hN = merge_fam_stat($hP, $h);
            @fams = sort {$hN->{$a}->[3] <=> $hN->{$b}->[3]} @fams;
            my $fam = $fams[0];

            $t->setElm($idxP, "end", $end);
            $t->setElm($idxP, "loc", locAry2Str([@$locP, @$loc]));
            $endP = $end;
            $locP = [@$locP, @$loc];

            my ($begQ, $endQ, $srdQ, $e, $src, $locQ) = @{$hN->{$fam}};
            $t->setElm($idxP, "begQ", $begQ);
            $t->setElm($idxP, "endQ", $endQ);
            $t->setElm($idxP, "locQ", locAry2Str($locQ));
            $t->setElm($idxP, "e", $e);
            $t->setElm($idxP, "source", $src);
            $t->setElm($idxP, "note", hash2note($hN));
            $hP = $hN;

            push @idxs_rm, $i;
            next;
        }
        $idxP = $i;
    }
    $t->delRows(\@idxs_rm);
    $log->info(sprintf "  %d removed / %d passed", scalar(@idxs_rm), $t->nofRow);
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub merge_sub_loc {
    my ($locI, $idxs_merge) = @_;
    $locI = [ sort {$a->[0] <=> $b->[0]} @$locI ];
    my ($idxs_nomerge) = posDiff([[0, @$locI-1]], $idxs_merge);
    my $locO = [];
    for (@$idxs_merge) {
        my ($idxB, $idxE) = @$_;
        my $b = min( map {$locI->[$_]->[0]} ($idxB..$idxE) );
        my $e = max( map {$locI->[$_]->[1]} ($idxB..$idxE) );
        push @$locO, [$b, $e];
    }
    for (@$idxs_nomerge) {
        my ($idxB, $idxE) = @$_;
        for my $idx ($idxB..$idxE) {
            push @$locO, $locI->[$idx];
        }
    }
    $locO = [ sort {$a->[0] <=> $b->[0]} @$locO ];
    return $locO;
}
sub merge_hits_within {
    my ($fi, $fo, $f_ref) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("merging within-group hits");

    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($chr, $beg, $end, $srd, $locS, $note) = map {$t->elm($i, $_)} qw/chr beg end srd loc note/;
        my $loc = locStr2Ary($locS);
        next if @$loc <= 1;
        $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
      
        my @idxpairs_merge;
        for my $i (1..@$loc-1) {
            my ($beg1, $end1) = @{$loc->[$i-1]};
            my ($beg2, $end2) = @{$loc->[$i]};
            my $loc_intv = [[$end1+1, $beg2-1]];
            my $intv = $beg2 - $end1 - 1;
            if($intv <= 60 && $intv % 3 == 0) {
                my $seqI = seqRet($loc_intv, $chr, $srd, $f_ref);
                my $prot = Bio::Seq->new(-seq=>$seqI)->translate()->seq;
                push @idxpairs_merge, [$i-1, $i] if $prot !~ /\*/;
            }
        }
        next unless @idxpairs_merge;
        
#    print "  $chr:$beg-$end\[$srd]\n";
        my $idxs_merge = posMerge(\@idxpairs_merge);
        my $locN = merge_sub_loc($loc, $idxs_merge);
        $t->setElm($i, "loc", locAry2Str($locN));

        my $h = note2hash($note);
        for my $fam (keys(%$h)) {
            my ($beg, $end, $srdQ, $e, $src, $locQo) = @{$h->{$fam}};
            $log->info("  $chr:$beg-$end\[$srd] [$note] hit != qry") if @$loc != @$locQo;
            $locQo = [ reverse @$locQo ] if is_opposite_strands($srd, $srdQ);
            $h->{$fam}->[-1] = merge_sub_loc($locQo, $idxs_merge);
        }
        my @fams = sort {$h->{$a}->[3] <=> $h->{$b}->[3]} keys(%$h);
        my $fam = $fams[0];
        my $locQ = $h->{$fam}->[-1];
        $t->setElm($i, "locQ", locAry2Str($locQ));
        $t->setElm($i, "note", hash2note($h));
    }
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}

sub remove_pseudo {
    my ($fi, $fo, $f_ref) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("removing pseudogenes");
    
    my $t = readTable(-in=>$fi, -header=>1);
    my @idxs_rm;
    for my $i (0..$t->nofRow-1) {
        my ($chr, $beg, $end, $srd, $locS) = map {$t->elm($i, $_)} qw/chr beg end srd loc/;
        my $loc = locStr2Ary($locS);
        my $seq_dna = seqRet($loc, $chr, $srd, $f_ref);
        my $seq = Bio::Seq->new(-seq=>$seq_dna)->translate()->seq;
        
        if($seq =~ /\*/) {
            $log->warn(sprintf("  %s:%s[%s] removed: premature stop codon", $chr, $locS, $srd));
            push @idxs_rm, $i;
        }
    }

    $t->delRows(\@idxs_rm);
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub trim_cds_boundary {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("trimming CDS boundary");
    my $len_min = 15 * 3;
  
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($chr, $beg, $end, $srd, $locS, $srdQ, $locQS) = map {$t->elm($i, $_)} qw/chr beg end srd loc srdQ locQ/;
        my $loc = locStr2Ary($locS);
        my $locQ = locStr2Ary($locQS);
        $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
        $loc = [ reverse @$loc ] if $srd eq "-";
        $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
        die "$chr:$beg-$end\[$srd] srdQ != +\n" unless $srdQ eq "+";
        if(@$loc == 1) {
            $loc->[0]->[1] -= 3*5 if locAryLen($loc) > $len_min;
            $loc->[0]->[0] += 3*5 if locAryLen($loc) > $len_min;
            next;
        }

        $log->info("  $chr:$beg-$end\[$srd] hit != qry") if @$loc != @$locQ;
        next if @$loc != @$locQ;
      
        for my $j (1..@$loc-1) {
            my ($locOvlp, $lenO) = posOvlp([$locQ->[$j-1]], [$locQ->[$j]]);
            $lenO = 3 * (int($lenO/3) + 1);
            my ($len1, $len2) = (locAryLen([$loc->[$j-1]]), locAryLen([$loc->[$j]]));
            my $lenT1 = max(min($len_min, $len1), $len1 - $lenO);
            my $lenT2 = max(min($len_min, $len2), $len2 - $lenO);
            if($srd eq "+") {
                $loc->[$j-1]->[1] = $loc->[$j-1]->[0] + $lenT1 - 1;
                $loc->[$j  ]->[0] = $loc->[$j  ]->[1] - $lenT2 + 1;
            } else {
                $loc->[$j-1]->[0] = $loc->[$j-1]->[1] - $lenT1 + 1;
                $loc->[$j  ]->[1] = $loc->[$j  ]->[0] + $lenT2 - 1;
            }
        }
        
        $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
        $t->setElm($i, 'loc', locAry2Str($loc));
        $t->setElm($i, 'beg', $loc->[0]->[0]);
        $t->setElm($i, 'end', $loc->[-1]->[1]);
    }

    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub output_hits {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("print hits in subgroup order");
    
    my $t = readTable(-in=>$fi, -header=>1);
    my $fams = $t->delCol("family"); 
    $t->addCol($fams, "family", 0);
    $t->sort("family", 1, 0, "chr", 1, 0, "beg", 0, 0);

    my @ids = map {sprintf "h%04d", $_} (1..$t->nofRow);
    $t->addCol(\@ids, "id", 0);
    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH $t->tsv(1);
    close FH;
}
  
sub hit2Gff {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo");
    print FH "##gff-version 3\n";

    my $flag_id = $t->hasCol("id");
    for my $i (0..$t->nofRow-1) {
        my ($fam, $chr, $srd, $locS, $e, $source) = map {$t->elm($i, $_)} qw/family chr srd loc e source/;
        my $loc = locStr2Ary($locS);
        my $id = $flag_id ? $t->elm($i, "id") : $i;
        for (@$loc) {
            my ($beg, $end) = @$_;
            print FH join("\t", $chr, $source, "match", $beg, $end, $e, $srd, ".", "ID=$id;Name=$id;Note=[$fam][$source][$e]")."\n";
        }
    }
    close FH;
}
sub align_core_seq {
    my ($fi, $dirO) = @_;
    make_path($dirO) unless -d $dirO;
    remove_tree($dirO, {keep_root => 1});

    my $t = readTable(-in=>$fi, -header=>1);
    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $strand, $locStr, $locQ, $e, $source, $seq) = $t->row($i);
        $h->{$fam} ||= [];
        my $seqObj = Bio::Seq->new(-id=>$id, -seq=>$seq);
        push @{$h->{$fam}}, $seqObj;
    }

    for my $fam (sort keys %$h) {
        my $seqs = $h->{$fam};
        next if @$seqs < 2;
        my $fo = "$dirO/$fam.aln";
        run_clustalo(-seqs=>$seqs, -out=>$fo);
    }
}
sub get_aln_score {
    my ($fi, $d_aln, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hits");
    $log->info("computing MSA scores");
    my $t = readTable(-in=>$fi, -header=>1);
    $t->addCol([('')x$t->nofRow], "aln_score", 9);
    
    my $f_bin = $ENV{"ClustalO"}."/bin/clustalo";
    my $ft1 = $ENV{"TMP_DIR"}."/aln_score1_".int(rand(1000)).".fa";
    my $ft2 = $ENV{"TMP_DIR"}."/aln_score2_".int(rand(1000)).".fa";
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $strand, $locStr, $locQ, $e, $aln_score, $source, $seq) = $t->row($i);
#    next unless $id == 1;
        my $f_aln = "$d_aln/$fam.aln";
        die "$f_aln is not there\n" unless -s $f_aln;

        writeFile($ft1, ">seq$id", $seq);
#    print "  pls ignore this warning: $id [$seq]\n" if $seq =~ /^[atcg]+$/i;
        runCmd("$f_bin --p1 $ft1 --p2 $f_aln --outfmt=fasta --force -o $ft2", 0);
        die "MSP $ft2 was not created\n" unless -s $ft2;
        my $score = aln_score_vector($ft2, "seq$id");
        $t->setElm($i, "aln_score", $score);
        printf "  %5d / %5d done...\r", $i+1, $t->nofRow;
    }
    print "\n";
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
    system("rm $ft1 $ft2");
}


sub pipe_hit {
    my ($dir, $f_hit, $p, $f_ref, $d_aln) = rearrange(['dir', 'in', 'p', 'ref', 'aln'], @_);
    make_path($dir) unless -d $dir;
    
    my $f01f = "$dir/01_hit_full.htb";
    my $f01 = "$dir/01_hit.htb";
    prepare_for_tiling($f_hit, $f01f, $f01);
    
    my $f02 = "$dir/02_tiled.htb";
    hit_tiling($f01, $f02, $f01f);
    my $f04 = "$dir/04_nr.htb";
    hit_noise_reduction($f02, $f04, $p);
    my $f05 = "$dir/05_resolve_rf.htb";
    hit_resolve_rf($f04, $f05);

    my $f11 = "$dir/11.tbl";
    reformat_hit_info($f05, $f11);
    my $f12 = "$dir/12.tbl";
    remove_noisy_hits($f11, $f12);
    my $f13 = "$dir/13_merged_multi.tbl";
    merge_hits_multi($f12, $f13);
    my $f17 = "$dir/17_merged_within.tbl";
    merge_hits_within($f13, $f17, $f_ref);
    
    my $f21 = "$dir/21.tbl";
    remove_pseudo($f17, $f21, $f_ref);
    my $f22 = "$dir/22_trimmed.tbl";
    trim_cds_boundary($f21, $f22);
    my $f29 = "$dir/29_hits.tbl";
    output_hits($f22, $f29);
    hit2Gff($f29, "$dir/29_hits.gff");
}




1;
__END__
