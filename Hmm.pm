package Hmm;
use strict;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Bio::SearchIO;
use Common;
use Seq;
use Gtb;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/parse_hmm_output parse_hmm_simple 
    score_hmm_by_hit
    hitTiling hitCollapse recoverCoordP get_longest_cds
    pipe_hmmsearch/;
@EXPORT_OK = qw//;

sub parse_aln {
    my ($alnQ, $alnT, $alnP) = @_;
    
    my @poss;
    my @locmap;
    my ($posA, $posQ, $posT) = (0, 0, 0);
    for my $i (1..length($alnQ)) {
        my $chQ = substr($alnQ, $i-1, 1);
        my $chT = substr($alnT, $i-1, 1);
        if($chQ !~ /[\.\-]/ || $chT !~ /[\.\-]/) {
            push @poss, $i;
            if($chQ =~ /[\.\-]/ && $chT !~ /[\.\-]/) {
                push @locmap, [++$posA, '', ++$posT];
            } elsif($chQ !~ /[\.\-]/ && $chT =~ /[\.\-]/) {
                push @locmap, [++$posA, ++$posQ, ''];
            } else {
                push @locmap, [++$posA, ++$posQ, ++$posT];
            }
        }
    }

    my ($lenA, $lenQ, $lenT) = ($posA, $posQ, $posT);
    my $gapQ = grep {$_->[1] eq ''} @locmap;
    my $gapT = grep {$_->[2] eq ''} @locmap;

    my $alnQ_n = join("", map {substr($alnQ, $_-1, 1)} @poss);
    my $alnT_n = join("", map {substr($alnT, $_-1, 1)} @poss);
    my $alnP_n = join("", map {substr($alnP, $_-1, 1)} @poss);

    my @idxs_ins = indexes {$_->[1] eq '' || $_->[2] eq ''} @locmap;
    my ($idxs_nogap) = posDiff([[0, @locmap-1]], [map {[$_, $_]} @idxs_ins]);
    my ($locA, $locQ, $locT) = ([], [], []);
    for (@$idxs_nogap) {
        my ($idxB, $idxE) = @$_;
        push @$locA, [$locmap[$idxB]->[0], $locmap[$idxE]->[0]];
        push @$locQ, [$locmap[$idxB]->[1], $locmap[$idxE]->[1]];
        push @$locT, [$locmap[$idxB]->[2], $locmap[$idxE]->[2]];
    }

    return ($alnQ_n, $alnT_n, $alnP_n, $lenA, $lenQ, $lenT, $gapQ, $gapT, $locA, $locQ, $locT);
}
sub parse_hmm_output {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("parsing HMM output");

    my $in = Bio::SearchIO->new(-file=>$fi, -format=>"hmmsearch3");
    my $t = Data::Table->new([], [qw/id idQ begQ endQ srdQ locQ idH begH endH srdH locH e score source locA alnQ alnH alnP/]);
    my $id = 1;
    while( my $res = $in->next_result ) {
        while( my $hit = $res->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                my $idQ = $res->query_name;
                my $idH = $hit->name;
                my ($begQ, $endQ) = ($hsp->start("query"), $hsp->end("query"));
                my ($begH, $endH) = ($hsp->start("hit"), $hsp->end("hit"));
                my ($alnQ, $alnH) = ($hsp->query_string, $hsp->hit_string);
                my $alnP = $hsp->{"PP_SEQ"};
                my ($e, $score) = ($hsp->{"IEVALUE"}, $hsp->score());
                $t->addRow([$id++, $idQ, $begQ, $endQ, "+", "", $idH, $begH, $endH, "+", "", $e, $score, '', '', $alnQ, $alnH, $alnP]);
            }
        }
    }

    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
    print FH $t->tsv(1);
    close FH;
}
sub get_aln_position {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("calculating alignment coordinates");
    
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $idQ, $begQ, $endQ, $srdQ, $locQS, $idH, $begH, $endH, $srdH, $locHS, $e, $score, $source, $locAS, $alnQ, $alnH, $alnP) = $t->row($i);
        
        my ($lenAQ, $lenAH, $lenAP) = map {length($_)} ($alnQ, $alnH, $alnP);
        $lenAQ == $lenAH || $log->error_die("$idQ $idH:$begH-$endH length not equal [alnQ:$lenAQ <> alnH:$lenAH");
        $lenAQ == $lenAP || $log->error_die("$idQ $idH:$begH-$endH length not equal [alnQ:$lenAQ <> alnH:$lenAP");

        my ($alnQ_n, $alnH_n, $alnP_n, $lenA, $lenQ, $lenH, $gapQ, $gapH, $locA, $locQ, $locH) = parse_aln($alnQ, $alnH, $alnP);
        $lenA == $lenQ + $gapQ || $log->error_die("$idQ $idH:$begH-$endH lenA[$lenA] != lenQ[$lenQ] + gapQ[$gapQ]\n$alnQ_n\n$alnH_n");
        $lenA == $lenH + $gapH || $log->error_die("$idQ $idH:$begH-$endH lenA[$lenA] != lenH[$lenH] + gapH[$gapH]\n$alnQ_n\n$alnH_n");
        $lenQ == locAryLen($locQ) + $gapH || $log->error_die("$idQ $idH:$begH-$endH lenQ[$lenQ] ".locAry2Str($locQ)." gapH[$gapH]\n$alnQ_n\n$alnH_n");
        $lenQ == $endQ - $begQ + 1 || $log->error_die("$idQ $idH:$begH-$endH lenQ[$lenQ] != $begQ-$endQ");
        $lenH == locAryLen($locH) + $gapQ || $log->error_die("$idQ $idH:$begH-$endH lenH[$lenH] ".locAry2Str($locH)." gapQ[$gapQ]\n$alnQ_n\n$alnH_n");
        $lenH == $endH - $begH + 1 || $log->error_die("$idQ $idH:$begH-$endH lenH[$lenH] != $begH-$endH");
        $t->setElm($i, "alnQ", $alnQ_n);
        $t->setElm($i, "alnH", $alnH_n);
        $t->setElm($i, "alnP", $alnP_n);

        $locQ = [ map {[$_->[0] + $begQ - 1, $_->[1] + $begQ - 1]} @$locQ ];
        $locH = [ map {[$_->[0] + $begH - 1, $_->[1] + $begH - 1]} @$locH ];
        $t->setElm($i, "locA", locAry2Str($locA));
        $t->setElm($i, "locQ", locAry2Str($locQ));
        $t->setElm($i, "locH", locAry2Str($locH));
    }
    
    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing\n");
    print FH $t->tsv(1);
    close FH;
}

sub score_homology_string {
    my ($str) = @_;
    my $score = 0;
    for my $i (0..length($str)-1) {
        my $ch = substr($str, $i, 1);
        if($ch eq " ") {
            $score -= 0.5;
        } elsif($ch eq "+") {
            $score += 1;
        } else {
            $score += 2;
        }
    }
    return $score;
}
sub score_probability_string {
    my ($alnP, $alnQ, $alnH) = @_;
    my $score = 0;
    for my $i (0..length($alnP)-1) {
        my ($ch, $chQ, $chH) = map {substr($_, $i, 1)} ($alnP, $alnQ, $alnH);
        if($chH =~ /[\.\-]/) {
        } elsif($ch eq "*") {
            $score += 10;
        } else {
            die "unknown pp[$ch] in $alnP\n" if $ch !~ /[0-9]/;
            $score += $ch;
        }
    }
    return $score;
}
sub parse_hmm_simple {
    my ($fi) = @_;
    my $h;
    my $in = Bio::SearchIO->new(-file=>$fi, -format=>"hmmsearch3");
    while( my $res = $in->next_result ) {
        while( my $hit = $res->next_hit ) {
            my ($id, $e, $score) = map {$hit->$_} qw/name significance score/;
            my @hsps;
            while( my $hsp = $hit->next_hsp ) {
                my ($begQ, $endQ) = ($hsp->start("query"), $hsp->end("query"));
                my ($begH, $endH) = ($hsp->start("hit"), $hsp->end("hit"));
                my ($e, $score_bit) = ($hsp->{"IEVALUE"}, $hsp->score());
                my ($alnQ, $alnH, $alnM, $alnP) = 
                    ($hsp->query_string, $hsp->hit_string, $hsp->homology_string, $hsp->{"PP_SEQ"});
                my $score_pp = score_probability_string($alnP, $alnQ, $alnH);
                push @hsps, [$begQ, $endQ, $e, $begH, $endH, $score_pp, $score_bit];
            }
            $h->{$id} = \@hsps if @hsps;
        }
    }
    return $h;
}
sub score_hmm_by_hit {
    my ($fi) = @_;
    my $h = parse_hmm_simple($fi);
    my $hs;
    for my $id (keys %$h) {
        my @hsps = @{$h->{$id}};
        my @locs = map {[$_->[0], $_->[1]]} @hsps;
        my @stats = map {$_->[2]} @hsps;
        my $ref = tiling(\@locs, \@stats, 1);
        my $score = 0;
        for (@$ref) {
            my ($beg, $end, $idx) = @$_;
            my ($begQ, $endQ, $e, $begH, $endH, $score_pp, $score_bit) = @{$h->{$id}->[$idx]};
            $score += ($end-$beg+1) * $score_pp / ($endQ-$begQ+1);
        }
        $hs->{$id} = $score;
    }
    return $hs;
}

sub coordTransform_aa2nt {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("transforming coordinates (Amino Acid -> Nucleotide)");
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $idQ, $begQ, $endQ, $srdQ, $locQS, $idH, $begH, $endH, $srdH, $locHS, $e, $score, $srcN, $locAS) = $t->row($i);
        my $locQL = locStr2Ary($locQS);
        my $locQG = [ map {[ 3*$_->[0]-2, 3*$_->[1] ]} @$locQL ];
        $locQG = [ sort {$a->[0] <=> $b->[0]} @$locQG ];
        ($begQ, $endQ) = ($locQG->[0]->[0], $locQG->[-1]->[1]);
        $t->setElm($i, "begQ", $begQ);
        $t->setElm($i, "endQ", $endQ);
        $t->setElm($i, "locQ", locAry2Str($locQG));
        
        my $locAL = locStr2Ary($locAS);
        my $locAG = [ map {[ 3*$_->[0]-2, 3*$_->[1] ]} @$locAL ];
        $t->setElm($i, "locA", locAry2Str($locAG));

        my $locHL = locStr2Ary($locHS);
        my $locHG = [ map {[ 3*$_->[0]-2, 3*$_->[1] ]} @$locHL ];
        $locHG = [ sort {$a->[0] <=> $b->[0]} @$locHG ];
        ($begH, $endH) = ($locHG->[0]->[0], $locHG->[-1]->[1]);
        $t->setElm($i, "begH", $begH);
        $t->setElm($i, "endH", $endH);
        $t->setElm($i, "locH", locAry2Str($locHG));
    }
    open(FH, ">", $fo) || die "Can't open file $fo: $!\n";
    print FH $t->tsv(1);
    close FH;
}

sub get_relative_loc {
    my ($locI, $srd) = @_;
    $locI = [ sort {$a->[0] <=> $b->[0]} @$locI ];
    my ($beg, $end) = ($locI->[0]->[0], $locI->[-1]->[1]);
    $locI = [ reverse @$locI ] if $srd =~ /^-1?$/;

    my $locO = [];
    my $len = 0;
    for (@$locI) {
        my ($bg, $eg) = @$_;
        my $len1 = $eg - $bg + 1;
        push @$locO, [$len+1, $len+$len1];
        $len += $len1;
    }
    return $locO;
}
sub trim_aln_loc {
    my ($locH, $begHn, $endHn, $srdH, $locQ, $srdQ, $locA, $alnQ, $alnH, $alnP) = @_;
    my $locHn = trimLoc($locH, $begHn, $endHn);
    $locHn = [ sort {$a->[0] <=> $b->[0]} @$locHn ];
    
    my $begQn = coordTransform($begHn, $locH, $srdH, $locQ, $srdQ);
    my $endQn = coordTransform($endHn, $locH, $srdH, $locQ, $srdQ);
    ($begQn, $endQn) = ($endQn, $begQn) if $begQn > $endQn;
    my $locQn = trimLoc($locQ, $begQn, $endQn);
    $locQn = [ sort {$a->[0] <=> $b->[0]} @$locQn ];
    
    my $begAn = coordTransform($begHn, $locH, $srdH, $locA, "+");
    my $endAn = coordTransform($endHn, $locH, $srdH, $locA, "+");
    ($begAn, $endAn) = ($endAn, $begAn) if $begAn > $endAn;
    my $locAn = trimLoc($locA, $begAn, $endAn);
    $locAn = [ map {[$_->[0] - $begAn + 1, $_->[1] - $begAn + 1]} @$locAn ];

    my ($ba, $ea) = ( ($begAn+2)/3, $endAn/3 );
    my $alnQn = substr($alnQ, $ba-1, $ea-$ba+1);
    my $alnHn = substr($alnH, $ba-1, $ea-$ba+1);
    my $alnPn = substr($alnP, $ba-1, $ea-$ba+1);
    return ($locQn, $locHn, $locAn, $alnQn, $alnHn, $alnPn);
}
sub recoverCoord {
    my ($fi, $fo, $f_ref) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("recovering to global coordinate");

    my $t = readTable(-in=>$fi, -header=>1);
    open(FH, ">", $fo) || die "Can't open file $fo: $!\n";
    print FH join("\t", $t->header)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $idQ, $begQ, $endQ, $srdQ, $locQS, $idH, $begH, $endH, $srdH, $locHS, $e, $score, $srcN, $locAS, $alnQ, $alnH, $alnP) = $t->row($i);
        $log->error_die("srd error: \n".join("\t", $t->row($i))."\n") if $srdQ ne "+" || $srdH ne "+";

        my ($seqid, $locHWS, $srd, $src) = reverse map {scalar reverse} split("_", reverse($idH), 4);
        $log->error_die("Unknown strand $srd\n") if $srd !~ /^[\+\-]$/;
        my $locQ = locStr2Ary($locQS);
        my $locA = locStr2Ary($locAS);
        my $locHL = locStr2Ary($locHS);
        
        my $locHWO = locStr2Ary($locHWS);
        $locHWO = [ sort {$a->[0] <=> $b->[0]} @$locHWO ];
        my $locHWI = get_relative_loc($locHWO, $srd);
        
        for my $j (0..@$locHWI-1) {
            my $locHW = [$locHWI->[$j]];
            my ($locHLn, $lenOvlp) = posOvlp($locHL, $locHW);
            next if $lenOvlp == 0;
            $locHLn = [ sort {$a->[0] <=> $b->[0]} @$locHLn ];
            my ($begHL, $endHL) = ($locHLn->[0]->[0], $locHLn->[-1]->[1]);
            my $locHGn = [];
            for (@$locHLn) {
                my ($begHL1, $endHL1) = @$_;
                my $begHG1 = coordTransform($begHL1, $locHWI, "+", $locHWO, $srd);
                my $endHG1 = coordTransform($endHL1, $locHWI, "+", $locHWO, $srd);
                ($begHG1, $endHG1) = ($endHG1, $begHG1) if $begHG1 > $endHG1;
                push @$locHGn, [$begHG1, $endHG1];
            }
            $locHGn = [ sort {$a->[0] <=> $b->[0]} @$locHGn ];
            ($begH, $endH) = ($locHGn->[0]->[0], $locHGn->[-1]->[1]);
            
            my ($locQn, $locHLn2, $locAn, $alnQn, $alnHn, $alnPn) = 
                trim_aln_loc($locHL, $begHL, $endHL, $srdH, $locQ, $srdQ, $locA, $alnQ, $alnH, $alnP);
            my ($begQn, $endQn) = ($locQn->[0]->[0], $locQn->[-1]->[1]);

            my $frame = ($srd eq "+") ? ($begH-1) % 3 : (seqLen($seqid, $f_ref) - $endH) % 3;
            my $idHn = join("_", $seqid, $srd, $frame);
            my $idn = sprintf("%s_%d", $id, $j+1);

            my ($locQSn, $locHSn, $locASn) = map {locAry2Str($_)} ($locQn, $locHGn, $locAn);
            print FH join("\t", $idn, $idQ, $begQn, $endQn, $srdQ, $locQSn, $idHn, $begH, $endH, $srd, $locHSn, $e, $score, $src, $locASn, $alnQn, $alnHn, $alnPn)."\n";
        }
    }
    $log->info(sprintf "  %d in total", $t->nofRow);
    close FH;
}

sub refine_hmm_hit {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Hmm");
    $log->info("refining HMM hits");
    
    my $t = readTable(-in=>$fi, -header=>1);
    my @idxs_rm;
    for my $i (0..$t->nofRow-1) {
        my ($id, $idQ, $begQ, $endQ, $srdQ, $locQS, $idH, $begH, $endH, $srdH, $locHS, $e, $score, $src, $locAS, $alnQ, $alnH, $alnP) = $t->row($i);
        my ($locQ, $locH, $locA) = map {locStr2Ary($_)} ($locQS, $locHS, $locAS);
        
        next if $src ne "x";
        my $seqP = $alnP;
        $seqP =~ s/[\-\.]//g;
        my @locs;
        while($seqP =~ /([9\*]+)/g) {
            my ($b, $e) = ($-[1]+1, $+[1]);
            my ($bl, $el) = (3*$b-2, 3*$e);
            my ($bg, $eg) = $srdH eq "+" ? ($begH+$bl-1, $begH+$el-1) : ($endH-$el+1, $endH-$bl+1);
            my ($locO, $lenO) = posOvlp($locH, [[$bg, $eg]]);
            push @locs, $locO->[0] if $lenO > 0;
        }
        @locs = sort {$a->[1]-$a->[0] <=> $b->[1]-$b->[0]} @locs;
        if(@locs == 0) {
            push @idxs_rm, $i;
            next;
        }

        my ($begHn, $endHn) = @{$locs[-1]};
        my ($locQn, $locHn, $locAn, $alnQn, $alnHn, $alnPn) = 
            trim_aln_loc($locH, $begHn, $endHn, $srdH, $locQ, $srdQ, $locA, $alnQ, $alnH, $alnP);
        my ($begQn, $endQn) = ($locQn->[0]->[0], $locQn->[-1]->[1]);
        die Dumper($begHn, $endHn, $locQ, $locH, $locA) if !@$locQn || !@$locHn || !@$locAn;
        
        $t->setElm($i, "begQ", $begQn);
        $t->setElm($i, "endQ", $endQn);
        $t->setElm($i, "begH", $begHn);
        $t->setElm($i, "endH", $endHn);
        $t->setElm($i, "locQ", locAry2Str($locQn));
        $t->setElm($i, "locH", locAry2Str($locHn));
        $t->setElm($i, "locA", locAry2Str($locAn));
        $t->setElm($i, "alnQ", $alnQn);
        $t->setElm($i, "alnH", $alnHn);
        $t->setElm($i, "alnP", $alnPn);
    }
    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
    print FH $t->tsv(1);
    close FH;
}


sub get_longest_cds {
    my ($loc, $srd) = @_;
    my @phases = getPhase($loc, $srd);

    my $idx_max = 0;
    my $len_max = 0;
    for my $i (0..(@$loc-1)) {
        my ($beg, $end) = @{$loc->[$i]};
        my $len = $end - $beg + 1;
        if($len > $len_max) { 
            $idx_max = $i;
            $len_max = $len;
        }
    }

    my $phase = $phases[$idx_max];
    my ($beg, $end) = @{$loc->[$idx_max]};
    if($srd eq "-") {
        $end -= $phase;
        my $mod = ($end - $beg + 1) % 3;
        $beg += $mod;
    } else {
        $beg += $phase;
        my $mod = ($end - $beg + 1) % 3;
        $end -= $mod;
    }
    return ($beg, $end);
}
sub coordinate_transform_cds {
    my ($locR, $srdR, $locO, $srdO) = @_;
    my $locI;
    my $srdI = "+";
    $locO = [ sort {$a->[0] <=> $b->[0]} @$locO ];
    $locO = [ reverse @$locO ] if $srdO =~ /^\-1?$/;
    my $cdsLen = 0;
    for (@$locO) {
        my $len = $_->[1] - $_->[0] + 1;
        push @$locI, [$cdsLen + 1, $cdsLen + $len];
        $cdsLen += $len;
    }
    my ($locR2, $oLen) = posOvlp($locR, $locI);
    my ($sPosAry1, $sLen1) = posDiff($locR, $locR2);
    die "not completely in CDS:\n".Dumper($locR).Dumper($locO) unless $sLen1 == 0;
    
    my $srdG = is_opposite_strands($srdI, $srdO) ? get_opposite_strand($srdR) : $srdR;
    my $locG;
    for (@$locR2) {
        my $pos1 = coordTransform($_->[0], $locI, $srdI, $locO, $srdO);
        my $pos2 = coordTransform($_->[1], $locI, $srdI, $locO, $srdO);
        push @$locG, [min($pos1, $pos2), max($pos1, $pos2)];
    }
    return ($locG, $srdG);
}

sub pipe_hmmsearch {
    my ($dir, $f_hmm, $f_tar, $f_ref, $f_gtb) = rearrange([qw/dir hmm target ref gtb/], @_);
    make_path($dir) unless -d $dir;
    my $log = Log::Log4perl->get_logger("Hmm");

    my $f_bin = $ENV{"HMMER"}."/bin/hmmsearch";
    $log->error_die("cannot execute $f_bin") unless -s $f_bin;

    my $f01 = "$dir/01_raw.txt";
    $log->info("running hmmsearch against $f_tar");
    runCmd("$f_bin --domE 10 --cpu 4 -o $dir/01_raw.txt --domtblout $dir/01_raw.tbl $f_hmm $f_tar", -1);

    my $f02 = "$dir/02.htb";
    parse_hmm_output($f01, $f02);
    my $f03 = "$dir/03_full.htb";
    get_aln_position($f02, $f03);
    my $f04 = "$dir/04_nt.htb";
    coordTransform_aa2nt($f03, $f04);
    my $f05 = "$dir/05_global.htb";
    recoverCoord($f04, $f05, $f_ref);
    my $f07 = "$dir/07_final.htb";
    refine_hmm_hit($f05, $f07);
}



1;
__END__
