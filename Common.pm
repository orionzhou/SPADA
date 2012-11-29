package Common;
use strict; 
use Bio::Seq;
use Data::Table;
use File::Basename;
use Data::Dumper;
use IPC::Open3; use Symbol qw/gensym/;
use Clone qw/clone/;
use Graph;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pretty tonum isnumber getDigits prettyStr extractAcc
    runCmd runCmd2 parse_gff_tags
    locStr2Ary locAry2Str locAryLen trimLoc cropLoc cropLoc_cds parse_old_loc_str
    posOvlp posCmp posMerge posSplit posDiff posMergeDeep
    rearrange group bsearch tiling readTable
    mergeArray aryCmp get_sliding_windows 
    get_opposite_strand is_opposite_strands coordTransform coordTransform_rough coordTransform_itv
    getIdxRange getPhase sample_serial
    rmRedPairs backOneLine scaleNumber writeFile/;
@EXPORT_OK = qw//;
sub rearrange {
    my( $order, @param ) = @_;
    return unless @param;
    my %param;
    if (ref $param[0] eq 'HASH') {
        %param = %{$param[0]};
    } else {
        return @param unless (defined($param[0]) && substr($param[0],0,1) eq '-');
        my $i;
        for ($i=0;$i<@param;$i+=2) {
            $param[$i]=~s/^\-//;  # get rid of initial - if present
            $param[$i]=~tr/a-z/A-Z/; # parameters are upper case
        }
        %param = @param; # convert into associative array
    }
    my (@return_array);
    local($^W) = 0;
    my ($key)='';
    foreach $key (@$order) {
        my ($value);
        if (ref($key) eq 'ARRAY') {
            foreach (@$key) {
            last if defined($value);
            $value = $param{uc($_)};
            delete $param{$_};
        }
        } else {
            $value = $param{uc($key)};
            delete $param{$key};
        }
        push(@return_array,$value);
    }
    push (@return_array,\%param) if %param;
    return @return_array;
}
sub pretty {
    my ($num, $sep) = @_;
    $sep = defined $sep ? $sep: ",";
    die "not a number: $num\n" if $num !~ /^\d+$/;
    $num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1$sep/g;
    return $num;
}
sub prettyStr {
    my ($str, $len, $sep) = @_;
    $len ||= 5;
    $sep ||= " ";
    my @ary;
    for(my $i = 0; $i*$len < length($str); $i++) {
        push @ary, substr($str, $i*$len, $len);
    }
    return join($sep, @ary);
}
sub tonum {
    my ($num) = @_;
    $num =~ s/[^\d]//g;
    return $num;
}
sub getDigits {
    my ($num) = @_;
    my $digit = 1;
    while(int($num/10) >= 1) {
        $num = int($num/10);
        $digit ++;
    }
    return $digit;
}
sub isnumber {
    my ($a) = @_;
    if( $a =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
        return 1;
    } else {
        return 0;
    }
}
sub runCmd {
    my ($cmd, $opt) = @_;
    $opt = defined($opt) ? $opt : 1;

    local *CATCHERR = IO::File->new_tmpfile;
    my $pid = open3(gensym, \*CATCHOUT, ">&CATCHERR", $cmd);
    if($opt == 1) {
        while(<CATCHOUT>) { print $_; }
    }
    waitpid($pid, 0);
    my $exit_status = $?;
    
    if($exit_status != 0 && $opt != -1) {
        seek CATCHERR, 0, 0;
        while(<CATCHERR>) { print $_; }
        die "!!!!! failed system call !!!!!\n$cmd\n";
    }
}
sub runCmd2 {
    my ($cmd) = @_;
    my @lines;
    open(my $out, "$cmd | ") or die "Failed: $!\n";
    while( <$out> ) {
        chomp;
        push @lines, $_;
    }
    return \@lines;
}
sub parse_gff_tags {
    my ($str) = @_;
    my @tagAry = split(";", $str);
    my $h;
    for (@tagAry) {
        my ($tag, $value) = split "=";
        die $str if $value eq "";
        $value =~ s/\=/\:/g;
        $value =~ s/\;/\|/g;
        $h->{$tag} = $value;
    }
    return $h;
}
sub extractAcc {
    my ($str) = @_;
    if($str =~ /\|?([A-Z]{2}\d{6}\.[0-9DF]{1,2})\|?/) {
        return $1;
    } elsif($str =~ /([A-Z]{2}\d{6})([^\.]|$)/) {
        return $1;
    } elsif($str =~ /^contig\_\d+$/) {
        return $str;
    } elsif($str =~ /^(fpc265)|(mtF83-89g3-P2)$/) {
        return $str;
    } else {
        die "cannot get a valid acc from $str\n";
    }
}
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
sub mergeArray {
    my @refAry = @_;
    my @rst;
    for my $ref (@refAry) {
        if(ref($ref) eq "ARRAY") {
            push @rst, mergeArray(@$ref);
        } elsif(!ref($ref)) {
            push @rst, $ref;
        } else{
            die("unknown type $ref\n");
        }
    }
    return @rst;
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
    return ($ref, $dLen);
}
sub tiling {
    #opt = 1(minimal) / 2(maximal)
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
sub aryCmp {
    my ($ary1, $ary2) = @_;
    my @ary1 = uniq(@$ary1);
    my @ary2 = uniq(@$ary2);
    my @ary_merged = uniq(@ary1, @ary2);
    my %hash1 = map {$_=>1} @ary1;
    my %hash2 = map {$_=>1} @ary2;
    my (@ary_share, @ary_1, @ary_2);
    for my $ele (@ary_merged) {
        if(exists $hash1{$ele} && exists $hash2{$ele}) {
            push @ary_share, $ele;
        } elsif(exists $hash1{$ele} && !exists $hash2{$ele}) {
            push @ary_1, $ele;
        } elsif(!exists $hash1{$ele} && exists $hash2{$ele}) {
            push @ary_2, $ele;
        } else {
            die "$ele not in ary1 nor ary2\n";
        }
    }
    return (\@ary_share, \@ary_1, \@ary_2);
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
sub readTable {
    my ($fi, $header, $skip) = rearrange(['in', 'header', 'skip'], @_);
    die "$fi is not there\n" unless -s $fi;
    $skip ||= 0;
    $header ||= 0;
    open(my $fh, "<$fi") or die "cannot open $fi\n";
    my $t;
    if($header ne 1 && $header ne 0) {
        $t = Data::Table::fromTSV($fh, 0, $header, {skip_lines=>$skip, skip_pattern=>'(^\s*#)|(^\s*$)'});
    } else {
        $t = Data::Table::fromTSV($fh, $header, {skip_lines=>$skip, skip_pattern=>'(^\s*#)|(^\s*$)'});
    }
#  printf "\t%s: cols[%d] rows[%d]\n", basename($fi), $t->nofCol, $t->nofRow;
    close $fh;
    return $t;
}
sub scaleNumber {
    my ($ary, $down, $up) = rearrange(["value", "down", "up"], @_);
    my $rst;
    ($down, $up) = $down<=$up ? ($down, $up) : ($up, $down);
    my $minX = min(@$ary);
    my $maxX = max(@$ary);
    $down ||= 0;
    $up ||= 100;
    my $factor = ($up-$down) / ($maxX-$minX);
    for my $x (@$ary) {
        my $y = $down + ($x-$minX) * $factor;
        push @$rst, $y;
    }
    return $rst;
}
sub backOneLine {
    my ($fH) = @_;
    my $char;
    my $flag_nonblank = 0;
    while(seek($fH, -1, 1)) {
        read($fH, $char, 1);
        last if ($char eq "\n" && $flag_nonblank == 1);
        $flag_nonblank = 1 if $char ne "\n";
        seek($fH, -1, 1);
    }
    return tell($fH);
}
sub get_sliding_windows {
    my ($beg, $end, $step, $size) = @_;
    my $n_win = ceil(($end-$beg-$size+1)/$step) + 1;
    $n_win = max(1, $n_win);

    my @wins;
    for my $i (0..$n_win-1) {
        my $begL = $step * $i + 1;
        my $endL = $step * $i + $size;
        $endL = min($endL, $end);
        push @wins, [$begL, $endL];
    }
    return \@wins;
}
sub get_opposite_strand {
    my ($srdI) = @_;
    return "-" if $srdI =~ /^[\+1]$/;
    return "+" if $srdI =~ /^\-1?$/;
    die "unknonw strand: $srdI\n";
}
sub is_opposite_strands {
    my ($srd1, $srd2) = @_;
    $srd1 = 1 if $srd1 eq "+";
    $srd1 = -1 if $srd1 eq "-";
    $srd2 = 1 if $srd2 eq "+";
    $srd2 = -1 if $srd2 eq "-";
    return 1 if $srd1 * $srd2 == -1;
    return 0 if $srd1 * $srd2 == 1;
    die "unknown strands: $srd1  $srd2\n";
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

sub rmRedPairs {
    my ($edgeRef) = @_;
    my $g = Graph::Undirected->new();
    $g->add_edges(map {@$_} @$edgeRef);
    my @pairs = ();
    my @es;
    my @vs = $g->vertices();
    my $cnt = 0;
    while( @pairs < @vs/2 ) {
        die Dumper($edgeRef) if ++$cnt > 10;
        @es = $g->edges();
#    print join("\t", map {join("-", @$_)} @es)."\n";
        @pairs = ();
        for my $e (@es) {
            my ($v1, $v2) = @$e;
            push @pairs, $e if $g->degree($v1) == 1 || $g->degree($v2) == 1;
        }
        if(@pairs < @vs/2) {
            my $v2 = first_value {$g->degree($_) == 2} @vs;
            $g->delete_edge(@{[$g->edges_at($v2)]->[0]});
        }
    }
    print @pairs." pairs obtained for @vs vertices\n" unless @vs/2 == @pairs;
    return \@pairs;
}
sub group {
    my ($ary) = @_;
    my $ref = {};
    for my $i (0..@$ary-1) {
        my $ele = $ary->[$i];
        if(exists $ref->{$ele}) {
            $ref->{$ele}->[1] ++;
        } else {
            $ref->{$ele} = [$i, 1];
        }
    }
    return $ref;
}
sub bsearch {
    my ($ary, $word) = @_;
    my $low = 0;
    my $high = @$ary - 1;
    my $try;
    while( $low <= $high ) {
        $try = int( ($low + $high) / 2 );
        die $word."\n".join("\t", @$ary)."\n" if $try < 0 || $try > @$ary - 1;
        $low  = $try + 1, next if $ary->[$try] < $word;
        $high = $try - 1, next if $ary->[$try] > $word;
        return $try;
    }
    return $try;
}
sub getIdxRange {
    my @a = @_;
    my @idx;
    my $str = join("", @a);
    my $i = 0;
    while($str =~ /0(1*)/g) {
        my $len = length($1);
        push @idx, [$i, $i+$len];
        $i += $len + 1;
    }
    return @idx;
}
sub getPhase {
    my ($loc, $srd) = @_;
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    $loc = [ reverse @$loc ] if $srd =~ /^\-1?$/;

    my @phases;
    my $len = 0;
    for (@$loc) {
        my ($beg, $end) = @$_;
        push @phases, (3-$len%3) % 3;
        $len += $end - $beg + 1;
    }
    return @phases;
}
sub sample_serial {
    my ($n, $m) = @_;
    my $inc = $n / $m;
    return map {1+$inc*$_-1} (0..$m-1);
}
sub parse_old_loc_str {
    my ($locS) = @_;
    my $srd;
    my $loc = [];
    if($locS =~ /^complement\(([\w\.\,]+)\)$/) {
        $srd = "-";
        $locS = $1;
    } else {
        $srd = "+";
    }
    while($locS =~ /(\d+)\.\.(\d+)/g) {
        push @$loc, [$1, $2];
    }
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    return ($loc, $srd);
}
sub writeFile {
    my ($fo, @strs) = @_;
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\n", @strs)."\n";
    close FH;
}



1;
__END__
