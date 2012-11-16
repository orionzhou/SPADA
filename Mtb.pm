package Mtb;
use strict;
use Data::Dumper;
use Bio::SearchIO;
use Common;
use Seq;
use Graph::Undirected;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/run_blast psl2Mtb bls2Mtb mtbFilter mtbExpand mtb2Gtb mtbRmDup/;
@EXPORT_OK = qw//;

sub run_blast {
    my ($f_qry, $fo, $param) = rearrange(["qry", "out", "param"], @_);
    die "$f_qry is not there\n" unless -s $f_qry;
    
    my ($db, $program, $e) = map {$param->{$_}} qw/db program e/;
    die "database not provided\n" unless $db;
    
    my @programs = qw/blastn blastp blastx tblastn tblastx/;
    my $idxProgram = first_index {$_ eq $program} @programs;
    die "program[$program] not supported\n" if $idxProgram == -1;

    $e ||= 0.01;
    runCmd("blastall -p $program -d $db -i $f_qry -o $fo -e $e -m 0 -F F -a 4", 1);
}

sub psl2Mtb {
    my ($fi, $fo) = @_;
#  my @cols = qw/matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts/;
    open(FHI, "<$fi");
    open(FHO, ">$fo");
    print FHO join("\t", qw/id idQ locQ srdQ gapQ idH locH srdH gapH matches alnLenQ lenQ/)."\n";
    my $id = 1;
    while(<FHI>) {
        chomp;
        my @ps = split " ";
        next unless @ps == 21;
        next if /^(psLayout)|(match)/;
        my ($match, $match_mis, $match_rep, $nN, $gapQ, $bpGapQ, $gapH, $bpGapH, $srdH, $idQ, $lenQ, $begQ, $endQ, $idH, $lenH, $begH, $endH, $n, $lens, $begsQ, $begsH) = @ps;
        my $srdQ = "+";
        my @begsQ = split(",", $begsQ);
        my @begsH = split(",", $begsH);
        my @lens = split(",", $lens);
        die "qry[$idQ] hit[$idH] not $n blocks: $begsQ\n" unless $n == @begsQ;
        die "qry[$idQ] hit[$idH] not $n blocks: $begsH\n" unless $n == @begsH;
        my ($locQ, $locH) = ([], []);
        for my $i (0..$n-1) {
            my ($begH, $len) = ($begsH[$i]+1, $lens[$i]);
            my $endH = $begH + $len - 1;
            
            my ($begQ, $endQ);
            if($srdH eq "+") {
                $begQ = $begsQ[$i] + 1;
                $endQ = $begQ + $len - 1;
            } else {
                die "unknown strand $srdH\n" unless $srdH eq "-";
                $endQ = $lenQ - $begsQ[$i];
                $begQ = $endQ - $len + 1;
            }
            push @$locQ, [$begQ, $endQ];
            push @$locH, [$begH, $endH];
        }
        $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
        $locH = [ reverse sort {$a->[0] <=> $b->[0]} @$locH ] if $srdH eq "-";;
        my ($locQStr, $locHStr) = map {locAry2Str($_)} ($locQ, $locH);
        my $alnLenQ = sum(@lens) + $bpGapQ;
        print FHO join("\t", $id++, $idQ, $locQStr, $srdQ, $bpGapQ, $idH, $locHStr, $srdH, $bpGapH, $match, $alnLenQ, $lenQ)."\n";
    }
    close FHI;
    close FHO;
}
sub bls2Mtb {
    my ($fi, $fo) = @_;
    my $in = Bio::SearchIO->new(-file=>$fi, -format=>'blast');
    open(FH, ">$fo");
    print FH join("\t", qw/id idQ locQ srdQ gapQ idH locH srdH gapH matches alnLenQ lenQ conserved e/)."\n";
    my $cnt = 0;
    while( my $rst = $in->next_result ) {
        my ($qId, $qLen) = ($rst->query_name, $rst->query_length);
        while( my $hit = $rst->next_hit ) {
            my ($hId, $hLen) = ($hit->name, $hit->length);
            while( my $hsp = $hit->next_hsp ) {
                my ($qS, $qE, $qL, $qG, $qStr, $hS, $hE, $hL, $hG, $hStr) =  
                    ($hsp->query->start, $hsp->query->end, $hsp->length('query'), $hsp->gaps('query'), 
                    $hsp->query->strand, $hsp->hit->start, $hsp->hit->end, $hsp->length('hit'), $hsp->gaps('hit'), 
                    $hsp->hit->strand);
                my $id = sprintf("%04d", ++$cnt);
                $hStr = $hStr==-1 ? "-" : "+";
                $qStr = $qStr==-1 ? "-" : "+";
                my $locQStr = "$qS-$qE"; 
                my $locHStr = "$hS-$hE"; 
                print FH join("\t", $id, $qId, $locQStr, $qStr, $qG, $hId, $locHStr, $hStr, $hG, $hsp->num_identical, $qL, $qLen, $hsp->num_conserved, $hsp->evalue)."\n";
            }
        }
    }
    close FH;
}

sub mtbFilter {
    my ($fi, $fo, $p) = @_;
    my ($lenCov, $pctCov, $pctIdty, $pctCnsv, $gap, $strand, $best) = 
        map {$p->{$_}} qw/lencov pctcov pctidty pctcnsv gap strand best/;
    $pctCov ||= 0.1;
    $lenCov ||= 10;
    $pctIdty ||= 0.3;
    $best ||= 0;
    $gap ||= 10_000;
    my $t = readTable(-in=>$fi, -header=>1);
    my @rows;
    for my $i (0..$t->nofRow-1) {
        my ($qLen, $qAlnLen, $nIdty, $qGap) = map {$t->elm($i, $_)} qw/lenQ alnLenQ matches gapQ/;
        my $flag = 1;
        $flag = 0 if $t->elm($i, "idQ") eq $t->elm($i, "idH");
        $flag = 0 if $qAlnLen < $lenCov;
        $flag = 0 if $qAlnLen / $qLen < $pctCov;
        $flag = 0 if $nIdty / $qAlnLen < $pctIdty;
        $flag = 0 if $strand && $t->elm($i, "srdH") ne $strand;
        $flag = 0 if $qGap > $gap;
        if($pctCnsv) {
            my $qCnsv = $t->elm($i, "conserved");
            $flag = 0 if $qCnsv / $qAlnLen < $pctCnsv;
        }
        push @rows, $i if $flag == 0;
    }
    $t->delRows(\@rows);
    my $qIds = $t->colRef("idQ");
    my $ref = group($qIds);
    if($best == 1) {
        my @rows2;
        for my $qId (keys %$ref) {
            my ($idx, $nofId) = @{$ref->{$qId}};
            next if $nofId == 1;
            my @idxs = ($idx..$idx+$nofId-1);
            my @tmp = map {$t->elm($_, "idQ")} @idxs;
            die "not 1 qId\n" unless (uniq(@tmp)) == 1;
            my %ref = map { $_ => $t->elm($_, 'matches') } @idxs;
            my $maxM = max( values %ref );
            for my $idx (@idxs) {
                push @rows2, $idx if $ref{$idx} < $maxM;
            }
            my @chrs = map {$t->elm($_, "idH")} @idxs;
            my @chrs2 = grep /^chr[1-8]$/, @chrs;
            if(@chrs2) {
                for my $idx (@idxs) {
#          push @rows2, $idx if $t->elm($idx, "hId") =~ /^chr[0TU]$/;
                }
            }
        }
        @rows2 = uniq(@rows2);
        $t->delRows(\@rows2);
    }
    printf "%d queries --> to %d hits\n", scalar(keys %$ref), $t->nofRow;
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub mtbRmDup {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    my @srdQs = uniq($ti->col("srdQ"));
    die "qry not all + strand\n" unless @srdQs == 1 && $srdQs[0] eq "+";
    $ti->sort("idQ", 1, 0, "idH", 1, 0);
    my $ref = group($ti->colRef("idQ"));
    my @qId_dup = grep {$ref->{$_}->[1] > 1} keys(%$ref);

    my $he;
    for my $idU (@qId_dup) {
        my ($idx, $cnt) = @{$ref->{$idU}};
        for my $i (0..$cnt-1) {
            my $idV = join("_", "@", map {$ti->elm($idx+$i, $_)} qw/srdH idH locH/);
            my $idE = join(" ", $idU, $idV);
            $he->{$idE} = [$idU, $idV, $idx+$i];
        }
    }
    my @idx_es_rm;
    my $g = Graph::Undirected->new();
    $g->add_edges(map {($_->[0], $_->[1])} values(%$he));
    for my $cc ($g->connected_components()) {
        my @vs = @$cc;
        my @vis = grep /^[^\@]/, @vs;
        my @vos = grep /^\@/, @vs;
        if(@vis != @vos) {
            printf "%d queries [%s] => %d hits [%s]\n", scalar(@vis), join(", ", @vis), scalar(@vos), join(", ", @vos);
        } else {
            my %idE_uniq = map {$vis[$_]." ".$vos[$_] => 1} (0..$#vis);
            for my $vi (@vis) {
                for my $vo (@vos) {
                    my $idE = "$vi $vo";
                    unless(exists $he->{$idE}) {
                        printf "%d queries [%s] => %d hits [%s]\n", scalar(@vis), join(", ", @vis), scalar(@vos), join(", ", @vos);
                        die "no edge: $vi => $vo\n";
                    }
                    push @idx_es_rm, $he->{$idE}->[2] unless exists $idE_uniq{$idE};
                }
            }
#      print join("\n", keys %idE_uniq)."\n\n";
        }
    }

    $ti->delRows(\@idx_es_rm);
    open(FH, ">$fo");
    print FH $ti->tsv(1);
    close FH;
}
sub mtbExpand {
    my ($fS, $fi, $fo) = rearrange(['fs', 'fi', 'fo'], @_);
    my $t = readTable(-in=>$fi, -header=>1);
    my $ref1 = readSeqInfo($fS, 1);
    my $ref2 = group($t->colRef("idQ"));
    $t->addCol([('') x $t->nofRow], "noteQ");
    for my $idQ (keys %$ref1) {
        my $noteQ = $ref1->{$idQ};
        if(exists $ref2->{$idQ}) {
            my ($idx, $n) = @{$ref2->{$idQ}};
            for my $i ($idx..$idx+$n-1) {
                $t->setElm($i, "noteQ", $noteQ);
            }
        } else {
            $t->addRow(["", $idQ, ("") x ($t->nofCol-3), $noteQ]);
        }
    }
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}

sub mtb2Gtb { # opt = 1 - multiple mappings allowed; 2 - multiple mappings NOT allowed
    my ($fi, $fo, $opt) = @_;
    $opt ||= 1;
    my $t = readTable(-in=>$fi, -header=>1);
  
    open(FH, ">$fo");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
    
    $t->sort('idQ', 1, 0);
    my $ref = group([$t->col("idQ")]);
    my $i = 0;
    for my $id (sort keys %$ref) {
        my ($idxS, $n) = @{$ref->{$id}};
        die "$id has $n mappings\n" if $n > 1 && $opt == 2;
        for my $j (0..$n-1) {
            my $idx = $idxS+$j;
            my ($mId, $qId, $qLoc, $qStr, $qGap, $chr, $hLoc, $hStr, $hGap, $matches, $qAlnLen, $qLen, $qDesc) = $t->row($idx);
            next if $chr eq "";
            my $idG = join("_", $id, $j+1);
            my $idM = join(".", $idG, 1);

            my $locC = locStr2Ary($hLoc);
            $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
            my ($beg, $end) = ($locC->[0]->[0], $locC->[-1]->[1]);
            $locC = [ reverse @$locC ] if $hStr =~ /^\-1?$/;
            my $locM = [[$beg, $end]];
            my ($locI) = posDiff($locM, $locC);
            my ($strC, $strI) = map {locAry2Str($_)} ($locC, $locI);
            
            my $cdsLen = 0;
            my @phases;
            for (@$locC) {
                my ($cdsBeg, $cdsEnd) = @$_;
                push @phases, (3 - $cdsLen % 3) % 3;
                $cdsLen += $cdsEnd - $cdsBeg + 1;
            }
            my $phase = join(",", @phases);

            $qDesc =~ s/\=/-/g;
            $qDesc =~ s/\;/./g;
            my $tag1 = ($j+1)."/$n";
            my $tag2 = "$matches/$qAlnLen/$qLen";
            my $note = "[$tag1][$tag2][$qDesc]";

            print FH join("\t", $idM, $idG, $chr, $beg, $end, $hStr, $strC, $strI, $strC, "", "", $phase, "", "", "gene", "mRNA", "", $note)."\n";
        }
    }
    close FH;
}



1;
__END__
