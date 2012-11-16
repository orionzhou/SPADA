package Gtb;
use strict;
use Common;
use Seq;
use Gene;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/readGtb gtbSum 
    gtb2Gff gtb2Bed gtb2Seq gtb2Tbl
    gtb_fix_phase pep_score_gtb/;
@EXPORT_OK = qw//;
=cut
print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note/)."\n";
my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = @$_;
=cut
sub readGtb {
    my ($fi, $opt) = rearrange(['in', 'opt'], @_);
    $opt ||= 1;
    die "fi doesn't exist: $fi\n" unless -s $fi;
    my $t = readTable(-in=>$fi, -header=>1);
    my $rst;
    if($opt == 1) {
        $rst = $t;
    } elsif($opt == 2) {
        $rst = { map {$t->elm($_, "id") => $t->rowRef($_)} (0..$t->nofRow-1) };
    } else {
        die "unsupported option: $opt\n";
    }
    return $rst;
}
sub gtbSum {
    my ($fi, $opt) = @_;
    $opt = 1 unless defined($opt);
    my $t = readTable(-in=>$fi, -header=>1);
    $t->sort("parent", 1, 0, "id", 1, 0);
    my @idGs = $t->col("parent");
    my $ref = group(\@idGs);
    @idGs = sort {$ref->{$a}->[0] <=> $ref->{$b}->[0]} keys %$ref;
    printf "%d genes : %d models\n", scalar(@idGs), $t->nofRow if $opt == 1;
    my @cnts = uniq(map {$_->[1]} values %$ref);
    for my $cnt (@cnts) {
        my @ids = grep {$ref->{$_}->[1] == $cnt} keys %$ref;
        printf "\t%5d x $cnt model(s)/gene [%s]\n", scalar(@ids), @ids>=0 ? $ids[0] : "" if $opt == 1;
    }
    return ($t, $ref, \@idGs);
}

sub gtb2Seq { # opt = 1 (add a 'seq' column at the end); 2 (write to fasta file)
    my ($fi, $fo, $f_seq, $opt) = rearrange(['in', 'out', 'seq', 'opt'], @_);
    $opt ||= 1;
    my $t = readTable(-in=>$fi, -header=>1);
    my $seqH;
    if($opt == 1) {
        my $idx_seq = first_index {$_ eq "seq"} $t->header;
        $t->addCol([("") x $t->nofRow], "seq") if $idx_seq < 0;
    } else {
        die "unknown opt: $opt\n" if $opt != 2;
        $seqH = Bio::SeqIO->new(-file=>">$fo", -format=>"fasta");
    }

    for my $i (0..$t->nofRow-1) {
        my ($id, $pa, $chr, $strand, $phase, $locS, $cat1, $cat2) = 
            map {$t->elm($i, $_)} qw/id parent chr strand phase locC cat1 cat2/;
        next if $cat2 ne "mRNA";
        die "no locCDS for $id\n" unless $locS;
        my $phase1 = [split(",", $phase)]->[0];

        my $loc = locStr2Ary($locS);
        my $seqStr = seqRet($loc, $chr, $strand, $f_seq);
        my $seq_cds = Bio::Seq->new(-id=>$id, -seq=>$seqStr);
        my $seq_pro = $seq_cds->translate(-frame=>$phase1);
        if($opt == 1) {
            $t->setElm($i, "seq", $seq_pro->seq);
        } else {
            $seqH->write_seq($seq_pro);
        }
        printf "  extracting sequence from Gtb... %5d out of %d done\r", $i+1, $t->nofRow;
    }
    print "\n";

    if($opt == 1) {
        open(FH, ">$fo");
        print FH $t->csv(1, {delimiter=>"\t"});
        close FH;
    }
}
sub gtb2Gff {
    my ($fi, $fo) = @_;
    my ($t, $ref, $idGs) = gtbSum($fi, 0);
    open(FH, ">$fo");
    print FH "##gff-version 3\n";
    my ($cntG, $cntR) = (1, 1);
    for my $idG (@$idGs) {
        my ($idxB, $cnt) = @{$ref->{$idG}};
        my $ts = $t->subTable([$idxB..$idxB+$cnt-1]);
        my $gene = Gene->new( -gtb=>$ts );
        print FH $gene->to_gff()."\n";
        for my $rna ($gene->get_rna()) {
            print FH $rna->to_gff()."\n";
            printf "  converting Gtb to Gff... ( %5d RNA | %5d gene ) done\r", $cntR++, $cntG;
        }
        $cntG ++; 
        print FH "\n";
    }
    print "\n";
    close FH;
}
sub gtb2Bed {
    my ($fi, $fo) = @_;
    my ($t, $ref, $idGs) = gtbSum($fi);
    open(FH, ">$fo");
    print FH "#track name=gene_models useScore=0\n";
    my ($cntG, $cntR) = (1, 1);
    for my $idG (@$idGs) {
        my ($idxB, $cnt) = @{$ref->{$idG}};
        my $ts = $t->subTable([$idxB..$idxB+$cnt-1]);
        my $gene = Gene->new( -gtb=>$ts );
        for my $rna ($gene->get_rna()) {
            my $idStr = $rna->id;
            $idStr .= $rna->note if $rna->note;
            printf "  converting Gtb to Bed... ( %5d RNA | %5d gene ) done\r", $cntR++, $cntG;
            print FH join("\t", $rna->seqid, $rna->beg-1, $rna->end, $idStr, 0, $rna->strand)."\n";
        }
        $cntG ++; 
    }
    print "\n";
    close FH;
}
sub gtb2Tbl {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    my @chrs = uniq($t->col("chr"));

    open(FH, ">$fo");
    print FH join("\t", qw/chr beg end type id/)."\n";
    my $h_type = {1=>"cds", 2=>"utr5", 3=>"utr3", 4=>"intron"};
    for my $chr (@chrs) {
        my $t2 = $t->match_pattern("\$_->[2] eq '$chr'");
        my @locs;
        my @stats;
        for my $i (0..$t2->nofRow-1) {
            my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $t2->row($i);
            $locC = locStr2Ary($locC);
            push @locs, @$locC;
            push @stats, ([1, $id]) x @$locC;
            
            if($loc5) {
                $loc5 = locStr2Ary($loc5);
                push @locs, @$loc5;
                push @stats, ([2, $id]) x @$loc5;
            }
            if($loc3) {
                $loc3 = locStr2Ary($loc3);
                push @locs, @$loc3;
                push @stats, ([3, $id]) x @$loc3;
            }
            if($locI) {
                $locI = locStr2Ary($locI);
                push @locs, @$locI;
                push @stats, ([4, $id]) x @$locI;
            }
        }
        my @stats_score = map {$_->[0]} @stats;
        my $ref = tiling(\@locs, \@stats_score, 1);
        for (@$ref) {
            my ($beg, $end, $idx) = @$_;
            my ($typeNum, $id) = @{$stats[$idx]};
            print FH join("\t", $chr, $beg, $end, $h_type->{$typeNum}, $id)."\n";
        }
    }
    close FH;
}
  

sub gtb_fix_phase {
    my ($fi, $fo, $f_ref) = @_;
    my $t = readTable(-in=>$fi, -header=>1);

    my $n_fixed = 0;
    for my $i (0..$t->nofRow-1) {
        my ($id, $pa, $chr, $beg, $end, $strand, $locE, $locI, $locC, $loc5, $loc3, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = $t->row($i);
        my $loc = locStr2Ary($locC);
        my $len = locAryLen($loc);
        my $seqStr = seqRet($loc, $chr, $strand, $f_ref);
        my $seq = Bio::Seq->new(-id=>"test", -seq=>$seqStr);
        my $frame = -1;
        for my $i (0..2) {
            my $prot = $seq->translate(-frame=>$i)->seq;
            if($prot =~ /^[A-Z]+\*?$/i) {
                $frame = $i;
                last;
            }
        }
        if($frame < 0) {
            printf "%s: cannot find frame\n", $id;
        } elsif($frame > 0) {
            $n_fixed ++;
            my @phases_old = split(",", $phase);
            my @phases_new = map {($frame + $_) % 3} @phases_old;
            $t->setElm($i, "phase", join(",", @phases_new));
        }
        printf "  fixing Gtb phases... ( %5d out of %5d done )\r", $i+1, $t->nofRow;
    }
    print "  \n$n_fixed non-0 frames fixed\n";
  
    open(FH, ">$fo") or die "cannot open $fo to write\n";
    print FH $t->csv(1, {delimiter=>"\t"});
    close FH;
}

sub pep_score_gtb {
    my ($fi, $fo) = @_;
    my $tg = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\t", qw/id codonStart codonStop preStop gap n_cds lenC lenI/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $chr, $strandG, $locCStr, $locIStr, $phaseG, $seq) = 
            map {$tg->elm($i, $_)} qw/id parent chr strand locC locI phase seq/;
        my $locCAry = locStr2Ary($locCStr);
        my $locIAry = locStr2Ary($locIStr);

        my ($codonStart, $codonStop, $preStop, $gap) = checkProtSeq($seq);
        
        my $n_cds = @$locCAry;
        my $lenC = locAryLen($locCAry);
        my $lenI = locAryLen($locIAry);
        print FH join("\t", $id, $codonStart, $codonStop, $preStop, $gap, $n_cds, $lenC, $lenI)."\n";
        printf "  assessing peptide scores: %5d / %5d done...\r", $i+1, $tg->nofRow;
    }
    print "\n";
    close FH;
}


1;
__END__
