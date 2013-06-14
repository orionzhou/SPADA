package ModelEval;
use strict;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Log::Log4perl;
use Data::Dumper;
use Common;
use Align;
use Seq;
use Hmm;
use Gtb;
use SignalP;
use CompareModel;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_model_evaluation filter_models/;
@EXPORT_OK = qw//;

sub get_stat_basic {
    my ($fi, $f_hit, $fo) = @_;
    
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("extracting basic stats");
   
    my $th = readTable(-in=>$f_hit, -header=>1);
    my $he = { map {$th->elm($_, "id") => $th->elm($_, "e")} (0..$th->nofRow-1) };

    my $t = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id parent family e seq/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $pa, $fam, $seq) = map {$t->elm($i, $_)} qw/id parent cat3 seq/;
        die "no E-value for $id:$pa\n" unless exists $he->{$pa};
        my $e = $he->{$pa};
        print FH join("\t", $id, $pa, $fam, $e, $seq)."\n";
    }
    close FH;
}
sub get_stat_aln {
    my ($fi, $d_aln, $f_sta, $fo) = @_;
    my $do = "$fo.dir";
   
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("computing MSA scores");
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $f_bin = $ENV{"ClustalO"}."/bin/clustalo";
    -s $f_bin || $log->error_die("cannot execute $f_bin");
    
    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id score/)."\n";
  
    my $ts = readTable(-in=>$f_sta, -header=>1);
    my $hs;
    for my $i (0..$ts->nofRow-1) {
        my ($fam, $nseq, $npos, $gap, $seq) = $ts->row($i);
        next if $gap == 1;
        my $f_aln = "$d_aln/$fam.aln";
        my $h_seq = read_aln_seq($f_aln);
        my @seqs = map { Bio::Seq->new(-id=>$_, -seq=>$h_seq->{$_}) } keys(%$h_seq);
        $hs->{$fam} = \@seqs;
    }
    
    my $t = readTable(-in=>$fi, -header=>1);
    my $ref = group($t->colRef("parent"));
    my $i = 1;
    
    my $f_fas = $ENV{"TMP_DIR"}."/aln_score_".int(rand(1000)).".fa";
    for my $pa (sort keys(%$ref)) {
        my ($idx, $cnt) = @{$ref->{$pa}};
        my $fam = $t->elm($idx, "cat3");
        my (@ids, @seqs);
        for my $i ($idx..$idx+$cnt-1) {
            my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
            push @ids, $id;
            push @seqs, Bio::Seq->new(-id=>$id, -seq=>$seq);
            die "family not consistent for gene[$pa]\n" unless $fam eq $fam2;
        }
        
        my $f_ao = "$do/$pa.fas";
        if(exists $hs->{$fam}) {
            writeSeq([@seqs, @{$hs->{$fam}}], $f_fas);
            runCmd("$f_bin -i $f_fas --outfmt=fasta --force -o $f_ao", 0);
        } else {
            my $f_ai = "$d_aln/$fam.aln";
            $log->error_die("alignment file [$f_ai] is not there") unless -s $f_ai;
            writeSeq(\@seqs, $f_fas);
            
            my $tag_input = @ids == 1 ? "--p1 $f_fas --p2 $f_ai" : "-i $f_fas --p1 $f_ai";
            runCmd("$f_bin $tag_input --outfmt=fasta --force -o $f_ao", 0);
        }

        my $h = aln_score_group($f_ao, \@ids);
        for my $id (@ids) {
            my $score = 0;
            $score = $h->{$id} if exists $h->{$id};
            print FH join("\t", $id, $score)."\n";
        }
        printf "  %5d / %5d done...\n", $i+1, scalar(keys(%$ref)) if ($i+1) % 1000 == 0;
    }
    runCmd("rm -rf $f_fas $do/*", 0);
    close FH;
}
sub get_stat_hmm {
    my ($fi, $d_hmm, $fo) = @_;
    my $do = "$fo.dir";

    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("computing HMM alignment scores");
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $f_bin = $ENV{"HMMER"}."/bin/hmmsearch";
    -s $f_bin || $log->error_die("cannot execute $f_bin");
    
    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id score/)."\n";
  
    my $f_fas = $ENV{"TMP_DIR"}."/hmmsearch_".int(rand(1000)).".fa";
    my $t = readTable(-in=>$fi, -header=>1);
    my $ref = group($t->colRef("parent"));
    my $i = 1;
    for my $pa (sort keys(%$ref)) {
        my ($idx, $cnt) = @{$ref->{$pa}};
        my $h_seq;
        my $fam = $t->elm($idx, "cat3");
        for my $i ($idx..$idx+$cnt-1) {
            my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
            $h_seq->{$id} = Bio::Seq->new(-id=>$id, -seq=>$seq);
            die "family not consistent for gene[$pa]\n" unless $fam eq $fam2;
        }

        my $f_hmm = "$d_hmm/$fam.hmm";
        $log->error_die("HMM file [$f_hmm] is not there") unless -s $f_hmm;
        writeSeq([values(%$h_seq)], $f_fas);
        
        my $f_ao = "$do/$pa.txt";
        runCmd("$f_bin -o $f_ao $f_hmm $f_fas", 0);
        my $h = score_hmm_by_hit($f_ao);
        for my $id (sort keys(%$h_seq)) {
            my ($e, $score) = (10, 0);
            $score = $h->{$id} if exists $h->{$id};
            print FH join("\t", $id, $score)."\n";
        }
        printf "  %5d / %5d done...\n", $i+1, scalar(keys(%$ref)) if ($i+1) % 1000 == 0;
    }
    runCmd("rm -rf $f_fas $do/*", 0);
    close FH;
}
sub merge_stats {
    my ($fi, $fo, $p) = @_;
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("merging stats");

    my $t = readTable(-in=>$fi, -header=>1);
    my ($eval_hmm, $fm) = @{$p->{"hmm"}};
    my ($eval_aln, $fa) = @{$p->{"aln"}};
    my ($eval_pep, $fp) = @{$p->{"pep"}};
    my ($eval_sp, $fs) = @{$p->{"sp"}};
    
    my @tag_sps = (0) x $t->nofRow;
    my @score_sps = (0) x $t->nofRow;
    if( $eval_sp ) {
        my $ts = readTable(-in=>$fs, -header=>1);
        @tag_sps = $ts->col("tag");
        @score_sps = $ts->col("score");
    }
    $t->addCol( \@tag_sps, "tag_sp" );
    $t->addCol( \@score_sps, "score_sp" );
   
    my @score_hmms = (0) x $t->nofRow;
    if( $eval_hmm ) {
        my $tm = readTable(-in=>$fm, -header=>1);
        @score_hmms = $tm->col("score");
    }
    $t->addCol( \@score_hmms, "score_hmm" );

    my @score_alns = (0) x $t->nofRow;
    if( $eval_aln ) {
        my $ta = readTable(-in=>$fa, -header=>1);
        @score_alns = $ta->col("score");
    }
    $t->addCol( \@score_alns, "score_aln" );

    my @n_cds = (0) x $t->nofRow;
    my @lenC = (0) x $t->nofRow;
    my @lenI = (0) x $t->nofRow;
    my @codonStart = (0) x $t->nofRow;
    my @codonStop = (0) x $t->nofRow;
    my @preStop = (0) x $t->nofRow;
    if( $eval_pep ) {
        my $tp = readTable(-in=>$fp, -header=>1);
        @n_cds = $tp->col("n_cds");
        @lenC = $tp->col("lenC");
        @lenI = $tp->col("lenI");
        @codonStart = $tp->col("codonStart");
        @codonStop = $tp->col("codonStop");
        @preStop = $tp->col("preStop");
    }
    $t->addCol( \@n_cds, "n_cds" );
    $t->addCol( \@lenC, "lenC" );
    $t->addCol( \@lenI, "lenI" );
    $t->addCol( \@codonStart, "codonStart" );
    $t->addCol( \@codonStop, "codonStop" );
    $t->addCol( \@preStop, "preStop" );

    $t->moveCol("seq", $t->nofCol-1);

    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
    print FH $t->tsv(1); 
    close FH;
}
sub pick_best_model {
    my ($f_gtb, $f_stat, $fo, $p) = @_;
    
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("picking best alternative models");

    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $ts = readTable(-in=>$f_stat, -header=>1);
    $log->error_die("not equal rows $f_gtb - $f_stat") unless $tg->nofRow == $ts->nofRow;

    my $h;
    for my $i (0..$ts->nofRow-1) {
        my ($id, $pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop) 
            = $ts->row($i);
        $h->{$pa} ||= [];
        push @{$h->{$pa}} , [$id, $tag_sp, $score_hmm+$score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop, $tg->rowRef($i)];
    }

    open(FH, ">$fo");
    print FH join("\t", $tg->header)."\n";
    my $cnt = 0;
    for my $pa (sort (keys(%$h))) {
        my @rows = sort {$a->[1]<=>$b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @{$h->{$pa}};
        my ($id, $tag_sp, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop, $row) = @{$rows[-1]};
        print FH join("\t", @$row)."\n";
        $cnt ++;
    }
    $log->info(sprintf "\t%5d / %5d picked", $cnt, $tg->nofRow);
    close FH;
}
sub remove_ovlp_models {
    my ($f_stat, $fi, $fo) = @_;

    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("removing overlapping models");

    my $ts = readTable(-in=>$f_stat, -header=>1);
    my $h;
    for my $i (0..$ts->nofRow-1) {
        my ($id, $e, $score_aln) = map {$ts->elm($i, $_)} qw/id e score_aln/;
        $h->{$id} = [$e, $score_aln];
    }

    my ($hl, $hs);
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (0..$ti->nofRow-1) {
        my ($id, $pa, $chr, $beg, $end) = map {$ti->elm($i, $_)} qw/id parent chr beg end/;
        my $loc = [$beg, $end, $i];
        $hl->{$chr} ||= [];
        push @{$hl->{$chr}}, $loc;
        die "no score_aln for $id\n" unless exists $h->{$id};
        $hs->{$i} = $h->{$id};
    }
    
    my @idxs_rm;
    for my $chr (keys %$hl) {
        my $ref = posMerge($hl->{$chr});
        for (@$ref) {
            my ($beg, $end, $idxs) = @$_;
            if(@$idxs > 1) {
                my @idxs = sort {$hs->{$a}->[0] <=> $hs->{$b}->[0] || $hs->{$b}->[1] <=> $hs->{$a}->[1]} @$idxs;
                push @idxs_rm, @idxs[1..$#idxs];
            }
        }
    }
    
    $ti->delRows(\@idxs_rm);
    $log->info(sprintf "\t%5d / %5d passed", $ti->nofRow, $ti->nofRow+@idxs_rm);

    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
    print FH $ti->tsv(1);
    close FH;
}
sub filter_models {
    my ($f_stat, $fi, $fo, $co_e, $co_aln, $opt_sp, $opt_codon, $opt_mt) = 
        rearrange([qw/stat in out e aln sp codon opt_mt/], @_);
    $co_e   ||= 10;
    $co_aln ||= -1000;
    $opt_sp ||= 0;
    $opt_codon ||= 0;
    $opt_mt ||= 0;

    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("final filter:");
    $log->info("\tSignal Peptide Present = $opt_sp");
    $log->info("\tComplete ORF = $opt_codon");
    $log->info("\tE value cutoff = $co_e");
    $log->info("\tAlignment score cutoff = $co_aln");
    $log->info("\t[optionl] Medicago evaluation filter = $opt_mt");
    
    my $ts = readTable(-in=>$f_stat, -header=>1);
    my $h;
    for my $i (0..$ts->nofRow-1) {
        my ($id, $pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop) 
            = $ts->row($i);
        $h->{$id} = [$tag_sp, $score_hmm, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop];
    }

    my $t = readTable(-in=>$fi, -header=>1);
    my @idxs_rm;
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam) = map {$t->elm($i, $_)} qw/id cat3/;
        die "no stat for $id\n" unless exists $h->{$id};
        my ($tag_sp, $score_hmm, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop) = @{$h->{$id}};
        if($opt_sp && $tag_sp == 0) {
            push @idxs_rm, $i;
        } elsif($opt_codon && (!$codonStart || !$codonStop || $preStop)) {
            push @idxs_rm, $i;
        } elsif($e > $co_e || $score_aln < $co_aln) {
            push @idxs_rm, $i;
        } elsif($opt_mt && $fam gt "CRP1530") {
            push @idxs_rm, $i;
        }
    }
    
    $t->delRows(\@idxs_rm);
    $log->info(sprintf "\t%5d / %5d passed", $t->nofRow, $t->nofRow+@idxs_rm);

    open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
    print FH $t->tsv(1);
    close FH;
}

sub crp_rename {
    my ($fi, $f_ref, $fo, $pre) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    $ti->sort("cat3", 1, 0, "chr", 1, 0, "beg", 0, 0);
    
    my @chrs = uniq($ti->col("chr"));
    my @len_digits = map {getDigits(seqLen($_, $f_ref) / 1000000)} @chrs;
    my $chr_digits = getDigits(scalar(grep /\d+/, @chrs));
    my $h = { map {$chrs[$_] => $len_digits[$_]} 0..$#chrs };

    for my $i (0..$ti->nofRow-1) {
        my ($chr, $beg, $fam) = map {$ti->elm($i, $_)} qw/chr beg cat3/;
        my $begStr = sprintf "%0".$h->{$chr}."d", $beg/1000000;
        my $chrStr = $chr;
        $chrStr =~ s/chr//i;
        $chrStr = sprintf "%0".$chr_digits."d", $chrStr if $chrStr =~ /^\d+$/;

        my $id = sprintf "$pre\_%s_chr%s_%sM", lc($fam), $chrStr, $begStr;
        $ti->setElm($i, "parent", $ti->elm($i, "id"));
        $ti->setElm($i, "id", $id);
    }
    
    my $ref = group($ti->colRef("id"));
    my $hd = { map { $_ => getDigits($ref->{$_}->[1]) } keys(%$ref) };
    my $hc;
    for my $i (0..$ti->nofRow-1) {
        my $id = $ti->elm($i, "id");
        $hc->{$id} ||= 0;
        my $cnt = ++$hc->{$id};
        
        $id = sprintf "%s_%0".$hd->{$id}."d", $id, $cnt;
        $ti->setElm($i, "id", $id);
    }
    open(FHO, ">$fo") or die "cannot open $fo for writing\n";
    print FHO $ti->tsv(1);
    close FHO;
}

sub gtb2Friendly {
    my ($fi, $fo, $f_stat) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    
    my $hs;
    my $ts = readTable(-in=>$f_stat, -header=>1);
    for my $i (0..$ts->nofRow-1) {
        my ($id, @stats) = map {$ts->elm($i, $_)} qw/id parent family e tag_sp score_sp score_hmm score_aln n_cds seq/;
        die "$id read twice in $f_stat\n" if exists $hs->{$id};
        $hs->{$id} = \@stats;
    }

    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\t", qw/id family chr beg end strand e score_sp score_hmm score_aln sequence/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $srd) = map {$t->elm($i, $_)} qw/id chr beg end strand/;
        my ($pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $seq) 
            = @{$hs->{$t->elm($i, "parent")}};
        print FH join("\t", $id, $fam, $chr, $beg, $end, $srd, $e, $score_sp, $score_hmm, $score_aln, $seq)."\n";
    }
    close FH;
}

sub align_by_group {
    my ($f_gtb, $f_hs, $dirO) = rearrange(['f_gtb', 'hmmstat', 'out'], @_);
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("making sub-family alignments");
    make_path($dirO) unless -d $dirO;
    remove_tree($dirO, {keep_root => 1});
    my $t = readTable(-in=>$f_gtb, -header=>1);

    my $ths = readTable(-in=>$f_hs, -header=>1);
    my $hs = { map {$ths->elm($_, "id") => $ths->elm($_, "consensus")} (0..$ths->nofRow-1) };

    my $h;
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
        $h->{$fam} ||= [];
        my $seqObj = Bio::Seq->new(-id=>$id, -seq=>$seq);
        push @{$h->{$fam}}, $seqObj;
    }

    my $i = 1;
    for my $fam (sort keys %$h) {
        my $seqs = $h->{$fam};
        die "no consensus for $fam\n" unless exists $hs->{"$fam"};
        push @$seqs, Bio::Seq->new(-id=>$fam, -seq=>$hs->{"$fam"});
        my $f_aln = "$dirO/$fam.aln";
        run_clustalo(-seqs=>$seqs, -out=>$f_aln);
        printf "  %5d / %5d done...\n", $i+1, scalar(keys(%$h)) if ($i+1) % 1000 == 0;
    }
}

sub pipe_model_evaluation {
    my ($dir, $f_hit, $f_gtb, $f_ref, $f_gtb_ref, $d_hmm, $d_aln, $f_sta) = 
        rearrange([qw/dir hit gtb_all ref gtb_ref d_hmm d_aln f_sta/], @_); 
    
    my $log = Log::Log4perl->get_logger("ModelEval");
    $log->info("#####  Stage 4 [Model Evaluation & Selection]  #####");
    make_path($dir) unless -d $dir;
   
    my ($eval_e, $eval_hmm, $eval_aln, $eval_pep) = (1) x 4;
    my $eval_sp = $ENV{"eval_sp"};

    my $f30 = "$dir/30_stat_basic.tbl";
    get_stat_basic($f_gtb, $f_hit, $f30);

    my $d31 = "$dir/31_stat";
    my $f31_02 = "$d31/02_hmm.tbl";
    my $f31_03 = "$d31/03_aln.tbl";
    my $f31_04 = "$d31/04_pep.tbl";
    my $f31_05 = "$d31/05_sigp.tbl";
    get_stat_hmm($f_gtb, $d_hmm, $f31_02) if $eval_hmm;
    get_stat_aln($f_gtb, $d_aln, $f_sta, $f31_03) if $eval_aln;
    pep_score_gtb($f_gtb, $f31_04) if $eval_pep;
    sigp_score_gtb($f_gtb, $f31_05) if $eval_sp;

    my $p = {
        "hmm" => [ $eval_hmm, $f31_02 ],
        "aln" => [ $eval_aln, $f31_03 ],
        "pep" => [ $eval_pep, $f31_04 ],
        "sp"  => [ $eval_sp,  $f31_05 ] };
    my $f41 = "$dir/41_stat.tbl";
    merge_stats($f30, $f41, $p);
    my $f51 = "$dir/51_best.gtb";
    pick_best_model($f_gtb, $f41, $f51, $p);
    my $f55 = "$dir/55_nonovlp.gtb";
    remove_ovlp_models($f41, $f51, $f55);
    my $f59 = "$dir/59.gtb";
    filter_models(-stat=>$f41, -in=>$f55, -out=>$f59, 
        -e=>$ENV{"evalue"}, -aln=>-1000, -sp=>$eval_sp, -codon=>1, -opt_mt=>0);
    
    my $f61 = "$dir/61_final.gtb";
    crp_rename($f59, $f_ref, $f61, "spada");
    gtb2Gff($f61, "$dir/61_final.gff");
    gtb2Friendly($f61, "$dir/61_final.tbl", $f41);
    
    my $f81 = "$dir/81_aln";
    align_by_group(-f_gtb=>$f61, -hmmstat=>$f_sta, -out=>$f81);
    my $f91 = "$dir/91_compare.tbl";
    if(defined($f_gtb_ref) && -s $f_gtb_ref) {
        compare_models($f61, $f_gtb_ref, $f91);
    }
}


1;
__END__
