package Model;
use strict;
use File::Basename;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Log::Log4perl;
use Data::Dumper;
use Common;
use Seq;
use Hmm;
use Align;
use Gtb;
use Genemodel;
use Software;
use SignalP;
use Spada;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_model pipe_model_postprocess
    filter_models/;
@EXPORT_OK = qw//;

sub remove_gap {
    my ($begr, $endr, $strand, $lenU, $lenH, $lenD, $seq) = @_;
#  print join("\t", $begr, $endr, $strand, $lenU, $lenH, $lenD, $seq)."\n";

    my ($beg2, $end2) = (1, $lenU+$lenH+$lenD);
    my ($lenU2, $lenD2) = ($lenU, $lenD);
    while($seq =~ /(^N+)|(N{90,})|(N+$)/ig) {
        my ($posl, $posr) = ($-[0] + 1, $+[0]);
        if($posl <= $lenU) {
            $posr = $posr<=$lenU ? $posr : $lenU;
            my $beg2_tmp = $posr + 1;
            my $lenU2_tmp = $lenU - $beg2_tmp + 1;
            if($lenU2_tmp < $lenU2) {
                $lenU2 = $lenU2_tmp;
                $beg2 = $beg2_tmp;
            }
        }
        if($posr >= $lenU+$lenH+1) {
            $posl = $posl>=$lenU+$lenH+1 ? $posl : $lenU+$lenH+1;
            my $end2_tmp = $posl - 1;
            my $lenD2_tmp = $end2_tmp - ($lenU+$lenH+1) + 1;
            if($lenD2_tmp < $lenD2) {
                $lenD2 = $lenD2_tmp;
                $end2 = $end2_tmp;
            }
        }
#    print join("\t", $posl, $posr, $beg2, $end2, $lenU2, $lenD2)."\n";
    }

    if($strand =~ /^\-1?$/) {
        $endr -= $beg2 - 1;
        $begr += ($lenU+$lenH+$lenD) - $end2;
    } else {
        $begr += $beg2 - 1;
        $endr -= ($lenU+$lenH+$lenD) - $end2;
    }
    my $seq2 = substr($seq, $beg2-1, $end2-$beg2+1);
    return ($begr, $endr, $lenU2, $lenD2, $seq2);
}
sub get_hit_seq {
    my ($fi, $fo, $f_ref) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("extracting hit sequence");
    my $ti = readTable(-in=>$fi, -header=>1);
    my $to = Data::Table->new([], [qw/id family chr beg end srd loc begr endr begL endL locL e seq_pro seq/], 0);

    $ti->sort("chr", 1, 0, "beg", 0, 0);
    for my $i (0..$ti->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locHS, $begQ, $endQ, $srdQ, $locQS, $e) = $ti->row($i);
#    next unless $id == 24;
        my $locG = locStr2Ary($locHS);
        my $seq_dna = seqRet($locG, $chr, $srd, $f_ref);
        my $seq_pro = Bio::Seq->new(-seq=>$seq_dna)->translate()->seq;

        my $len_hit = $end - $beg + 1;
        $log->warn("hit[$id] too long: $len_hit bp") if $len_hit > 10000;
        my $len_chr = seqLen($chr, $f_ref);
        
        my ($len_le, $len_ri, $len_up, $len_dw);
        my ($len_up_max, $len_dw_max) = (2500, 1500);
        my ($intv_hit_up, $intv_hit_dw) = (2500, 1500);
        if($srd =~ /^[\+1]$/) {
            $intv_hit_up = $beg - $ti->elm($i-1, "end") - 1 if ($i > 0 && $chr eq $ti->elm($i-1, "chr"));
            $len_le = min($beg-1, $intv_hit_up, $len_up_max);
            $len_ri = min($len_chr-$end, $intv_hit_dw, $len_dw_max);
            ($len_up, $len_dw) = ($len_le, $len_ri);
        } elsif($srd =~ /^\-1?$/) {
            $intv_hit_up = $ti->elm($i+1, "beg") - $end - 1 if ($i < $ti->nofRow-1 && $chr eq $ti->elm($i+1, "chr"));
            $len_le = min($beg-1, $intv_hit_dw, $len_dw_max);
            $len_ri = min($len_chr-$end, $intv_hit_up, $len_up_max);
            ($len_up, $len_dw) = ($len_ri, $len_le);
        } else {
            die "unknown strand: $srd\n";
        }
        my $begr = $beg - $len_le;
        my $endr = $end + $len_ri;
        die join("\t", $id, $chr, $begr, $endr, $intv_hit_up, $intv_hit_dw)."\n" if $begr >= $endr;
        
        my $seq = seqRet([[$begr, $endr]], $chr, $srd, $f_ref);

        ($begr, $endr, $len_up, $len_dw, $seq) = remove_gap($begr, $endr, $srd, $len_up, $len_hit, $len_dw, $seq);
#    my $locr2 = Bio::Location::Simple->new(-seq_id=>$chr, -start=>$begr, -end=>$endr, -strand=>$strand);
#    my $seq2 = seqRet($locr2, $f_ref);
#    die "$seq\n$seq2\n" unless $seq eq $seq2;

        my $locL = $srd =~ /^\-1?$/ ? [ map {[$endr-$_->[1]+1, $endr-$_->[0]+1]} @$locG ]
            : [ map {[$_->[0]-$begr+1, $_->[1]-$begr+1]} @$locG ];
        $locL = [ sort {$a->[0] <=> $b->[0]} @$locL ];
        my ($begL, $endL) = ($locL->[0]->[0], $locL->[-1]->[1]);
        my $locLS = locAry2Str($locL);
        $to->addRow([$id, $fam, $chr, $beg, $end, $srd, $locHS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq]);
    }

    open(FH, ">$fo") or die "cannot open $fo\n";
    $to->sort("id", 1, 0);
    print FH $to->tsv(1);
    close FH;
}
sub prefilter_hits {
    my ($fi, $fo, $co_e) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("prefiltering hits by E-value");

    my @idxs_rm;
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (0..$ti->nofRow-1) {
        my $e = $ti->elm($i, "e");
        push @idxs_rm, $i if $e > $co_e;
    }
    
    $log->info(sprintf "  %d out of %d passed", $ti->nofRow-@idxs_rm, $ti->nofRow);
    $ti->delRows(\@idxs_rm);
    open(FH, ">$fo") or die "cannot write to $fo\n";
    print FH $ti->tsv(1);
    close FH;
}
sub pipe_model_prepare {
    my ($dir, $f_hit, $f_ref) = rearrange(['dir', 'hit', 'ref'], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("preparing to run");
    make_path($dir) unless -d $dir;
    my $f01 = "$dir/01_hit_seq.tbl";
    get_hit_seq($f_hit, $f01, $f_ref);
    my $f05 = "$dir/05_hits.tbl";
    prefilter_hits($f01, $f05, 10);
#  prefilter_hits($f01, $f05, $ENV{'evalue'});
}

sub pipe_model_run {
    my ($dir, $f_ref, $soft) = rearrange(['dir', 'ref', 'soft'], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("____working on $soft pipeline___");
    my $f_hit = "$dir/05_hits.tbl";
    
    my $d15 = "$dir/15";
    my $p = {$soft => "$d15/11.gtb"};
    if($soft eq "SPADA") {
        my $d12 = "$dir/12_augustus";
        pipe_augustus($f_hit, $d12);
        my $d13 = "$dir/13_spada";
        pipe_spada($f_hit, $d13);
        $p = {"augustus"=>"$d12/11.gtb", "spada"=>"$d13/11.gtb"};
    } elsif($soft eq "Augustus") {
        pipe_augustus_simple($f_hit, $d15);
    } elsif($soft eq "GeneMark") {
        pipe_genemark($f_hit, $d15);
    } elsif($soft eq "GlimmerHMM") {
        pipe_glimmerhmm($f_hit, $d15);
    } elsif($soft eq "GeneID") {
        pipe_geneid($f_hit, $d15);
    } else {
        die "unsupported program: $soft\n";
    }
}

sub collect_models { 
    my ($f_hit, $fo, $f_ref, $p) = rearrange(['hit', 'out', 'ref', 'p'], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("collecting all prediction results");
    my $hm;
    for my $source (keys(%$p)) {
        my $f_gtb = $p->{$source};
        next unless -s $f_gtb;
        my $tg = readTable(-in=>$f_gtb, -header=>1);
        for my $i (0..$tg->nofRow-1) {
            my ($id) = $tg->elm($i, "parent");
            $tg->setElm($i, "source", $source);
            $hm->{$id} ||= [];
            my $row = $tg->rowRef($i);
            push @{$hm->{$id}}, [@$row[0..18]];
        }
    }
    
    my $th = readTable(-in=>$f_hit, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note seq/)."\n";
    for my $i (0..$th->nofRow-1) {
        my ($id, $fam) = $th->row($i);
        next if( !exists $hm->{$id} );
        for my $row (@{$hm->{$id}}) {
            $row->[14] = 'gene';
            $row->[15] = 'mRNA';
            $row->[16] = $fam;
            print FH join("\t", @$row)."\n";
        }
    }
    close FH;
}
sub refine_incomplete_models {
    my ($f_hit, $fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("refining incomplete models");
    my $th = readTable(-in=>$f_hit, -header=>1);
    my $h = { map {$th->elm($_, "id") => $th->rowRef($_)} (0..$th->nofRow-1) };

    my $tg = readTable(-in=>$fi, -header=>1);
    my @idxs_rm;
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $locStr, $phase, $seq_gene) = map {$tg->elm($i, $_)} qw/id parent locC phase seq/;
        my $loc = locStr2Ary($locStr);
        $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
        my @phases = split(",", $phase);

        $log->error_die("no hit info for $pa") unless exists $h->{$pa};
        my ($id_hit, $fam, $chr, $begG, $endG, $strand, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq)
            = @{$h->{$pa}};
        
        my ($codonStart, $codonStop, $preStop, $gap) = checkProtSeq($seq_gene);
        push @idxs_rm, $i if $preStop;
        next if $codonStart && $codonStop;
        
        unless($codonStart) {
            my ($beg, $end, $phase) = (@{$loc->[0]}, $phases[0]);

            my $codon_up = min(20, int(($beg-1)/3));
            my $begR = $beg + $phase - 3*$codon_up;
            
            my $seq_dna = substr($seq, $begR-1, $end-$begR+1);
            my $seq_pro = Bio::Seq->new(-seq=>$seq_dna)->translate->seq;
            if($seq_pro =~ /(M[^\*]+)$/) {
                my $codon_offset = $-[0];
                my $begN = $begR + $codon_offset * 3;
                $loc->[0]->[0] = $begN;
                $phases[0] = 0;
            }
        }
        unless($codonStop) {
            my ($beg, $end, $phase) = (@{$loc->[-1]}, $phases[-1]);
            
            my $codon_dw = min(200, int((length($seq)-$end)/3));
            my $endR = $end + 3*$codon_dw;

            my $seq_dna = substr($seq, $beg+$phase-1, $endR-$beg-$phase+1);
            die "$id\n" unless $seq_dna;
            my $seq_pro = Bio::Seq->new(-seq=>$seq_dna)->translate->seq;
            if($seq_pro =~ /^([^\*]+\*)/) {
                my $codon_offset = $+[0];
                my $endN = $beg+$phase + $codon_offset * 3 - 1;
                $loc->[-1]->[-1] = $endN;
                $phases[-1] = $phase;
            }
        }
        $tg->setElm($i, "locC", locAry2Str($loc));
        $tg->setElm($i, "phase", join(",", @phases));
        my $seq_new = Bio::Seq->new(-seq=>getSubSeq($seq, $loc))->translate->seq;
        $tg->setElm($i, "seq", $seq_new);
    }

    open(FH, ">$fo");
    $tg->delRows(\@idxs_rm);
    print FH $tg->tsv(1);
    close FH;
}
sub merge_redundant_models {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("merging redundant models");
    my $ti = readTable(-in=>$fi, -header=>1);
    
    my $h;
    for my $i (0..$ti->nofRow-1) {
        my ($id) = $ti->elm($i, "parent");
        $h->{$id} ||= [];
        my $row = $ti->rowRef($i);
        push @{$h->{$id}}, [@$row[0..18]];
    }
  
    my $cnt = 0;
    open(FH, ">$fo");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note seq/)."\n";
    for my $id (sort keys %$h) {
        my @rows = @{$h->{$id}};
        my @locCs;
        for (@rows) {
            my $locAry = locStr2Ary($_->[8]);
            $locAry = [ sort {$a->[0] <=> $b->[0]} @$locAry ];
            push @locCs, join(",", map {$_->[0]."-".$_->[1]} @$locAry);
        }
        my @locCs_uniq = (uniq(@locCs));
        my $n = min($#locCs_uniq+1, 99);
        $log->warn("$id has >99 models: first 99 kept") if @locCs_uniq > 99;
        for my $i (0..$n-1) {
            my $locC = $locCs_uniq[$i];
            my @idxs = indexes {$_ eq $locC} @locCs;
            my @sources = map {$rows[$_]->[12]} @idxs;
            my $source = join(" ", uniq(@sources));
            my $row = @rows[$idxs[0]];
            $row->[0] = sprintf "%s.%02d", $id, ($i+1);
            $row->[12] = $source;
            $row->[17] = $source;
            print FH join("\t", @$row)."\n";
        }
        $cnt += $n;
    }
    close FH;
    $log->info("  $cnt in total");
}
sub locRel2Abs {
    my ($beg, $end, $strand, $loc_rel_ary, $strand_rel) = @_;

    my $strand_abs = is_opposite_strands($strand, $strand_rel) ? "-" : "+";
    my $loc_abs_ary;
    if($strand_abs eq "-") {
        $loc_abs_ary = [ map {[$end - $_->[1] + 1, $end - $_->[0] + 1]} @$loc_rel_ary ];
    } else {
        $loc_abs_ary = [ map {[$beg + $_->[0] - 1, $beg + $_->[1] - 1]} @$loc_rel_ary ];
    }
    return ($loc_abs_ary, $strand_abs);
}
sub recover_global_coordinate {
    my ($f_hit, $fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("converting coordinates to global genomic positions");
    my $th = readTable(-in=>$f_hit, -header=>1);
    my $h = { map {$th->elm($_, "id") => $th->rowRef($_)} (0..$th->nofRow-1) };

    my $tg = readTable(-in=>$fi, -header=>1);
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $loc_cds_rel_str, $seq_pro) = map {$tg->elm($i, $_)} qw/id parent locC seq/;
        my $loc_cds_rel_ary = locStr2Ary($loc_cds_rel_str);

        $log->error_die("no hit info for $pa") unless exists $h->{$pa};
        my ($id_hit, $fam, $chr, $beg, $end, $strand, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq)
            = @{$h->{$pa}};

        my ($locCAry, $strand_abs) = locRel2Abs($begr, $endr, $strand, $loc_cds_rel_ary, "+");
        $log->error_die("strand inconsistent: $strand_abs, $strand") unless $strand_abs eq $strand;
        $locCAry = [ sort {$a->[0] <=> $b->[0]} @$locCAry ];
        my ($begM, $endM) = ($locCAry->[0]->[0], $locCAry->[-1]->[-1]);
        my ($locIAry) = posDiff([[$begM, $endM]], $locCAry);

        my $locCStr = locAry2Str($locCAry);
        my $locIStr = locAry2Str($locIAry);
        my $phase = join(",", getPhase($locCAry, $strand));

        $tg->setElm($i, "chr", $chr);
        $tg->setElm($i, "beg", $begM);
        $tg->setElm($i, "end", $endM);
        $tg->setElm($i, "strand", $strand_abs);
        $tg->setElm($i, "locE", $locCStr);
        $tg->setElm($i, "locI", $locIStr);
        $tg->setElm($i, "locC", $locCStr);
        $tg->setElm($i, "phase", $phase);
    }
    open(FH, ">$fo");
    print FH $tg->tsv(1);
    close FH;
}
sub remove_incompatible_models {
    my ($f_hit, $fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("removing incompatible models");
    my $th = readTable(-in=>$f_hit, -header=>1);
    my $h;
    for my $i (0..$th->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $th->row($i);
        my $loc = locStr2Ary($locS);
        my $phase = join(",", (0) x @$loc);
        $h->{$id} = [$loc, $srd, $phase];
    }

    my @idxs_rm;
    my $ti = readTable(-in=>$fi, -header=>1);
    for my $i (0..$ti->nofRow-1) {
        my ($id, $pa, $chr, $srdG, $locGS, $locIS, $phaseG) = map {$ti->elm($i, $_)} qw/id parent chr strand locC locI phase/;
        my $locG = locStr2Ary($locGS);
        my $locI = locStr2Ary($locIS);
      
        my $tag_rm = 0;
        my ($locH, $srdH, $phaseH) = @{$h->{$pa}};
        if($srdG eq $srdH) {
            my ($tag, $lenO, $len1, $len2) = compare_2_model($locH, $phaseH, $locG, $phaseG, $srdG);
            $tag_rm = 1 if $tag > 2 || $len1 > 30;
        }
        for (@$locI) {
            my $lenI = $_->[1] - $_->[0] + 1;
            $tag_rm = 1 if $lenI > 1500;
        }
        push @idxs_rm, $i if $tag_rm;
    }
    $log->info(sprintf "  %d removed", scalar(@idxs_rm));
    
    $ti->delRows(\@idxs_rm);
    open(FH, ">$fo");
    print FH $ti->tsv(1);
    close FH;
}
sub pipe_model_check {
    my ($dir, $f_ref, $soft) = rearrange(['dir', 'ref', 'soft'], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("checking & refining models");
    my $f_hit = "$dir/05_hits.tbl";
  
    my $p = {$soft=>"$dir/15/11.gtb"};
    $p = {"augustus"=>"$dir/12_augustus/11.gtb", "spada"=>"$dir/13_spada/11.gtb"} if $soft eq "SPADA";
    
    my $f21 = "$dir/21_rel.gtb";
    collect_models(-hit=>$f_hit, -out=>$f21, -ref=>$f_ref, -p=>$p);
    my $f22 = "$dir/22_refined.gtb";
    refine_incomplete_models($f_hit, $f21, $f22);
    my $f23 = "$dir/23_merged.gtb";
    merge_redundant_models($f22, $f23);
    my $f24 = "$dir/24_abs.gtb";
    recover_global_coordinate($f_hit, $f23, $f24);
    my $f26 = "$dir/26_all.gtb";
    remove_incompatible_models($f_hit, $f24, $f26);
    gtb2Gff($f26, "$dir/26_all.gff");
}

sub get_aln_score {
    my ($fi, $d_aln, $f_sta, $do, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("computing MSA scores");
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $f_bin = $ENV{"ClustalO"}."/bin/clustalo";
    $log->error_die("cannot execute $f_bin") unless -s $f_bin;
    
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
        printf "  %5d / %5d done\r", $i++, scalar(keys(%$ref));
    }
    print "\n";
    system("rm -rf $f_fas");
    close FH;
}
sub get_hmm_score {
    my ($fi, $d_hmm, $do, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("computing HMM alignment scores");
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $f_bin = $ENV{"HMMER"}."/bin/hmmsearch";
    $log->error_die("cannot execute $f_bin") unless -s $f_bin;
    
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
        printf "  %5d / %5d done\r", $i++, scalar(keys(%$ref));
    }
    print "\n";
    system("rm -rf $f_fas");
    close FH;
}
sub merge_stats {
    my ($f_gtb, $fe, $fm, $fa, $fs, $fp, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("merging stats");
    my $tg = readTable(-in=>$f_gtb, -header=>1);
    my $te = readTable(-in=>$fe, -header=>1);
    my $tm = readTable(-in=>$fm, -header=>1);
    my $ta = readTable(-in=>$fa, -header=>1);
    my $ts = readTable(-in=>$fs, -header=>1);
    my $tp = readTable(-in=>$fp, -header=>1);

    my $he = { map {$te->elm($_, "id") => $te->elm($_, "e")} (0..$te->nofRow-1) };

    open(FH, ">$fo");
    print FH join("\t", qw/id parent fam e tag_sp score_sp score_hmm score_aln n_cds lenC lenI codonStart codonStop preStop seq/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $pa, $fam, $seq) = map {$tg->elm($i, $_)} qw/id parent cat3 seq/;
        my ($idM, $score_hmm) = map {$tm->elm($i, $_)} qw/id score/;
        my ($idA, $score_aln) = map {$ta->elm($i, $_)} qw/id score/;
        my ($idS, $tag_sp, $score_sp) = map {$ts->elm($i, $_)} qw/id tag score/;
        my ($idP, $codonStart, $codonStop, $preStop, $gap, $n_cds, $lenC, $lenI) = $tp->row($i);
        die "id conflict: $id != $idM [hmm]\n" unless $id eq $idM;
        die "id conflict: $id != $idA [aln]\n" unless $id eq $idA;
        die "id conflict: $id != $idS [sigp]\n" unless $id eq $idS;
        die "id conflict: $id != $idP [pep]\n" unless $id eq $idP;
        $score_sp = 0 unless $tag_sp;
        my $e = $he->{$pa};
        
        print FH join("\t", $id, $pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop, $seq)."\n";
    }
    close FH;
}
sub pick_best_model {
    my ($f_gtb, $f_stat, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("picking best models & removing incomplete models");
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
    my $i = 0;
    for my $pa (sort (keys(%$h))) {
        my @rows = sort {$a->[1]<=>$b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @{$h->{$pa}};
        my ($id, $tag_sp, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop, $row) = @{$rows[-1]};
        next if $tag_sp == 0 || !$codonStart || !$codonStop || $preStop;
        print FH join("\t", @$row)."\n";
        $i ++;
    }
    $log->info("  $i picked");
    close FH;
}
sub remove_ovlp_models {
    my ($f_stat, $fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Model");
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
    $log->info(sprintf "  %d passed", $ti->nofRow);
    open(FH, ">$fo");
    print FH $ti->tsv(1);
    close FH;
}
sub filter_models {
    my ($f_stat, $fi, $fo, $co_e, $co_aln, $opt_mt) = @_;
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("final e-value filter");
    
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
        if($e > $co_e || $score_aln < $co_aln) {
            push @idxs_rm, $i;
        } elsif($opt_mt && $ENV{"SPADA_ORG"} eq "Mtruncatula" && $fam gt "CRP1530") {
            push @idxs_rm, $i;
        }
    }
    
    $t->delRows(\@idxs_rm);
    $log->info(sprintf "  %d passed", $t->nofRow);
    open(FH, ">$fo");
    print FH $t->tsv(1);
    close FH;
}
sub output_models {
    my ($fi, $fo, $f_stat) = @_;
    my $t = readTable(-in=>$fi, -header=>1);
    
    my $hs;
    my $ts = readTable(-in=>$f_stat, -header=>1);
    for my $i (0..$ts->nofRow-1) {
        my ($id, @stats) = map {$ts->elm($i, $_)} qw/id parent fam e tag_sp score_sp score_hmm score_aln n_cds seq/;
        die "$id read twice in $f_stat\n" if exists $hs->{$id};
        $hs->{$id} = \@stats;
    }

    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\t", qw/id family chr start end strand e score_sp score_hmm score_aln sequence/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $chr, $beg, $end, $srd) = map {$t->elm($i, $_)} qw/id chr beg end strand/;
        my ($pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $seq) = @{$hs->{$id}};
        print FH join("\t", $id, $fam, $chr, $beg, $end, $srd, $e, $score_sp, $score_hmm, $score_aln, $seq)."\n";
    }
    close FH;
}
sub pipe_model_pick {
    my ($dir, $f_ref, $f_gtb, $d_hmm, $d_aln, $f_sta) = rearrange([qw/dir ref gtb_all d_hmm d_aln f_sta/], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("evaluating & picking models");

    my $f_hit = "$dir/05_hits.tbl";

    my $d31 = "$dir/31_stat";
    my $d31_01 = "$d31/01_hmm_aln";
    my $f31_02 = "$d31/02_hmm_score.tbl";
    get_hmm_score($f_gtb, $d_hmm, $d31_01, $f31_02);
    my $d31_11 = "$d31/11_aln";
    my $f31_12 = "$d31/12_aln_score.tbl";
    get_aln_score($f_gtb, $d_aln, $f_sta, $d31_11, $f31_12);
    my $f31_21 = "$d31/21_sigp_score.tbl";
    sigp_score_gtb($f_gtb, $f31_21);
    my $f31_31 = "$d31/31_pep_score.tbl";
    pep_score_gtb($f_gtb, $f31_31);
    my $f41 = "$dir/41_stat.tbl";
    merge_stats($f_gtb, $f_hit, $f31_02, $f31_12, $f31_21, $f31_31, $f41);
    my $f51 = "$dir/51_best.gtb";
    pick_best_model($f_gtb, $f41, $f51);
    my $f55 = "$dir/55_nonovlp.gtb";
    remove_ovlp_models($f41, $f51, $f55);
    my $f61 = "$dir/61_final.gtb";
    filter_models($f41, $f55, $f61, $ENV{"evalue"}, -1000);
    gtb2Gff($f61, "$dir/61_final.gff");
    output_models($f61, "$dir/61_final.tbl", $f41);
}

sub pipe_model {
    my ($dir, $f_hit, $f_ref, $d_hmm, $d_aln, $f_sta, $soft) = rearrange([qw/dir hit ref d_hmm d_aln f_sta soft/], @_); 
    pipe_model_prepare(-dir=>$dir, -hit=>$f_hit, -ref=>$f_ref);
    pipe_model_run(-dir=>$dir, -ref=>$f_ref, -soft=>$soft);
    my $f29 = "$dir/26_all.gtb";
    pipe_model_check(-dir=>$dir, -ref=>$f_ref, -soft=>$soft);
    my $f61 = "$dir/61_final.gtb";
    pipe_model_pick(-dir=>$dir, -ref=>$f_ref, -gtb_all=>$f29, -d_hmm=>$d_hmm, -d_aln=>$d_aln, -f_sta=>$f_sta);
}

sub align_by_group {
    my ($f_gtb, $f_hs, $dirO) = rearrange(['f_gtb', 'hmmstat', 'out'], @_);
    my $log = Log::Log4perl->get_logger("Model");
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
        printf "  %5d / %5d done...\r", $i++, scalar(keys(%$h));
    }
    print "\n";
}
sub pipe_model_postprocess {
    my ($dir, $f_sta, $f_gtb) = rearrange([qw/dir f_sta f_gtb/], @_);
    my $f_gtb_picked = "$dir/61_final.gtb";
    my $f81 = "$dir/81_aln";
    align_by_group(-f_gtb=>$f_gtb_picked, -hmmstat=>$f_sta, -out=>$f81);
    my $f91 = "$dir/91_compare.tbl";
    if(-s $f_gtb) {
        compare_models($f_gtb_picked, $f_gtb, $f91);
    }
}



1;
__END__
