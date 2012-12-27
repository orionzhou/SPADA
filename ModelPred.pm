package ModelPred;
use strict;
use File::Path qw/make_path remove_tree/;
use Log::Log4perl;
use Data::Dumper;
use Common;
use Seq;
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
@EXPORT = qw/pipe_model_prediction/;
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
    $log->info("preparing hit sequences");
    make_path($dir) unless -d $dir;
    my $f01 = "$dir/01_hit_seq.tbl";
    get_hit_seq($f_hit, $f01, $f_ref);
    my $f05 = "$dir/05_hits.tbl";
    prefilter_hits($f01, $f05, 10);
#  prefilter_hits($f01, $f05, $ENV{'evalue'});
}

sub pipe_model_run {
    my ($dir, $f_hit, $f_ref, $soft) = rearrange(['dir', 'hit', 'ref', 'soft'], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("### working on $soft pipeline ###");
    
    if($soft eq "Augustus_evidence") {
        pipe_augustus($f_hit, $dir);
    } elsif($soft eq "SPADA") {
        pipe_spada($f_hit, $dir);
    } elsif($soft eq "Augustus_de_novo") {
        pipe_augustus_simple($f_hit, $dir);
    } elsif($soft eq "GeneMark") {
        pipe_genemark($f_hit, $dir);
    } elsif($soft eq "GlimmerHMM") {
        pipe_glimmerhmm($f_hit, $dir);
    } elsif($soft eq "GeneID") {
        pipe_geneid($f_hit, $dir);
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
    my ($dir, $f_hit, $f_ref, $p) = rearrange([qw/dir hit ref p/], @_);
    my $log = Log::Log4perl->get_logger("Model");
    $log->info("checking & refining models");
  
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

sub pipe_model_prediction {
    my ($dir, $f_hit, $f_ref) = rearrange([qw/dir hit ref/], @_); 
    my $log = Log::Log4perl->get_logger("ModelPrediction");
    $log->info("#####  Stage 3 [Model Prediction]  #####");

    my $f05 = "$dir/05_hits.tbl";
    pipe_model_prepare(-dir=>$dir, -hit=>$f_hit, -ref=>$f_ref);
    
    my $p;
    for my $method (keys %{$ENV{"method"}}) {
        next if $ENV{"method"}->{$method} != 1;
        my $dirs = "$dir/$method";
        pipe_model_run(-dir=>$dirs, -hit=>$f05, -ref=>$f_ref, -soft=>$method);
        $p->{$method} = "$dirs/11.gtb";
    }
    my $f26 = "$dir/26_all.gtb";
    pipe_model_check(-dir=>$dir, -hit=>$f05, -ref=>$f_ref, -p=>$p);
}

1;
__END__
