package ModelEval;
use strict;
use Cwd qw/getcwd abs_path/;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Log::Log4perl;
use Data::Dumper;
use Common;
use Location;
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
  open(my $fh, ">$fo") || die "cannot open $fo for writing\n";
  print $fh join("\t", qw/id parent family e seq/)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($id, $par, $fam, $seq) = map {$t->elm($i, $_)} qw/id par cat3 seq/;
    die "no E-value for $id:$par\n" unless exists $he->{$par};
    my $e = $he->{$par};
    print $fh join("\t", $id, $par, $fam, $e, $seq)."\n";
  }
  close $fh;
}
sub get_stat_hmm {
  my ($fi, $d_hmm, $fo) = @_;
  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("computing HMM alignment scores");
 
  my $cwd = abs_path(getcwd());
  ($fi, $fo, $d_hmm) = map {abs_path $_} ($fi, $fo, $d_hmm);
  my $dir = "$fo.dir";
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  
  my $f_bin = $ENV{"HMMER"}."/bin/hmmsearch";
  -x $f_bin || $log->error_die("cannot execute $f_bin");
  
  -d "01_seq" || make_path("01_seq");
  
  my $t = readTable(-in => $fi, -header => 1);
  my $ref = group($t->colRef("par"));
  
  open(my $fhc, ">03.cmds") or die "cannot write 03.cmds\n";
  for my $par (sort keys(%$ref)) {
    my ($idx, $cnt) = @{$ref->{$par}};
    my $h_seq;
    my $fam = $t->elm($idx, "cat3");
    for my $i ($idx..$idx+$cnt-1) {
      my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
      $h_seq->{$id} = Bio::Seq->new(-id=>$id, -seq=>$seq);
      die "family not consistent for gene[$par]\n" unless $fam eq $fam2;
    }

    my $f_hmm = "$d_hmm/$fam.hmm";
    -s $f_hmm || $log->error_die("HMM file [$f_hmm] not there");
    
    writeSeq([values(%$h_seq)], "01_seq/$par");
    print $fhc "$f_bin -o $par $f_hmm 01_seq/$par\n";
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");

  open(my $fho, ">$fo") || die "cannot write $fo\n";
  print $fho join("\t", qw/id score/)."\n";
  for my $par (sort keys(%$ref)) {
    my ($idx, $cnt) = @{$ref->{$par}};
    my $h = score_hmm_by_hit($par);
    for my $i ($idx..$idx+$cnt-1) {
      my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
      my $score = exists $h->{$id} ? $h->{$id} : 0;
      print $fho join("\t", $id, $score)."\n";
    }
  }
  close $fho;
  
  chdir $cwd || die "cannot chdir to $cwd\n";
  runCmd("rm -rf $dir", 0);
}
sub get_stat_aln {
  my ($fi, $d_aln, $f_sta, $fo) = @_;
  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("computing MSA scores");
  
  my $cwd = abs_path(getcwd());
  ($fi, $fo, $d_aln, $f_sta) = map {abs_path $_} ($fi, $fo, $d_aln, $f_sta);
  my $dir = "$fo.dir";
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
 
  my $hs;
  my $ts = readTable(-in => $f_sta, -header => 1);
  for my $i (0..$ts->lastRow) {
    my ($fam, $nseq, $npos, $gap, $seq) = $ts->row($i);
    next if $gap == 1;
    my $f_aln = "$d_aln/$fam.aln";
    my $h_seq = read_aln_seq($f_aln);
    my @seqs = map { Bio::Seq->new(-id=>$_, -seq=>$h_seq->{$_}) } keys(%$h_seq);
    $hs->{$fam} = \@seqs;
  }
  
  my $f_bin = $ENV{"ClustalO"}."/bin/clustalo";
  -x $f_bin || $log->error_die("cannot execute $f_bin");
  
  -d "01_seq" || make_path("01_seq");
  
  my $t = readTable(-in => $fi, -header => 1);
  my $ref = group($t->colRef("par"));
  
  open(my $fhc, ">03.cmds") or die "cannot write 03.cmds\n";
  for my $par (sort keys(%$ref)) {
    my ($idx, $cnt) = @{$ref->{$par}};
    my $fam = $t->elm($idx, "cat3");
    my (@ids, @seqs);
    for my $i ($idx..$idx+$cnt-1) {
      my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
      push @ids, $id;
      push @seqs, Bio::Seq->new(-id=>$id, -seq=>$seq);
      die "family not consistent for gene[$par]\n" unless $fam eq $fam2;
    }
    
    if(exists $hs->{$fam}) {
      writeSeq([@seqs, @{$hs->{$fam}}], "01_seq/$par");
      print $fhc "$f_bin -i 01_seq/$par --outfmt=fasta --force -o $par\n";
    } else {
      my $f_ai = "$d_aln/$fam.aln";
      -s $f_ai || $log->error_die("aln [$f_ai] not there");
      writeSeq(\@seqs, "01_seq/$par");
      
      my $tag_input = @ids == 1 ? 
        "--p1 01_seq/$par --p2 $f_ai" : "-i 01_seq/$par --p1 $f_ai";
      print $fhc "$f_bin $tag_input --outfmt=fasta --force -o $par\n";
    }
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");

  open(my $fho, ">$fo") || die "cannot write $fo\n";
  print $fho join("\t", qw/id score/)."\n";
  for my $par (sort keys(%$ref)) {
    my ($idx, $cnt) = @{$ref->{$par}};
    my (@ids, @seqs);
    for my $i ($idx..$idx+$cnt-1) {
      my ($id, $fam2, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
      push @ids, $id;
    }
    my $h = aln_score_group($par, \@ids);
    for my $id (@ids) {
      my $score = exists $h->{$id} ? $h->{$id} : 0;
      print $fho join("\t", $id, $score)."\n";
    }
  }
  close $fho;
  
  chdir $cwd || die "cannot chdir to $cwd\n";
  runCmd("rm -rf $dir", 0);
}
sub get_stat_pep {
  my ($fi, $fo) = @_;
  
  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("assessing peptide scores");
  
  my $tg = readTable(-in=>$fi, -header=>1);
  open(my $fho, ">$fo") or die "cannot open $fo for writing\n";
  print $fho join("\t", qw/id codonStart codonStop preStop gap n_cds lenC lenI/)."\n";
  for my $i (0..$tg->nofRow-1) {
    my ($id, $par, $chr, $srdG, $locCS, $locIS, $phaseG, $seq) = 
      map {$tg->elm($i, $_)} qw/id par chr srd cloc iloc phase seq/;
    my $locC = locStr2Ary($locCS);
    my $locI = locStr2Ary($locIS);

    my ($codonStart, $codonStop, $preStop, $gap) = checkProtSeq($seq);
    
    my $n_cds = @$locC;
    my $lenC = locAryLen($locC);
    my $lenI = locAryLen($locI);
    print $fho join("\t", $id, $codonStart, $codonStop, $preStop, $gap, $n_cds, $lenC, $lenI)."\n";
  }
  close $fho;
}
sub get_stat_sigp {
  my ($fi, $fo) = @_;
  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("assessing signalp scores");
  
  my $cwd = abs_path(getcwd());
  ($fi, $fo) = map {abs_path $_} ($fi, $fo);
  my $dir = "$fo.dir";
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  
  -d "01_seq" || make_path("01_seq");
  
  my $f_bin = $ENV{"SignalP"}."/signalp";
  -x $f_bin || die "$f_bin not there\n";
 
  my $tg = readTable(-in => $fi, -header => 1);
  open(my $fhc, ">03.cmds") or die "cannot write 03.cmds\n";
  for my $i (0..$tg->nofRow-1) {
    my ($id, $seq) = map {$tg->elm($i, $_)} qw/id seq/;
    $seq =~ s/\*$//;
    writeFile("01_seq/$id", ">tmp", $seq);
    print $fhc "perl $f_bin -t euk -s notm 01_seq/$id > $id\n";
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");
  
  open(my $fho, ">$fo") || die "cannot write $fo\n";
  print $fho join("\t", qw/id tag score pos/)."\n";
  for my $i (0..$tg->nofRow-1) {
    my ($id, $seq) = map {$tg->elm($i, $_)} qw/id seq/;
    my ($tag, $d, $pos) = parse_sigp($id);
    print $fho join("\t", $id, $tag, $d, $pos)."\n";
  }
  close $fho;
  
  chdir $cwd || die "cannot chdir to $cwd\n";
  runCmd("rm -rf $dir", 0);
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

  my $tg = readTable(-in => $f_gtb, -header => 1);
  my $ts = readTable(-in => $f_stat, -header => 1);
  $log->error_die("not equal rows $f_gtb - $f_stat") unless $tg->nofRow == $ts->nofRow;

  my $h;
  for my $i (0..$ts->nofRow-1) {
    my ($id, $par, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop) 
      = $ts->row($i);
    $h->{$par} ||= [];
    push @{$h->{$par}} , [$id, $tag_sp, $score_hmm+$score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop, $tg->rowRef($i)];
  }

  open(my $fh, ">$fo") || die "cannot write $fo\n";
  print $fh join("\t", $tg->header)."\n";
  my $cnt = 0;
  for my $par (sort (keys(%$h))) {
    my @rows = sort {$a->[1]<=>$b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @{$h->{$par}};
    my ($id, $tag_sp, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop, $row) = @{$rows[-1]};
    print $fh join("\t", @$row)."\n";
    $cnt ++;
  }
  $log->info(sprintf "\t%5d / %5d picked", $cnt, $tg->nofRow);
  close $fh;
}
sub remove_ovlp_models {
  my ($f_stat, $fi, $fo) = @_;

  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("removing overlapping models");

  my $ts = readTable(-in => $f_stat, -header => 1);
  my $h;
  for my $i (0..$ts->nofRow-1) {
    my ($id, $e, $score_aln) = map {$ts->elm($i, $_)} qw/id e score_aln/;
    $h->{$id} = [$e, $score_aln];
  }

  my ($hl, $hs);
  my $ti = readTable(-in => $fi, -header => 1);
  for my $i (0..$ti->nofRow-1) {
    my ($id, $par, $chr, $beg, $end) = map {$ti->elm($i, $_)} qw/id par chr beg end/;
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
  my ($f_stat, $fi, $fo, $co_e, $co_aln, $opt_sp, $opt_codon) = 
    rearrange([qw/stat in out e aln sp codon/], @_);
  $co_e   ||= 10;
  $co_aln ||= -1000;
  $opt_sp ||= 0;
  $opt_codon ||= 0;

  my $log = Log::Log4perl->get_logger("ModelEval");
  $log->info("final filter:");
  $log->info("\tSignal Peptide Present = $opt_sp");
  $log->info("\tComplete ORF = $opt_codon");
  $log->info("\tE value cutoff = $co_e");
  $log->info("\tAlignment score cutoff = $co_aln");
  
  my $ts = readTable(-in => $f_stat, -header => 1);
  my $h;
  for my $i (0..$ts->nofRow-1) {
    my ($id, $par, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $lenC, $lenI, $codonStart, $codonStop, $preStop) 
      = $ts->row($i);
    $h->{$id} = [$tag_sp, $score_hmm, $score_aln, $score_sp, $e, $codonStart, $codonStop, $preStop];
  }

  my $t = readTable(-in => $fi, -header => 1);
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
    }
  }
  
  $t->delRows(\@idxs_rm);
  $log->info(sprintf "\t%5d / %5d passed", $t->nofRow, $t->nofRow+@idxs_rm);

  open(FH, ">$fo") || $log->error_die("cannot open $fo for writing");
  print FH $t->tsv(1);
  close FH;
}

sub crp_rename {
  my ($fi, $f_ref, $fo) = @_;
  my $ti = readTable(-in => $fi, -header => 1);
  $ti->sort("cat3", 1, 0, "chr", 1, 0, "beg", 0, 0);
  
  my @chrs = uniq($ti->col("chr"));
  my @len_digits = map {getDigits(seqLen($_, $f_ref) / 1000000)} @chrs;
#  my $chr_digits = getDigits(scalar(grep /\d+/, @chrs));
  my $h = { map {$chrs[$_] => $len_digits[$_]} 0..$#chrs };

  for my $i (0..$ti->nofRow-1) {
    my ($chr, $beg, $fam) = map {$ti->elm($i, $_)} qw/chr beg cat3/;
    my $begStr = sprintf "%0".$h->{$chr}."d", $beg/1000000;
#    my $chrStr = $chr;
#    $chrStr =~ s/chr//i;
#    $chrStr = sprintf "%0".$chr_digits."d", $chrStr if $chrStr =~ /^\d+$/;

    my $id = sprintf "%s_%s_%sM", lc($fam), $chr, $begStr;
    $ti->setElm($i, "par", $ti->elm($i, "id"));
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
  my $t = readTable(-in => $fi, -header => 1);
  
  my $hs;
  my $ts = readTable(-in=>$f_stat, -header=>1);
  for my $i (0..$ts->nofRow-1) {
    my ($id, @stats) = map {$ts->elm($i, $_)} qw/id parent family e tag_sp score_sp score_hmm score_aln n_cds seq/;
    die "$id read twice in $f_stat\n" if exists $hs->{$id};
    $hs->{$id} = \@stats;
  }

  open(FH, ">$fo") or die "cannot open $fo for writing\n";
  print FH join("\t", qw/id family chr beg end srd e score_sp score_hmm score_aln sequence/)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($id, $chr, $beg, $end, $srd) = map {$t->elm($i, $_)} qw/id chr beg end srd/;
    my ($pa, $fam, $e, $tag_sp, $score_sp, $score_hmm, $score_aln, $n_cds, $seq) 
      = @{$hs->{$t->elm($i, "par")}};
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
  my $t = readTable(-in => $f_gtb, -header => 1);

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
  
  my $cwd = abs_path(getcwd());
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
 
  my ($eval_e, $eval_hmm, $eval_aln, $eval_pep) = (1) x 4;
  my $eval_sp = $ENV{"eval_sp"};

  get_stat_basic($f_gtb, $f_hit, '30_stat_basic.tbl');
  $eval_hmm && get_stat_hmm($f_gtb, $d_hmm, '32_hmm.tbl');
  $eval_aln && get_stat_aln($f_gtb, $d_aln, $f_sta, '33_aln.tbl');
  $eval_pep && get_stat_pep($f_gtb, '34_pep.tbl');
  $eval_sp  && get_stat_sigp($f_gtb, '35_sigp.tbl');
  my $p = {
    "hmm" => [ $eval_hmm, '32_hmm.tbl' ],
    "aln" => [ $eval_aln, '33_aln.tbl' ],
    "pep" => [ $eval_pep, '34_pep.tbl' ],
    "sp"  => [ $eval_sp,  '35_sigp.tbl' ] };
  merge_stats('30_stat_basic.tbl', '41_stat.tbl', $p);
  pick_best_model($f_gtb, '41_stat.tbl', '51_best.gtb', $p);
  remove_ovlp_models('41_stat.tbl', '51_best.gtb', '55_nonovlp.gtb');
  filter_models(-stat=>'41_stat.tbl', -in=>'55_nonovlp.gtb', -out=>'59.gtb', 
    -e=>$ENV{"evalue"}, -aln=>-1000, -sp=>$eval_sp, -codon=>1);
  
  crp_rename('59.gtb', $f_ref, '61_final.gtb');
  runCmd("gtb2gff.pl -i 61_final.gtb -o 61_final.gff", 1);
  gtb2Friendly('61_final.gtb', "61_final.tbl", '41_stat.tbl');
  
  align_by_group(-f_gtb=>'61_final.gtb', -hmmstat=>$f_sta, -out=>'81_aln');
  if(defined($f_gtb_ref) && -s $f_gtb_ref) {
    runCmd("gtb.cmp.pl -q 61_final.gtb -t $f_gtb_ref -o 91_compare.tbl", 1);
  }

  chdir $cwd || die "cannot chdir to $cwd\n";
}


1;
__END__
