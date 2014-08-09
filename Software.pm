package Software;
use strict;
use Cwd qw/getcwd abs_path/;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Common;
use Location; 
use Seq;
use Gtb;
use Align;
use Data::Dumper;
use Log::Log4perl;
use List::Util qw/min max sum/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_model_run/;
@EXPORT_OK = qw//;

sub pipe_model_run {
  my ($dir, $f_hit, $f_ref, $soft) = 
    rearrange(['dir', 'hit', 'ref', 'soft'], @_);
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("### working on $soft pipeline ###");
  
  if($soft eq "GeneWise_SplicePredictor") {
    pipe_genewise_splicepredictor($f_hit, $dir);
  } else {
    pipe_software($f_hit, $dir, $soft);
  }
}
sub pipe_software {
  my ($fi, $dir, $opt) = @_;
  my $cwd = abs_path(getcwd());
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $cwd\n";

  if($opt eq "Augustus_evidence") {
    run_augustus_batch($fi, "02_raw");
  } elsif($opt eq "Augustus_de_novl") {
    run_augustus_batch_simple($fi, "02_raw");
  } elsif($opt eq "Genemark") {
    run_genemark_batch($fi, "02_raw");
  } elsif($opt eq "Glimmerhmm") {
    run_glimmerhmm_batch($fi, "02_raw");
  } elsif($opt eq "GeneID") {
    run_geneid_batch($fi, "02_raw");
  } else {
    die "unknown opt: $opt\n";
  }
  sum_gff_batch($fi, "02_raw", "11.gtb");

  chdir $cwd || die "cannot chdir to $cwd\n";
}
sub pipe_genewise_splicepredictor {
  my ($fi, $dir) = @_;
  my $cwd = abs_path(getcwd());
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  runCmd("rm -rf $dir/*", 0);

  run_genewise_batch($fi, "02_raw");
  build_model($fi, "05.tbl", "02_raw");
  output_gtbx("05.tbl", "11.gtb");
  
  chdir $cwd || die "cannot chdir to $cwd\n";
}

sub run_augustus_batch {
  my ($fi, $dir) = @_;
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running Augustus(evidence mode)");
  
  my $cwd = abs_path(getcwd());
  $fi = abs_path($fi);
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  
  -d "01_seq" || make_path("01_seq");
  -d "02_hint" || make_path("02_hint");

  my $f_bin = "$ENV{'Augustus'}/bin/augustus";
  -x $f_bin || die "$f_bin not there\n";
  my $d_cfg = "$ENV{'Augustus'}/config";
  -d $d_cfg || die "$d_cfg not there\n";
  my $species = 'arabidopsis';
  if($ENV{'SPADA_ORG'} eq "Athaliana") {
    -d "$d_cfg/species/arabidopsis" || die "no augustus athal cfg\n";
  } elsif($ENV{'SPADA_ORG'} eq "Mtruncatula") {
    $species = "medicago" if -d "$d_cfg/species/medicago";
  } elsif($ENV{'SPADA_ORG'} eq "Osativa") {
    $species = "osativa" if -d "$d_cfg/species/osativa";
  }
  $log->info("using $species matrix");

  open(my $fhc, ">03.cmds") or die "cannot write 03.cmds\n";  
  my $t = readTable(-in=>$fi, -header=>1);
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 2576;
    writeFile("01_seq/$id", ">tmp", $seq);
    
    my $locL = locStr2Ary($locLS);
    my @hints_hit = map {join("\t", "tmp", "pz", "CDSpart", @$_, ".", "+", ".", "source=M")} @$locL;
    writeFile("02_hint/$id", @hints_hit);
  
    print $fhc "$f_bin --AUGUSTUS_CONFIG_PATH=$d_cfg" . 
      " --species=$species --gff3=on --strand=forward" .
      " --noInFrameStop=true --hintsfile=02_hint/$id 01_seq/$id > $id\n";
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");
  runCmd("rm -rf 01_seq 02_hint", 0);

  chdir $cwd || die "cannot chdir to $cwd\n";
}
sub run_augustus_batch_simple {
  my ($fi, $dir) = @_;
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running Augustus(evidence mode)");

  my $cwd = abs_path(getcwd());
  $fi = abs_path($fi);
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  
  -d "01_seq" || make_path("01_seq");

  my $f_bin = "$ENV{'Augustus'}/bin/augustus";
  -x $f_bin || die "$f_bin not there\n";
  my $d_cfg = "$ENV{'Augustus'}/config";
  -d $d_cfg || die "$d_cfg not there\n";
  my $species = 'arabidopsis';
  if($ENV{'SPADA_ORG'} eq "Athaliana") {
    -d "$d_cfg/species/arabidopsis" || die "no augustus athal cfg\n";
  } elsif($ENV{'SPADA_ORG'} eq "Mtruncatula") {
    $species = "medicago" if -d "$d_cfg/species/medicago";
  } elsif($ENV{'SPADA_ORG'} eq "Osativa") {
    $species = "osativa" if -d "$d_cfg/species/osativa";
  }
  $log->info("using $species matrix");

  open(my $fhc, ">03.cmds") or die "cannot write 03.cmds\n";  
  my $t = readTable(-in=>$fi, -header=>1);
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
    writeFile("01_seq/$id", ">tmp", $seq);
    
    print $fhc "$f_bin --AUGUSTUS_CONFIG_PATH=$d_cfg" . 
      " --species=$species --gff3=on --strand=forward" .
      " --noInFrameStop=true 01_seq/$id > $id\n";
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");
  runCmd("rm -rf 01_seq", 0);

  chdir $cwd || die "cannot chdir to $cwd\n";
}
sub run_genewise_batch {
  my ($fi, $dir) = @_;
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running genewise");
  
  my $cwd = abs_path(getcwd());
  $fi = abs_path($fi);
  -d $dir || make_path($dir);
  chdir $dir || die "cannot chdir to $dir\n";
  
  -d "01_dna" || make_path("01_dna");
  -d "02_pro" || make_path("02_pro");

  my $f_bin = $ENV{"GeneWise"}."/bin/genewise";
  $ENV{"WISECONFIGDIR"} = $ENV{"GeneWise"}."/wisecfg";
  -x $f_bin || die "$f_bin not there\n";
  
  open(my $fhc, ">03.cmds") || die "cannot write 03.cmds\n";
  my $t = readTable(-in => $fi, -header => 1);
  for my $i (0..$t->lastRow) {
    my ($id, $fam, $chr, $begG, $endG, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
    my $loc = locStr2Ary($locLS);
    next if @$loc == 1;
    writeFile("01_dna/$id", ">seq", $seq);
    writeFile("02_pro/$id", ">pro", $seqP);
    print $fhc "$f_bin -gff 02_pro/$id 01_dna/$id > $id.txt\n";
  }
  close $fhc;
  runCmd("parallel -j $ENV{'threads'} --no-notice < 03.cmds");
  runCmd("rm -rf 01_dna 02_pro", 0);

  chdir $cwd || die "cannot chdir to $cwd\n";
}
sub run_genemark_batch {
  my ($fi, $dir) = @_;
  make_path($dir) unless -d $dir;
  remove_tree($dir, {keep_root=>1});
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running GeneMark");
  
  my $f_bin = $ENV{"GeneMark"}."/gmhmme3";
  $log->error_die("$f_bin not there") unless -s $f_bin;

  my $f_mod;
  if($dir =~ /Athaliana/) {
    $f_mod = "a_thaliana.mod";
  } elsif($dir =~ /Mtruncatula/) {
    $f_mod = "m_truncatula.mod";
  } elsif($dir =~ /Osativa/) {
    $f_mod = "o_sativa.mod";
  } else {
    $log->warn("No proper MOD file to use -> using a_thaliana.mod instead");
    $f_mod = "a_thaliana.mod";
  }
  $f_mod = $ENV{"GeneMark"}."/".$f_mod;

  my $t = readTable(-in=>$fi, -header=>1);
  my $f_fas = $ENV{"TMP_DIR"}."/genemark_".int(rand(1000)).".fa";
  
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 1;
    writeFile($f_fas, ">tmp", $seq);
    my $fo = "$dir/$id";
    runCmd("$f_bin -m $f_mod -f gff3 -o $fo $f_fas", 0);
    printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
  }
  system("rm $f_fas"); 
}
sub run_glimmerhmm_batch {
  my ($fi, $dir) = @_;
  make_path($dir) unless -d $dir;
  remove_tree($dir, {keep_root=>1});
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running GlimmerHMM");
  
  my $f_bin = $ENV{"GlimmerHMM"}."/bin/glimmerhmm";
  $log->error_die("$f_bin not there") unless -s $f_bin;

  my $d_train;
  if($dir =~ /Athaliana/) {
    $d_train = "arabidopsis";
  } elsif($dir =~ /Osativa/) {
    $d_train = "rice";
  } else {
    $log->warn("no training directory available -> using arabidopsis instead");
    $d_train = "arabidopsis";
  }
  $d_train = $ENV{"GlimmerHMM"}."/trained_dir/$d_train";

  my $t = readTable(-in=>$fi, -header=>1);
  my $f_fas = $ENV{"TMP_DIR"}."/glimmerhmm_".int(rand(1000)).".fa";
  
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 1;
    writeFile($f_fas, ">tmp", $seq);
    my $fo = "$dir/$id";
    runCmd("$f_bin $f_fas $d_train -g > $fo", 0);
    printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
  }
  system("rm $f_fas"); 
}
sub run_geneid_batch {
  my ($fi, $dir) = @_;
  make_path($dir) unless -d $dir;
  remove_tree($dir, {keep_root=>1});
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("running geneid");

  my $f_bin = $ENV{"GeneID"}."/bin/geneid";
  $log->error_die("$f_bin not there") unless -s $f_bin;

  my $t = readTable(-in=>$fi, -header=>1);
  my $f_fas = $ENV{"TMP_DIR"}."/geneid_".int(rand(1000)).".fa";

  my $f_param;
  if($dir =~ /Athaliana/) {
    $f_param = "arabidopsis";
  } elsif($dir =~ /Osativa/) {
    $f_param = "rice";
  } else {
    $log->warn("no param file available -> using arabidopsis instead");
    $f_param = "arabidopsis";
  }
  $f_param = $ENV{"GeneID"}."/param/$f_param.param";

  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 1;
    writeFile($f_fas, ">tmp", $seq);
    my $fo = "$dir/$id";
    runCmd("$f_bin -3 -W -P $f_param $f_fas > $fo", 0);
    printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
  }
  system("rm $f_fas"); 
}
sub sum_gff {
  my ($fi, $locH) = @_;
  my $h = {};
  open(IN, "<$fi") or die "Failed: $!\n";
  while(<IN>) {
    chomp;
    my @ps = split "\t";
    next if @ps < 9 || $ps[6] ne "+" || $ps[2] ne "CDS";
    my $note = $ps[8];
    $note .= ";" if $note !~ /\;$/;
    if($note =~ /Parent=([\w\.\_\-]+)\;/) {
      my $id = $1;
      $h->{$id} ||= [];
      push @{$h->{$id}}, \@ps;
    } else {
      die "cannot get mRNA ID from\n\t".join("\t", @ps)."\n";
    }
  }
  return undef if keys(%$h) == 0;
  
  my @stats;
  for my $id (keys(%$h)) {
    my $row = $h->{$id};
    $row = [ sort {$a->[3] <=> $b->[3]} @$row ];
    $row->[0]->[3] += $row->[0]->[7];
    $row->[0]->[7] = 0;

    my $len = sum( map {$_->[4] - $_->[3] + 1} @$row );
    $row->[-1]->[4] -= $len % 3;

    my @locs = map {[$_->[3], $_->[4]]} @$row;
    my @phases = map {$_->[7]} @$row;
    my ($locO, $lenO) = posOvlp(\@locs, $locH);
    push @stats, [\@locs, $lenO];
  }
  @stats = sort {$a->[1] <=> $b->[1]} @stats;
  return undef if $stats[-1]->[1] == 0;
  return $stats[-1]->[0];
}
sub sum_gff_batch {
  my ($fi, $dir, $fo) = @_;
  my $log = Log::Log4perl->get_logger("Software");
  $log->info("collecting prediction results");
  my $t = readTable(-in=>$fi, -header=>1);
  
  open(my $fho, ">$fo") || die "cannot write to $fo\n";
  print $fho join("\t", @HEAD_GTBX)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
    my $locL = locStr2Ary($locLS);

    my $fm = "$dir/$id";
    -s $fm || next;
    my $locC = sum_gff($fm, $locL);
    next unless $locC;

    my $seq_cds = getSubSeq($seq, $locC);
    my $seq_pro = Bio::Seq->new(-seq=>$seq_cds)->translate()->seq;
    my $locCStr = locAry2Str($locC);
    my $phase = join(",", @{getPhase($locC, "+")});
    
    print $fho join("\t", "$id.1", $id, ("") x 3, "+", ("") x 2, $locCStr, "", "", $phase, "", ("") x 5, $seq_pro)."\n";
  }
  close $fho;
  runCmd("rm -rf $dir", 0);
}
sub sum_augustus1 {
  my ($fi) = @_;
  my $h;
  open(IN, "<$fi") or die "Failed: $!\n";
  my $gene;
  while(<IN>) {
    chomp;
    if(/^\# start gene (\w+)/) {
      $gene = $1;
      $h->{$gene} ||= [];
    } elsif(/^\# end gene (\w+)/) {
      $gene = "";
    } elsif(/^\# hint groups fully obeyed:\s*(\d+)/) {
      $h->{$gene}->[0] = $1;
    }
    next if /^\#/;
    my @ps = split "\t";
    next unless @ps >= 9;
    if($ps[2] eq "CDS" && $ps[8] =~ /t1$/) {
      $h->{$gene}->[1] ||= [];
      push @{$h->{$gene}->[1]}, \@ps;
    }
  }
  my @genes = grep {$h->{$_}->[0] > 0} (keys(%$h));
  if(@genes > 0) {
    @genes = sort {$h->{$a}->[0] <=> $h->{$b}->[0]} @genes;
    my $row = $h->{$genes[-1]}->[1];
    die Dumper($row)." not on + strand\n" unless $row->[0]->[6] eq "+";
    
    $row = [ sort {$a->[3] <=> $b->[3]} @$row ];
    $row->[0]->[3] += $row->[0]->[7];
    $row->[0]->[7] = 0;

    my $len = sum( map {$_->[4] - $_->[3] + 1} @$row );
    $row->[-1]->[4] -= $len % 3;

    my @locs = map {[$_->[3], $_->[4]]} @$row;
    my @phases = map {$_->[7]} @$row;
    return (\@locs, \@phases);
  } else {
    return undef;
  }
}
sub crop_cds { # cut model from specified start codon
  my ($loc, $pos) = @_;
  my $locN = [];
  my $len_crop = 0;
  for (sort {$a->[0] <=> $b->[0]} @$loc) {
    my ($beg, $end) = @$_;
    if($pos > $beg) {
      if($pos <= $end) {
        $len_crop += $pos-1 - $beg + 1;
        push @$locN, [$pos, $end];
      } else {
        $len_crop += $end - $beg + 1;
      }
    } else {
      push @$locN, [$beg, $end];
    }
  }
  return ($locN, $len_crop);
}

sub sum_genewise1 {
  my ($fi) = @_;
  open(FHT, "<$fi") || die "cannot open $fi for reading\n";
  my $loc = [];
  while(<FHT>) {
    chomp;
    my @ps = split "\t";
    if(@ps == 9 && $ps[1] eq "GeneWise" && $ps[2] eq "cds") {
      push @$loc, [$ps[3], $ps[4]];
    } 
  }
  close FHT;
  return $loc;
}
sub pick_valid_splice_sites {
  my ($stat, $beg, $end, $seq) = @_;
  
  my @poss_d = ();
  my @poss_a = ();
  for (@$stat) {
    my ($head, $type, $id, $pos, $c, $U, $s, $p, $rho, $gamma, $parse, $q1, $q2, $q3, $q4, $q5, $q) = @$_;
    $q =~ s/[^\d]//g;
    $pos = int($pos);
    push @poss_d, $pos if($type eq "DONOR");
    push @poss_a, $pos if($type eq "ACPTR");
  }

  my @pairs;
  for my $pos_d (@poss_d) {
    my @tmp = grep {$_ > $pos_d} @poss_a;
    for my $pos_a (@tmp) {
      my $tag = 1;
      my $len_padding = $pos_d-$beg + $end-$pos_a;
      $tag = 0 if $len_padding % 3 != 0;
      if($len_padding >= 3) {
        my $dna = getSubSeq($seq, [[$beg, $pos_d-1], [$pos_a+1, $end]]); 
        my $prot = Bio::Seq->new(-seq=>$dna)->translate()->seq;
        $tag = 0 if $prot =~ /\*/;
      }
      push @pairs, [$pos_d, $pos_a] if $tag == 1;
    }
  }
  return @pairs;
}
sub get_splice_sites {
  my ($seq, $beg, $end) = @_;

  my $f_bin = $ENV{"SplicePredictor"}."/bin/SplicePredictor";
  -x $f_bin || die "$f_bin not there\n";
  my $f_fas = $ENV{"TMP_DIR"}."/sppr_".int(rand(1000)).".fas";
  writeFile($f_fas, ">tmp", $seq);

  my $cmd = "$f_bin -s Arabidopsis -c -99.9 -p 5 -a $beg -b $end -L $f_fas";
  my $lines = runCmd($cmd, 2);

  my @stats = map { [split("\t", $_)] } @$lines;
  @stats = grep {defined($_->[1]) && $_->[1] =~ /^(ACPTR)|(DONOR)$/} @stats;
  my @stats_f = grep {$_->[4] >=0} @stats;

  my @pairs = pick_valid_splice_sites(\@stats_f, $beg, $end, $seq);
  if(@pairs == 0) {
    @pairs = pick_valid_splice_sites(\@stats, $beg, $end, $seq);
  }
  system("rm $f_fas");
  return @pairs;
}
sub build_model {
  my ($fi, $fo, $dirG) = @_;
  my $log = Log::Log4perl->get_logger("Spada");
  $log->info("building exon models");
  my $t = readTable(-in=>$fi, -header=>1);
  
  open(FH, ">$fo") || die "cannot open $fo for writing\n";
  print FH join("\t", qw/id parent family beg end strand n_cds locC seq/)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($id, $fam, $chr, $begG, $endG, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
    my $loc = locStr2Ary($locLS);
    $loc = [sort {$a->[0] <=> $b->[0]} @$loc];
    my $pos_end = get_stop_codon($seq, $endL+1);
    my @pos_begs = get_start_codons($seq, $begL+2);
    my $n_cds = @$loc;

    my @locCs;
    if($n_cds == 1) {
      for my $pos_beg (@pos_begs) {
        push @locCs, [[$pos_beg, $pos_end]];
      }
    } else {
      my $locC = sum_genewise1("$dirG/$id.txt");
#      die Dumper($locC, $pos_end, \@pos_begs) if $id == 200;
      if(@$locC == $n_cds) {
        $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
        next if @pos_begs == 0 || $locC->[0]->[0] < max(@pos_begs) || $locC->[-1]->[1] > $pos_end;
        for my $pos_beg (@pos_begs) {
          my $locT = [ map {[$_->[0], $_->[1]]} @$locC ];
          $locT->[0]->[0] = $pos_beg;
          $locT->[-1]->[1] = $pos_end;
          push @locCs, $locT;
        }
      } elsif($n_cds == 2) {
        my @posIs = get_splice_sites($seq, $loc->[0]->[1]+1, $loc->[1]->[0]-1);
        for (@posIs) {
          my ($endE2, $begE1) = ($_->[0]-1, $_->[1]+1);
          for my $pos_beg (@pos_begs) {
            push @locCs, [[$pos_beg, $endE2], [$begE1, $pos_end]];
          }
        }
      }
    }
    
    for my $i (1..@locCs) {
      my $locC = $locCs[$i-1];
      $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
      my ($beg, $end) = ($locC->[0]->[0], $locC->[-1]->[1]);
      my $seq = getSubSeq($seq, $locC);
      my $prot = Bio::Seq->new(-seq=>$seq)->translate()->seq;
      print FH join("\t", "$id.$i", $id, $fam, $beg, $end, "+", $n_cds, locAry2Str($locC), $prot)."\n";
    }
  }
  close FH;
}
sub output_gtbx {
  my ($fi, $fo) = @_;
  my $t = readTable(-in=>$fi, -header=>1);
  
  open(my $fh, ">$fo") || die "cannot write $fo\n";
  print $fh join("\t", @HEAD_GTBX)."\n";
  for my $i (0..$t->nofRow-1) {
    my ($id, $par, $fam, $beg, $end, $strand, $n_cds, $locCStr, $seq) = $t->row($i);
    my $locC = locStr2Ary($locCStr);
    my $phase = join(",", @{getPhase($locC, "+")});
    print $fh join("\t", $id, $par, ("") x 3, "+", ("") x 2, $locCStr, "", "", $phase, "", ("") x 5, $seq)."\n";
  }
  close $fh;
}





1;
__END__
