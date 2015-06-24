package Align;
use strict;
use Bio::AlignIO;
use Bio::Matrix::IO;
use Common;
use Location;
use Seq;
use Data::Dumper; 
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/read_aln_ids aln_fmt_convert
  run_clustalw run_tcoffee run_clustalo run_pal2nal
  run_pw_aln
  aln_score_pw aln_score_vector aln_score_group aln_score_matrix/;
@EXPORT_OK = qw//;

sub read_aln_ids {
    my ($fi, $format) = @_;
    $format ||= "clustalw";
    my $ai = Bio::AlignIO->new(-file=>"<$fi", -format=>$format);
    my @ids;
    while(my $aln = $ai->next_aln()) {
        for my $seq ($aln->each_seq()) {
            my $id = $seq->id;
            push @ids, $id;
        }
    }
    return @ids;
}
sub aln_fmt_convert {
  my ($fi, $fo, $if, $of) = @_;
  my $ai = Bio::AlignIO->new(-file=>"<$fi", -format=>$if);
  my $ao = Bio::AlignIO->new(-file=>">$fo", -format=>$of);
  while(my $aln = $ai->next_aln()) {
    $ao->write_aln($aln);
  }
}
sub run_clustalw {
    my ($seqAry, $fOut, $type) = rearrange(['seqs', 'out', 'type'], @_);
    $type = uc($type);
    my $f_bin = $ENV{"ClustalW"}."/bin/clustalw2";
    my $f1 = $ENV{"TMP_DIR"}."/clustalw_".int(rand(1000)).".fa";
    writeSeq($seqAry, $f1);
    my $cmd = qq/$f_bin -INFILE=$f1 -ALIGN -OUTFILE=$fOut -OUTORDER=INPUT -TYPE=$type/;
    my $matrixP ||= "BLOSUM";
    my $matrixD ||= "CLUSTALW";
    $cmd .= " -PWMATRIX=$matrixP -MATRIX=$matrixP" if $type eq "PROTEIN";
    $cmd .= " -PWDNAMATRIX=$matrixD -DNAMATRIX=$matrixD" if $type eq "DNA";
    runCmd($cmd, 0);
    my $fTree = $f1;
    $fTree =~ s/\.\w+$/\.dnd/i;
    system("rm $fTree $f1"); 
}
sub run_clustalo {
  my ($seqs, $fi, $fo, $f_dis, $f_hmm) = rearrange(['seqs', 'in', 'out', 'dis', 'f_hmm'], @_);
  
  my $mode = defined($fi) ? 1 : 2;
  if($mode == 2) {
    if( defined $ENV{"TMP_DIR"} ) {
      $fi = $ENV{"TMP_DIR"}."/clustalo_".int(rand(1000)).".fa";
    } else {
      $fi = "$fo.tmp";
    }
    writeSeq($seqs, $fi);
  }
  my $f_bin = $ENV{"ClustalO"} ? $ENV{"ClustalO"}."/bin/clustalo" : "clustalo";

  my $cmd = qq/$f_bin -i $fi -o $fo --outfmt=clu --force/;
  $cmd .= " --hmm-in=$f_hmm" if $f_hmm;
  $cmd .= " --full --distmat-out=$f_dis" if $f_dis;
  runCmd($cmd, 0);
  system("rm $fi") if $mode == 2; 
}
sub run_pal2nal {
    my ($fIn1, $fIn2, $fOut) = rearrange(['in1', 'in2', 'out'], @_);
    die "$fIn1 is not there\n" unless -s $fIn1;
    die "$fIn2 is not there\n" unless -s $fIn2;
    my $cmd = qq/pal2nal.pl $fIn2 $fIn1 > $fOut/;
    runCmd($cmd, 0);
}
sub run_tcoffee {
    my ($seqAry, $fo) = rearrange(['seqs', 'out'], @_);
    my $f1 = $ENV{"TMP_DIR"}."/tcoffee_".int(rand(1000)).".fa";
    writeSeq(-seqs=>$seqAry, -out=>$f1);
    my $cmd = qq/t_coffee $f1 -outfile $fo -gapopen -100 -gapext -5/;
    runCmd($cmd, 0);
    system("rm $f1 $fo.html"); 
    system("rm $fo.dnd") if -f "$fo.dnd";
    system("rm tcoffee.dnd") if -f "tcoffee.dnd";
}

sub run_pw_aln {
  my ($seq1, $seq2, $prog) = @_;
  $prog ||= "water";
  my $fi1 = "s1.".int(rand(1000)).".fa";
  my $fi2 = "s2.".int(rand(1000)).".fa";
  my $fo = "o".int(rand(1000));
  writeFile($fi1, ">seq1", $seq1);
  writeFile($fi2, ">seq2", $seq2);
  runCmd("$prog $fi1 $fi2 -gapopen 10 -gapextend 0.5 -aformat fasta -outfile $fo");
  
  my $seqHI = Bio::SeqIO->new(-file => "<$fo", -format => 'fasta');
  my $seqa1 = $seqHI->next_seq()->seq;
  my $seqa2 = $seqHI->next_seq()->seq;
  my $alen = length($seqa1);
  
  my @locs = ([1, $alen, 0]);
  while($seqa1 =~ /([\-\_ ]+)/ig) {
    push @locs, [$-[0] + 1, $+[0], 1];
  }
  while($seqa2 =~ /([\-\_ ]+)/ig) {
    push @locs, [$-[0] + 1, $+[0], 2];
  }
  my $ref = posSplit(\@locs);
  my (@loc1, @loc2);
  my ($p1, $p2) = (1, 1);
  for (@$ref) {
    my ($beg, $end, $idxs) = @$_;
    my $len = $end - $beg + 1;
    my @nidxs = grep {$_ != 0} @$idxs;
    @nidxs <= 1 || die "double gaps: $seqa1\n$seqa2\n";
    if(@nidxs == 0) {
      push @loc1, [$p1, $p1 + $len - 1];
      push @loc2, [$p2, $p2 + $len - 1];
      $p1 += $len;
      $p2 += $len;
    } elsif($nidxs[0] == 1) {
      $p2 += $len;
    } else {
      $p1 += $len;
    }
  }
  
  my ($ilen, $mlen) = (0, 0);
  my ($gap1, $gap2) = (0, 0);
  for my $i (0..$alen-1) {
    my ($ch1, $ch2) = map {substr($_, $i, 1)} ($seqa1, $seqa2);
    if($ch1 =~ /[\-\_ ]/ || $ch2 =~ /[\-\_ ]/) {
      $gap1 ++ if $ch1 =~ /[\-\_ ]/;
      $gap2 ++ if $ch2 =~ /[\-\_ ]/;
    } elsif($ch1 eq $ch2) {
      $ilen ++;
    } else {
      $mlen ++;
    }
  }
  $alen == $ilen + $mlen + $gap1 + $gap2 
    || die "aln parse err\n$seqa1\n$seqa2\n" ;
  system("rm $fi1 $fi2 $fo");
  return (\@loc1, \@loc2, $alen, $ilen, $mlen, $gap1, $gap2);
}

sub run_dotplot {
    my ($seqAry, $fo, $ws) = rearrange(['seqs', 'out', "ws"], @_);
    $ws ||= 10;
    my $fi1 = $ENV{"TMP_DIR"}."/dottup_1_".int(rand(1000)).".fa";
    my $fi2 = $ENV{"TMP_DIR"}."/dottup_2_".int(rand(1000)).".fa";
    writeSeq(-seqs=>$seqAry->[0], -out=>$fi1);
    writeSeq(-seqs=>$seqAry->[1], -out=>$fi2);
    my ($sid1, $sid2) = map {$_->id} @$seqAry;
    my $cmd = qq/dottup $fi1 $fi2 -graph png -word $ws/;
    runCmd($cmd, 0);
    system("rm $fi1 $fi2");
    if( -s "dottup.1.png" ) {
        system("mv dottup.1.png $fo");
    } else {
        die "dotplot failed\n";
    }
}

sub aln_score_pw {
    my ($seq1, $seq2) = @_;
    die "seq length not equal\n$seq1\n$seq2\n" unless length($seq1) == length($seq2);
    ($seq1, $seq2) = (uc($seq1), uc($seq2));
    $seq1 =~ s/\./\-/g;
    $seq2 =~ s/\./\-/g;
    
    my ($gap_open, $gap_ext) = (-8, -1);
    
    my $parser = Bio::Matrix::IO->new(-format=>'scoring', -file=>$ENV{"BLOSUM80"});
    my $matrix = $parser->next_matrix();

    my $score = 0;
    my $tag_gap_open = 1;
    for my $i (0..length($seq1)-1) {
        my $aa1 = substr($seq1, $i, 1);
        my $aa2 = substr($seq2, $i, 1);
        if($aa1 eq "-" && $aa2 eq "-") {
        } elsif($aa1 ne "-" && $aa2 ne "-") {
            $score += $matrix->get_entry($aa1, $aa2);
            $tag_gap_open = 1;
        } else {
            if($tag_gap_open == 1) {
                $score += $gap_open;
                $tag_gap_open = 0;
            } else {
                $score += $gap_ext;
            }
        }
    }
    return $score;
}
sub aln_score_vector {
    my ($fi, $id) = @_;
    my $h_seq = readSeq($fi, 2);
    die "$id not in $fi\n" unless exists $h_seq->{$id};
    my $seq1 = $h_seq->{$id};
    my @scores;
    for my $id2 (keys(%$h_seq)) {
        next if $id2 eq $id;
        my $seq2 = $h_seq->{$id2};
        push @scores, aln_score_pw($seq1, $seq2);
    }
    my $score = sprintf "%d", sum(@scores) / @scores;
    return $score;
}
sub aln_score_group {
  my ($fi, $idsQ, $idsT) = @_;
  my $h_seq = readSeq($fi, 2);

  my @idsQ_n = grep { !exists($h_seq->{$_}) } @$idsQ;
  die "idsQ not all found\n".Dumper($idsQ).Dumper($h_seq) if @idsQ_n;
  my $h = { map {$_ => 0} @$idsQ };
  $idsT ||= [ grep { !exists($h->{$_}) } keys(%$h_seq) ];

  for my $idQ (@$idsQ) {
    my @scores;
    for my $idT (@$idsT) {
      my $score = aln_score_pw($h_seq->{$idQ}, $h_seq->{$idT});
      push @scores, $score;
    }
    $h->{$idQ} = sprintf "%.01f", sum(@scores) / @scores;
  }
  return $h;
}
sub aln_score_matrix {
    my ($fi, $ids, $hi) = @_;

    my $h_seq = readSeq($fi, 2);
    my @seqs = map {$h_seq->{$_}} @$ids;
    my $nseq = @seqs;
    
    my ($hs, $m) = ({}, []);
    for my $i (0..$nseq-1) {
        $m->[$i] = [];
        my @js;
        if(defined($hi)) {
            @js = @{$hi->[$i]};
        } else {
            push @js, (0..$i-1) if $i >= 1;
            push @js, ($i+1, $nseq-1) if $i <= $nseq-2;
        }
        for my $j (@js) {
            if(exists($hs->{"$i.$j"})) {
                push @{$m->[$i]}, $hs->{"$i.$j"};
            } else {
                my $score = aln_score_pw($seqs[$i], $seqs[$j]);
                push @{$m->[$i]}, $score;
                $hs->{"$j.$i"} = $score;
            }
        }
    }
    return $m;
}


1;
__END__
