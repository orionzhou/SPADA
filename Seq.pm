package Seq;
use strict;
use Common; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/readSeq readSeqInfo writeSeq writeSeqInOrder read_aln_seq
  revcom translate translate6 seqret_simple
  get_start_codons get_stop_codon
  seqRet seqLen getLongestOrf getSubSeq slideSeq
  validate_dna_seq checkDnaSeq checkProtSeq
  seqCompare translate_smart/;
@EXPORT_OK = qw//;

my (%genetic_code) = (
  'TCA' => 'S',    # Serine
  'TCC' => 'S',    # Serine
  'TCG' => 'S',    # Serine
  'TCT' => 'S',    # Serine
  'TTC' => 'F',    # Phenylalanine
  'TTT' => 'F',    # Phenylalanine
  'TTA' => 'L',    # Leucine
  'TTG' => 'L',    # Leucine
  'TAC' => 'Y',    # Tyrosine
  'TAT' => 'Y',    # Tyrosine
  'TAA' => '*',    # Stop
  'TAG' => '*',    # Stop
  'TGC' => 'C',    # Cysteine
  'TGT' => 'C',    # Cysteine
  'TGA' => '*',    # Stop
  'TGG' => 'W',    # Tryptophan
  'CTA' => 'L',    # Leucine
  'CTC' => 'L',    # Leucine
  'CTG' => 'L',    # Leucine
  'CTT' => 'L',    # Leucine
  'CCA' => 'P',    # Proline
  'CCC' => 'P',    # Proline
  'CCG' => 'P',    # Proline
  'CCT' => 'P',    # Proline
  'CAC' => 'H',    # Histidine
  'CAT' => 'H',    # Histidine
  'CAA' => 'Q',    # Glutamine
  'CAG' => 'Q',    # Glutamine
  'CGA' => 'R',    # Arginine
  'CGC' => 'R',    # Arginine
  'CGG' => 'R',    # Arginine
  'CGT' => 'R',    # Arginine
  'ATA' => 'I',    # Isoleucine
  'ATC' => 'I',    # Isoleucine
  'ATT' => 'I',    # Isoleucine
  'ATG' => 'M',    # Methionine
  'ACA' => 'T',    # Threonine
  'ACC' => 'T',    # Threonine
  'ACG' => 'T',    # Threonine
  'ACT' => 'T',    # Threonine
  'AAC' => 'N',    # Asparagine
  'AAT' => 'N',    # Asparagine
  'AAA' => 'K',    # Lysine
  'AAG' => 'K',    # Lysine
  'AGC' => 'S',    # Serine
  'AGT' => 'S',    # Serine
  'AGA' => 'R',    # Arginine
  'AGG' => 'R',    # Arginine
  'GTA' => 'V',    # Valine
  'GTC' => 'V',    # Valine
  'GTG' => 'V',    # Valine
  'GTT' => 'V',    # Valine
  'GCA' => 'A',    # Alanine
  'GCC' => 'A',    # Alanine
  'GCG' => 'A',    # Alanine
  'GCT' => 'A',    # Alanine
  'GAC' => 'D',    # Aspartic Acid
  'GAT' => 'D',    # Aspartic Acid
  'GAA' => 'E',    # Glutamic Acid
  'GAG' => 'E',    # Glutamic Acid
  'GGA' => 'G',    # Glycine
  'GGC' => 'G',    # Glycine
  'GGG' => 'G',    # Glycine
  'GGT' => 'G',    # Glycine
  'NAA' => 'X',    
  'NAT' => 'X',    
  'NAC' => 'X',    
  'NAG' => 'X',    
  'NTA' => 'X',    
  'NTT' => 'X',    
  'NTC' => 'X',    
  'NTG' => 'X',    
  'NGA' => 'X',    
  'NGT' => 'X',    
  'NGC' => 'X',    
  'NGG' => 'X',    
  'NCA' => 'X',    
  'NCT' => 'X',    
  'NCC' => 'X',    
  'NCG' => 'X',    
  'ANA' => 'X',    
  'ANT' => 'X',    
  'ANG' => 'X',    
  'ANC' => 'X',    
  'TNA' => 'X',    
  'TNT' => 'X',    
  'TNG' => 'X',    
  'TNC' => 'X',    
  'GNA' => 'X',    
  'GNT' => 'X',    
  'GNG' => 'X',   
  'GNC' => 'X',    
  'CNA' => 'X',  
  'CNT' => 'X',   
  'CNG' => 'X',    
  'CNC' => 'X',    
  'AAN' => 'X',
  'ATN' => 'X', # in 3 of 4 codes for I otherwise M-start
  'AGN' => 'X',
  'ACN' => 'T',
  'TAN' => 'X',
  'TTN' => 'X',
  'TGN' => 'X',
  'TCN' => 'S',
  'GAN' => 'X',
  'GTN' => 'V',
  'GGN' => 'G',
  'GCN' => 'A',
  'CAN' => 'X',
  'CTN' => 'L',
  'CGN' => 'R',
  'CCN' => 'P',
  'ANN' => 'X',	
  'TNN' => 'X',	
  'GNN' => 'X',
  'CNN' => 'X',
  'NAN' => 'X',
  'NTN' => 'X',
  'NGN' => 'X',
  'NCN' => 'X',
  'NNA' => 'X',
  'NNT' => 'X',
  'NNG' => 'X',
  'NNC' => 'X',
  'NNN' => 'X',
);

sub revcom {
  my($dna) = @_;
#  $dna =~ s/[^ATCGN]/N/ig;
  # First reverse the sequence
  my $revcom = reverse($dna);
  # Next, complement the sequence, dealing with upper and lower case
  # A->T, T->A, C->G, G->C
  $revcom =~ tr/ACGTacgt/TGCAtgca/;
  return $revcom;
}
sub translate {
  my ($dna) = @_;
  my $seqlen = length($dna);
  my $n_codon = int($seqlen / 3);
  my $pro = '';
  for(my $i=0; $i < $n_codon ; $i ++) {
    my $codon = uc( substr($dna, $i*3, 3) );
    my $aa = exists $genetic_code{$codon} ? $genetic_code{$codon} : "X";
    $pro .= $aa;
  }
  return $pro;
}
sub translate6 {
  my ($fi, $fo, $sep) = @_;
  my $log = Log::Log4perl->get_logger("Seq");
  $log->info("translating sequences in 6 reading frames");
  $sep ||= "|";
  
  my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
  my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');
  while(my $seqO = $seqHI->next_seq()) {
    my ($id, $seq) = ($seqO->id, $seqO->seq);
    my $len = length($seq);

    for my $i (0..2) {
      my $n_codon = int( ($len - $i) / 3 );
      my $lenC = $n_codon * 3;
      my ($beg, $end) = ($i+1, $i+$lenC);
      my $dna = substr($seq, $beg-1, $lenC);
      my $pro = translate($dna);
      my $sid = join($sep, $id, $beg, $end, "+");
      $seqHO->write_seq( Bio::Seq->new(-id=>$sid, -seq=>$pro) );
    }
    $seq = revcom($seq);
    for my $i (0..2) {
      my $n_codon = int( ($len - $i) / 3 );
      my $lenC = $n_codon * 3;
      my ($beg, $end) = ($i+1, $i+$lenC);
      my $dna = substr($seq, $beg-1, $lenC);
      my $pro = translate($dna);
      my ($begO, $endO) = ($len-$end+1, $len-$beg+1);
      my $sid = join($sep, $id, $begO, $endO, "-");
      $seqHO->write_seq( Bio::Seq->new(-id=>$sid, -seq=>$pro) );
    }
  }
  $seqHI->close();
  $seqHO->close();
}

sub readSeq {
  my ($fi, $opt) = @_;
  $opt ||= 1;
  my $h = {};
  my $seqs = [];
  my $seqH = Bio::SeqIO->new(-file=>$fi);
  while( my $seq = $seqH->next_seq ) {
    $h->{$seq->id} = $seq->seq;
    push @$seqs, $seq;
  }
  if($opt == 1) {
    return $seqs;
  } else {
    die "unknown opt: $opt\n" unless $opt == 2;
    return $h;
  }
}
sub readSeqInfo {
  my ($fi, $opt, $inf) = @_;
  die "$fi is not there\n" unless -s $fi;
  $opt ||= 1;
  $inf ||= "fasta";
  my $ref = {};
  my $seqH = Bio::SeqIO->new(-file => $fi, -format => $inf);
  while(my $seq = $seqH->next_seq()) {
    my $id = $seq->id;
    die "$id appeared >1 times\n" if exists $ref->{$id};
    if($opt == 1) {
      $ref->{$id} = $seq->description;
    } elsif($opt == 2) {
      $ref->{$id} = [$seq->length, $seq->description];
    } else {
      die "unsupported option: $opt\n";
    }
  }
  return $ref;
}
sub read_aln_seq {
  my ($fi) = @_;
  my $ai = Bio::AlignIO->new(-file=>"<$fi");
  my $h;
  while(my $aln = $ai->next_aln()) {
    for my $seq ($aln->each_seq) {
      $h->{$seq->id} = $seq->seq;
    }
  }
  return $h;
}

sub writeSeq {
  my ($seqAry, $fOut, $format) = @_;
  $format ||= "fasta";
  $seqAry = ref($seqAry) eq "ARRAY" ? $seqAry : [$seqAry];
  my $seqH = Bio::SeqIO->new(-file=>">$fOut", -format=>$format);
  for my $seq (@$seqAry) {
    next unless $seq;
    $seqH->write_seq($seq);
  }
  $seqH->close();
}
sub writeSeqInOrder {
  my ($fo, $seqHash, $ids) = rearrange(['out', 'seqs', 'ids'], @_);
  my $seqH = Bio::SeqIO->new(-file => ">$fo", -format => "fasta");
  for my $id (@$ids) {
    my $seq = Bio::Seq->new(-id=>$id, -seq=>$seqHash->{$id});
    $seqH->write_seq($seq);
  }
}

sub seqret_simple {
  my ($db, $chr, $beg, $end, $srd) = @_;
  my $seqStr = $db->seq($chr, $beg, $end);
  die "cannot find seq: $chr:$beg-$end\n" unless $seqStr;
  $seqStr = Bio::Seq->new(-seq=>$seqStr)->revcom->seq if $srd =~ /^\-1?$/;
  return $seqStr;
}
sub seqRet {
  my ($loc, $seqid, $srd, $f_seq) = @_;
  die "no sequence file: $f_seq\n" unless -s $f_seq;
  my $db = Bio::DB::Fasta->new($f_seq);
  
  my @seqStrs;
  $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
  $loc = [ reverse @$loc ] if $srd =~ /^\-1?$/;
  for (@$loc) {
    my ($beg, $end) = @$_;
    my $seqStr = $db->seq($seqid, $beg, $end);
    die "cannot find seq: $seqid:$beg-$end in $f_seq\n" unless $seqStr;
    $seqStr = Bio::Seq->new(-seq=>$seqStr)->revcom->seq if $srd =~ /^\-1?$/;
    push @seqStrs, $seqStr;
  }
  my $seqStr = join("", @seqStrs);
  return $seqStr;
}
sub seqLen {
  my ($seqid, $f_seq) = @_;
  die "no sequence file: $f_seq\n" unless -s $f_seq;
  my $db = Bio::DB::Fasta->new($f_seq);
  return $db->length($seqid);
}

sub seqCompare {
  my ($qry, $tgt) = @_;
  my ($lenQ, $lenT) = (length($qry), length($tgt));
  $lenQ == $lenT || die "unequal seq length: $lenQ <> $lenT\n$qry\n$tgt\n";
  my ($mat, $mis, $qN, $tN) = (0, 0, 0, 0);
  for my $j (0..$lenQ-1) {
    my $chQ = uc(substr($qry, $j, 1));
    my $chT = uc(substr($tgt, $j, 1));
    if($chQ eq "N") {
      $qN ++;
    } elsif($chT eq "N") {
      $tN ++;
    } elsif($chQ eq $chT) {
      $mat ++;
    } else {
      $mis ++;
    }
  }
  return ($mat, $mis, $qN, $tN);
}

sub get_start_codons {
  my ($seq, $end) = @_;
  my $codons = min(60, int($end/3));
  my $begT = $end - $codons*3 + 1;
  
  my $seq2 = substr($seq, $begT-1, $end-$begT+1);
  my $pro = Bio::Seq->new(-seq=>$seq2)->translate->seq;

  if($pro =~ /\*(\w+)$/) {
    $begT = $end - length($1)*3 + 1;
    $pro = substr($pro, $codons-length($1));
  }

  my @begs;
  while($pro =~ /(M)/ig) {
    push @begs, $begT + $-[1]*3;
  }
  return @begs;
}
sub get_stop_codon {
  my ($seq, $beg) = @_;
  my $endT = min(length($seq), $beg + 60*3 - 1);
  my $end = $endT;
  if($endT-$beg+1 >= 2) {
    my $seq2 = substr($seq, $beg-1, $endT-$beg+1);
    my $pro = Bio::Seq->new(-seq=>$seq2)->translate->seq;
    if($pro =~ /^(\w{0,})\*/) {
        $end = $beg + length($1)*3 + 2;
    }
  }
  return $end;
}
sub getLongestOrf {
  my ($loc, $seqid, $srd, $f_ref) = @_;
  $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
  my ($beg, $end) = ($loc->[0]->[0], $loc->[-1]->[1]);
  $loc = [ reverse @$loc ] if $srd eq "-";
  
  my $seqStr1D = seqRet($loc, $seqid, $srd, $f_ref);
  my $seqStr1P = Bio::Seq->new(-seq=>$seqStr1D)->translate->seq;
  my ($aaL, $aaR) = (100, 100);
  my $aaLM = int( ($beg - 1) / 3 );
  my $seqLen = seqLen($seqid, $f_ref);
  my $aaRM = int( ($seqLen - $end) / 3 );
  $aaL = min($aaL, $aaLM);
  $aaR = min($aaR, $aaRM);
  my ($naL, $naR) = (3*$aaL, 3*$aaR);
  my $loc2 = $seqid.":".($beg-$naL)."..".($end+$naR)."[$srd]";
  my $seqStr2D = seqRet($loc2, $seqid, $srd, $f_ref);
  my $seqStr2P = Bio::Seq->new(-seq=>$seqStr2D)->translate->seq;
  if( $seqStr2P =~ /([\*M])([A-LN-Za-z]*)\Q$seqStr1P\E/i ) {
    my $aaU = length($2);
    $aaL = $aaU if $srd eq "+";
    $aaR = $aaU if $srd eq "-";
    if($1 eq "M") {
      $aaL ++;
      $aaR ++;
    }
  }
  if( $seqStr2P =~ /\Q$seqStr1P\E([A-Za-z]*)\*/ ) {
    my $aaD = length($1);
    $aaL = $aaD + 1 if $srd eq "-";
    $aaR = $aaD + 1 if $srd eq "+";
  }
  ($naL, $naR) = (3*$aaL, 3*$aaR);
  my $loc3 = $seqid.":".($beg-$naL)."..".($end+$naR)."[$srd]";
  my $seqStr3D = seqRet($loc3, $seqid, $srd, $f_ref);
  my $lenDna = $end - $beg + 1;
  $seqStr3D .= "N" if $lenDna % 3 == 2;
  my $lenPro = length($seqStr3D) / 3;
  die "protein length [$lenPro]" unless $lenPro =~ /^\d+$/;
  my ($startD, $endD, $startA, $endA);
  ($startA, $endA) = ($aaL+1, $aaL+$lenPro);
  ($startD, $endD) = ($naL+1, $naL+$lenDna);
  if($srd eq "-") {
    ($startA, $endA) = ($aaR+1, $aaR+$lenPro);
    ($startD, $endD) = ($naR+1, $naR+$lenDna);
  }
  my $relLoc = "left[$aaL] right[$aaR]";
  return ($loc3, $seqStr3D, $relLoc);
}


sub checkProtSeq {
  my ($seqPro) = @_;
  my ($codonStart, $codonStop, $preStop, $gap) = (0) x 4;
  $codonStart = 1 if $seqPro =~ /^M/i;
  $codonStop = 1 if $seqPro =~ /\*$/i;
  $preStop = 1 if $seqPro !~ /^[A-Z]+\*?$/i;
  $gap = 1 if $seqPro =~ /X{10,}/i;
  return ($codonStart, $codonStop, $preStop, $gap);
}
sub getSubSeq {
  my ($seq, $loc) = @_;
  my $seq2 = "";
  for (sort {$a->[0] <=> $b->[0]} @$loc) {
    my ($beg, $end) = @$_;
    $seq2 .= substr($seq, $beg-1, $end-$beg+1);
  }
  return $seq2;
}
sub slideSeq {
  my ($seqI, $size, $step) = rearrange(['seq', 'size', 'step'], @_);
  $size ||= 2000;
  $step ||= 1000;
  my $id = $seqI->id;
  my $seqLen = $seqI->length;
  my $i = 0;
  my $wins = get_sliding_windows(1, $seqLen, $step, $size);
  return sub {
    return undef if $i >= @$wins;
    my ($beg, $end) = @{$wins->[$i++]};
    my $seqStr = $seqI->subseq($beg, $end);
    my $id2 = "$id\_$beg\_$end";
    return Bio::Seq->new(-id=>$id2, -seq=>$seqStr);
  }
}

sub translate_smart {
  my ($seq) = @_;
  my $seqCds = $seq->seq;
  my $red = length($seqCds) % 3;
  my $seqPro;
  my $tag = -1;
  my $seqCdsT;
  for my $i (0..2) {
    $seqCdsT = substr($seqCds, $i, length($seqCds)-$red);
    $seqPro = Bio::Seq->new(-seq=>$seqCdsT)->translate->seq;
    if( $seqPro !~ /\*\w+/ ) {
      $tag = $i;
      last;
    }
  }
  if($tag==-1) {
    die "cannot determine from start of end for smart translate: ".$seq->id."\n";
  }
  return $seqPro;
}




1;
__END__

