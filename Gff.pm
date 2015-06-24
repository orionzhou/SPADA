package Gff;
use strict;
use Common;
use Seq;
use Rna;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/parse_gene_gff
    gffMerge/;
@EXPORT_OK = qw//;
our $h_gff = {
  'mRNA'   => [qw/exon intron CDS five_prime_UTR three_prime_UTR/],
  'rRNA'   => ['exon'],
  'tRNA'   => ['exon'],
  'miRNA'  => ['exon'],
  'ncRNA'  => ['exon'],
  'snRNA'  => ['exon'],
  'snoRNA' => ['exon'],
};

sub parse_gene_gff {
  my ($fhi) = @_;
  my $header = [qw/chr source type beg end score strand phase tag/];
  my $t = readTable(-inh => $fhi, -header => $header);
  my @bins;
  my @idxs = indexes {$_ =~ /gene/i} $t->col("type");
  for my $i (0..$#idxs) {
    my $beg = $idxs[$i];
    my $end = ($i < $#idxs) ? $idxs[$i+1]-1 : $t->nofRow-1;
    push @bins, [$beg, $end];
  }
  my $i = 0;
  return sub {
    if($i <= $#bins) {
      my ($ib, $ie) = @{$bins[$i]};
      my @idxs_r = grep {exists $h_gff->{$t->elm($_, "type")}} ($ib..$ie);
      my @rnas;
      for my $j (0..$#idxs_r) {
        my $jb = $idxs_r[$j];
        my $je = ($j < $#idxs_r) ? $idxs_r[$j+1]-1 : $ie;
        my @stats = map {parse_gff_line($t->row($_))} ($jb..$je);
        push @rnas, Rna->new(-gff => \@stats);
      }
      $i ++;
      return [reverse @rnas];
    } else { 
      return undef;
    }
  }
}
sub parse_gff_line {
  my ($seqid, $src, $type, $beg, $end, $score, $srd, $phase, $tag) = @_;
  my $ht = parse_gff_tags($tag);
  my $id = exists $ht->{'ID'} ? $ht->{'ID'} : '';
  my $par = exists $ht->{'Parent'} ? $ht->{'Parent'} : '';
  my $note = exists $ht->{'Note'} ? $ht->{'Note'} : '';
  return [$id, $par, $type, $seqid, $src, $beg, $end, $score, $srd, $phase, $note];
}

sub gffMerge {
  my ($fis, $fo) = rearrange(['in', 'out'], @_);
  my $fho = new IO::File $fo, "w";
  my (@fSeq, @fGff);
  for my $fi (@$fis) {
    if($fi =~ /\.(fa|fasta|fas)$/i) {
      push @fSeq, $fi;
    } elsif($fi =~ /\.(gff|gff3|gff2)/i) {
      push @fGff, $fi;
    } else {
      die("unknown format: $fi\n");
    }
  }
  print $fho "##gff-version 3\n";
  for my $fi (@fGff) {
    my $fhi = new IO::File $fi, "r";
    while( <$fhi> ) {
      chomp;
      next if /^\#/;
      print $fho $_."\n";
    }
  }
  print $fho "##FASTA\n" if @fSeq > 0;
  for my $fi (@fSeq) {
    my $seqIH = Bio::SeqIO->new(-file=>$fi, -format=>'fasta');
    my $seqOH = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
    while( my $seq = $seqIH->next_seq ) {
      $seqOH->write_seq($seq);
    }
  }
}


1;
__END__
