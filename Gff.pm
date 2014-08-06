package Gff;
use strict;
use Common;
use Seq;
use Rna;
use Gene;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/parse_gff
    gffMerge makeAsblTbl makeAsblGff/;
@EXPORT_OK = qw//;
sub parse_gff {
  my ($fhi) = @_;
  my $header = [qw/chr source type beg end score strand phase tag/];
  my $t = readTable(-inh=>$fhi, -header=>$header);
  my @bins;
  my @idxs = indexes {exists $h_gff->{$_}} $t->col("type");
  for my $i (0..$#idxs) {
    my $beg = $idxs[$i];
    my $end = ($i < $#idxs) ? $idxs[$i+1]-1 : $t->nofRow-1;
    push @bins, [$beg..$end];
  }
  my $i = 0;
  return sub {
    if($i <= $#bins) {
      my $t_gene = $t->subTable($bins[$i++]);
      return Gene->new(-gff=>$t_gene);
    } else { 
      return undef;
    }
  }
}

sub makeAsblTbl {
  my ($fi, $fo) = @_;
  my $fhi = new IO::File $fi, "r";
  my $fho = new IO::File $fo, "w";
  my $phaseDict = {D=>1, U=>1, A=>2, F=>3};
  my $hLen;
  while(<$fhi>) {
    chomp;
    next unless $_;
    my @ps = split "\t";
    my ($chrC, $begC, $endC, $tag, $chrB, $begB, $endB, $strand) = @ps[0..2,4..8];
    $chrC =~ s/^(\d+)$/chr$1/;
    $strand = $strand eq "-" ? -1 : 1;
    my ($type, $score) = ("") x 2;

    $hLen->{$chrC} = 0 unless exists $hLen->{$chrC};
    $hLen->{$chrC} = max($endC, $hLen->{$chrC});
    
    if($tag eq "N") {
      next;
    } else {
      unless( exists $phaseDict->{$tag} ) {
        print "unknown phase $tag => phase 1\n";
        $score = 1;
      } else {
        $score = $phaseDict->{$tag};
      }
      $type = $chrB =~ /^contig/ ? "contig" : "BAC";
    }
    print $fho join("\t", $chrC, $begC, $endC, $strand, $chrB, $begB, $endB, $type, $score, "")."\n";
  }
  for my $chr (keys %$hLen) {
    my $chrLen = $hLen->{$chr};
    print $fho join("\t", $chr, 1, $chrLen, 1, $chr, 1, $chrLen, 'chromosome', '', '')."\n";
  }
}
sub makeAsblGff {
  my ($fi, $fo1, $fo2) = @_;
  my ($cHash, $i) = ({}, 0);
  my $fhi = new IO::File $fi, "r";
  my $fho1 = new IO::File $fo1, "w";
  my $fho2 = new IO::File $fo2, "w";
  print $fho1 "##gff-version 3\n";
  while(<$fhi>) {
    chomp;
    my ($chrC, $begC, $endC, $strand, $chrB, $begB, $endB, $type, $score, $note) = split "\t";
    my $fe = Bio::SeqFeature::Generic->new(-seq_id=>$chrC, -start=>$begC, -end=>$endC, 
      -strand=>$strand, -primary_tag=>$type);
    my $id = sprintf "%05d", ++$i;
    $fe->add_tag_value("ID", $id);
    $fe->add_tag_value("Name", $chrB);
    $fe->score($score) if $score;
    $fe->add_tag_value("Note", "Phase $score") if $score;
    
    print $fho1 join("\n", fe2GffLine($fe))."\n";
    
    $strand = $strand == -1 ? "-" : "+";
    print $fho2 join("\t", $chrC, $begC, $endC, $chrB, $score, $strand)."\n";
  }
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
