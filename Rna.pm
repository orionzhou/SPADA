package Rna;
use strict; 
use Common; 
use Location; 
use Seq;
use Data::Dumper;
use List::Util qw/min max sum/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw//;

sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, $class;
  $self->set_attributes(@args);
  return $self;
}
sub set_attributes {
  my ($self, @args) = @_;
  my ($id, $pa, $beg, $end, $seqid, $strand) = 
    rearrange([qw/id parent beg end seqid strand/], @args);
  defined $id     && $self->id($id);
  defined $pa     && $self->parent($pa);
  defined $seqid  && $self->seqid($seqid);
  defined $beg    && $self->beg($beg);
  defined $end    && $self->end($end);
  defined $strand && $self->strand($strand);

  my ($eloc, $iloc, $cloc, $floc, $tloc, $phase) = 
    rearrange([qw/exon intron cds utr5 utr3 phase/], @args);
  defined $eloc  && $self->exon($eloc);
  defined $iloc  && $self->intron($iloc);
  defined $cloc  && $self->cds($cloc);
  defined $floc  && $self->utr5($floc);
  defined $tloc  && $self->utr3($tloc);
  defined $phase && $self->phase($phase);

  my ($src, $conf, $cat1, $cat2, $cat3, $note) = 
    rearrange([qw/source conf cat1 cat2 cat3 note/], @args);
  defined $src   && $self->source($src);
  defined $conf  && $self->conf($conf);
  defined $cat1  && $self->cat1($cat1);
  defined $cat3  && $self->cat3($cat3);
  defined $cat3  && $self->cat3($cat3);
  defined $note  && $self->note($note);
  
  my ($t_gff, $t_gtb) = rearrange(['gff', 'gtb'], @args);
  defined $t_gff  && $self->from_gff($t_gff);
  defined $t_gtb  && $self->from_gtb($t_gtb);
}
sub id {
  my $obj = shift;
  return $obj->{'_id'} = shift if @_;
  return $obj->{'_id'};
}
sub parent {
  my $obj = shift;
  return $obj->{'_parent'} = shift if @_;
  return $obj->{'_parent'};
}
sub seqid {
  my $obj = shift;
  return $obj->{'_seqid'} = shift if @_;
  return $obj->{'_seqid'};
}
sub beg {
  my $obj = shift;
  return $obj->{'_beg'} = shift if @_;
  return $obj->{'_beg'};
}
sub end {
  my $obj = shift;
  return $obj->{'_end'} = shift if @_;
  return $obj->{'_end'};
}
sub strand {
  my ($obj, $srd) = @_;
  if(defined($srd)) {
    $srd = "+" if $srd =~ /^\+?1$/;
    $srd = "-" if $srd eq "-1";
    $srd =~ /^[\+\-]$/ || die "unknown strand: $srd\n";
    $obj->{'_strand'} = $srd;
  } else {
    return $obj->{'_strand'};
  }
}

sub exon {
  my $obj = shift;
  return $obj->{'_exon'} = shift if @_;
  return $obj->{'_exon'};
}
sub intron {
  my $obj = shift;
  return $obj->{'_intron'} = shift if @_;
  return $obj->{'_intron'};
}
sub cds {
  my $obj = shift;
  return $obj->{'_cds'} = shift if @_;
  return $obj->{'_cds'};
}
sub utr5 {
  my $obj = shift;
  return $obj->{'_utr5'} = shift if @_;
  return $obj->{'_utr5'};
}
sub utr3 {
  my $obj = shift;
  return $obj->{'_utr3'} = shift if @_;
  return $obj->{'_utr3'};
}
sub phase {
  my $obj = shift;
  return $obj->{'_phase'} = shift if @_;
  return $obj->{'_phase'};
}

sub source {
  my $obj = shift;
  return $obj->{'_source'} = shift if @_;
  return $obj->{'_source'};
}
sub conf {
  my $obj = shift;
  return $obj->{'_conf'} = shift if @_;
  return $obj->{'_conf'};
}
sub cat1 {
  my $obj = shift;
  return $obj->{'_cat1'} = shift if @_;
  return $obj->{'_cat1'};
}
sub cat2 {
  my $obj = shift;
  return $obj->{'_cat2'} = shift if @_;
  return $obj->{'_cat2'};
}
sub cat3 {
  my $obj = shift;
  return $obj->{'_cat3'} = shift if @_;
  return $obj->{'_cat3'};
}
sub note {
  my $obj = shift;
  return $obj->{'_note'} = shift if @_;
  return $obj->{'_note'};
}

sub check_mRNA { # infer exon/intron/utr/cds
  my $self = shift;

  my $srd = $self->strand;
  my $locR = [ [1, $self->end-$self->beg+1] ];
  my ($locE, $locI, $locC, $loc5, $loc3) = 
    map {$self->$_} qw/exon intron cds utr5 utr3/;
  $locC = [ @$locE ] if @$locC == 0;
  
  my $cdsB = min( map {$_->[0]} @$locC );
  my $cdsE = max( map {$_->[1]} @$locC );
  if( @$locE && @$locC ) { #infer utr5/3 from exon+cds
    ($loc5, $loc3) = ([], []);
    my ($locU, $lenU) = posDiff($locE, $locC);
    for (@$locU) {
      my ($b, $e) = @$_;
      if( $e < $cdsB ) {
        push @$loc5, [$b, $e];
      } elsif( $b > $cdsE ) {
        push @$loc3, [$b, $e];
      } else {
        die "\nUTR [ $b - $e ] not outside CDS [ $cdsB - $cdsE ]\n";
      }
    }
  } elsif( !@$locE ) {
    my @locC2 = map {[$_->[0], $_->[1]]} @$locC;
    $loc5 = [ grep { $_->[1] < $cdsB } @$loc5 ];
    $loc3 = [ grep { $_->[0] > $cdsE } @$loc3 ];
    $locE = posMerge([@locC2, @$loc5, @$loc3]);
    $locE = [ map {[$_->[0], $_->[1]]} @$locE ];
  }
  ($locI) = posDiff($locR, $locE);
  
  $self->exon($locE);
  $self->intron($locI);
  $self->cds($locC);
  $self->utr5($loc5);
  $self->utr3($loc3);
  defined $self->phase || $self->phase(getPhase($locC));
}
sub check_phase { # check & fix phases
  my $self = shift;
  my $fs = shift;
  
  my $rloc = $self->cds;
  my $loc = $self->strand eq "-" ? 
    [ map {[$self->end-$_->[1]+1, $self->end-$_->[0]+1]} @$rloc ] :
    [ map {[$self->beg+$_->[0]-1, $self->beg+$_->[1]-1]} @$rloc ]; 
  my $seqStr = seqRet($loc, $self->seqid, $self->strand, $fs);
  my $seq = Bio::Seq->new(-id=>"test", -seq=>$seqStr);
  my $frame = -1;
  for my $i (0..2) {
    my $prot = $seq->translate(-frame => $i)->seq;
    if($prot =~ /^[A-Z]+\*?$/i) {
      $frame = $i;
      last;
    }
  }
  
  $frame >= 0 || die $self->id.": no proper frame[$frame]\n$seqStr\n";

  my $flag = 0;
  if($frame > 0) {
    $flag = 1;
    ($rloc) = cropLoc($rloc, $frame, 1);
  }
  my $res = locAryLen($rloc) % 3;
  if($res > 0) {
    ($rloc) = cropLoc($rloc, $res, 2);
    $rloc = [ sort {$a->[0] <=> $b->[0]} @$rloc ];
  }
#  die Dumper($rloc) if $self->id eq "Medtr1047te0010.1";
  $self->cds($rloc);
  $self->phase(getPhase($rloc));
  $self->check_mRNA();
  return $flag;
}

sub from_gff {
  my ($self, $rows) = @_;

  $self->id($rows->[0]->[0]);
  $self->parent($rows->[0]->[1]);
  $self->cat1($rows->[0]->[2]);
  $self->seqid($rows->[0]->[3]);
  $self->source($rows->[0]->[4]);
  $self->beg($rows->[0]->[5]);
  $self->end($rows->[0]->[6]);
  $self->strand($rows->[0]->[8]);
  $self->note($rows->[0]->[10]);
#  die join("\n", map {join("\t", @$_)} @$rows)."\n";
  
  my ($locE, $locI, $locC, $loc5, $loc3) = ([], [], [], [], []);
  for my $i (1..@$rows-1) {
    my ($id, $par, $type, $seqid, $src, $beg, $end, $score, $srd, $phase, $note) = @{$rows->[$i]};
    my $errstr = join("\t", $par, $type, $seqid, $srd, $beg, $end)."\n";
    $seqid eq $self->seqid || die "seqid error: [".$self->seqid."] $errstr\n";
    $srd eq $self->strand || die "srd error: [".$self->strand."] $errstr\n";
    $beg >= $self->beg || die "out of boundary beg[".$self->beg."]: $errstr\n";
    $end <= $self->end || die "out of boundary end[".$self->end."]: $errstr\n";
    my ($rb, $re) = $self->strand eq "+" ? 
      ($beg - $self->beg + 1, $end - $self->beg + 1) :
      ($self->end - $end + 1, $self->end - $beg + 1);
    if($type eq 'exon') {
      push @$locE, [$rb, $re];
    } elsif($type eq 'intron') {
      push @$locI, [$rb, $re];
    } elsif($type eq 'CDS') {
      push @$locC, [$rb, $re, $phase];
    } elsif($type eq 'five_prime_UTR') {
      push @$loc5, [$rb, $re];
    } elsif($type eq 'three_prime_UTR') {
      push @$loc3, [$rb, $re];
    } else {
      die "unsupported type: $type for Rna\n";
    }
  }
  $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
  $self->exon($locE);
  $self->intron($locI);
  $self->cds($locC);
  $self->utr5($loc5);
  $self->utr3($loc3);
  my @phases = map {$_->[2]} @$locC;
  @phases = grep {defined($_) && $_=~/^[012]$/} @phases;
  $self->phase(\@phases) if @phases > 0 && @phases == @$locC;
  
  $self->check_mRNA() if $self->cat1 eq "mRNA";
}
sub to_gff {
  my $self = shift;
  my $flag_utr = shift;
  $flag_utr = 1 unless defined($flag_utr);
  my @rows;
  
  my $seqid = $self->seqid;
  my $srd = $self->strand;
  my @tags = ("ID=".$self->id, "Parent=".$self->parent);
  push @tags, "Note=".$self->note if $self->note;
  push @rows, [$seqid, $self->source, $self->cat1, $self->beg, $self->end, '.', $srd, '.', join(";", @tags)];
  
  my @types = qw/exon CDS five_prime_UTR three_prime_UTR/;
  my @locs = ($self->exon, $self->cds, $self->utr5, $self->utr3);
  for my $i (0..$#types) {
    my $type = $types[$i];
    next if $flag_utr != 1 && $type =~ /prime\_UTR$/;
    my $loc = $locs[$i];
    for my $j (0..@$loc-1) {
      my ($rb, $re) = @{$loc->[$j]};
      my $phase = $type eq 'CDS' ? $self->phase->[$j] : '.';
      my ($beg, $end) = $srd eq "-" ? 
        ($self->end-$re+1, $self->end-$rb+1) : 
        ($self->beg+$rb-1, $self->beg+$re-1);
      push @rows, [$seqid, '.', $type, $beg, $end, '.', $srd, $phase, 'Parent='.$self->id];
    }
  }
  return join("\n", map {join("\t", @$_)} @rows);
}


sub from_gtb {
  my ($self, $row) = @_;
  my ($id, $pa, $seqid, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, 
    $phase, $src, $conf, $cat1, $cat2, $cat3, $note) = @$row;
  my $locE = locStr2Ary($locES);
  my $locI = locStr2Ary($locIS);
  my $locC = locStr2Ary($locCS);
  my $loc5 = locStr2Ary($loc5S);
  my $loc3 = locStr2Ary($loc3S);
  $self->id($id);
  $self->parent($pa);
  $self->seqid($seqid);
  $self->beg($beg);
  $self->end($end);
  $self->strand($srd);
  $self->source($src);
  $self->conf($conf);
  $self->cat1($cat1);
  $self->cat2($cat2);
  $self->cat3($cat3);
  $self->note($note);

  my @phases = split(",", $phase);
  if(@phases >= 1 && @$locC == @phases) {
    $self->phase(\@phases);
  }

  $self->exon($locE);
  $self->intron($locI);
  $self->cds($locC);
  $self->utr5($loc5);
  $self->utr3($loc3);
  $self->check_mRNA() if $cat1 eq "mRNA";
}
sub to_gtb {
  my ($self) = @_;
  my ($locE, $locI, $locC, $loc5, $loc3) = 
    map {locAry2Str($self->$_)} qw/exon intron cds utr5 utr3/;
  my $phase = defined $self->phase ? join(",", @{$self->phase}) : ".";
  my $src  = defined $self->source ? $self->source : "";
  my $conf = defined $self->conf ? $self->conf : "";
  my $cat1 = defined $self->cat1 ? $self->cat1 : "";
  my $cat2 = defined $self->cat2 ? $self->cat2 : "";
  my $cat3 = defined $self->cat3 ? $self->cat3 : "";
  my $note = defined $self->note ? $self->note : "";
  my $row = [$self->id, $self->parent, $self->seqid, $self->beg, $self->end, $self->strand,
    $locE, $locI, $locC, $loc5, $loc3, $phase,
    $src, $conf, $cat1, $cat2, $cat3, $note];
  return join("\t", @$row);
}


1;
__END__
