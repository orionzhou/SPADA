package Gene;
use strict; 
use Data::Dumper;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use Common; 
use Rna; 
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/$h_gff/;

sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, $class;
  $self->set_attributes(@args);
  return $self;
}
sub set_attributes {
  my ($self, @args) = @_;
  my ($id, $beg, $end, $seqid, $strand, $type, $source, $score, $note) = 
    rearrange([qw/id beg end seqid strand type source score note/], @args);
  defined $id     && $self->id($id);
  defined $beg    && $self->beg($beg);
  defined $end    && $self->end($end);
  defined $seqid  && $self->seqid($seqid);
  defined $strand && $self->strand($strand);
  defined $type   && $self->type($type);
  defined $source && $self->source($source);
  defined $score  && $self->score($score);
  defined $note   && $self->note($note);
  my ($t_gff, $t_gtb) = rearrange(['gff', 'gtb'], @args);
  defined $t_gff  && $self->from_gff($t_gff);
  defined $t_gtb  && $self->from_gtb($t_gtb);
}
sub id {
  my $obj = shift;
  return $obj->{'_id'} = shift if @_;
  return $obj->{'_id'} || '';
}
sub beg {
  my $obj = shift;
  return $obj->{'_beg'} = shift if @_;
  return $obj->{'_beg'} || '';
}
sub end {
  my $obj = shift;
  return $obj->{'_end'} = shift if @_;
  return $obj->{'_end'} || '';
}
sub seqid {
  my $obj = shift;
  return $obj->{'_seqid'} = shift if @_;
  return $obj->{'_seqid'} || '';
}
sub strand {
  my $obj = shift;
  return $obj->{'_strand'} = shift if @_;
  return $obj->{'_strand'} || '';
}
sub type {
  my $obj = shift;
  return $obj->{'_type'} = shift if @_;
  return $obj->{'_type'} || '';
}
sub source {
  my $obj = shift;
  return $obj->{'_source'} = shift if @_;
  return $obj->{'_source'} || '.';
}
sub note {
  my $obj = shift;
  return $obj->{'_note'} = shift if @_;
  return $obj->{'_note'} || '';
}



1;
__END__
