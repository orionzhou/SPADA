package Gene;
use strict; 
use Common; 
use Data::Dumper;
use Rna; 
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
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

our $h_gff = {
    'gene' => {
        'mRNA'   => ['exon', 'intron', 'CDS', 'five_prime_UTR', 'three_prime_UTR'],
        'rRNA'   => ['exon'],
        'tRNA'   => ['exon'],
        'miRNA'  => ['exon'],
        'ncRNA'  => ['exon'],
        'snRNA'  => ['exon'],
        'snoRNA' => ['exon'],
    },
    'transposable_element_gene' => {
        'mRNA'   => ['exon', 'intron', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
    },
    'pseudogene' => {
        'mRNA'   => ['exon', 'intron', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
    }
};

sub from_gff {
    my ($self, $t) = @_;
    
    my @stats;
    for my $i (0..$t->nofRow-1) {
        my ($seqid, $source, $type, $beg, $end, $score, $strand, $phase, $tag) = $t->row($i);
        my $ht = parse_gff_tags($tag);
        my $id = exists $ht->{'ID'} ? $ht->{'ID'} : '';
        my $pa = exists $ht->{'Parent'} ? $ht->{'Parent'} : '';
        my $note = exists $ht->{'Note'} ? $ht->{'Note'} : '';
        push @stats, [$id, $pa, $type, $seqid, $source, $beg, $end, $score, $strand, $phase, $note];
    }
    my $statG = $stats[0];
    $self->id($statG->[0]);
    $self->type($statG->[2]);
    $self->seqid($statG->[3]);
    $self->source($statG->[4]);
    $self->beg($statG->[5]);
    $self->end($statG->[6]);
    $self->strand($statG->[8]);
    $self->note($statG->[10]);
        
    my ($idG, $typeG) = @{$stats[0]}[0,2];
    die "type1[$typeG] not supported\n" unless exists $h_gff->{$typeG};
    my @idxs_typeR;
    my $hr = $h_gff->{$typeG};
    for my $i (1..$#stats) {
        my ($id, $pa, $typeR) = @{$stats[$i]};
        if(exists($hr->{$typeR})) {
            die "Parent of $id is $pa [ not $idG ]\n" unless $pa eq $idG;
            push @idxs_typeR, $i;
        }
    }
    
    my $h_idx;
    for my $i (0..$#idxs_typeR) {
        my $idx_typeR = $idxs_typeR[$i];
        my ($idR, $typeR) = @{$stats[$idx_typeR]}[0,2];
        my @typesC = @{$hr->{$typeR}};
      
        my @idxs_typeC = indexes {$_->[1] eq $idR} @stats;
#    my $idx_typeC_beg = $idxs_typeR[$i] + 1;
#    my $idx_typeC_end = ($i == $#idxs_typeR) ? $#stats : $idxs_typeR[$i+1]-1;
#    for my $j ($idx_typeC_beg..$idx_typeC_end) {
        for my $j (@idxs_typeC) {
            my ($id, $pa, $typeC) = @{$stats[$j]};
            my $idx_tmp = first_index {$_ eq $typeC} @typesC;
            $h_idx->{$idx_typeR} ||= [];
            push @{$h_idx->{$idx_typeR}}, $j if $idx_tmp > -1;
        }
    }

    for my $idx (sort(keys(%$h_idx))) {
        my @idxs = @{$h_idx->{$idx}};
        $self->add_rna( Rna->new(-gff=>[@stats[$idx, @idxs]], -cat1=>$self->type) );
    }
}
sub to_gff {
    my $self = shift;
    my $seqid = $self->seqid;
    my $strand = $self->strand;
    my @tags = ("ID=".$self->id);
    push @tags, "Note=".$self->note if $self->note;
    return join("\t", $seqid, $self->source, $self->type, $self->beg, $self->end, '.', $strand, '.', join(";", @tags));
}

sub from_gtb {
    my ($self, $t) = @_;
    for my $i (0..$t->nofRow-1) {
        $self->add_rna( Rna->new(-gtb=>$t->rowRef($i)) );
    }
}

sub add_rna {
    my ($self, $rna) = @_;
    $self->{'_rna'} ||= [];
    push @{$self->{'_rna'}}, $rna;
    
    $self->id     || $self->id($rna->parent);
    die "gene id conflict: ".$self->id." <> ".$rna->parent."\n" unless $self->id eq $rna->parent;
    $self->type   || $self->type($rna->cat1);
    die "gene [".$self->id."] type conflict: ".$self->type." <> ".$rna->cat1."\n" unless $self->type eq $rna->cat1;
    $self->seqid  || $self->seqid($rna->seqid);
    die "seqid conflict: ".$self->seqid." <> ".$rna->seqid."\n" unless $self->seqid eq $self->seqid;
    $self->strand || $self->strand($rna->strand);
    die "strand conflict: ".$self->strand." <> ".$rna->strand."\n" unless $self->strand eq $self->strand;
    $self->beg    || $self->beg($rna->beg);
    $self->end    || $self->end($rna->end);
    $self->beg > $rna->beg && $self->beg($rna->beg);
    $self->end < $rna->end && $self->end($rna->end);
}
sub get_rna { 
    my $self = shift;
    die Dumper($self) unless $self->{'_rna'};
    return @{$self->{'_rna'}};
}



1;
__END__
