package Rna;
use strict; 
use Common; 
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
    my ($id, $pa, $beg, $end, $seqid, $strand, $type, $source, $note) = 
        rearrange([qw/id parent beg end seqid strand type source note/], @args);
    defined $id     && $self->id($id);
    defined $pa     && $self->parent($pa);
    defined $beg    && $self->beg($beg);
    defined $end    && $self->end($end);
    defined $seqid  && $self->seqid($seqid);
    defined $strand && $self->strand($strand);
    defined $type   && $self->type($type);
    defined $source && $self->source($source);
    defined $note   && $self->note($note);

    my ($exon, $intron, $cds, $utr5, $utr3) = rearrange([qw/exon intron cds utr5 utr3/], @args);
    defined $exon   && $self->exon($exon);
    defined $intron && $self->intron($intron);
    defined $cds    && $self->cds($cds);
    defined $utr5   && $self->utr5($utr5);
    defined $utr3   && $self->utr3($utr3);

    my ($conf, $cat1, $cat2, $cat3) = rearrange([qw/conf cat1 cat2 cat3/], @args);
    defined $conf   && $self->conf($conf);
    defined $cat1   && $self->cat1($cat1);
    defined $cat3   && $self->cat3($cat3);
    defined $cat3   && $self->cat3($cat3);
    
    my ($t_gff, $t_gtb) = rearrange(['gff', 'gtb'], @args);
    defined $t_gff  && $self->from_gff($t_gff);
    defined $t_gtb  && $self->from_gtb($t_gtb);
}
sub id {
    my $obj = shift;
    return $obj->{'_id'} = shift if @_;
    return $obj->{'_id'} || '';
}
sub parent {
    my $obj = shift;
    return $obj->{'_parent'} = shift if @_;
    return $obj->{'_parent'} || '';
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
    my ($obj, $strand) = @_;
    if(defined($strand)) {
        $strand = "+" if $strand =~ /^\+?1$/;
        $strand = "-" if $strand eq '-1';
        die "unknown strand: $strand\n" unless $strand =~ /^[\+\-]$/;
        $obj->{'_strand'} = $strand;
    } else {
        return $obj->{'_strand'} || '.';
    }
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

sub exon {
    my $obj = shift;
    return $obj->{'_exon'} = shift if @_;
    return $obj->{'_exon'} || '';
}
sub intron {
    my $obj = shift;
    return $obj->{'_intron'} = shift if @_;
    return $obj->{'_intron'} || '';
}
sub cds {
    my $obj = shift;
    return $obj->{'_cds'} = shift if @_;
    return $obj->{'_cds'} || '';
}
sub utr5 {
    my $obj = shift;
    return $obj->{'_utr5'} = shift if @_;
    return $obj->{'_utr5'} || '';
}
sub utr3 {
    my $obj = shift;
    return $obj->{'_utr3'} = shift if @_;
    return $obj->{'_utr3'} || '';
}

sub conf {
    my $obj = shift;
    return $obj->{'_conf'} = shift if @_;
    return $obj->{'_conf'} || '';
}
sub cat1 {
    my $obj = shift;
    return $obj->{'_cat1'} = shift if @_;
    return $obj->{'_cat1'} || '';
}
sub cat2 {
    my $obj = shift;
    return $obj->{'_cat2'} = shift if @_;
    return $obj->{'_cat2'} || '';
}
sub cat3 {
    my $obj = shift;
    return $obj->{'_cat3'} = shift if @_;
    return $obj->{'_cat3'} || '';
}

sub check_mRNA { # infer exon/intron/utr/cds, sort cds location, check phase and length
    my $self = shift;

    my $strand = $self->strand;
    my $locR = [ [$self->beg, $self->end] ];
    my ($locE, $locI, $locC, $loc5, $loc3) = map {$self->$_} qw/exon intron cds utr5 utr3/;
    $locC = [ @$locE ] if @$locC == 0;
    $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
    $locC = [ reverse @$locC ] if $strand eq "-";

    if(defined($locC->[0]->[2]) && $locC->[0]->[2] != 0) {
        my $opt_crop = $strand eq "-" ? 2 : 1;
        my $loc5x;
        ($locC, $loc5x) = cropLoc($locC, $locC->[0]->[2], $opt_crop);
        push @$loc5, @$loc5x;
    }
    if(locAryLen($locC) % 3 != 0) {
        my $res = locAryLen($locC) % 3;
        my $opt_crop = $strand eq "-" ? 1 : 2;
        my $loc3x;
        ($locC, $loc3x) = cropLoc($locC, $res, $opt_crop);
        push @$loc3, @$loc3x;
    }

    $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
    $locC = [ reverse @$locC ] if $strand eq "-";
    my $len = 0;
    for (@$locC) {
        my $lenT = $_->[1] - $_->[0] + 1;
        $_->[2] = ( 3 - $len % 3 ) % 3;
        $len += $lenT;
    }

    my $cdsB = min( map {$_->[0]} @$locC );
    my $cdsE = max( map {$_->[1]} @$locC );
    if( @$locE && @$locC ) { #infer utr5/3 from exon+cds
        ($loc5, $loc3) = ([], []);
        my ($locU, $lenU) = posDiff($locE, $locC);
        for (@$locU) {
            my ($b, $e) = @$_;
            if( ($strand eq "+" && $e<$cdsB) || ($strand eq "-" && $b>$cdsE) ) {
                push @$loc5, [$b, $e];
            } elsif( ($strand eq "+" && $b>$cdsE) || ($strand eq "-" && $e<$cdsB) ) {
                push @$loc3, [$b, $e];
            } else {
                die "\nUTR [ $b - $e ] not outside CDS [ $cdsB - $cdsE ]\n";
            }
        }
    } elsif( !@$locE ) {
        my @locC2 = map {[$_->[0], $_->[1]]} @$locC;
        $loc5 = [ grep { ($strand eq "+" && $_->[1] < $cdsB) || ($strand eq "-" && $_->[0] > $cdsE) } @$loc5 ];
        $loc3 = [ grep { ($strand eq "+" && $_->[0] > $cdsE) || ($strand eq "-" && $_->[1] < $cdsB) } @$loc3 ];
        $locE = posMerge([@locC2, @$loc5, @$loc3]);
        $locE = [ map {[$_->[0], $_->[1]]} @$locE ];
    }
    ($locI) = posDiff($locR, $locE);

    $self->exon($locE);
    $self->intron($locI);
    $self->cds($locC);
    $self->utr5($loc5);
    $self->utr3($loc3);
}

sub from_gff {
    my ($self, $rows) = @_;

    $self->id($rows->[0]->[0]);
    $self->parent($rows->[0]->[1]);
    $self->type($rows->[0]->[2]);
    $self->seqid($rows->[0]->[3]);
    $self->source($rows->[0]->[4]);
    $self->beg($rows->[0]->[5]);
    $self->end($rows->[0]->[6]);
    $self->strand($rows->[0]->[8]);
    $self->note($rows->[0]->[10]);
    $self->cat2($self->type);
#  die join("\n", map {join("\t", @$_)} @$rows)."\n";

    my $strand = $self->strand;
    my ($locR, $locE, $locI, $locC, $loc5, $loc3) = ([[$self->beg, $self->end]], [], [], [], [], []);
    for my $i (1..@$rows-1) {
        my ($id, $pa, $type, $seqid, $source, $beg, $end, $score, $strand, $phase, $note) = @{$rows->[$i]};
        die "seqid conflict: [".$self->seqid."] ".join("\t", $pa, $type, $seqid, $strand, $beg, $end)."\n" if $seqid ne $self->seqid;
        die "strand conflict: [".$self->strand."] ".join("\t", $pa, $type, $seqid, $strand, $beg, $end)."\n" if $strand ne "." && $strand ne $self->strand;
        die "out of boundary beg[".$self->beg."]: ".join("\t", $pa, $type, $seqid, $strand, $beg, $end)."\n" if $beg < $self->beg;
        die "out of boundary end[".$self->end."]: ".join("\t", $pa, $type, $seqid, $strand, $beg, $end)."\n" if $end > $self->end;
        if($type eq 'exon') {
            push @$locE, [$beg, $end];
        } elsif($type eq 'intron') {
            push @$locI, [$beg, $end];
        } elsif($type eq 'CDS') {
            push @$locC, $phase=~/^[012]$/ ? [$beg, $end, $phase] : [$beg, $end];
        } elsif($type eq 'five_prime_UTR') {
            push @$loc5, [$beg, $end];
        } elsif($type eq 'three_prime_UTR') {
            push @$loc3, [$beg, $end];
        } else {
            die "unsupported type: $type for Rna\n";
        }
    }
    $self->exon($locE);
    $self->intron($locI);
    $self->cds($locC);
    $self->utr5($loc5);
    $self->utr3($loc3);
    $self->check_mRNA() if $self->type eq "mRNA";
}
sub to_gff {
    my $self = shift;
    my $flag_utr = shift;
    $flag_utr = 1 unless defined($flag_utr);
    my @rows;
    
    my $seqid = $self->seqid;
    my $strand = $self->strand;
    my @tags = ("ID=".$self->id, "Parent=".$self->parent);
    push @tags, "Note=".$self->note if $self->note;
    push @rows, [$seqid, $self->source, $self->type, $self->beg, $self->end, '.', $strand, '.', join(";", @tags)];

    for (@{$self->exon}) { push @rows, [$seqid, '.', 'exon', $_->[0], $_->[1], '.', $strand, '.', 'Parent='.$self->id]; }
    for (@{$self->cds}) { push @rows, [$seqid, '.', 'CDS', $_->[0], $_->[1], '.', $strand, $_->[2], 'Parent='.$self->id]; }
    if($flag_utr eq 1) {
        for (@{$self->utr5}) { push @rows, [$seqid, '.', 'five_prime_UTR', $_->[0], $_->[1], '.', $strand, '.', 'Parent='.$self->id]; }
        for (@{$self->utr3}) { push @rows, [$seqid, '.', 'three_prime_UTR', $_->[0], $_->[1], '.', $strand, '.', 'Parent='.$self->id]; }
    }
    return join("\n", map {join("\t", @$_)} @rows);
}


sub from_gtb {
    my ($self, $row) = @_;
    my ($id, $pa, $seqid, $beg, $end, $strand, $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, $source, $conf, $cat1, $cat2, $cat3, $note) = @$row;
    $self->id($id);
    $self->parent($pa);
    $self->type($cat2);
    $self->seqid($seqid);
    $self->source($source);
    $self->beg($beg);
    $self->end($end);
    $self->strand($strand);
    $self->note($note);
    $self->cat1($cat1);
    $self->cat2($self->type);
    my $locE = locStr2Ary($locES);
    my $locI = locStr2Ary($locIS);
    my $locC = locStr2Ary($locCS);
    my $loc5 = locStr2Ary($loc5S);
    my $loc3 = locStr2Ary($loc3S);

    $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
    $locC = [ reverse @$locC ] if $strand eq "-";
    my @phases = split(",", $phase);
    if(@phases >= 1 && @$locC == @phases) {
        for my $i (0..$#phases) { $locC->[$i]->[2] = $phases[$i]; }
    }

    $self->exon($locE);
    $self->intron($locI);
    $self->cds($locC);
    $self->utr5($loc5);
    $self->utr3($loc3);
    $self->check_mRNA() if $self->type eq "mRNA";
}
sub to_gtb {
    my ($self) = @_;
    my ($locE, $locI, $locC, $loc5, $loc3) = map {locAry2Str($self->$_)} qw/exon intron cds utr5 utr3/;
    my $phase = join(",", map {$_->[2]} @{$self->cds});
    my $row = [$self->id, $self->parent, $self->seqid, $self->beg, $self->end, $self->strand,
        $locE, $locI, $locC, $loc5, $loc3, $phase, $self->source,
        $self->conf, $self->cat1, $self->cat2, $self->cat3, $self->note];
    return join("\t", @$row);
}


1;
__END__
