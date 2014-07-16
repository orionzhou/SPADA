package Gtb;
use strict;
use Common;
use Location; 
use Seq;
use Gene;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/@HEAD_GTB @HEAD_GTBX
    readGtb gtbSum gtb2gff gtb2bed gtb2gtbx
    get_ovlp_gtb/;
@EXPORT_OK = qw//;

our @HEAD_GTB  = qw/id par chr beg end srd locE locI locC loc5 loc3 phase src conf cat1 cat2 cat3 note/;
our @HEAD_GTBX = qw/id par chr beg end srd locE locI locC loc5 loc3 phase src conf cat1 cat2 cat3 note seq/;
  my ($id, $par, $chr, $beg, $end, $srd, 
    $locES, $locIS, $locCS, $loc5S, $loc3S, $phase, 
    $src, $conf, $cat1, $cat2, $cat3, $note) = ();
sub readGtb {
    my ($fi, $opt) = rearrange(['in', 'opt'], @_);
    $opt ||= 1;
    die "fi doesn't exist: $fi\n" unless -s $fi;
    my $t = readTable(-in=>$fi, -header=>1);
    my $rst;
    if($opt == 1) {
        $rst = $t;
    } elsif($opt == 2) {
        $rst = { map {$t->elm($_, "id") => $t->rowRef($_)} (0..$t->nofRow-1) };
    } else {
        die "unsupported option: $opt\n";
    }
    return $rst;
}
sub gtbSum {
    my ($fhi, $opt) = @_;
    $opt = 1 unless defined($opt);
    my $t = readTable(-inh=>$fhi, -header=>1);
    $t->sort("par", 1, 0, "id", 1, 0);
    my @idGs = $t->col("par");
    my $ref = group(\@idGs);
    @idGs = sort {$ref->{$a}->[0] <=> $ref->{$b}->[0]} keys %$ref;
    printf "%d genes : %d models\n", scalar(@idGs), $t->nofRow if $opt == 1;
    my @cnts = uniq(map {$_->[1]} values %$ref);
    for my $cnt (@cnts) {
        my @ids = grep {$ref->{$_}->[1] == $cnt} keys %$ref;
        printf "\t%5d x $cnt model(s)/gene [%s]\n", scalar(@ids), @ids>=0 ? $ids[0] : "" if $opt == 1;
    }
    return ($t, $ref, \@idGs);
}

sub gtb2gtbx {
    my ($fhi, $fho, $fs) = @_;
    my $t = readTable(-inh=>$fhi, -header=>1);
    
    $t->addCol([("") x $t->nofRow], "seq") if $t->colIndex("seq") < 0;
    for my $i (0..$t->nofRow-1) {
        my ($id, $par, $chr, $beg, $end, $srd, $phase, $locS, $cat1, $cat2) = 
            map {$t->elm($i, $_)} qw/id par chr beg end srd phase locC cat1 cat2/;
        next if $cat2 ne "mRNA";
        die "no locCDS for $id\n" unless $locS;
        my $phase1 = [split(",", $phase)]->[0];

        my $rloc = locStr2Ary($locS);
        my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
            [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
        my $seqStr = seqRet($loc, $chr, $srd, $fs);
        my $seq_cds = Bio::Seq->new(-id=>$id, -seq=>$seqStr);
        my $seq_pro = $seq_cds->translate(-frame=>$phase1);
        $t->setElm($i, "seq", $seq_pro->seq);
        printf "%5d | %5d\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    print $fho $t->tsv(1);
    close $fho;
}
sub gtb2gff {
    my ($fhi, $fho) = @_;
    my ($t, $ref, $idGs) = gtbSum($fhi, 0);
    my ($cntG, $cntR) = (1, 1);
    for my $idG (@$idGs) {
        my ($idxB, $cnt) = @{$ref->{$idG}};
        my $ts = $t->subTable([$idxB..$idxB+$cnt-1]);
        my $gene = Gene->new( -gtb=>$ts );
        print $fho $gene->to_gff()."\n";
        for my $rna ($gene->get_rna()) {
            print $fho $rna->to_gff()."\n";
            printf "  Gtb -> Gff %5d RNA | %5d gene...\n", $cntR, $cntG if $cntR % 1000 == 0;
            $cntR ++;
        }
        $cntG ++; 
    }
}

sub get_ovlp_gtb { # opt_srd = 1 (srd sensitive); 2 (srd Insensitive)
    my ($locQ, $chr, $srd, $t, $opt_srd) = rearrange(['loc', 'chr', 'srd', 'tgt', 'opt'], @_);
    $opt_srd ||= 1;
    die "no strand provided\n" if $opt_srd == 1 && !defined($srd);
    $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
    my ($beg, $end) = ($locQ->[0]->[0], $locQ->[-1]->[1]);

    my $t2 = $t->match_pattern_hash("\$_{'chr'} eq '$chr' && \$_{'beg'} <= $end && \$_{'end'} >= $beg");
    if($opt_srd == 1) {
        $t2 = $t2->match_pattern_hash("\$_{'srd'} eq '$srd'");
    }
    return $t2;
}


1;
__END__
