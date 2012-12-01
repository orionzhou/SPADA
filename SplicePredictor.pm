package SplicePredictor;
use strict;
use Data::Dumper;
use Common;
use Seq;
use Mtb;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/get_splice_sites/;
@EXPORT_OK = qw//;

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
    die "$f_bin not there\n" unless -s $f_bin;
    my $f_fas = $ENV{"TMP_DIR"}."/splice_predictor_".int(rand(1000)).".fa";
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

sub search_splice_sites {
    my ($f_seq, $f_sp, $fo) = @_;
    my $ts = readTable(-in=>$f_seq, -header=>1);
    my $hs = { map {$ts->elm($_, "id") => [$ts->elm($_, "len"), $ts->elm($_, "seq")]} (0..$ts->nofRow-1) };
    
    my $t = readTable(-in=>$f_sp, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id begE1 endE1 begE2/)."\n";
    printf "searching for compatible splice sites... might take long time: %4d in total\n", $t->nofRow;
    for my $i (0..$t->nofRow-1) {
        my ($id, $begE1, $score, $seq_sp) = $t->row($i);
#    next if $id != 3;
        die "no seq info for $id\n" unless exists $hs->{$id};
        my ($lenStr, $seq) = @{$hs->{$id}};
        
        my $len_sp = min(15, length($seq_sp));
        my $begIp = $begE1 + $len_sp*3;
        my ($len_up, $len, $len_dw) = split(/\+/, $lenStr);
        my $endIp = $len_up;

        my @locIs = get_splice_sites($seq, $begIp, $endIp);
        for (@locIs) {
            my ($begI, $endI) = @$_;
            my ($endE1, $begE2) = ($begI-1, $endI+1);
            
            my $locC = [[$begE1, $endE1], [$begE2, $endIp]];
            my $seq_cds = getSubSeq($seq, $locC);
            my $seq_pro = Bio::Seq->new(-seq=>$seq_cds)->translate->seq;
            next if $seq_pro =~ /\*/;

            print FH join("\t", $id, $begE1, $endE1, $begE2)."\n";
        }
        printf "%4d done...\n", $i+1 if ($i+1) % 100 == 0;
    }
    close FH; 
}


1;
__END__

