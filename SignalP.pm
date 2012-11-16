package SignalP;
use strict;
use Data::Dumper;
use Common;
use Seq;
use Mtb;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/run_sigp pipe_search_sigp search_sigp sigp_score_gtb/;
@EXPORT_OK = qw//;

sub run_sigp_hmm {
    my ($seq) = @_;
    die "no SignalP 3.0\n" unless $ENV{"SIGNALP-3.0"};

    $seq =~ s/\*$//;
    my $f_fas = $ENV{"TMP_DIR"}."/signalp_".int(rand(1000)).".fa";
    writeFile($f_fas, ">tmp", $seq);
    my $f_bin = $ENV{"SIGNALP-3.0"}."/signalp";
    die "$f_bin not there\n" unless -s $f_bin;

    open(OUT, "perl $f_bin -t euk -m hmm $f_fas |") or die "failed: $!\n";
    my @stats;
    my $flag = 0;
    while( <OUT> ) {
        chomp;
        my $line = $_;
        if(/^\#\s*pos/) {
            $flag = 1;
            next;
        } elsif(/\>tmp/) {
            $flag = 0;
        }
        next if $flag == 0;
        
        my @ps = split(/\s+/, $line);
        @ps = grep {$_ ne ""} @ps;
        next if @ps < 7;
        push @stats, \@ps;
    }
    my $seq2 = join("", map {$_->[1]} @stats);
    die "parsing error: $seq2\n$seq\n" unless $seq eq $seq2;
    system("rm /tmp/sigp.fa");
  
    my $co_prob = 0.7;
    my @cs = map {$_->[3]} @stats;
    my @cs_sp = grep {$_ > $co_prob} @cs[0..14];
    
#  print join("\t", @cs)."\n" if $seq =~ /^METHVLSRIFLLVLCIYSLKT/;
    my ($tag, $score, $len) = (0, "", 0);
    if(@cs_sp > 10) {
        my @idxs = indexes {$_ > $co_prob} @cs;
        $tag = 1;
        $len = $idxs[$#idxs] + 1;
        $score = sprintf "%.03f", sum(@cs[0..$len-1]) / $len;
    }
    return ($tag, $score, $len);
}
sub run_sigp {
    my ($seq) = @_;
    $seq =~ s/\*$//;
    my $f_fas = $ENV{"TMP_DIR"}."/signalp_".int(rand(1000)).".fa";
    writeFile($f_fas, ">tmp", $seq);
    my $f_bin = $ENV{"SignalP"}."/signalp";
    die "$f_bin not there\n" unless -s $f_bin;
    
    my ($prob1, $prob2, $site) = ("") x 3;
    open(OUT, "perl $f_bin -t euk -s notm $f_fas |") or die "failed: $!\n";
    my ($tag, $d, $pos) = (0, "", "");
    while( <OUT> ) {
        chomp;
        next unless /^tmp/;
        my ($id, $Cmax, $posC, $Ymax, $posY, $Smax, $posS, $Smean, $D, $sp, $Dmaxcut, $network) = split " ";
        ($tag, $pos, $d) = (1, $posY, $D) if $sp eq "Y";
    }
    system("rm $f_fas");
    return ($tag, $d, $pos);
}

sub get_sigp_candidates {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);

    open(FH, ">$fo");
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $strand, $locHStr, $begr, $endr, $lenStr, $e, $seq) = $t->row($i);
        my ($len_up, $len_hit, $len_dw) = split(/\+/, $lenStr);
        my $seq_up = substr($seq, 0, $len_up);
        while($seq_up =~ /ATG/ig) {
            my $pos = $-[0];
            my $seq_test = substr($seq, $pos, 90);
            my $prot = Bio::Seq->new(-seq=>$seq_test)->translate(-frame=>0)->seq;
            $prot =~ s/\*.*//g;
            next if length($prot) < 18;
            print FH ">$id\_$pos\n";
            print FH "$prot\n";
        }
    }
    close FH;
}
sub get_sigp_final {
    my ($fi, $fo) = @_;
    my $ti = readTable(-in=>$fi, -header=>1);
    my $h;
    open(FH, ">$fo");
    print FH join("\t", qw/id pos score/)."\n";
    for my $i (0..$ti->nofRow-1) {
        my ($idStr, $e) = map {$ti->elm($i, $_)} qw/qId e/;
        next if exists $h->{$idStr};
        $h->{$idStr} = $e;
        my ($id, $pos) = split("_", $idStr);
        if($e <= 1) {
            print FH join("\t", $id, $pos, $e)."\n";
        }
    }
    close FH;
}
sub pipe_search_sigp {
    my ($fi, $dir, $f_sp) = @_;
    make_path($dir) unless -d $dir;
    my $f01 = "$dir/01_candidates.fa";
    get_sigp_candidates($fi, $f01);
    my $f02 = "$dir/02.bls";
#  runCmd("blastall -p blastp -d $f_sp -i $f01 -o $f02 -m 0 -F F", 1);
    my $f03 = "$dir/03.mtb";
#  bls2Mtb($f02, $f03);
    my $f04 = "$dir/04_filtered.mtb";
#  mtbFilter($f03, $f04, {pctidty=>0.3, best=>1});
    my $f05 = "$dir/05_sigp.tbl";
#  get_sigp_final($f04, $f05);
}


sub search_sigp {
    my ($fi, $fo) = @_;
    my $t = readTable(-in=>$fi, -header=>1);

    open(FH, ">$fo");
    print FH join("\t", qw/id pos score sp/)."\n";
    printf "  searching sig-pep (might take long time)... (total: %4d)\n", $t->nofRow;
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $strand, $locHStr, $begr, $endr, $lenStr, $e, $seq) = $t->row($i);
        my ($len_up, $len_hit, $len_dw) = split(/\+/, $lenStr);
#    next if $id != 1;
        my $seq_up = substr($seq, 0, $len_up);
        while($seq_up =~ /ATG/ig) {
            my $pos = $-[0];
            my $seq_test = substr($seq, $pos, 120);

            my $prot = Bio::Seq->new(-seq=>$seq_test)->translate(-frame=>0)->seq;
            $prot =~ s/\*.*//g;
            next if length($prot) < 15;

            my ($tag, $score) = run_sigp_hmm($prot);
#      my ($tag, $score) = run_sigp($prot);
            next if $tag == 0;
            print FH join("\t", $id, $pos+1, $score, $prot)."\n";
        }
        printf "%4d done...\r", $i+1 if ($i+1) % 100 == 0;
    }
    print "\n";
    close FH;
}
sub sigp_score_gtb {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("SignalP");
    $log->info("assessing signalp scores");
    my $tg = readTable(-in=>$fi, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id tag score pos/)."\n";
    for my $i (0..$tg->nofRow-1) {
        my ($id, $seq) = map {$tg->elm($i, $_)} qw/id seq/;
#    $seq = substr($seq, 0, 20);
#    next unless $id eq "AT1G70250.1";
        my ($tag, $d, $pos) = run_sigp($seq);
        print FH join("\t", $id, $tag, $d, $pos)."\n";
        printf "  %5d / %5d done...\r", $i+1, $tg->nofRow;
    }
    print "\n";
    close FH;
}



1;
__END__
