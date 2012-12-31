package Spada;
use strict;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper;
use Common;
use Seq;
use Align;
use Log::Log4perl;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_genewise_splicepredictor/;
@EXPORT_OK = qw//;

sub run_genewise_batch {
    my ($fi, $dirO) = @_;
    my $log = Log::Log4perl->get_logger("Spada");
    $log->info("running genewise");
    make_path($dirO) unless -d $dirO;
    remove_tree($dirO, {keep_root => 1});
    
    my $f_bin = $ENV{"GeneWise"}."/bin/genewise";
    $ENV{"WISECONFIGDIR"} = $ENV{"GeneWise"}."/wisecfg";
    die "$f_bin not there\n" unless -s $f_bin;
    
    my $f_dna = $ENV{"TMP_DIR"}."/genewise_dna_".int(rand(1000)).".fa";
    my $f_pro = $ENV{"TMP_DIR"}."/genewise_pro_".int(rand(1000)).".fa";
    
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $begG, $endG, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
        my $loc = locStr2Ary($locLS);
        next if @$loc == 1;
        writeFile($f_dna, ">seq", $seq);
        writeFile($f_pro, ">pro", $seqP);
        my $fo = "$dirO/$id.txt";
        runCmd("$f_bin -gff $f_pro $f_dna > $fo", 0);
        printf "  %5d / %5d done\r", ($i+1), $t->nofRow;
    }
    print "\n";
    system("rm $f_dna $f_pro");
}
sub sum_genewise1 {
    my ($fi) = @_;
    open(FHT, "<$fi") || die "cannot open $fi for reading\n";
    my $loc = [];
    while(<FHT>) {
        chomp;
        my @ps = split "\t";
        if(@ps == 9 && $ps[1] eq "GeneWise" && $ps[2] eq "cds") {
            push @$loc, [$ps[3], $ps[4]];
        } 
    }
    close FHT;
    return $loc;
}

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
sub build_model {
    my ($fi, $fo, $dirG) = @_;
    my $log = Log::Log4perl->get_logger("Spada");
    $log->info("building exon models");
    my $t = readTable(-in=>$fi, -header=>1);
    
    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id parent family beg end strand n_cds locC seq/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $begG, $endG, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
        my $loc = locStr2Ary($locLS);
        $loc = [sort {$a->[0] <=> $b->[0]} @$loc];
        my $pos_end = get_stop_codon($seq, $endL+1);
        my @pos_begs = get_start_codons($seq, $begL+2);
        my $n_cds = @$loc;

        my @locCs;
        if($n_cds == 1) {
            for my $pos_beg (@pos_begs) {
                push @locCs, [[$pos_beg, $pos_end]];
            }
        } else {
            my $locC = sum_genewise1("$dirG/$id.txt");
#      die Dumper($locC, $pos_end, \@pos_begs) if $id == 200;
            if(@$locC == $n_cds) {
                $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
                next if @pos_begs == 0 || $locC->[0]->[0] < max(@pos_begs) || $locC->[-1]->[1] > $pos_end;
                for my $pos_beg (@pos_begs) {
                    my $locT = [ map {[$_->[0], $_->[1]]} @$locC ];
                    $locT->[0]->[0] = $pos_beg;
                    $locT->[-1]->[1] = $pos_end;
                    push @locCs, $locT;
                }
            } elsif($n_cds == 2) {
                my @posIs = get_splice_sites($seq, $loc->[0]->[1]+1, $loc->[1]->[0]-1);
                for (@posIs) {
                    my ($endE2, $begE1) = ($_->[0]-1, $_->[1]+1);
                    for my $pos_beg (@pos_begs) {
                        push @locCs, [[$pos_beg, $endE2], [$begE1, $pos_end]];
                    }
                }
            }
        }
        
        for my $i (1..@locCs) {
            my $locC = $locCs[$i-1];
            $locC = [ sort {$a->[0] <=> $b->[0]} @$locC ];
            my ($beg, $end) = ($locC->[0]->[0], $locC->[-1]->[1]);
            my $seq = getSubSeq($seq, $locC);
            my $prot = Bio::Seq->new(-seq=>$seq)->translate()->seq;
            print FH join("\t", "$id.$i", $id, $fam, $beg, $end, "+", $n_cds, locAry2Str($locC), $prot)."\n";
        }
    }
    close FH;
}
sub output_gtb {
    my ($fi, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Spada");
    $log->info("writing output in Gtb format");
    my $t = readTable(-in=>$fi, -header=>1);
    
    open(FH, ">$fo");
    print FH join("\t", qw/id parent chr beg end strand locE locI locC loc5 loc3 phase source conf cat1 cat2 cat3 note seq/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $pa, $fam, $beg, $end, $strand, $n_cds, $locCStr, $seq) = $t->row($i);
        my $locC = locStr2Ary($locCStr);
        my $phase = join(",", getPhase($locC, "+"));
        print FH join("\t", $id, $pa, ("") x 3, "+", ("") x 2, $locCStr, "", "", $phase, "", ("") x 5, $seq)."\n";
    }
    close FH;
}
sub pipe_genewise_splicepredictor {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;

    my $d02 = "$dir/02_genewise";
    run_genewise_batch($fi, $d02);
    my $f05 = "$dir/05.tbl";
    build_model($fi, $f05, $d02);
    my $f11 = "$dir/11.gtb";
    output_gtb($f05, $f11);
}



1;
__END__
