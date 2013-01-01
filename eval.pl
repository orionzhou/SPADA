#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/dirname/;
BEGIN { unshift @INC, dirname(abs_path($0)); }
use Pod::Usage;
use Getopt::Long;
use Log::Log4perl;
use Time::HiRes qw/gettimeofday tv_interval/;
use File::Path qw/make_path remove_tree/;

use Common;
use ConfigSetup;
use ModelEval;

my ($f_cfg, $dir);
GetOptions(
    'config|cfg|c=s'    => \$f_cfg, 
    'directory|dir|d=s' => \$dir, 
) || pod2usage(2);
pod2usage("$0: argument required [--cfg]") if ! defined $f_cfg;
pod2usage("cfg file not there: $f_cfg") if ! -s $f_cfg;

config_setup($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, $org, $cutoff_e);

$dir = $ENV{"SPADA_OUT_DIR"};

my $f_log = sprintf "$dir/log.eval.%02d%02d%02d%02d.txt", (localtime(time))[4]+1, (localtime(time))[3,2,1];
my $log_conf = qq/
    log4perl.category                  = INFO, Logfile, Screen

    log4perl.appender.Logfile          = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename = $f_log
    log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.ConversionPattern = %d{HH:mm:ss} %F{1} %L> %m %n

    log4perl.appender.Screen           = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr    = 0
    log4perl.appender.Screen.layout    = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Screen.layout.ConversionPattern = [%d{HH:mm:ss}] %m %n
/;
Log::Log4perl->init(\$log_conf);

my $log = Log::Log4perl->get_logger("main");

my $dp = $ENV{"SPADA_HMM_DIR"};
my $dp_aln = "$dp/12_aln_trim";
my $dp_hmm = "$dp/15_hmm";
my $fp_sta = "$dp/16_stat.tbl";
my $fp_hmm = "$dp/21_all.hmm";

my $d01 = "$dir/01_preprocessing";
my $f01_01 = "$d01/01_refseq.fa";
my $f01_61 = "$d01/61_gene.gtb";

my $d21 = "$dir/21_model_prediction";
my $f21 = "$d21/26_all.gtb";

my $p = {
    All => "Augustus_evidence GeneWise_SplicePredictor Augustus_de_novo GeneMark GlimmerHMM GeneID",
    SPADA => "Augustus_evidence GeneWise_SplicePredictor",
    Augustus_evidence => "Augustus_evidence",
    GeneWise_SplicePredictor => "GeneWise_SplicePredictor",
    Augustus_de_novo => "Augustus_de_novo",
    GeneMark => "GeneMark",
    GlimmerHMM => "GlimmerHMM",
    GeneID => "GeneID"
};
my $d41 = "$dir/41_perf_eval";
split_pred_sets($d41, $f21, $p);
sub split_pred_sets {
    my ($dir, $fi, $p) = @_;
    make_path($dir) unless -d $dir;

    my @softs = split(" ", $p->{"All"});
    my $ref_idx = { map {$_=>[]} @softs };
    my $t = readTable(-in=>$f21, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my @srcs = split(" ", $t->elm($i, "source"));
        for my $src (@srcs) {
            push @{$ref_idx->{$src}}, $i;
        }
    }

    for my $soft (@softs) {
        my $fo = "$dir/$soft.gtb";
        my $to = $t->subTable($ref_idx->{$soft}, [$t->header]);
        open(FH, ">$fo") or die "cannot open $fo for writing\n";
        print FH $to->tsv(1);
        close FH;
    }
}

my $f_gtb_gs = "$d41/01_model.gtb";

for my $soft (keys %{$ENV{"method"}}) {
}
sub model_eval {
    my ($f_eval, $f_ref_gtb, $fo) = @_;
    my $tv = readTable(-in=>$f_eval, -header=>1);
    my $tr = readTable(-in=>$f_ref_gtb, -header=>1);

    my $hg;
    for my $i (0..$tr->nofRow-1) {
        my ($id, $locS) = map {$tr->elm($i, $_)} qw/id locC/;
        my $loc = locStr2Ary($locS);
        my $len = locAryLen($loc);
        my $n_cds = @$loc;
        $hg->{$id} = [0, $len, $n_cds];
    }

    open(FH, ">$fo") || die "cannot open $fo for writing\n";
    print FH join("\t", qw/id tag gene lenTP lenFP lenFN exonTP exonFP exonFN/)."\n";
    for my $i (0..$tv->nofRow-1) {
        my ($id, $gene, $tag, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN) = $tv->row($i);
        $tag = 5 if $tag == 7 || $tag == 8 || ($tag==2 && $lenFP+$lenFN>30);

        print FH join("\t", $id, $tag, $gene, $lenTP, $lenFP, $lenFN, $exonTP, $exonFP, $exonFN)."\n";
        $hg->{$gene}->[0] ++ if $gene ne "";
    }
    for my $gene (keys(%$hg)) {
        my ($cnt, $len, $n_cds) = @{$hg->{$gene}};
        if($cnt == 0) {
            print FH join("\t", '', 10, $gene, 0, 0, $len, 0, 0, $n_cds)."\n";
        } elsif($cnt > 1) {
            print "  $gene hit $cnt times\n";
        }
    }
    close FH;
}
sub get_sn_sp {
        my ($fi) = @_;
        my $t = readTable(-in=>$fi, -header=>1);
        my $sn_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFN")) );
        my $sp_nt = sum($t->col("lenTP")) / ( sum($t->col("lenTP")) + sum($t->col("lenFP")) );
        my $sn_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFN")) );
        my $sp_ex = sum($t->col("exonTP")) / ( sum($t->col("exonTP")) + sum($t->col("exonFP")) );
        return ($sn_nt, $sp_nt, $sn_ex, $sp_ex);
}


sub align_for_visual_inspection {
    my ($fg1, $fg2, $fs, $dirO) = rearrange(['gtb1', 'gtb2', 'stat', 'out'], @_);
    make_path($dirO) unless -d $dirO;
    system("rm -rf $dirO/*");

    my $ts = readTable(-in=>$fs, -header=>1);
    my $hg1 = readGtb(-in=>$fg1, -opt=>2);
    my $hg2 = readGtb(-in=>$fg2, -opt=>2);
    die "not equal rows btw $fg1 and $fs\n" unless scalar(keys(%$hg1)) == $ts->nofRow;

    my $h;
    for my $i (0..$ts->nofRow-1) {
        my ($id1, $id2, $tag) = $ts->row($i);
        my ($fam, $seq1) = @{$hg1->{$id1}}[16,18];
        $h->{$fam} ||= [];
        push @{$h->{$fam}}, Bio::Seq->new(-id=>$id1, -seq=>$seq1);
        if($tag > 1 && $tag < 9) {
            my $seq2 = $hg2->{$id2}->[18];
            push @{$h->{$fam}}, Bio::Seq->new(-id=>"$tag\_$id2", -seq=>$seq2);
        }
    }

    for my $fam (sort keys %$h) {
        my $seqs = $h->{$fam};
        next if @$seqs == 1;
        my $f_aln = "$dirO/$fam.aln";
        run_clustalo(-seqs=>$seqs, -out=>$f_aln);
    }
}
sub pipe_eval {
    my ($dir, $f_gtb_ref, $p, $org) = @_;
    
    my $fo = "$dir/51_stat.tbl";
    open(FHO, ">$fo") || die "cannot open $fo for writing\n";
    print FHO join("\t", qw/org soft e score sn_nt sp_nt sn_exon sp_exon/)."\n";
    
    my $f_gtb_qry = $ENV{"TMP_DIR"}."/qry.gtb";
    my $f_eval = $ENV{"TMP_DIR"}."/eval.tbl";
    my $f_stat = $ENV{"TMP_DIR"}."/stat.tbl";
    
    my @es = (1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 0.01, 0.05, 0.1, 0.5, 1);
    my $score = -1000;
    my $opt_mt = 1 if $org eq "Mtruncatula";
    for my $soft (keys(%$p)) {
        print "Evaluating performace for $soft\n";
        my $dirI = $p->{$soft};
        my $fi_sta = "$dirI/41_stat.tbl";
        my $fi_gtb = "$dirI/55_nonovlp.gtb";
        die "$fi_gtb is not there\n" unless -s $fi_gtb;
        for my $e (@es) {
            filter_models($fi_sta, $fi_gtb, $f_gtb_qry, $e, $score, $opt_mt);
            compare_models($f_gtb_qry, $f_gtb_ref, $f_eval);
            model_eval($f_eval, $f_gtb_ref, $f_stat);
            my ($sn_nt, $sp_nt, $sn_ex, $sp_ex) = get_sn_sp($f_stat);
            print FHO join("\t", $org, $soft, $e, $score, $sn_nt, $sp_nt, $sn_ex, $sp_ex)."\n";

            system("cp $f_stat $dir/21_stat_$soft\_$e.tbl") if $soft eq "SPADA" && $e == 0.001;
        }
    }
    close FHO;
    system("rm $f_gtb_qry $f_eval $f_stat");
}


