k#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/dirname/;
BEGIN { unshift @INC, dirname(abs_path($0)); }
use Pod::Usage;
use Getopt::Long;
use Log::Log4perl;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;

use Common;
use ConfigSetup;
use ModelEval;
use CompareModel;

my ($f_cfg, $dir, $dir_hmm);
GetOptions(
    'config|cfg|c=s'    => \$f_cfg, 
    'directory|dir|d=s' => \$dir, 
    'profile|hmm|h=s'   => \$dir_hmm, 
) || pod2usage(2);
pod2usage("$0: argument required [--cfg]") if ! defined $f_cfg;
pod2usage("cfg file not there: $f_cfg") if ! -s $f_cfg;

config_setup($f_cfg, $dir, $dir_hmm);

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
my $f21_05 = "$d21/05_hits.tbl";
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
my $f_gtb_gs = "$d41/01_model.gtb";
my $d41_10 = "$d41/10_raw";
#split_pred_sets($d41_10, $f21, $p);

#complete_model_eval($p, $d41, $f21_05, $f01_01, $f_gtb_gs, $dp_hmm, $dp_aln, $fp_sta);

my $f41 = "$d41/51_stat.tbl";
perf_eval($d41, $f_gtb_gs, $f41, $p);

sub split_pred_sets {
    my ($dir, $fi, $p) = @_;
    make_path($dir) unless -d $dir;
    remove_tree($dir, {keep_root=>1});

    my @softs = split(" ", $p->{"All"});
    my $ref_idx = { map {$_=>[]} @softs };
    my @ary_src;
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my @srcs = split(" ", $t->elm($i, "source"));
        for my $src (@srcs) {
            push @{$ref_idx->{$src}}, $i;
        }
        push @ary_src, {map {$_=>1} @srcs};
    }
    
    for my $key (keys %$p) {
        my @softs = split(" ", $p->{$key});
        my $fo = "$dir/$key.gtb";

        my $h;
        for my $soft (@softs) {
            my @idxs = @{$ref_idx->{$soft}};
            for my $idx (@idxs) {
                $h->{$idx} ||= [];
                push @{$h->{$idx}}, $soft;
            }
        }

        my $to = $t->subTable([0..$t->nofRow-1], [$t->header]);
        for my $idx (keys %$h) {
            $to->setElm( $idx, "source", join(" ", @{$h->{$idx}}) );
        }
        $to = $to->subTable([sort {$a<=>$b} keys(%$h)], [$to->header]);

        open(FH, ">$fo") or die "cannot open $fo for writing\n";
        print FH $to->tsv(1);
        close FH;
    }
}

sub complete_model_eval {
    my ($p, $dir, $f_hit, $f_ref, $f_gtb_ref, $d_hmm, $d_aln, $f_sta) = @_;
    my $log = Log::Log4perl->get_logger("main");
    for my $key (keys(%$p)) {
        $log->info("==========  Evaluating $key  ==========");
        my $f_gtb_all = "$dir/10_raw/$key.gtb";
        my $do = "$dir/$key";
        pipe_model_evaluation(-dir=>$do, -hit=>$f_hit, -gtb_all=>$f_gtb_all, -ref=>$f_ref, -gtb_ref=>$f_gtb_ref, -d_hmm=>$d_hmm, -d_aln=>$d_aln, -f_sta=>$f_sta);
    }
}

sub perf_eval {
    my ($dir, $f_gtb_gs, $fo, $p) = @_;
    
    my $log = Log::Log4perl->get_logger("Main");
    
    my @es = (1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 0.01, 0.05, 0.1, 0.5, 1);
    my $score = -1000;
    
    open(FHO, ">$fo") || die "cannot open $fo for writing\n";
    print FHO join("\t", qw/soft e score sn_nt sp_nt sn_exon sp_exon/)."\n";
    
    my $f_gtb_qry = $ENV{"TMP_DIR"}."/qry.gtb";
    my $f_eval = $ENV{"TMP_DIR"}."/eval.tbl";
    my $f_stat = $ENV{"TMP_DIR"}."/stat.tbl";
    
    for my $key (keys %$p) {
        my $fs = "$dir/$key/41_stat.tbl";
        my $fg = "$dir/$key/59.gtb";
        die "$fg is not there\n" unless -s $fg;
        for my $e (@es) {
            filter_models(-stat=>$fs, -in=>$fg, -out=>$f_gtb_qry, -e=>$e, -aln=>$score, -sp=>1, -codon=>1);
            compare_models($f_gtb_qry, $f_gtb_gs, $f_eval);
            model_eval($f_eval, $f_gtb_gs, $f_stat);
            my ($sn_nt, $sp_nt, $sn_ex, $sp_ex) = get_sn_sp($f_stat);
            print FHO join("\t", $key, $e, $score, $sn_nt, $sp_nt, $sn_ex, $sp_ex)."\n";

            system("cp $f_stat $dir/21_stat_$key\_$e.tbl") if $key eq "SPADA" && $e == 0.001;
        }
    }
    close FHO;
    system("rm $f_gtb_qry $f_eval $f_stat");
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



