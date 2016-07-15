#!/usr/bin/perl -w
#
# POD documentation
#-------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  spada.pl - Small Peptide Alignment Detection Application

=head1 SYNOPSIS
  
  spada.pl [-help] [options...]

  Options:
    -h (--help)    brief help message
    -c (--cfg)     config file
                     defaut: 'cfg.txt'
    -d (--dir)     SPADA output directory
                     defaut: 'spada.crp.test'
    -p (--hmm)     directory with profile alignments and HMM files
                     defaut: 'hmm.crp'
    -f (--fas)     genome sequence (FASTA) file
                     defaut: 'test/01_refseq.fas'
    -g (--gff)     gene annotation file (GFF3 format, not required)
                     defaut: ''
    -o (--org)     organism: this tells GeneMark / GlimmerHMM / GeneID 
                     which *.mod / training dir / *.param file to use
                     supported: 'Athaliana', 'Mtruncatula', 'Osativa'
                     defaut: 'Athaliana'
    -s (--sp)      if set, predictions without signal peptide are filtered
                     default: TRUE
    -e (--evalue)  E-value threshold
                     default: 0.05
    -t (--threads) threads (precessors) to use
                     default: 1
    -m (--method)  gene prediction programs to run (semicolon-seperated)
                   supported:
                     'Augustus_evidence', 'GeneWise_SplicePredictor', 
                     'Augustus_de_novo', 'GeneMark', 'GlimmerHMM', 'GeneID'
                   default: 'Augustus_evidence;GeneWise_SplicePredictor'

=cut
  
#### END of POD documentation.
#-------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Pod::Usage;
use Getopt::Long;
use Log::Log4perl;
use Time::HiRes qw/gettimeofday tv_interval/;
use File::Path qw/make_path remove_tree/;

use ConfigSetup;
use PrepareGenome;
use MotifMining;
use ModelPred;
use ModelEval;

my $help_flag;
my $f_cfg    = "cfg.txt";
my $dir      = "spada.crp.test";
my $dir_hmm  = "hmm.crp";
my $f_fas    = "test/01_refseq.fas";
my $f_gff    = "";
my $org      = "Athaliana";
my $sp       = 1;
my $e        = 0.05;
my $methods  = "Augustus_evidence;GeneWise_SplicePredictor";
my $ncpu     = 1;

GetOptions(
  "help|h"     => \$help_flag,
  'cfg|c=s'    => \$f_cfg, 
  'dir|d=s'    => \$dir, 
  'hmm|p=s'    => \$dir_hmm, 
  'fas|f=s'    => \$f_fas, 
  'gff|g=s'    => \$f_gff,
  'org|o=s'    => \$org, 
  'sp|s'       => \$sp,
  'evalue|e=f' => \$e,
  'method|m=s' => \$methods,
  'threads|t=i'=> \$ncpu,
) || pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !defined($f_cfg) || ! -s $f_cfg;

config_setup($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, 
  $org, $sp, $e, $methods, $ncpu);

$dir = $ENV{"SPADA_OUT_DIR"};
my $t0 = [gettimeofday];

my $f_log = sprintf "$dir/log.%02d%02d%02d%02d.txt", (localtime(time))[4]+1, (localtime(time))[3,2,1];
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

$log->info("##########  Starting pipeline  ##########");

# Pre-processing
my $d01 = "$dir/01_preprocessing";
my $f01_01 = "$d01/01_refseq.fas";
my $f01_61 = "$d01/61_gene.gtb";
my $f01_12 = "$d01/12_orf_genome.fas";
my $f01_71 = "$d01/71_orf_proteome.fas";

# Motif Mining
my $d11 = "$dir/11_motif_mining";
my $f11 = "$d11/21_hits/29_hits.tbl";

# Model Prediction
my $d21 = "$dir/21_model_prediction";
my $f21_05 = "$d21/05_hits.tbl";
my $f21 = "$d21/30_all.gtb";

# Model Evaluation & Selection 
my $d31 = "$dir/31_model_evaluation";

#pipe_pre_processing($d01);
#pipe_motif_mining(-dir=>$d11, -hmm=>$fp_hmm, -orf_g=>$f01_12, -orf_p=>$f01_71, -ref=>$f01_01, -gtb=>$f01_61);
#pipe_model_prediction(-dir=>$d21, -hit=>$f11, -ref=>$f01_01);
pipe_model_evaluation(-dir=>$d31, -hit=>$f21_05, -gtb_all=>$f21, -ref=>$f01_01, -gtb_ref=>$f01_61, -d_hmm=>$dp_hmm, -d_aln=>$dp_aln, -f_sta=>$fp_sta);
runCmd("ln -sf $d31/61_final.tbl $dir/61_final.tbl", 0);
runCmd("ln -sf $d31/61_final.gtb $dir/61_final.gtb", 0);

$log->info("##########  Pipeline successfully completed  ##########");
$log->info(sprintf("time elapsed: %.01f min", tv_interval($t0, [gettimeofday]) / 60));


__END__


