#!/usr/bin/perl -w
#
# POD documentation
#-------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  spada.pl - Secreted Peptide Alignment Detection Algorithm

=head1 SYNOPSIS
  
  spada.pl [-help] <-cfg config-file> [options...]

  Options:
      -help   brief help message
      -cfg    config file
      -dir    SPADA output directory
      -hmm    directory containing profile alignments and HMM files
      -fas    genome sequence file (FASTA format)
      -gff    gene annotation file (GFF3 format, optional)
      -org    organism to run
      -evalue E-value threshold

=head1 DESCRIPTION

  This program identify and predict the structure of cysteine-rich peptides in plant genomes

=cut
  
#### END of POD documentation.
#-------------------------------------------------------------------------

use strict;
use Cwd qw/abs_path/;
use File::Basename qw/dirname/;
BEGIN { unshift @INC, dirname(abs_path($0)); }
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

my ($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, $org, $cutoff_e);
GetOptions(
    'config|cfg|c=s'    => \$f_cfg, 
    'directory|dir|d=s' => \$dir, 
    'profile|hmm|h=s'   => \$dir_hmm, 
    'fas|f=s'           => \$f_fas, 
    'gff|g=s'           => \$f_gff,
    'organism|org|o=s'  => \$org, 
    'evalue|e=f'        => \$cutoff_e,
) || pod2usage(2);
pod2usage("$0: argument required [--cfg]") if ! defined $f_cfg;
pod2usage("cfg file not there: $f_cfg") if ! -s $f_cfg;

config_setup($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, $org, $cutoff_e);

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
my $f01_01 = "$d01/01_refseq.fa";
my $f01_61 = "$d01/61_gene.gtb";
my $f01_12 = "$d01/12_orf_genome.fa";
my $f01_71 = "$d01/71_orf_proteome.fa";
pipe_pre_processing($d01);

# Motif Mining
my $d11 = "$dir/11_motif_mining";
my $f11 = "$d11/21_hits/29_hits.tbl";
pipe_motif_mining(-dir=>$d11, -hmm=>$fp_hmm, -orf_g=>$f01_12, -orf_p=>$f01_71, -ref=>$f01_01, -gtb=>$f01_61);

# Model Prediction
my $d21 = "$dir/21_model_prediction";
my $f21_05 = "$d21/05_hits.tbl";
my $f21 = "$d21/26_all.gtb";
pipe_model_prediction(-dir=>$d21, -hit=>$f11, -ref=>$f01_01);

# Model Evaluation & Selection 
my $d31 = "$dir/31_model_evaluation";
pipe_model_evaluation(-dir=>$d31, -hit=>$f21_05, -gtb_all=>$f21, 
    -ref=>$f01_01, -gtb_ref=>$f01_61, -d_hmm=>$dp_hmm, -d_aln=>$dp_aln, -f_sta=>$fp_sta);

$log->info("##########  Pipeline successfully completed  ##########");
$log->info(sprintf("time elapsed: %.01f min", tv_interval($t0, [gettimeofday]) / 60));


__END__


