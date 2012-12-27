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

use ConfigSetup;
use Model;

my ($f_cfg, $org, $f_fas, $f_gff) = ('') x 4;
my $cutoff_e = ""; 
GetOptions(
    'config|cfg|c=s'    => \$f_cfg, 
    'organism|org|o=s'  => \$org, 
    'fasta|fas|f=s'     => \$f_fas, 
    'gff|g=s'           => \$f_gff,
    'evalue|e=f'        => \$cutoff_e
) || pod2usage(2);

config_setup($f_cfg, $org, $f_fas, $f_gff, $cutoff_e);

my $dir = $ENV{"SPADA_DATA"}."/".$org;

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

my $dp = $ENV{"SPADA_PROFILE"};
my $dp_aln = "$dp/12_aln_trim";
my $dp_hmm = "$dp/15_hmm";
my $fp_sta = "$dp/16_stat.tbl";
my $fp_hmm = "$dp/21_all.hmm";

my $d01 = "$dir/01_preprocessing";
my $f01_01 = "$d01/01_refseq.fa";
my $f01_61 = "$d01/61_gene.gtb";

my $d21 = "$dir/21_hits";
my $f21 = "$d21/29_hits.tbl";

my $d31 = "$dir/31_model_SPADA";
my $d34 = "$dir/34_model_Augustus";
my $d35 = "$dir/35_model_GeneMark";
my $d36 = "$dir/36_model_GlimmerHMM";
my $d37 = "$dir/37_model_GeneID";
for my $soft (keys(%$p)) {
    next if $soft eq "SPADA";
    next unless $soft eq "GeneID";
    my $dir = $p->{$soft};
#  pipe_model(-dir=>$d37, -hit=>$f21, -ref=>$f01_01, -d_hmm=>$dp_hmm, -d_aln=>$dp_aln, -f_sta=>$fp_sta, -soft=>$soft);
}

my $d41 = "$dir/41_eval";
my $f_gtb_ref = "$d41/01_model.gtb";
#pipe_eval($d41, $f_gtb_ref, $p, $org);


