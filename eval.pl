#!/usr/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Path qw/make_path/;
BEGIN {
        my ($f_cfg, $org) = ('') x 2;
        GetOptions(
                'config|cfg|c=s'    => \$f_cfg, 
                'organism|org|o=s'  => \$org, 
        ) || pod2usage(2);
    
        open(FH, "<$f_cfg") || die "config file $f_cfg is not there\n";
        while(<FH>) {
                chomp;
                next unless $_;
                next if /^\#/;
                $_ =~ s/\s//g;
                my ($k, $v) = split "=";
                while( $v =~ /\$\{(\w+)\}/g ) {
                        die "no env variable named $1\n" unless exists $ENV{$1};
                        my $rep = $ENV{$1};
                        $v =~ s/\$\{$1\}/$rep/;
                }
                if($k eq "SPADA_ORG" && $org) {
                        $ENV{$k} = $org;
                } else {
                        $ENV{$k} = $v;
                }
        }

        my @keys = qw/SPADA_SOURCE SPADA_DATA SPADA_ORG
                ClustalO SignalP HMMER Augustus GeneMark GlimmerHMM GeneID/;
        for my $key (@keys) {
                exists $ENV{$key} || die "$key not defined\n";
        }
        $ENV{"TMP_DIR"} = $ENV{"SPADA_DATA"}."/".$ENV{"SPADA_ORG"};
        make_path($ENV{"TMP_DIR"});
        push @INC, $ENV{"SPADA_SOURCE"};
        $ENV{'PATH'} = join(":", $ENV{"SPADA_SOURCE"}, $ENV{'PATH'});
}
use Common;
use Model;
use Eval;

my $org = $ENV{"SPADA_ORG"};
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
$log->info("##########  working on $org  ##########");

my $dp = $ENV{"SPADA_PROFILE"};
my $dp_aln = "$dp/12_aln_trim";
my $dp_hmm = "$dp/15_hmm";
my $fp_sta = "$dp/16_stat.tbl";

my $d01 = "$dir/01_genome";
my $f01_01 = "$d01/01_refseq.fa";
my $f01_61 = "$d01/61_gene.gtb";

my $d21 = "$dir/21_hits";
my $f21 = "$d21/29_hits.tbl";

my $d31 = "$dir/31_model_SPADA";
my $d34 = "$dir/34_model_Augustus";
my $d35 = "$dir/35_model_GeneMark";
my $d36 = "$dir/36_model_GlimmerHMM";
my $d37 = "$dir/37_model_GeneID";
my $p = { 'SPADA'=>$d31, 'Augustus'=>$d34, 'GeneMark'=>$d35, 'GlimmerHMM'=>$d36, 'GeneID'=>$d37 };
for my $soft (keys(%$p)) {
        next if $soft eq "SPADA";
        next unless $soft eq "GeneID";
        my $dir = $p->{$soft};
#  pipe_model(-dir=>$d37, -hit=>$f21, -ref=>$f01_01, -d_hmm=>$dp_hmm, -d_aln=>$dp_aln, -f_sta=>$fp_sta, -soft=>$soft);
}

my $d41 = "$dir/41_eval";
my $f_gtb_ref = "$d41/01_model.gtb";
pipe_eval($d41, $f_gtb_ref, $p, $org);


