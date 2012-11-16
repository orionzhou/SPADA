#!/usr/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use Log::Log4perl;
use Time::HiRes qw/gettimeofday tv_interval/;
use File::Path qw/make_path/;
BEGIN {
    my ($f_cfg, $org, $f_fas, $f_gff) = ('') x 4;
    my $cutoff_e = ""; 
    GetOptions(
        'config|cfg|c=s'    => \$f_cfg, 
        'organism|org|o=s'  => \$org, 
        'fasta|fas|f=s'     => \$f_fas, 
        'gff|g=s'           => \$f_gff,
        'evalue|e=f'        => \$cutoff_e
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
        } elsif($k eq "SPADA_FAS" && $f_fas) {
            $ENV{$k} = $f_fas;
        } elsif($k eq "SPADA_GFF" && $f_gff) {
            $ENV{$k} = $f_gff;
        } elsif($k eq "evalue" && $cutoff_e) {
            $ENV{$k} = $cutoff_e;
        } else {
            $ENV{$k} = $v;
        }
    }

    my @keys = qw/SPADA_SOURCE SPADA_PROFILE SPADA_DATA SPADA_ORG SPADA_FAS
        ClustalO GeneWise SplicePredictor SignalP HMMER Augustus/;
    for my $key (@keys) {
        exists $ENV{$key} || die "$key not defined\n";
    }
    $ENV{"TMP_DIR"} = $ENV{"SPADA_DATA"}."/".$ENV{"SPADA_ORG"};
    make_path($ENV{"TMP_DIR"});
    push @INC, $ENV{"SPADA_SOURCE"};
    $ENV{'PATH'} = join(":", $ENV{"SPADA_SOURCE"}, $ENV{'PATH'});
}
use Common;
use PrepareGenome;
use Hmm;
use Hits;
use Model;
use Spada;

my $org = $ENV{"SPADA_ORG"};
my $dir = $ENV{"SPADA_DATA"}."/".$org;
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

my $dp = $ENV{"SPADA_PROFILE"};
my $dp_aln = "$dp/12_aln_trim";
my $dp_hmm = "$dp/15_hmm";
my $fp_sta = "$dp/16_stat.tbl";
my $fp_hmm = "$dp/21_all.hmm";

$log->info("##########  working on $org  ##########");

my $d01 = "$dir/01_genome";
my $f01_01 = "$d01/01_refseq.fa";
my $f01_61 = "$d01/61_gene.gtb";
my $f01_12 = "$d01/12_orf_genome.fa";
my $f01_71 = "$d01/71_orf_protein.fa";
pipe_prepare_genome($d01, $org);

my $d11 = "$dir/11_hmmsearch_x";
my $f11 = "$d11/07_final.htb";
pipe_hmmsearch(-dir=>$d11, -hmm=>$fp_hmm, -target=>$f01_12, -ref=>$f01_01, -gtb=>$f01_61);
my $d12 = "$dir/12_hmmsearch_p";
my $f12 = "$d12/07_final.htb";
pipe_hmmsearch(-dir=>$d12, -hmm=>$fp_hmm, -target=>$f01_71, -ref=>$f01_01, -gtb=>$f01_61) if -s $f01_71;
my $d21 = "$dir/21_hits";
my $f21 = "$d21/29_hits.tbl";
pipe_hit(-dir=>$d21, -in=>[$f11, $f12], -ref=>$f01_01, -aln=>$dp_aln, -p=>{min_e=>5, min_len=>30});

my $d31 = "$dir/31_model_SPADA";
my $soft = "SPADA";
pipe_model(-dir=>$d31, -hit=>$f21, -ref=>$f01_01, -d_hmm=>$dp_hmm, -d_aln=>$dp_aln, -f_sta=>$fp_sta, -soft=>'SPADA');

pipe_model_postprocess(-dir=>$d31, -f_sta=>$fp_sta, -f_gtb=>$f01_61);

$log->info(sprintf("time elapsed: %.01f min", tv_interval($t0, [gettimeofday]) / 60));


__END__


