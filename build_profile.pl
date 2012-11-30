#!/usr/bin/perl -w
use strict; 
use Pod::Usage;
use Getopt::Long;
use File::Path qw/make_path remove_tree/;
BEGIN {
        my ($f_cfg) = ('') x 1;
        GetOptions(
                'config|cfg|c=s'    => \$f_cfg, 
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
                $ENV{$k} = $v;
        }
}
use Data::Dumper;
use Common; 
use Seq;
use Align;
use List::Util qw/min max sum/;

sub get_subgroups {
    my ($dir, $fo) = @_;
    print "Extracting gene family IDs\n";
    
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\t", qw/family/)."\n";

    opendir(DH, $dir) or die "cannot open $dir: $!\n";
    for my $fname (sort readdir(DH)) {
        next if $fname =~ /^\./;
        my @ids = read_aln_ids( "$dir/$fname" );
        my $fam = $fname;
        $fam =~ s/\.[^.]+$//;
        print FH join("\t", $fam)."\n"; 
    }
    closedir DH;
    close FH;
}

sub trim_aln {
    my ($f_fam, $di, $do) = @_;
    print "trimming alignments\n";
    
    make_path($do) unless -d $do;
    system("rm -rf $do/*");
    
    my $f_bin = $ENV{"trimAl"}."/bin/trimal";
    die("trimAl not installed: $f_bin not there") unless -s $f_bin;

    my $t = readTable(-in=>$f_fam, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $fi = "$di/$fam.aln";
        my $fo = "$do/$fam.aln";
#    my @ids = read_aln_ids($fi);
        runCmd("$f_bin -in $fi -out $fo -gappyout", 0);
#    runCmd("trimal -in $ft1 -out $ft2 -resoverlap 0.75 -seqoverlap 80", 0);
        die "not work for $fam\n" unless -s $fo;
    }
}
sub check_hmm {
    my ($fi) = @_;
    my $f_bin = $ENV{"HMMER"}."/bin/hmmemit";
    die("cannot execute hmmemit: $f_bin not there") unless -s $f_bin;
    
    my $tag = 0;
    if( -s $fi && open(JJ, "$f_bin $fi |") ) {
        my @lines;
        while ( <JJ> ){
            chomp;
            push @lines, $_;
        }
        $tag = 1 if $lines[0] =~ /^\>/;
    }
    return $tag;
}
sub aln2hmm {
    my ($f_fam, $di, $do) = @_;
    print "converting alignments to HMMs\n";
    
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $f_bin = $ENV{"HMMER"}."/bin/hmmbuild";
    die("cannot execute hmmbuild: $f_bin not there") unless -s $f_bin;
    
    my $t = readTable(-in=>$f_fam, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $fi = "$di/$fam.aln";
        my $fo = "$do/$fam.hmm";
        my $ft1 = "$do/$fam.selex";
        aln_fmt_convert($fi, $ft1, 'clustalw', 'selex');
        die "$ft1 is empty\n" unless -s $ft1;
        
        my $ft2 = "$do/$fam.sum";
        while( !check_hmm($fo) ) {
            runCmd("$f_bin --informat selex -o $ft2 $fo $ft1", -1);
        }
        die "$fo is empty\n" unless -s $fo;
        system("rm $ft1 $ft2");
    }
}

sub check_gap {
    my ($fi) = @_;
    my $ai = Bio::AlignIO->new(-file=>"<$fi");
    my $gap = 0;
    while(my $aln = $ai->next_aln()) {
        for my $seq ($aln->each_seq()) {
            if($seq->seq =~ /[\-\.]/) {
                $gap = 1;
                last;
            }
        }
    }
    return $gap;
}
sub get_hmm_stat {
    my ($f_fam, $d_aln, $d_hmm, $fo) = @_;
    print "extracting MSP statistics\n";
    
    my $f_bin1 = $ENV{"HMMER"}."/bin/hmmstat";
    my $f_bin2 = $ENV{"HMMER"}."/bin/hmmemit";
    die("cannot execute hmmstat: $f_bin1 not there") unless -s $f_bin1;
    die("cannot execute hmmemit: $f_bin2 not there") unless -s $f_bin2;

    my $t = readTable(-in=>$f_fam, -header=>1);
    open(FH, ">$fo");
    print FH join("\t", qw/id nseq length gap consensus/)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $f_hmm = "$d_hmm/$fam.hmm";
        my $f_aln = "$d_aln/$fam.aln";
        
        my $lines = runCmd("$f_bin1 $f_hmm", 2);
        die "cannot get stat for $fam from $f_hmm\n" unless $lines->[-1] =~ /^\s*\d+\s+$fam/;
        my @ps = split " ", $lines->[-1];
        my ($nseq, $len) = @ps[3,5];
        
        $lines = runCmd("$f_bin2 -c $f_hmm", 2);
        die "cannot get con seq for $fam from $f_hmm\n" unless $lines->[0] =~ /^\>$fam/;
        my $seq = join("", @$lines[1..@$lines-1]);

        my $gap = check_gap($f_aln);
        print FH join("\t", $fam, $nseq, $len, $gap, $seq)."\n";
    }
    close FH;
}
sub split_hmm_single {
    my ($fh, $fs, $dir) = @_;
    make_path($dir) unless -d $dir;
    remove_tree($dir, {keep_root => 1});
    
    my $ts = readTable(-in=>$fs, -header=>1);
    for my $i (0..$ts->nofRow-1) {
        my ($id, $nseq, $len) = $ts->row($i);
        my $fo = "$dir/$id.hmm";
        system("hmmfetch $fh $id > $fo");
    }
}

my $dir = $ENV{"SPADA_PROFILE"};
die "cannot open $dir for writing\n" unless -d $dir;
my $d11 = "$dir/11_aln";
die "cannot open $d11 for reading\n" unless -d $d11;
my $f03 = "$dir/03_fam.tbl";
get_subgroups($d11, $f03);
my $d12 = "$dir/12_aln_trim";
trim_aln($f03, $d11, $d12);
my $d15 = "$dir/15_hmm";
aln2hmm($f03, $d12, $d15);
my $f16 = "$dir/16_stat.tbl";
get_hmm_stat($f03, $d12, $d15, $f16);
my $f21 = "$dir/21_all.hmm";
runCmd("cat $d15/*.hmm > $f21");


__END__

