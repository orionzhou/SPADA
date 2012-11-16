#!/usr/bin/perl -w
use strict;
use File::Path qw/make_path remove_tree/;
use Common;
use Seq;
use Align;

my $dir = "/project/youngn/zhoup/Data/misc3/spada_profile";
my $f01 = "$dir/01_fam.tbl";
my $d_sup = "$dir/crp_sup";
my $d02 = "$dir/02_seqs";
my $d04 = "$dir/04_aln_core";
my $d05 = "$dir/05_aln_core_hmm";
#copy_core_profile($f01, $d_sup, $d04, $d05);
my $d09 = "$dir/09_seqs_selected";
#pick_selected_seqs($f01, $d02, $d04, $d09);
my $d11 = "$dir/11_aln";
#align_seqs($f01, $d09, $d11);

sub copy_core_profile {
    my ($f_fam, $d_sup, $do1, $do2) = @_;
    my $t = readTable(-in=>$f_fam, -header=>1);
    make_path($do1) unless -d $do1;
    make_path($do2) unless -d $do2;
    remove_tree($do1, {keep_root => 1});
    remove_tree($do2, {keep_root => 1});

    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $fi1 = "$d_sup/$fam/$fam.trim.msf";
        die "$fi1 is not there\n" unless -s $fi1;
        my $fo1 = "$do1/$fam.aln";
        aln_fmt_convert($fi1, $fo1, "msf", "clustalw");
        my $fi2 = "$d_sup/$fam/$fam.trim.hmm";
        die "$fi2 is not there\n" unless -s $fi2;
        my $fo2 = "$do2/$fam.hmm";
        system("cp $fi2 $fo2");
    }
}
sub pick_selected_seqs {
    my ($f_fam, $ds, $di, $do) = @_;
    my $t = readTable(-in=>$f_fam, -header=>1);
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $fs = "$ds/$fam.fas";
        my $fi = "$di/$fam.aln";
        my $hs = readSeq($fs, 2);
        my @ids = read_aln_ids($fi);
        my @seqs;
        for my $id (@ids) {
            next unless exists $hs->{$id};
            my $seq = $hs->{$id};
#      $seq =~ s/\*//g;
            push @seqs, Bio::Seq->new(-id=>$id, -seq=>$seq);
        }
        die "$fi no seq\n" if @seqs <= 1;
        my $fo = "$do/$fam.fas";
        writeSeq(\@seqs, $fo);
    }
}
sub align_seqs {
    my ($fi, $di, $do) = @_;
    print "making MSPs\n";
    make_path($do) unless -d $do;
    remove_tree($do, {keep_root => 1});
    
    my $t = readTable(-in=>$fi, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($fam) = $t->row($i);
        my $fi = "$di/$fam.fas";
        die("$fi is not there") unless -s $fi;
        
        my $h_seq = readSeq($fi, 2);
        my @seqs;
        for my $id (sort(keys(%$h_seq))) {
            my $seq = $h_seq->{$id};
            $seq =~ s/\*//g;
            push @seqs, Bio::Seq->new(-id=>$id, -seq=>$seq);
        }
        my $fo = "$do/$fam.aln";
        run_clustalo(-seqs=>\@seqs, -out=>$fo);
    }
}


