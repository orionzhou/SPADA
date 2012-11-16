#!/usr/bin/perl -w
use strict; 
use Data::Dumper;
use Common; 
use Seq;
use Align;
use List::Util qw/min max sum/;

my $dir = "/project/youngn/zhoup/Data/misc3/spada";
my @orgs = ("Athaliana", "Mtruncatula", "Osativa");
my @subgroups = ("CRP1280", "CRP1300", "CRP1510");
for my $subgroup (@subgroups) {
        my @seqs;
        for my $org (@orgs) {
                my $f_gtb = "$dir/$org/31_model_SPADA/61_final.gtb";
                my $t = readTable(-in=>$f_gtb, -header=>1);
                for my $i (0..$t->nofRow-1) {
                        my ($id, $fam, $seq) = map {$t->elm($i, $_)} qw/id cat3 seq/;
                        next unless $fam eq $subgroup;
                        my $id2 = sprintf "%s_%s", substr($org, 0, 2), $id;
                        my $seqObj = Bio::Seq->new(-id=>$id2, -seq=>$seq);
                        push @seqs, $seqObj;
                }
        }
        my $f_seq = "$dir/cs_aln/$subgroup.fas";
#  writeSeq(\@seqs, $f_seq);
        my $f_aln = "$dir/cs_aln/$subgroup.aln";
        run_clustalo(-in=>$f_seq, -out=>$f_aln);
}


