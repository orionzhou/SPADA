package Eval;
use strict;
use File::Path qw/make_path/;
use Common;
use Genemodel;
use Data::Dumper;
use Gtb;
use Align;
use ModelEval;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/pipe_eval model_eval/;

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


sub align_for_eval {
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

                        system("cp $f_stat $dir/21_stat_$soft\_$e.tbl") if($soft eq "SPADA" && $e == 0.001);
                }
        }
        close FHO;
        system("rm $f_gtb_qry $f_eval $f_stat");
}



1;
__END__


