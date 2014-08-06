#!/usr/bin/perl -w
use strict;
use Log::Log4perl;
use File::Path qw/make_path remove_tree/;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use Common;
use Hmm;
use Hits;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_motif_mining/;
@EXPORT_OK = qw//;

sub pipe_motif_mining {
  my ($dir, $f_hmm, $f_orf_g, $f_orf_p, $f_ref, $f_gtb) =
    rearrange([qw/dir hmm orf_g orf_p ref gtb/], @_);

  my $log = Log::Log4perl->get_logger("MotifMining");
  $log->info("#####  Stage 2 [Motif Mining]  #####");

  my $d11 = "$dir/11_hmmsearch_x";
  my $f11 = "$d11/07_final.htb";
  pipe_hmmsearch(-dir=>$d11, -hmm=>$f_hmm, -target=>$f_orf_g, -ref=>$f_ref, -gtb=>$f_gtb);
  my $d12 = "$dir/12_hmmsearch_p";
  my $f12 = "$d12/07_final.htb";
  pipe_hmmsearch(-dir=>$d12, -hmm=>$f_hmm, -target=>$f_orf_p, -ref=>$f_ref, -gtb=>$f_gtb) if defined($f_orf_p) && -s $f_orf_p;
  
  my $d21 = "$dir/21_hits";
  pipe_hit(-dir=>$d21, -in=>[$f11, $f12], -ref=>$f_ref, -p=>{min_e=>5, min_len=>30});
}


1;
__END__
