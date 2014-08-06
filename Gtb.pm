package Gtb;
use strict;
use Common;
use Location; 
use Seq;
use Gene;
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/@HEAD_GTB @HEAD_GTBX
  read_gtb read_gtb_hash get_ovlp_gtb/;
@EXPORT_OK = qw//;

our @HEAD_GTB  = qw/id par chr beg end srd eloc iloc cloc floc tloc phase src conf cat1 cat2 cat3 note/;

our @HEAD_GTBX = qw/id par chr beg end srd eloc iloc cloc floc tloc phase src conf cat1 cat2 cat3 note seq/;

my ($id, $par, $chr, $beg, $end, $srd, 
  $elocS, $ilocS, $clocS, $flocS, $tlocS, $phase, 
  $src, $conf, $cat1, $cat2, $cat3, $note) = ();

sub read_gtb_hash {
  my ($fi) = @_;
  my $t = readTable(-in => $fi, -header => 1);
  return { map {$t->elm($_, "id") => $t->rowRef($_)} (0..$t->nofRow-1) };
}

sub get_ovlp_gtb { # opt_srd = 1 (srd sensitive); 2 (srd Insensitive)
  my ($locQ, $chr, $srd, $t, $opt_srd) = rearrange(['loc', 'chr', 'srd', 'tgt', 'opt'], @_);
  $opt_srd ||= 1;
  die "no strand provided\n" if $opt_srd == 1 && !defined($srd);
  $locQ = [ sort {$a->[0] <=> $b->[0]} @$locQ ];
  my ($beg, $end) = ($locQ->[0]->[0], $locQ->[-1]->[1]);

  my $t2 = $t->match_pattern_hash("\$_{'chr'} eq '$chr' && \$_{'beg'} <= $end && \$_{'end'} >= $beg");
  if($opt_srd == 1) {
    $t2 = $t2->match_pattern_hash("\$_{'srd'} eq '$srd'");
  }
  return $t2;
}


1;
__END__
