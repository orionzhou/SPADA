#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gtb2gtbx.pl - convert a Gtb file to Gtbx format

=head1 SYNOPSIS
  
  gtb2gtbx.pl [-help] [-i input-file] [-d refseq-fasta] [-o output-file]

  Options:
      -h (--help)   brief help message
      -i (--in)     input file
      -o (--out)    output file
      -d (--db)     genome sequence database file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Location;
use Seq;
use Gtb;

my ($fi, $fo, $db) = ('') x 3;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "db|d=s"   => \$db,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$db;

my ($fhi, $fho);
if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot open $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $t = readTable(-inh => $fhi, -header => 1);
$t->addCol([("") x $t->nofRow], "seq") if $t->colIndex("seq") < 0;

my $cntp = 0;
my @idxs;
for my $i (0..$t->lastRow) {
  my ($id, $par, $chr, $beg, $end, $srd, $phase, $locS, $cat1, $cat2) = 
    map {$t->elm($i, $_)} qw/id par chr beg end srd phase cloc cat1 cat2/;
  $cat1 eq "mRNA" || next;
  $locS || die "no CDS for $id\n";
  my $phase1 = [split(",", $phase)]->[0];
 
  my $rloc = locStr2Ary($locS);
  my $loc = $srd eq "-" ? [map {[$end-$_->[1]+1, $end-$_->[0]+1]} @$rloc] : 
    [map {[$beg+$_->[0]-1, $beg+$_->[1]-1]} @$rloc];
  my $seqStr = seqRet($loc, $chr, $srd, $db);
  my $seq_cds = Bio::Seq->new(-id=>$id, -seq=>$seqStr);
  my $seq_pro = $seq_cds->translate(-frame=>$phase1);
  my $seqpro = $seq_pro->seq;
  
  $cntp ++ if $seqpro =~ /\*[\w\*]/;
  $seqpro =~ s/\*$//;
  $seqpro =~ s/\*/X/;
  $t->setElm($i, "seq", $seqpro);
  my $lent = length($seqpro);

  $seqpro =~ s/X//g;
  my $lenp = length($seqpro);
  push @idxs, $i if $lenp < 10 || $lenp/$lent < 0.5;
}
$t->delRows(\@idxs);
print $fho $t->tsv(1);
close $fho;

printf "%d w. pre-mature stop codons\n", $cntp;
printf "%d short (<10 aa) and incomplete (>50%% 'X') peptides removed\n", scalar(@idxs);


__END__
