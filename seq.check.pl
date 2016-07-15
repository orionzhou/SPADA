#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seq.check.pl - replace non-ATCGN characters with 'N's

=head1 SYNOPSIS
  
  seq.check.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (default: stdin)
    -o (--out)    output file (default: stdout)

=head1 DESCRIPTION

  This program reads any number of fasta entries from files or stdin, and 
  replace non-ATCGN characters with 'N's

=cut
  
#### END of POD documentation.
#--------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#-------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $seqHI = Bio::SeqIO->new(-fh => $fhi, -format => 'fasta');
my $seqHO = Bio::SeqIO->new(-fh => $fho, -format => 'fasta');
while(my $seqO = $seqHI->next_seq()) {
  my ($id, $seq) = ($seqO->id, $seqO->seq);
  $seq =~ s/[^ATCGN]/N/ig;
  $seqHO->write_seq(Bio::Seq->new(-id => $id, -seq => $seq));
}
$seqHI->close();
$seqHO->close();

exit 0;



