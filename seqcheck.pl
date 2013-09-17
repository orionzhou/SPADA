#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqcheck.pl - replace non-ATCGN characters in the input sequence file with 'N's

=head1 SYNOPSIS
  
  seqcheck.pl [-help] [-in input-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file (could be 'stdin')
      -out    output file (could be 'stdout')

=head1 DESCRIPTION

  This program reads any number of fasta entries from files or stdin, and 
  replace non-ATCGN characters with 'N's

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=back
  
=head1 BUGS
  
=head1 REFERENCES
  
=head1 VERSION
  
  0.1
  
=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $fhi;
my $fho;
my $seq_id = '';
my $seq_desc = '';
my $seq_id_old ='';
my $seq_desc_old = '';
my $wait_4_1st_entry = 'true';
my $seq = '';
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

if ($fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if ($fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHI = Bio::SeqIO->new(-fh=>$fhi, -format=>'fasta');
my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
while(my $seqO = $seqHI->next_seq()) {
    my ($id, $seq) = ($seqO->id, $seqO->seq);
    $seq =~ s/[^ATCGN]/N/ig;
    $seqHO->write_seq(Bio::Seq->new(-id=>$id, -seq=>$seq));
}
$seqHI->close();
$seqHO->close();

exit 0;



