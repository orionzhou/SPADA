#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  validate_dna_seq.pl - replace non-ATCGN characters in the input sequence file with 'N's

=head1 SYNOPSIS
  
  validate_dna_seq.pl [-help] [-out output-file] <input-file|stdin|->

  Options:
      -help   brief help message
      -out    output file, instead of stdout

=head1 DESCRIPTION

  This program reads any number of fasta entries from files or stdin, and 
  replace non-ATCGN characters with 'N's

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<input-file>

  To read fasta entries from stdin, the user should specify 'stdin' or '-'
  in place of the filename argument. For example, one 
  could run this script in a pipe as follows:
      gzip -cd seq.fa.gz | validate_dna_seq.pl - | gzip -9 > seq2.fa.gz

=item B<output-file>

  To write to stdout, the user could either specify 'stdout' or simply leave this
  augument empty.

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

my $fi = '';
my $fo = '';
my $fhi;
my $fho;
my $seq_id = '';
my $seq_desc = '';
my $seq_id_old ='';
my $seq_desc_old = '';
my $wait_4_1st_entry = 'true';
my $seq = '';
my $help_flag;
my %options = (
                              "help|h" => \$help_flag,
                              "out|o=s" => \$fo,
	      );

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;

($fi)= @ARGV;
if(!$fi) {
    pod2usage(2);
} elsif ($fi eq '-' || $fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

if(!$fo || $fo eq "stdout") {
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



