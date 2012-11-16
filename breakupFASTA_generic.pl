#!/usr/bin/perl -w
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#           Copyright 2002, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#  Kevin Silverstein              11 March 2004
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  breakupFASTA_generic.pl - breakup FASTA sequences into smaller overlapping chunks

=head1 SYNOPSIS
  
  breakupFASTA_generic.pl [-help] [-ovl nres] [-size nres] <fasta-file+|stdin>

  Options:
      -help     brief help message
      -ovl      overlap in adjacent segments of the sequence
      -size     size of each segment

=head1 DESCRIPTION

  This program reads any number of fasta entries from files or stdin, and 
  outputs the same sequences in overlapping chunks of an arbitrary size.
  Unlike breakupFASTA.pl, this program has no size limit on the size of the
  input sequences.

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<-ovl nres>

  Specify the number of residues that each block overlaps.
  Default = 1000.

=item B<-size nres>

  Specify the number of residues in each sequence block to output.
  Default = 2000.

=item B<fasta-files>

  To read fasta entries from stdin, the user should specify 'stdin'
  in place of the filename argument. For example, one 
  could uncompress a file and translate it as follows:
        gzip -cd ex.fas.gz | breakupFASTA.pl stdin | gzip -9 > ex2.fas.gz

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

my $fhandle;
my $infile = '';
my $seq_id = '';
my $seq_desc = '';
my $old_seq_id ='';
my $old_seq_desc = '';
my $wait_4_1st_entry = 'true';
my $seq = '';
my $curr_pos_in_seq = 0;
my $buffer_size;
my $increment;
my $seg_ovl = 1000;
my $seg_size = 2000;
my $help_flag;
my $tab_output = '';
my %options = (
                              "help|h" => \$help_flag,
                              "ovl|o=s" => \$seg_ovl,
                              "size|s=s" => \$seg_size
	      );

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage({-message => "Overlap ($seg_ovl) cannot exceed Sequence segment size ($seg_size).\n", -exitval => 2}) if ($seg_ovl >= $seg_size);

# We will read in one line at at time and process a portion of the sequence
# when we exceed a defined buffer size.  So we will set the buffer size to be
# the size of our segment.
$buffer_size = $seg_size;
$increment = $seg_size - $seg_ovl;

($infile)= @ARGV;
if ($infile eq 'stdin') {
    $fhandle = \*STDIN;
}
elsif ($infile =~ /\.g?z/i) {
    open ($fhandle, "gzip -cd |") || die "Can't uncompress file $infile: $!\n";
}
else {
    open ($fhandle, $infile) || die "Can't open file $infile: $!\n";
}

while (<$fhandle>) {
    chomp;
    if (/^>\s*(\S+)(.*)/) {
        $old_seq_id = $seq_id;
        $old_seq_desc = $seq_desc;
        ($seq_id, $seq_desc) = ($1, $2);
	
        if ($wait_4_1st_entry eq 'true') {
            $wait_4_1st_entry = 'false';
        }
        else { 
            &dump_remainder($old_seq_id, $old_seq_desc, $curr_pos_in_seq, $seq);
            $seq='';
        }

        $curr_pos_in_seq = 0;
    }
    elsif (!(/^\s*$/)) {
        $seq .= $_;
#    print "DEBUG $seq_id $curr_pos_in_seq: $seq\n";
        while (length($seq) >= $buffer_size) {
	my $segment_defline = $seq_id . "_" . ($curr_pos_in_seq+1) .
	                 $seq_desc;
	&print_sequence($segment_defline, substr($seq, 0, $seg_size), 60);
	$seq = substr($seq, $increment);
	$curr_pos_in_seq += $increment;
        }
    }
}

#print last record
&dump_remainder($seq_id, $seq_desc, $curr_pos_in_seq, $seq);

close ($fhandle);

exit 0;

###############################################################################
#subroutines
#############################################################################

sub print_sequence {

        my ($defline, $sequence, $length) = @_;
        my $seqlen = length($sequence);

        print ">$defline\n";
        # Print sequence in lines of $length
        for ( my $pos = 0 ; $pos < $seqlen ; $pos += $length ) {
                print substr($sequence, $pos, $length), "\n";
        }
}

sub dump_remainder {

        my ($seq_id, $desc, $pos, $seq) = @_;

        my $segment_defline = $seq_id . "_" . ($pos+1) . $desc;
        &print_sequence($segment_defline, $seq, 60);

        # Note that in cases where the remaining sequence is less than the
        # size of our segments, we still may have a little bit more to output
        # because we have *overlapping* windows.  For example, if there were
        # 87 bases left and our segment size was 100 with overlap of 50, we
        # would still need to print out the last 87 - (100 - 50) = 37 bp.
        if (length($seq) > $increment) { 
            $pos += $increment;
            $segment_defline = $seq_id . "_" . ($pos+1) . $desc;
            &print_sequence($segment_defline, substr($seq, $increment), 60);
        }
}