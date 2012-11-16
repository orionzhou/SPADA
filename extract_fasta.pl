#!/usr/bin/perl -w
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#           Copyright 2006, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  Kevin Silverstein
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN

=head1 NAME

    extract_fasta.pl - extract or exclude a subset of fasta entries

=head1 SYNOPSIS

    extract_fasta.pl [-help] [-debug] [-col column] [-delim delimiter] 
                                                                        [-exclude] [-ids id_file] [-kwsep kw_regexp] [-multim] 
                                                                        [-nocase] [-outf outfile] [-range locol:hicol] 
                                                                        [-separate] [-verbose] [-xactid] <fasta-repository>+

        Options:
                -help       brief help message
                -col        column in the id_list spreadsheet to scan for IDs
                -delim      delimiter used to separate columns in the id_list spreadsheet
                -exclude    flag to exclude (prune out) seqs with IDs in the list
                -ids        a tab-delim spreadsheet with IDs to extract/exclude
                -kwsep      an alternative set of token separators for deflines
                -multim     flag to re-use terms/ids in the list for multiple matches
                -nocase     flag to eliminate case-sensitivity when matching IDs
                -outf       filename for the output fasta file
                -range      columns in the id_list spreadsheet to scan for seq ranges
                -separate   flag to put each output fasta sequence in its own file
                -verbose    output basic status messages
                -xactid     flag to disable defline tokenization and match IDs exactly

=head1 DESCRIPTION

    This flexible program extracts (or excludes) selected sequences
    from a master fasta sequence repository.  By default, the first
    occurrence of any fasta entry is extracted that has each specified 
    ID anywhere in its defline.  With the appropriate command-line
    options, you can exclude (i.e., prune out ) sequences on the list 
    instead.  You also have the flexibility to extract or exclude
    multiple fasta entries with the same ID or keyword burried in
    their deflines (e.g., all seqs with the keyword 'hydrogenase').
    Or you can force exact matches of IDs instead of scanning all
    tokens in a defline.

=head1 OPTIONS

=over 6

=item B<-help>

    Print a usage summary.

=item B<-debug>

    Print debugging information

=item B<-col column>

    Specifies which column in the input ID list to use for input
    IDs.  Columns follow the perl convention and thus start with
    0.  By default, column 0 is used.

=item B<-delim delimiter>

    Specifies the delimiter that separates columns in the input ID 
    list to use for input IDs.  
    By default, the delimiter is "\t".

=item B<-exclude>

    Sets the behavior to exclude or prune out all fasta entries
    on the list (i.e., print out all fasta entries *except* those
    listed).  By default, this flag is *not* set and only the 
    sequences on the list are extracted and printed out.

=item B<-ids id_file>

    Mandatory file specifying the IDs or terms to be extracted or 
    excluded.  The file could be as simple as a single ID per line,
    e.g.:

                TC10211
                TC11422

    or

                monooxygenase
                oxidoreductase

=item B<-kwsep kw_regexp>

    A regular expression to use to tokenize fasta definition lines
    when searching for matches to IDs or terms on the list.  By
    default, the following regular expression is used:
    '[\s*\b\|\,\.\;\(\)\:]+'

    This will separate the defline:

    >gi|85666109|ref|NC_001133.6| Saccharomyces cerevisiae chromosome I

    Into the terms:

    gi
    85666109
    ref
    NC_001133
    6
    Saccharomyces
    cerevisiae
    chromosome
    I

    each of which will be scanned against the input id/term list for
    possible matches.

    If you want to preserve the standard GenBank Accession.version
    as a single term, eliminate the '\.' from the separator list
    like this:

    -kwsep '[\s*\b\|\,\;\(\)\:]+'

    If you wish to keep entire IDs unaltered, it is better to use
    the -xactid option which overrides the default separator.

=item B<-multim>

    Enable multiple fasta entry matches for each input ID/term. This
    overrides the default behavior which removes each ID/term on the list
    after the first matching fasta/entry is found.

=item B<-nocase>

    Enable case-insensitivity when matching IDs/terms in the list to
    IDs/terms in the defline.  Default behavior when this flag is not
    set is to perform case-sensitive matches.

=item B<-outf outfile>

    Filename for the fasta file of extracted sequences.
    Default is STDOUT if no file is specified.

=item B<-range locol:hicol>

    Specifies which columns in the input ID list to use for sequence range
    low and high values, to allow subsequences to be output.  Columns 
    follow the perl convention and thus start with 0.  By default, 
    entire sequences are extracted (not seq ranges) if this option is
    not set.

=item B<-separate>

    Flag to output every fasta sequence in a file of its own in the
    current directory with filename = 'identifier.seq'.
    If a filename is also specified by the -outf option, then the
    string entered there will be used as a suffix.

=item B<-xactid>

    Flag to match IDs/terms exactly to the first non-whitespace term
    on the defline.  This overrides the default behavior to tokenize
    the defline and scan each token for matches.

    For example, using the switch -xactid will enable one to match
    the list entry 'gi|85666109|ref|NC_001133.6|' to the defline:

    >gi|85666109|ref|NC_001133.6| Saccharomyces cerevisiae chromosome I

=item B<fasta-repository>

    Any fasta-formated file(s) or STDIN stream to use as the sequence 
    repository from which to extract or exclude sequences.  Extracted
    sequences will be sent to STDOUT.

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
use Bio::SeqIO;

my $help_flag;
my $debug;
my $verbose;
my $id_col = 0;
my $delim = '\t';
my $range;
my $low_col;
my $high_col;
my $exclude_flag;
my $id_file;
my $kw_regexp = '[\s*\b\|\,\.\;\(\)\:]+';
my $multimatch_flag;
my $case_insensitive;
my $outfile;
my $out_suffix = ".seq";
my $sep_file_flag;
my $exactid_flag;
my %options = (
	       "help|h"     => \$help_flag,
	       "debug|d"    => \$debug,
	       "column|c=s" => \$id_col,
	       "delim=s"    => \$delim,
	       "exclude|e"  => \$exclude_flag,
	       "ids|i=s"    => \$id_file,
	       "kwsep|k=s"  => \$kw_regexp,
	       "multim|m"   => \$multimatch_flag,
	       "nocase|n"   => \$case_insensitive,
	       "outf|o=s"   => \$outfile,
	       "range=s"    => \$range,
	       "separate"   => \$sep_file_flag,
	       "verbose|v"  => \$verbose,
	       "xactid|x"   => \$exactid_flag
	       );

#
# VARIABLES
#

my %saved_ids = ();  # Hash storing IDs to find in deflines

#----------------------------------- MAIN -----------------------------------#
#
# Check the command-line options
#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;

pod2usage({-message => "A file with the desired search IDs/terms must be specified with the -ids option.\n", -exitval => 2}) if (!$id_file);

if ($sep_file_flag && $outfile) { # both were specified
        $out_suffix = $outfile;
}

if ($range) {
        pod2usage({-message => "-range option must be followed by positive integer values separated by a colon (low-column:high-column).\n", -exitval => 2}) if (!($range =~ /^\d+:\d+/));
        ($low_col, $high_col) = $range =~ /^(\d+):(\d+)/;
}

#
# 1. Read in IDs/terms from the ID file
#
if ($verbose || $debug) { warn "\nReading list of IDs from $id_file.\n";}

open (FILE, $id_file) || 
        die "Can't open file listing search IDs/terms $id_file: $!\n";
while (<FILE>) {

        if (/^\s*$/) { next; }  # skip blank lines

        if (/^\#/) { next; } # skip comment lines

        # get rid of leading whitespace which messes up column identification
        $_ =~ s/^\s*//;

        chomp;
        my @fields = split /$delim/;
        my $new_id = &clean_field($fields[$id_col]);

        if ($case_insensitive) {
                $new_id = uc($new_id);  # switch all IDs to upper case before storing
        }

        if ($range) {
#    my ($low, $high) = (&clean_field($fields[$low_col]), $fields[$high_col]));
                my $low = &min(&clean_field($fields[$low_col]), 
		   &clean_field($fields[$high_col]));
                my $high = &max(&clean_field($fields[$low_col]),
		    &clean_field($fields[$high_col]));
                $saved_ids{$new_id}{"$low:$high"} = 0;
                if ($verbose || $debug) {warn "Storing $new_id ($low:$high) in the hash.\n";}
        }
        else {
                $saved_ids{$new_id} = 0;
                if ($verbose || $debug) {warn "Storing $new_id in the hash.\n";}
        }

}
close(FILE);

#
# 2. Setup the sequence output file/stream for FASTA writing
#
my $seqout;
if (!$sep_file_flag) {
        if ($outfile) {
                $seqout = Bio::SeqIO->new(-file => "> $outfile", -format => 'Fasta');
        }
        else { # no file specified, so we will use STDOUT
                $seqout = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
        }
}

#
# 3. Now walk through the FASTA source input file(s)/stream and pluck out
#    the desired seqs.
#
my $seqin; # Bio::SeqIO Fasta Input stream

if ($#ARGV < 0) { # No fasta given, expect it from STDIN
        if ($verbose || $debug) { warn "\nReading input FASTA seqs from STDIN\n";}
        $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'Fasta');
        &output_selected_seqs($seqin);
}
else {
        foreach my $fasta_file (@ARGV) {

                if ($verbose || $debug) { 
                        warn "\nReading input FASTA seqs from $fasta_file\n";
                }

                my $fhandle;

                if ($fasta_file =~ /\.(gz|Z)$/) { # compressed file
                        # Couldn't get the following to work, so I opened my own file handle!
                        # $seqin = Bio::SeqIO(-file => "gzip -cd $fasta_file |", 
                        #     		  -format => 'Fasta');
                        open ($fhandle, "gzip -cd $fasta_file |") || 
	die "Can't uncompress $fasta_file: $!\n";
                }
                else {
                        # Couldn't get the above to work, so I open my own file handles!
                        # $seqin = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');
                        open ($fhandle, "$fasta_file") || 
	die "Can't uncompress $fasta_file: $!\n";
                }

                $seqin = Bio::SeqIO->new(-fh => \*$fhandle, -format => 'Fasta');

                &output_selected_seqs($seqin);
                close(FILE);
        }
}

#
# 3. Finally, print out the IDs/terms that weren't found (if any).
#
if ((scalar keys %saved_ids) > 0) {

        if ($multimatch_flag) {
                if ($range) {
                        warn "\nHere are the number of sequence ranges found with each ID/term:\n";
                }
                else {
                        warn "\nHere are the number of sequences found with each ID/term:\n";
                }
        }
        else {
                warn "\nWARNING: FASTA entries with the following IDs/terms could not be found:\n";
        }

        my @missing_items = ();

        foreach my $seqid (sort keys %saved_ids) {
                if ($range) {
                        my $found_missed_item = 'FALSE';
                        foreach my $range_pair (sort keys %{$saved_ids{$seqid}}) {
	warn "$seqid\t$range_pair\t$saved_ids{$seqid}{$range_pair}\n";
	if ($saved_ids{$seqid}{$range_pair} == 0) {
	  $found_missed_item = 'TRUE';
	}
                        }
                        if ($found_missed_item eq 'TRUE') {
	push @missing_items, $seqid;
                        }
                }
                else {
                        warn "$seqid\t$saved_ids{$seqid}\n";
                        if ($saved_ids{$seqid} == 0) {
	push @missing_items, $seqid;
                        }
                }
                
        }

        if (($multimatch_flag) && (scalar(@missing_items) > 0)) {
                warn "\nWARNING: FASTA entries with the following IDs/terms could not be found:\n";
                foreach my $item (@missing_items) {
                        warn "$item\n";
                }
        }
                
}
elsif ($verbose) {
        warn "\nProgram has successfully terminated.\n";
}

exit 0;

#------------------------------ SUBROUTINES --------------------------------#

sub clean_field {
        # clean of leading/trailing whitespace or surrounding "-chars from excel
        my ($input) = @_;

        $input =~ s/\"//g;
        $input =~ s/^\s+//;
        $input =~ s/\s+$//;

        return $input;
}

sub min {
        my @vals = @_;

        my $minval = $vals[0];

        for (my $i=0; $i<=$#vals; $i++) {
                if ($vals[$i] < $minval) {
                        $minval = $vals[$i];
                }
        }
        return $minval;
}

sub max {
        my @vals = @_;

        my $maxval = $vals[0];

        for (my $i=0; $i<=$#vals; $i++) {
                if ($vals[$i] > $maxval) {
                        $maxval = $vals[$i];
                }
        }
        return $maxval;
}

sub output_selected_seqs {
        my ($seqin) = @_;

        while (my $seqobj = $seqin->next_seq() ) {

                my $seq_matches = 'false';

                my $defline = $seqobj->display_id();

                my @range_list; # only used if a subseq range is desired.

                if (!$exactid_flag) {  # add on the rest of the defline if inexact
                        if ($seqobj->desc()) {
	$defline .= " " . $seqobj->desc();
                        }
                }

                if ($case_insensitive) {
                        $defline = uc($defline); # switch to upper case to match stored IDs
                }

                if ($debug) {warn "\nScanning $defline for a match.\n";}  ## DEBUG

                my @possible_ids = ();
                if ($exactid_flag) {
                        # store defline (ID) exactly as it appears
                        push @possible_ids, $defline;
                }
                else {
                        # store all word-terms in the defline
                        push @possible_ids, split /$kw_regexp/, $defline;
                }

                # Now search among the terms found in the def-line for a recognized id
                foreach my $seq_id (@possible_ids) {
                        if ($debug) {warn "Checking $seq_id\n";} ## DEBUG

                        if (defined($saved_ids{$seq_id})) {  # we found a match!
	$seq_matches = 'true';

	if ($range) {
	  @range_list = sort keys %{$saved_ids{$seq_id}};
	}

	if ($debug) {print STDERR "RANGE BEFORE: ", join(" ", @range_list), "\n";}

	if ($multimatch_flag) { 
	  if ($range) {
	    foreach my $range_pair (@range_list) {
	      $saved_ids{$seq_id}{$range_pair}++;
	    }
	  }
	  else {
	    $saved_ids{$seq_id}++;
	  }
	}
	else {
	  # only delete key if we don't want to match again!
	  delete($saved_ids{$seq_id});
	}

	if ($verbose || $debug) {
	  if ($exclude_flag) {
	    warn "Excluding seq matching $seq_id\n";
	  }
	  else {
	    warn "Outputting seq matching $seq_id\n";
	  }
	}

	last;  # don't bother to scan the remaining tokens for this defline
                        }
                }

                if (($exclude_flag && ($seq_matches eq 'false')) ||
	(!$exclude_flag && ($seq_matches eq 'true'))) {

                        if ($debug) {warn "Outputting sequence.\n";} ## DEBUG

                        if ($sep_file_flag) {
	my $fname = $seqobj->display_id() . $out_suffix;
	$seqout = Bio::SeqIO->new(-file => "> $fname", -format => 'Fasta');
                        }

                        if ($range) {

	if ($debug) {print STDERR "RANGE AFTER: ", join(" ", @range_list), "\n";}

	foreach my $range_pair (@range_list) {
	  my ($low, $high) = $range_pair =~ /^(\d+):(\d+)/;
	  my $new_id = $seqobj->display_id . "_$low" . "_$high";

	  if ($debug) {warn "Creating subseq $new_id\n";}

	  my $subseqobj = Bio::Seq->new( -display_id => $new_id,
					 -desc => $seqobj->desc(),
					 -seq => $seqobj->subseq($low, $high));
	  $seqout->write_seq($subseqobj);
	}
                        }
                        else {
	$seqout->write_seq($seqobj);
                        }
                }
        }
}
