#!/usr/bin/perl
#
#------------------------------------------------------------------------------
#                      University of Minnesota
#            Academic Health Center Computational Biology Centers
#           Copyright 2001, Regents of the University of Minnesota
#------------------------------------------------------------------------------
#
# File:
#  pruneFASTA.pl <FASTA-Repository> <column> <exclude_flag> <filenames>+
#
# Description:
#  A quick hack script to extract a selected set of FASTA entries from 
#  a large FASTA repository file <FASTA-Repository> from the files <filenames>.
#  The <column> of the tab-separated <filename> will be searched for
#  identifiers, and the FASTA entries with this identifier will be extracted 
#  to standard out.
#  if <exclude_flag> is T, the FASTA entries for identifiers found in 
#  <filename> will be excluded from the output.  Otherwise, the identifiers
#  will serve as a list of entries to *include*.
#
# Author:
#  Kevin Silverstein          (ksilvers@cbs.umn.edu)  2001-10-26
#

$DEBUG = 2; #### 0-off, 1-low, 2-medium, 3-high, levels of debug messages

($fasta_source, $column, $exclude_flag, @input_files) = @ARGV;

# First make a hash to hold all seq_ids we want
%target_seq_hash = ();

# March through the input files and populate the seq_id hash
foreach $fname (@input_files) {
        if ($DEBUG >= 2) {print STDERR "Parsing $fname for target seq_ids...\n";} #### DEBUG
        open(INFILE, $fname) || die "Can't open $fname: $!\n";

        while (<INFILE>) {
	$input_line = $_;
	@fields = split /\t/;
	$new_id = $fields[$column];
	if ($DEBUG >= 3) {print STDERR "Encountered $new_id.\n";} #### DEBUG
	$new_id =~ s/\"//g;
	$new_id =~ s/^\s+//;
	$new_id =~ s/\s+$//;
	if ($DEBUG >= 2) {print STDERR "Storing $new_id in the hash.\n";} #### DEBUG
	$target_seq_hash{$new_id} = "";
        }

        close(INFILE);
}

# Now walk through the FASTA source file and pluck out the ones in the hash.
if ($DEBUG >= 2) {print STDERR "Opening the source FASTA file $fasta_source...\n"} #### DEBUG

if ($fasta_source =~ /\.(gz|Z)$/) {
        open(FASTAIN, "gzip -cd $fasta_source |") || die "Can't uncompress and read $fasta_source: $!\n";
}
else {
        open(FASTAIN, $fasta_source) || die "Can't open $fasta_source: $!\n";
}

$printStatus = "true";
while (<FASTAIN>) {
        if (/^>\s*/) {
	# store all word-terms in the def-line
	@possible_ids = split /[\s*\b\|\,\;]+/, $';
#	@possible_ids = split /[\s*\b\|\,\.\;]+/, $';
#	@possible_ids = split /\s+/, $';

	# If no id is actually found in the whole def-line, the following
	# defaults will remain in place!
	if ($exclude_flag =~ /T/i) {
	    $printStatus = "true";
	}
	else {
	    $printStatus = "false";
	}

	# Now search among the terms found in the def-line for a recognized id
	foreach $seq_id (@possible_ids) {

	    if ($DEBUG >= 3) {print STDERR "Encountered $seq_id.\n";} #### DEBUG

	    if (defined($target_seq_hash{$seq_id})) {
		delete($target_seq_hash{$seq_id});
		if ($exclude_flag =~ /T/i) {
		    $printStatus = "false";
		    if ($DEBUG >= 2) {print STDERR "Excluding $seq_id\n";} #### DEBUG
		}
		else {
		    $printStatus = "true";
		    if ($DEBUG >= 2) {print STDERR "Outputting $seq_id\n";} #### DEBUG
		    print;
		}
	    }
	}

	if (($exclude_flag =~ /T/i) && ($printStatus eq 'true')) {
	    print;
	}
        }
        elsif ($printStatus eq "true") {
	print;
        }
}
close(FASTAIN);

# Now see if any identifiers have been missed.

if ((scalar keys %target_seq_hash) > 0) {
        print STDERR "WARNING: the following FASTA entries could not be found:\n";
        foreach $seq_id (sort keys %target_seq_hash) {
	print STDERR "$seq_id\n";
        }
}
elsif ($DEBUG >= 2){print STDERR "Program has successfully terminated.\n";} #### DEBUG

exit 0;

	    
