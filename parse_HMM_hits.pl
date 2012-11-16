#!/usr/local/bin/perl -w
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#           Copyright 2006, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#  Kevin Silverstein
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  parse_HMM_hits.pl - parse HMMer output to extract out CRP hits

=head1 SYNOPSIS
  
  parse_HMM_hits.pl [-help] [-aln] [-E evalue-cutoff]
                                      [-hitpat hit-regex] [-querypat query-regex] 
                                      [-top ] <HMMoutput>

  Options:
      -help     Brief help message
      -aln      Report alignment
      -E        Evalue cutoff above which hits are ignored
      -hitpat   Regular expr pattern for core ID extraction from the hit
      -querypat Regular expr pattern for core ID extraction from the query
      -top      Report top-hit only

=head1 DESCRIPTION

  This program reads the output of hmmpfam or hmmsearch and outputs 
  a tab-delim spreadsheet of top hits with the following information:

  Query-ID Hit-ID E-value Query-Desc Hit-Desc

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<-aln>
  
  Print out the alignment for each hit.
  By default, no alignment is printed.

=item B<-E evalue-cutoff>

  If set, only output hits that are above the specified cutoff.
  By default, all hits are reported.

=item B<-hitpat hit-regex>

  A standard regular expression to use to extract a sub-portion 
  of each hit ID.  For example, 
        -hitpat '^(CRP\d{4}).trim$'
  will extract CRP3670 from the hit CRP3670.trim.
  By default, the entire hit ID will be reported.

=item B<-querypat query-regex>

  A standard regular expression to use to extract a sub-portion
  of each query ID.  For example, 
        -querypat '^(LOC_Os\d{2}g\d+)\S+$'
  will extract LOC_Os01g01010 from the query ID 
  LOC_Os01g01010.1|11971.m06748|protein
  By default, the entire query ID will be reported.

=item B<-top>

  Only print the top hit.
  By default, information on all hits satisfying the E-value cutoff
  (if used) will be reported.

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
use Bio::SearchIO;

my $help_flag;
my $debug_flag;
my $align_flag;
my $evalue_cutoff;
my $hit_regex;
my $query_regex;
my $top_hit_only;
my %options = (
                              "help|h"     => \$help_flag,
                              "debug|d"    => \$debug_flag,
	       "aln|a"      => \$align_flag,
	       "E=s"        => \$evalue_cutoff,
	       "hitpat=s"   => \$hit_regex,
	       "querypat=s" => \$query_regex,
	       "top|t"      => \$top_hit_only
	      );

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;

my $hmmoutfile = shift @ARGV;

my $search_data = new Bio::SearchIO(-format => 'hmmer',
				    -file   => $hmmoutfile);

# Print header info
print "Query ID\tHit ID\tE-value\tQuery Desc\tHit Desc\n";

my $query_id;
my $query_desc;
my $hit_id;
my $hit_desc;
my $evalue;

# Loop through each Query result
# Note: $result is a Bio::Search::Result::HMMERResult object
while (my $result = $search_data->next_result) {

        $query_id = $result->query_name();
        if ($query_regex) {
	$query_id =~ s/$query_regex/$1/;
        }
        $query_desc = $result->query_description();
        if ($debug_flag) { warn "Query: $query_id $query_desc\n";}  ## DEBUG

        #
        # Loop through each hit to the query
        #
        # NOTE: bioPerl returns hits in the order of the *alignments*
        #  returned by HMMer (which are pretty random), NOT in the ranked
        #  order of descreasing significance!!!  So we will have to be
        #  very careful with 'top-hit' requests!

        my $best_evalue = 100;
        my $top_hit_result = '';
        while (my $hit = $result->next_hit ) {

	$hit_id = $hit->name();
	if ($hit_regex) {
	    $hit_id =~ s/$hit_regex/$1/;
	}
	$hit_desc = $hit->description();
	if ($debug_flag) { warn "Hit: $hit_id $hit_desc\n";}  ## DEBUG

	$evalue = $hit->significance();
	if (!$evalue_cutoff || ($evalue < $evalue_cutoff)) {
	    my $output_line = "$query_id\t$hit_id\t$evalue\t$query_desc\t$hit_desc\n";

	    # If we are only looking for the top hit, 
	    # check to make sure this hit has a better
	    # evalue than what's on record.  If so,
	    # store the info in $top_hit_result
	    if ($top_hit_only) {
		if ($evalue < $best_evalue) {
		    $top_hit_result = $output_line;
		    $best_evalue = $evalue;
		}
		else {
		    next;
		}
	    }
	    else {
		print $output_line;
	    }

	    if ($align_flag) {

		while (my $hsp = $hit->next_hsp) {
		    # print alignment the hard way since bioperl doesn't
		    # support obtaining a SimpleAlign object for HMMER.

		    my $align_output = &format_align($hsp);
		    if ($top_hit_only) {
			$top_hit_result .= $align_output;
		    }
		    else {
			print $align_output;
		    }
		}

	    }
	}

        }
        if ($top_hit_only && ($top_hit_result ne '')) {
	print $top_hit_result;
        }
}

exit 0;

sub format_align {
        my $hsp = shift @_;

        my $align_str = sprintf "%10d ", $hsp->start('query');

        $align_str .= $hsp->query_string . " " . $hsp->end('query') . "\n";
        $align_str .= (" " x 11) . $hsp->homology_string . "\n";
        $align_str .= sprintf "%10d ", $hsp->start('hit');
        $align_str .= $hsp->hit_string . " " . $hsp->end('hit') . "\n";

        return $align_str;
}
