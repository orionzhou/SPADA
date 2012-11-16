#!/usr/bin/perl
#
#------------------------------------------------------------------------------
#                      University of Minnesota
#           Copyright 2002, Regents of the University of Minnesota
#------------------------------------------------------------------------------
#
# File:
#  verifyHMMmatch.pl <reference-fasta> <HMMsearch_output>+
#
# Description:
#  HMMer has an undesirable feature when searching for motifs in translated
#  nucleotide sequences -- it completely ignores '*' characters.  Thus, 
#  apparent hits (e.g., vs. Arabidopsis) may span stop codons.  This hack
#  script analyzes all one-line alignments in the HMMsearch output and
#  sees if the search string appears WITHOUT STOP CODONS in the intended
#  translated target.
#

use strict;

my %regex_hash = ();
my $fasta_file = shift @ARGV;
my $hmmout_file;

my $header;
my $seq = "";
my $seq_count = 0;
my $verbose = 0;
my $debug = 0;
my $seq_id = "";

foreach $hmmout_file (@ARGV) {
        open (FILE, $hmmout_file) || die ("Can't open $hmmout_file: $!\n");
        while (<FILE>) {
	if (/^(\S+): domain \d+ of \d+, from \d+ to \d+: score/) {
	    $seq_id = $1;
	}
	if (/^\s*\S+\s+\d+\s+([A-Za-z\-]+)\s+\d+\s*$/) {
	    my $motif_regexp = $1;
	    $motif_regexp =~ s/\-//g;
	    if ($debug) {print STDERR "Recording $seq_id -> $motif_regexp\n";}
	    if (defined($regex_hash{$seq_id})) {
		$regex_hash{$seq_id} .= "|$motif_regexp";
	    }
	    else {
		$regex_hash{$seq_id} = $motif_regexp;
	    }
	    $seq_id = "";
	}
        }
        close(FILE);
}

open (FILE, $fasta_file) || die ("Can't open $fasta_file: $!\n");
while (<FILE>) {
        if (/^>/) {
            if ($seq_count > 0) {
	&processSeq($header, $seq, \%regex_hash);
            }
            chomp;
            $header = $_;
            $seq = "";
            $seq_count++;
        }
        else {
            $seq .= $_;
            chomp $seq;
        }
}
close(FILE);

# process the last one.
&processSeq($header, $seq, \%regex_hash);

exit 0;

sub processSeq {
    my ($header, $seq, $motif_href) = @_;

    my (@orfs, $orf, $orf_size);
    my $start_pos = 0;
    my ($seqid) = $header =~ /^>\s*(\S+)/;

    if (defined($$motif_href{$seqid})) {
        my $motif = $$motif_href{$seqid};
        if ($verbose) {print STDERR "Searching for motif $motif in $header\n";}

        if ($seq =~ /\*/) {
            @orfs = split /\*/, $seq;
        }
        else {
            push @orfs, $seq;
        }

        foreach $orf (@orfs) {
            $orf_size = length($orf);
            if ($orf =~ /$motif/i) {
                print "$header";
                if ($header =~ /RF ([-]?\d)/) {
                    my $frame = $1;
                    my @match_bounds = &get_nt_match_bounds($seq, $orf, $start_pos, 
						$motif, $frame);
                    print " ($match_bounds[0],$match_bounds[1])";
                }
                print "\n$orf\n";
            }
            $start_pos += $orf_size+1;
        }
    }
}

sub get_nt_match_bounds {
    my ($full_seq, $orf, $start_pos, $motif, $frame) = @_;

    if ($orf =~ /($motif)/i) {
        my $start = length($`)+1+$start_pos;  # amino acid start position
        my $end = $start + length($1) - 1;        # amino acid end position
#    print STDERR "Frame: $frame, orf start: $start_pos, aa start: $start, aa end: $end\n";
        if ($frame > 0) {
            return(3*($start-1)+$frame-1, 3*$end+$frame-1);
        }
        else {
            my $prot_seq_len = length($full_seq);
            return(3*($prot_seq_len-$start+1)-$frame, 
	     3*($prot_seq_len-$end)-$frame);
        }
    }
    else {
        if ($debug) { print STDERR "Can't find motif $motif in ORF $orf\n";}
        return(-1,-1)
    }
}

