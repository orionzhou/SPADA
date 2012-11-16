#!/usr/bin/perl -w
#
# Usage: screen_out_known_hits.pl <master-hit-file> <genome-hits-file>
#
# Description:
#    Quick hack script to screen out genomic hits in the parse_HMM_hits.pl
# generated output file <genome-hits-file> that are already in
# our <master-hit-file> which was generated with the script
# create_master_rice_list.pl.
# 
# This script also screens out hits to
#   o CRP2830 with E-value > 1e-05
#   o regions with an internal stop codon with E-value > 0.001
#
# Note: the code as written is only applicable to searches done on
# the 2000 bp pieces of the genomic sequence!

use strict;

my %chr_locs;  # $chr_hits{$chr}{$pos} = "";
my ($master_hit_file, $genome_hits_file) = @ARGV;
my $frag_spacing = 1000;
my $win_size = 2000;

#
# 1. Load in the known sequence locations from the master hit file
#
open (FILE, $master_hit_file) || die "Can't open $master_hit_file: $!\n";
while (<FILE>) {
        if (/^Pub LOC/) { next;}  # skip header line

        my ($id, $grp, $evalue, $chr, $orient, $start, $stop, @rest) = split /\t/;
        my $rounded_start = (int($start/$frag_spacing)*$frag_spacing)+1;
        my $rounded_stop = (int($stop/$frag_spacing)*$frag_spacing)+1;

        # mark the bin plus surrounding bins as seen
        # Note that some windows may have multiple entries whose
        # storage is attempted in duplicate, and hence
        # we need a hash to store this info.
        # print STDERR "Marking Chr $chr pos $rounded_start\n"; #### DEBUG
        $chr_locs{$chr}{$rounded_start}{$orient}{"$start:$stop"} = '';
        $chr_locs{$chr}{$rounded_start-$frag_spacing}{$orient}{"$start:$stop"} = '';
        $chr_locs{$chr}{$rounded_start+$frag_spacing}{$orient}{"$start:$stop"} = '';

        # print STDERR "Marking Chr $chr pos $rounded_stop\n"; #### DEBUG
        $chr_locs{$chr}{$rounded_stop}{$orient}{"$start:$stop"} = '';
        $chr_locs{$chr}{$rounded_stop-$frag_spacing}{$orient}{"$start:$stop"} = '';
        $chr_locs{$chr}{$rounded_stop+$frag_spacing}{$orient}{"$start:$stop"} = '';
}
close(FILE);

#
# 2. Open the parsed HMM search output file and print out any novel hits
#    (Assuming that distinct real hits will not overlap unless they are
#    in opposite orientations!)

my $printOn = 'true';
open (FILE, $genome_hits_file) || die "Can't open $genome_hits_file: $!\n";
while (<FILE>) {
        if (/^chr0*(\d+)\S*_(\d+)\s+(CRP\d+)\s+(\S+)\s+RF ([-]?\d)/) {
	my ($chr, $window_pos, $grp, $evalue, $frame) = ($1, $2, $3, $4, $5);

	# NOTE: CRP2830 has a lot of spurious hits especially with low signif
	if (($grp eq 'CRP2830') && ($evalue >= 1e-05)) {
	    $printOn = 'false';
	    next;
	}

	my $orient = ($frame < 0) ? '-' : '+';
	my $old_line = $_;

	# get the next (alignment query) line and check the start/stop pos
	$_ = <FILE>;
	if (/^\s*(\d+)\s+\S+\s+(\d+)\s*$/) {
	    my ($rel_aa_start, $rel_aa_end) = ($1, $2);
	    my $nt_start = &get_nt_pos($window_pos, $orient, $rel_aa_start);
	    my $nt_end = &get_nt_pos($window_pos, $orient, $rel_aa_end);
	    my $hit_mid_pos = 0.5*($nt_start + $nt_end);

#	    warn "NT CONVERSION: ($nt_start, $nt_end) MID: $hit_mid_pos\n";
#
#	    if (defined($chr_locs{$chr}{$window_pos}{$orient})) {
#		warn "DEFINED ($chr, $window_pos, $orient): " .
#		    join (" ", keys %{$chr_locs{$chr}{$window_pos}{$orient}}) 
#		    . "\n";
#	    }  ## DEBUG
	    
	    # if we see a known overlapping gene in the same window
	    # shut the print mode off until we encounter an unknown gene
	    if (defined($chr_locs{$chr}{$window_pos}{$orient}) &&
		&overlaps($hit_mid_pos, 
			  keys %{$chr_locs{$chr}{$window_pos}{$orient}})) {
		$printOn = 'false';
	    }
	    else {

		# NOTE: lets skip pseudogenes (those with stop codon,
		# 'X', in the middle of the hit) with poor evalues too
		# or we'll be here all year!!
		if ((/X/) && ($evalue >= 0.001)) {
		    $printOn = 'false';
		}
		else {
		    $printOn = 'true';
		    # Don't forget that we read an extra line to get to
		    # the position info, so print that old line!
		    chomp $old_line;
		    print "$old_line\t$nt_start\t$nt_end\n";
		}
	    }
	}
        }
        # print anything as long as printOn is true!
        if ($printOn eq 'true') { print;}
}
close (FILE);

exit 0;

sub get_nt_pos {
        my ($window_start, $orient, $rel_aa_pos) = @_;

        if ($orient eq '+') {
	return ($window_start + 3*$rel_aa_pos);
        }
        else {
	return ($window_start + $win_size - 3*$rel_aa_pos);
        }
}

sub overlaps {
        my ($pos, @range_list) = @_;

        foreach my $gene_coords (@range_list) {
	my ($start, $stop) = split /:/, $gene_coords;

	if ((($pos >= $start) && ($pos <= $stop)) ||  # within + orient gene
	    (($pos >= $stop) && ($pos <= $start))) {  # within - orient gene
	    return 1;
	}
        }

        # if we got here without returning, then we never found a set of
        # gene coordinates that overlapped our input position.
        return 0;
}
