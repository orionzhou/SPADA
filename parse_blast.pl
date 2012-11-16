#!/usr/bin/perl -w
#
# Usage: parse_blast.pl <blast-output>+
#
# Description:
#  Quick script to use Bioperl to parse blast data for all
# the information we will need to remove redundancy
#
# NOTE1: Hits against oneself are removed!

use strict;
use Bio::SearchIO;

# print headers
print "Query ID\tHit ID\tPcnt ID\tPcnt Ovlp\tQuery Cov\tHit Cov\tQuery Len\tHit Len\tHsp Len\tMatches\n";

foreach my $file (@ARGV) {
        my $in = new Bio::SearchIO(-format => 'blast',
			       -file   => $file);

        while ( my $result = $in->next_result ) {

	my $query_name = $result->query_name;
	my $query_length = $result->query_length;

	while ( my $hit = $result->next_hit ) {

	    my $nmatches = $hit->matches('id');
	    my $hit_name = $hit->name;
	    my $hit_length = $hit->length;
#	    my $pcnt_id = &round_pcnt($hit->frac_identical);
	    my $pcnt_id = &round_pcnt($nmatches/
				      ($hit->length_aln+$hit->gaps));
	    my $query_cov = &round_pcnt($hit->frac_aligned_query);
	    my $hit_cov = &round_pcnt($hit->frac_aligned_hit);

	    while ( my $hsp = $hit->next_hsp ) {
		
		my $hsp_length = $hsp->hsp_length;
		my $pcnt_ovlp = &round_pcnt($hsp_length/
					    ($hsp_length +
					     $hit->num_unaligned_query +
					     $hit->num_unaligned_sbjct));

		if ($query_name ne $hit_name) {
		    print "$query_name\t$hit_name\t$pcnt_id\t$pcnt_ovlp\t" .
			"$query_cov\t$hit_cov\t$query_length\t$hit_length\t" .
			"$hsp_length\t$nmatches\n";
		}
	    }
	}
        }
}

exit 0;

sub round_pcnt {
        my $frac = shift @_;

        return int(100*$frac + 0.5);
}
