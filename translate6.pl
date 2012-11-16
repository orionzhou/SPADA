#!/usr/bin/perl -w
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#           Copyright 2002, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#  Kevin Silverstein
#  Michelle Graham
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  translate6.pl - translate a nucleotide FASTA file in all six frames.

=head1 SYNOPSIS
  
  translate6.pl [-help] [-tab] <fasta-file|stdin|->

  Options:
      -help     brief help message
      -tab      output results in tabular format (seqID,frame,sequence)

=head1 DESCRIPTION

  This program reads any number of nucleotide fasta entries from files or 
  stdin, and outputs the protein translations in all six frames using the
  Standard genetic code.

=head1 OPTIONS

=over 6
  
=item B<-help>
  
  Print a usage summary.

=item B<-tab>

  Output the translations in tabular format (seqID, frame, sequence).
  By default, the output appears in standard FASTA format.

=item B<fasta-files>

  To read fasta entries from stdin, the user should specify either 
  'stdin' or '-' in place of the filename argument. For example, one 
  could uncompress a file and translate it as follows:
        gzip -cd ex.fas.gz | translate6.pl - | gzip -9 > ex.aa.fas.gz

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
my $help_flag;
my $tab_output = '';
my %options = (
    "help|h" => \$help_flag,
    "tab|t"  => \$tab_output
);

# global genetic code
my(%genetic_code) = (
        
        'TCA' => 'S',    # Serine
        'TCC' => 'S',    # Serine
        'TCG' => 'S',    # Serine
        'TCT' => 'S',    # Serine
        'TTC' => 'F',    # Phenylalanine
        'TTT' => 'F',    # Phenylalanine
        'TTA' => 'L',    # Leucine
        'TTG' => 'L',    # Leucine
        'TAC' => 'Y',    # Tyrosine
        'TAT' => 'Y',    # Tyrosine
        'TAA' => '*',    # Stop
        'TAG' => '*',    # Stop
        'TGC' => 'C',    # Cysteine
        'TGT' => 'C',    # Cysteine
        'TGA' => '*',    # Stop
        'TGG' => 'W',    # Tryptophan
        'CTA' => 'L',    # Leucine
        'CTC' => 'L',    # Leucine
        'CTG' => 'L',    # Leucine
        'CTT' => 'L',    # Leucine
        'CCA' => 'P',    # Proline
        'CCC' => 'P',    # Proline
        'CCG' => 'P',    # Proline
        'CCT' => 'P',    # Proline
        'CAC' => 'H',    # Histidine
        'CAT' => 'H',    # Histidine
        'CAA' => 'Q',    # Glutamine
        'CAG' => 'Q',    # Glutamine
        'CGA' => 'R',    # Arginine
        'CGC' => 'R',    # Arginine
        'CGG' => 'R',    # Arginine
        'CGT' => 'R',    # Arginine
        'ATA' => 'I',    # Isoleucine
        'ATC' => 'I',    # Isoleucine
        'ATT' => 'I',    # Isoleucine
        'ATG' => 'M',    # Methionine
        'ACA' => 'T',    # Threonine
        'ACC' => 'T',    # Threonine
        'ACG' => 'T',    # Threonine
        'ACT' => 'T',    # Threonine
        'AAC' => 'N',    # Asparagine
        'AAT' => 'N',    # Asparagine
        'AAA' => 'K',    # Lysine
        'AAG' => 'K',    # Lysine
        'AGC' => 'S',    # Serine
        'AGT' => 'S',    # Serine
        'AGA' => 'R',    # Arginine
        'AGG' => 'R',    # Arginine
        'GTA' => 'V',    # Valine
        'GTC' => 'V',    # Valine
        'GTG' => 'V',    # Valine
        'GTT' => 'V',    # Valine
        'GCA' => 'A',    # Alanine
        'GCC' => 'A',    # Alanine
        'GCG' => 'A',    # Alanine
        'GCT' => 'A',    # Alanine
        'GAC' => 'D',    # Aspartic Acid
        'GAT' => 'D',    # Aspartic Acid
        'GAA' => 'E',    # Glutamic Acid
        'GAG' => 'E',    # Glutamic Acid
        'GGA' => 'G',    # Glycine
        'GGC' => 'G',    # Glycine
        'GGG' => 'G',    # Glycine
        'GGT' => 'G',    # Glycine
        'NAA' => 'X',    
        'NAT' => 'X',    
        'NAC' => 'X',    
        'NAG' => 'X',    
        'NTA' => 'X',    
        'NTT' => 'X',    
        'NTC' => 'X',    
        'NTG' => 'X',    
        'NGA' => 'X',    
        'NGT' => 'X',    
        'NGC' => 'X',    
        'NGG' => 'X',    
        'NCA' => 'X',    
        'NCT' => 'X',    
        'NCC' => 'X',    
        'NCG' => 'X',    
        'ANA' => 'X',    
        'ANT' => 'X',    
        'ANG' => 'X',    
        'ANC' => 'X',    
        'TNA' => 'X',    
        'TNT' => 'X',    
        'TNG' => 'X',    
        'TNC' => 'X',    
        'GNA' => 'X',    
        'GNT' => 'X',    
        'GNG' => 'X',   
        'GNC' => 'X',    
        'CNA' => 'X',  
        'CNT' => 'X',   
        'CNG' => 'X',    
        'CNC' => 'X',    
        'AAN' => 'X',
        'ATN' => 'I', # in 3 of 4 codes for I otherwise M-start
        'AGN' => 'X',
        'ACN' => 'T',
        'TAN' => 'X',
        'TTN' => 'X',
        'TGN' => 'X',
        'TCN' => 'S',
        'GAN' => 'X',
        'GTN' => 'V',
        'GGN' => 'G',
        'GCN' => 'A',
        'CAN' => 'X',
        'CTN' => 'L',
        'CGN' => 'R',
        'CCN' => 'P',
        'ANN' => 'X',	
        'TNN' => 'X',	
        'GNN' => 'X',
        'CNN' => 'X',
        'NAN' => 'X',
        'NTN' => 'X',
        'NGN' => 'X',
        'NCN' => 'X',
        'NNA' => 'X',
        'NNT' => 'X',
        'NNG' => 'X',
        'NNC' => 'X',
        'NNN' => 'X',
);


#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag;


($infile)= @ARGV;
if (!$infile) {
    pod2usage(2);
} elsif (($infile eq '-') || ($infile eq 'stdin')) {
        $fhandle = \*STDIN;
} else {
        open ($fhandle, $infile) || die "Cannot open file $infile \n";
}

while (<$fhandle>) {
        if (/^>\s*(\S+)(.*)/) {
	$old_seq_id = $seq_id;
                $old_seq_desc = $seq_desc;
	($seq_id, $seq_desc) = ($1, $2);
	
	if ($wait_4_1st_entry eq 'true') {
	    $wait_4_1st_entry = 'false';
	}
	else { 
	    &print_3_frames($old_seq_id, "+", $old_seq_desc, $seq, $tab_output);
	    &print_3_frames($old_seq_id, "-", $old_seq_desc, &revcom($seq),
                                                        $tab_output);
	    $seq='';
	}    
        }
        elsif (!(/^\s*$/i)) {
	chomp;
	$seq .= $_;
        }
}
#print last record
&print_3_frames($seq_id, "+", $seq_desc, $seq, $tab_output);
&print_3_frames($seq_id, "-", $seq_desc, &revcom($seq),
                                $tab_output);

close ($fhandle);

exit 0;

###############################################################################
#subroutines
#############################################################################

sub dna2peptide {

        my($dna) = @_;
        my $seqlen = length($dna);

        # Initialize variables
        my $protein = '';

        # Translate each three-base codon to an amino acid, and append to a protein 
        for(my $i=0; $i < ($seqlen - 2) ; $i += 3) {
                $protein .= codon2aa( substr($dna,$i,3) );
        }

        return $protein;
}


sub codon2aa {
        my($codon) = @_;

        $codon = uc $codon;

        if(exists $genetic_code{$codon}) {
                return $genetic_code{$codon};
        }
        else { 
#	print STDERR "Bad codon \"$codon\"!!\n";
                return 'X'
        }
}


sub translate_frame {

        my($dna, $frame) = @_;

        die "unknown frame: $frame\n" unless $frame =~ /^[012]$/;
        my $protein;
        my $len = length($dna);

        # Finally, calculate and return the translation
        if ($len > 2) {
	return dna2peptide ( substr( $dna, $frame ) );
        }
        else {
	return '';
        }
}


sub print_sequence {

        my($sequence, $length) = @_;
        my $seqlen = length($sequence);

        # Print sequence in lines of $length
        for ( my $pos = 0 ; $pos < $seqlen ; $pos += $length ) {
                print substr($sequence, $pos, $length), "\n";
        }
}

sub print_3_frames {

        my($seq_id, $sense, $seq_desc, $dna, $tab_flag) = @_;
        my $line_width = 60;
        my $frame;
        
        if ($dna eq '') {
                print STDERR "No sequence for $seq_id.  Skipping...\n";
            }
        else {
            foreach $frame (0 .. 2) {
	if ($tab_flag) {
	  print "$seq_id\t$sense$frame\t", &translate_frame($dna,$frame), "\n";
	}
	else { # output in FASTA format
	  print ">$seq_id\_$sense\_$frame $seq_desc\n";
	  print_sequence(translate_frame($dna, $frame), $line_width);
	}
            }
        }
}    

sub revcom {

        my($dna) = @_;

        # First reverse the sequence
        my $revcom = reverse($dna);

        # Next, complement the sequence, dealing with upper and lower case
        # A->T, T->A, C->G, G->C
        $revcom =~ tr/ACGTacgt/TGCAtgca/;

        return $revcom;
}
