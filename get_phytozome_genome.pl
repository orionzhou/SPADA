#!/usr/bin/perl -w
=pod BEGIN
  
=head1 NAME
  
  get_phytozome_genome.pl - download genome information (sequence & annotation) from Phytozome 

=head1 SYNOPSIS
  
  get_phytozome_genome.pl [-help] -g <genome> -o <out directory>

  Options:
      -help    brief help message
      -genome  genome to download 
      -outdir  directory to put genome files in 

=head1 DESCRIPTION

  This program downloads genome files (sequence & gene annotation) from Phytozome, 
  put them in a local directory.

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw/make_path/;

my $phytozome_genome_id = {
    'Alyrata' => 107,
    'Athaliana' => 167,
    'Cpapaya' => 113,
    'Gmax' => 109,
    'Mtruncatula' => 135,
    'Osativa' => 193,
    'Sbicolor' => 79,
    'Vvinifera' => 145,
    'Zmays' => 181,
};
my $help_flag;
my ($genome, $f_seq, $f_gff);
my %options = (
    "help|h" => \$help_flag,
    "genome|g=s" => \$genome,
    "seq|s=s" => \$f_seq,
    "gff|f=s" => \$f_gff
);

#----------------------------------- MAIN -----------------------------------#
GetOptions(%options) or pod2usage(2);
pod2usage(1) if $help_flag || !$genome;
$f_seq ||= "01_refseq.fa";
$f_gff ||= "51_gene.gff";
pod2usage({-message => "Genome ($genome) not supported yet\n", -exitval => 2}) if ! exists $phytozome_genome_id->{$genome};

my $id = $phytozome_genome_id->{$genome};
download_genome($genome, $id, $f_seq, $f_gff);

sub download_genome {
    my ($genome, $id, $f_seq, $f_gff) = @_;
    my $url_pat1 = "ftp://ftp.jgi-psf.org/pub/JGI_data/phytozome/v8.0/%s/assembly/%s_%d.fa.gz";
    my $url_pat2 = "ftp://ftp.jgi-psf.org/pub/JGI_data/phytozome/v8.0/%s/annotation/%s_%d_gene.gff3.gz";
    my $dir_seq = dirname($f_seq);
    my $dir_gff = dirname($f_gff);
    make_path($dir_seq) unless -d $dir_seq;
    make_path($dir_gff) unless -d $dir_gff;
    my $url_fas = sprintf $url_pat1, $genome, $genome, $id;
    my $url_gff = sprintf $url_pat2, $genome, $genome, $id;
    system("wget $url_fas -O $f_seq.gz");
    system("gunzip -f $f_seq.gz");
    system("wget $url_gff -O $f_gff.gz");
    system("gunzip -f $f_gff.gz");
}



