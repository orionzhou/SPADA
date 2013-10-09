#!/usr/bin/perl -w
#
# POD documentation
#------------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  gff2gtb.pl - convert a Gff file to Gtb file and validates

=head1 SYNOPSIS
  
  gff2gtb.pl [-help] [-in input-file] [-seq refseq-file] [-out output-file]

  Options:
      -help   brief help message
      -in     input file
      -out    output file
      -seq    refseq file

=head1 DESCRIPTION

  This program converts an input Gff file to an output Gtb file

=cut
  
#### END of POD documentation.
#-----------------------------------------------------------------------------
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Common;
use Gff;
use Gtb;

my ($fi, $fo, $fs) = ('') x 3;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "in|i=s"  => \$fi,
    "out|o=s" => \$fo,
    "seq|s=s" => \$fs,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo || !$fs;

gff2gtb($fi, "$fo.tmp");
gtb_fix_phase("$fo.tmp", $fo, $fs);
runCmd("rm -f $fo.tmp", 0);

__END__
