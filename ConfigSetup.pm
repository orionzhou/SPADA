package ConfigSetup;
use strict; 
use File::Path qw/make_path remove_tree/;
use File::Basename; 
use Cwd qw/abs_path/;
use Common; 
use Data::Dumper;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/config_setup config_setup_simple/; 

sub read_cfg_hash {
  my ($f_cfg) = @_;
  my $h;
  open(FH, "<$f_cfg") || die "config file $f_cfg is not there\n";
  while(<FH>) {
    chomp;
    next unless $_;
    next if /^\#/;
    $_ =~ s/\s//g;
    my ($k, $v) = split "=";
    while( $v =~ /\$\{(\w+)\}/g ) {
      die "no env variable named $1\n" unless exists $ENV{$1};
      my $rep = $ENV{$1};
      $v =~ s/\$\{$1\}/$rep/;
    }
    $h->{$k} = $v;
  }
  return $h;
}
sub config_setup_simple {
  my ($f_cfg, $dir_hmm) = @_;
 
  print "=====  setting up environment variables  =====\n";
  my $h = read_cfg_hash($f_cfg);
  for my $k (keys %$h) { $ENV{$k} = $h->{$k}; }

  $ENV{"SPADA_HMM_DIR"} = $dir_hmm if defined $dir_hmm;
  make_path($ENV{"SPADA_HMM_DIR"}) if ! -d $ENV{"SPADA_HMM_DIR"};
  $ENV{"TMP_DIR"} = $ENV{"SPADA_HMM_DIR"};

  my @keys = qw/ClustalO trimAl HMMER/;
  for my $key (@keys) {
    exists $ENV{$key} || die "$key not defined\n";
  }
} 

sub config_setup {
  my ($f_cfg, $dir, $dir_hmm, $f_fas, $f_gff, 
    $org, $sp, $e, $methods, $ncpu) = @_;
 
  print "=====  setting up environment variables  =====\n";
  my $h = read_cfg_hash($f_cfg);
  for my $k (keys %$h) { $ENV{$k} = $h->{$k}; }

  $ENV{"SPADA_SRC_DIR"} = dirname(__FILE__);
  $ENV{"SPADA_OUT_DIR"} = abs_path($dir);
  $ENV{"SPADA_HMM_DIR"} = abs_path($dir_hmm);
  $ENV{"SPADA_FAS"} = abs_path($f_fas);
  $ENV{"SPADA_GFF"} = abs_path($f_gff);
  $ENV{"SPADA_ORG"} = $org;
  $ENV{"evalue"} = $e;
  $ENV{"eval_sp"} = $sp;
  $ENV{"spada_methods"} = $methods;
  $ENV{"threads"} = $ncpu;

  my @keys = qw/ClustalO SignalP HMMER/;
  for my $key (@keys) {
    exists $ENV{$key} || die "$key not defined\n";
  }

  $dir = $ENV{"SPADA_OUT_DIR"};
  $dir_hmm = $ENV{"SPADA_HMM_DIR"};
  $f_fas = $ENV{"SPADA_FAS"};
  $f_gff = $ENV{"SPADA_GFF"};
  
  printf "  using %s matrix\n", $ENV{"SPADA_ORG"};

  # check HMM / profile directory
  -s "$dir_hmm/21_all.hmm" || die "$dir_hmm/21_all.hmm is not there\n";

  # check availability of called programs
  my $bin_rel_path = {
      "Augustus" => "bin/augustus",
      "GeneWise" => "bin/genewise",
      "SplicePredictor" => "bin/SplicePredictor",
      "GeneMark" => "gmhmme3",
      "GlimmerHMM" => "bin/glimmerhmm",
      "GeneID" => "bin/geneid"
  };

  my @method_ary = split(";", $ENV{"spada_methods"});
  for my $method (@method_ary) {
    my @softs;
    if($method eq "Augustus_evidence") {
      @softs = ("Augustus");
    } elsif($method eq "Augustus_de_novo") {
      @softs = ("Augustus");
    } elsif($method eq "GeneWise_SplicePredictor") {
      @softs = ("GeneWise", "SplicePredictor");
    } elsif($method eq "GeneMark") {
      @softs = ("GeneMark");
    } elsif($method eq "GlimmerHMM") {
      @softs = ("GlimmerHMM");
    } elsif($method eq "GeneID") {
      @softs = ("GeneID");
    }
    $ENV{$method} = join(";", @softs);
    for my $soft (@softs) {
      exists $bin_rel_path->{$soft} || die "no rel path for $soft\n";
      exists $ENV{$soft} || die "no dir-env path for $soft\n";
      my $fb = $ENV{$soft}."/".$bin_rel_path->{$soft};
      -x $fb || die "$soft: $fb is not there\n";
      my $softr = $soft."_bin";
      $ENV{$softr} = $bin_rel_path->{$soft};
    }
    printf "\twill run %s: %s\n", $method, $ENV{$method};
  }

  
  -d $ENV{"SPADA_OUT_DIR"} || make_path($ENV{"SPADA_OUT_DIR"});
  $ENV{"TMP_DIR"} = $ENV{"SPADA_OUT_DIR"};
  
  $ENV{"BLOSUM62"} = $ENV{"SPADA_SRC_DIR"}."/BLOSUM62";
  $ENV{"BLOSUM80"} = $ENV{"SPADA_SRC_DIR"}."/BLOSUM80";
  push @INC, $ENV{"SPADA_SRC_DIR"};
  $ENV{'PATH'} = join(":", $ENV{"SPADA_SRC_DIR"}, $ENV{'PATH'});
}

1;
__END__
