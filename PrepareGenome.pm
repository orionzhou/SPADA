package PrepareGenome;
use strict; 
use File::Compare;
use File::Path qw/make_path remove_tree/;
use Common; 
use Data::Dumper;
use Seq;
use Gff;
use Gtb;
use Log::Log4perl;
use List::Util qw/min max sum/;

use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT_OK = qw//;
@EXPORT = qw/pipe_prepare_genome/; 

sub get_genome_orf {
    my ($fi, $fo, $cutoff_missing) = @_;
    my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
    my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');
    $cutoff_missing ||= 0.5;
    
    my $log = Log::Log4perl->get_logger("PrepareGenome");
    $log->info("extracting ORFs from translated genomic sequence");
    while(my $seq = $seqHI->next_seq()) {
        my $seqStr = $seq->seq;
        my ($id, $beg, $end, $srd) = ($seq->id, 1, 3*length($seqStr), "+");
        if($seq->id =~ /^(\w+)\_(\d+)\_(\d+)\_([\+\-])$/) {
            ($id, $beg, $end, $srd) = ($1, $2, $3, $4);
        }
        while( $seqStr =~ /([^\*]{15,})/g ) {
            my ($begL, $endL) = ($-[1]+1, $+[1]);
            my ($begG, $endG);
            if($srd eq "-") {
                $begG = $end - $endL*3 + 1;
                $endG = $end - ($begL*3-2) + 1;
            } else {
                $begG = $beg + ($begL*3-2) - 1;
                $endG = $beg + $endL*3 - 1;
            }
            my $n_x =()= $1 =~ /X/gi;
            if($n_x / ($endL-$begL+1) <= $cutoff_missing) {
                my $seqObj = Bio::Seq->new(-id=>join("_", $id, "$begG-$endG", $srd, "x"), -seq=>$1);
                $seqHO->write_seq($seqObj);
            }
        }
    }
    $seqHI->close();
    $seqHO->close();
}
sub get_protein_orf {
    my ($f_gtb, $fo, $f_seq) = @_;
    my $t = readTable(-in=>$f_gtb, -header=>1);
    my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');
    
    my $log = Log::Log4perl->get_logger("PrepareGenome");
    $log->info("extracting ORFs from predicted protein sequence");

    my $h = {};
    for my $i (0..$t->nofRow-1) {
        my ($idM, $idG, $chr, $srd, $phaseS, $locS, $cat1, $cat2) = 
            map {$t->elm($i, $_)} qw/id parent chr strand phase locC cat1 cat2/;
        next if $cat2 ne "mRNA";
        die "no locCDS for $idM\n" unless $locS;
        my $loc = locStr2Ary($locS);

        $loc = cropLoc_cds($loc, $srd, $phaseS);
        my $id = join("_", $chr, locAry2Str($loc), $srd, "p");
        next if exists $h->{$id};
        $h->{$id} = 1;

        my $seqStr = seqRet($loc, $chr, $srd, $f_seq);
        my $seq_cds = Bio::Seq->new(-id=>$id, -seq=>$seqStr);
        my $seq_pro = $seq_cds->translate();
        $seqHO->write_seq($seq_pro);
        printf "  %5d / %5d done\r", $i+1, $t->nofRow;
    }
    print "\n";
    $seqHO->close();
}

sub pipe_prepare_genome {
    my ($dir, $genome) = @_;
    make_path($dir) unless -d $dir;

    my $log = Log::Log4perl->get_logger("PrepareGenome");
    $log->error_die("Genome Seq file not there: $ENV{'SPADA_FAS'}") unless( -s $ENV{"SPADA_FAS"} );
    $log->info("___Preparing Genome Files___");

    my $f01 = "$dir/01_refseq.fa";
    if(-s $f01 && $f01 eq $ENV{"SPADA_FAS"}) {
        $log->info("Sequence FAS already in data directory");
    } else {
        $log->info("Copying Sequence FAS to data directory");
        system("cp -f $ENV{'SPADA_FAS'} $f01");
    }
    system("rm -rf $f01.index") if -s "$f01.index";

    my $f11 = "$dir/11_refseq_trans6.fa";
    translate6($f01, $f11);
    my $f12 = "$dir/12_orf_genome.fa";
    get_genome_orf($f11, $f12);

    my $f71 = "$dir/71_orf_protein.fa";
    if(!exists $ENV{"SPADA_GFF"}) {
        $log->info("Annotation GFF not defined");
        $log->warn("Proceeding without GFF");
    } elsif(! -s $ENV{"SPADA_GFF"}) {
        $log->warn("Annotation GFF not there: $ENV{'SPADA_GFF'}");
        $log->warn("Proceeding without GFF");
    } else {
        my $f51 = "$dir/51_gene.gff";
        if(-s $f51 && $f51 eq $ENV{"SPADA_GFF"}) {
            $log->info("Annotation GFF already in data directory");
        } else {
            $log->info("copying Annotation GFF to data directory");
            system("cp -f $ENV{'SPADA_GFF'} $f51");
        }
        my $f61 = "$dir/61_gene.gtb";
        gff2Gtb($f51, $f61);
        my $f62 = "$dir/62_gene.gff";
        gtb2Gff($f61, $f62);
        get_protein_orf($f61, $f71, $f01);
    }
}




1;
__END__
