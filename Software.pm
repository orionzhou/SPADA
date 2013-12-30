package Software;
use strict;
use File::Path qw/make_path remove_tree/;
use Bio::Seq;
use Bio::SeqIO;
use Common;
use Location; 
use Seq;
use Gtb;
use Data::Dumper;
use Log::Log4perl;
use List::Util qw/min max sum/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pipe_augustus pipe_augustus_simple pipe_genemark pipe_glimmerhmm pipe_geneid
    /;
@EXPORT_OK = qw//;

sub sum_gff {
    my ($fi, $locH) = @_;
    my $h = {};
    open(IN, "<$fi") or die "Failed: $!\n";
    while(<IN>) {
        chomp;
        my @ps = split "\t";
        next if @ps < 9 || $ps[6] ne "+" || $ps[2] ne "CDS";
        my $note = $ps[8];
        $note .= ";" if $note !~ /\;$/;
        if($note =~ /Parent=([\w\.\_\-]+)\;/) {
            my $id = $1;
            $h->{$id} ||= [];
            push @{$h->{$id}}, \@ps;
        } else {
            die "cannot get mRNA ID from\n\t".join("\t", @ps)."\n";
        }
    }
    return undef if keys(%$h) == 0;
    
    my @stats;
    for my $id (keys(%$h)) {
        my $row = $h->{$id};
        $row = [ sort {$a->[3] <=> $b->[3]} @$row ];
        $row->[0]->[3] += $row->[0]->[7];
        $row->[0]->[7] = 0;

        my $len = sum( map {$_->[4] - $_->[3] + 1} @$row );
        $row->[-1]->[4] -= $len % 3;

        my @locs = map {[$_->[3], $_->[4]]} @$row;
        my @phases = map {$_->[7]} @$row;
        my ($locO, $lenO) = posOvlp(\@locs, $locH);
        push @stats, [\@locs, $lenO];
    }
    @stats = sort {$a->[1] <=> $b->[1]} @stats;
    return undef if $stats[-1]->[1] == 0;
    return $stats[-1]->[0];
}
sub sum_gff_batch {
    my ($fi, $dir, $fo) = @_;
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("collecting prediction results");
    my $t = readTable(-in=>$fi, -header=>1);
    
    open(my $fho, ">$fo") || die "cannot write to $fo\n";
    print $fho join("\t", @HEAD_GTBX)."\n";
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seqP, $seq) = $t->row($i);
        my $locL = locStr2Ary($locLS);
    
        my $fm = "$dir/$id";
        next unless -s $fm;
        my $locC = sum_gff($fm, $locL);
        next unless $locC;

        my $seq_cds = getSubSeq($seq, $locC);
        my $seq_pro = Bio::Seq->new(-seq=>$seq_cds)->translate()->seq;
        my $locCStr = locAry2Str($locC);
        my $phase = join(",", @{getPhase($locC, "+")});
        
        print $fho join("\t", "$id.1", $id, ("") x 3, "+", ("") x 2, $locCStr, "", "", $phase, "", ("") x 5, $seq_pro)."\n";
    }
    close $fho;
    runCmd("rm -rf $dir/*", 0);
}

sub run_augustus_batch {
    my ($fi, $dirO) = @_;
    make_path($dirO) unless -d $dirO;
    remove_tree($dirO, {keep_root => 1});
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("running Augustus(evidence mode)");

    my $soft = $ENV{"Augustus"};
    my $f_bin = "$soft/bin/augustus";
    die "$f_bin not there\n" unless -s $f_bin;
    my $d_cfg = "$soft/config";
    my $t = readTable(-in=>$fi, -header=>1);
    my $f_fas = $ENV{"TMP_DIR"}."/augustus_seq_".int(rand(1000)).".fa";
    my $f_hin = $ENV{"TMP_DIR"}."/augustus_hint_".int(rand(1000)).".gff";
    
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 2576;
        writeFile($f_fas, ">tmp", $seq);
        
        my $locL = locStr2Ary($locLS);
        my @hints_hit = map {join("\t", "tmp", "pz", "CDSpart", @$_, ".", "+", ".", "source=M")} @$locL;
        writeFile($f_hin, @hints_hit);
      
        my $fo = "$dirO/$id";
        my $cmd = "$f_bin --AUGUSTUS_CONFIG_PATH=$d_cfg --species=arabidopsis --hintsfile=$f_hin --gff3=on --strand=forward --noInFrameStop=true $f_fas > $fo";
        runCmd($cmd, 0);
        printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    system("rm $f_fas $f_hin");
}
sub sum_augustus1 {
    my ($fi) = @_;
    my $h;
    open(IN, "<$fi") or die "Failed: $!\n";
    my $gene;
    while(<IN>) {
        chomp;
        if(/^\# start gene (\w+)/) {
            $gene = $1;
            $h->{$gene} ||= [];
        } elsif(/^\# end gene (\w+)/) {
            $gene = "";
        } elsif(/^\# hint groups fully obeyed:\s*(\d+)/) {
            $h->{$gene}->[0] = $1;
        }
        next if /^\#/;
        my @ps = split "\t";
        next unless @ps >= 9;
        if($ps[2] eq "CDS" && $ps[8] =~ /t1$/) {
            $h->{$gene}->[1] ||= [];
            push @{$h->{$gene}->[1]}, \@ps;
        }
    }
    my @genes = grep {$h->{$_}->[0] > 0} (keys(%$h));
    if(@genes > 0) {
        @genes = sort {$h->{$a}->[0] <=> $h->{$b}->[0]} @genes;
        my $row = $h->{$genes[-1]}->[1];
        die Dumper($row)." not on + strand\n" unless $row->[0]->[6] eq "+";
        
        $row = [ sort {$a->[3] <=> $b->[3]} @$row ];
        $row->[0]->[3] += $row->[0]->[7];
        $row->[0]->[7] = 0;

        my $len = sum( map {$_->[4] - $_->[3] + 1} @$row );
        $row->[-1]->[4] -= $len % 3;

        my @locs = map {[$_->[3], $_->[4]]} @$row;
        my @phases = map {$_->[7]} @$row;
        return (\@locs, \@phases);
    } else {
        return undef;
    }
}
sub crop_cds { # cut model from specified start codon
    my ($loc, $pos) = @_;
    my $locN = [];
    my $len_crop = 0;
    for (sort {$a->[0] <=> $b->[0]} @$loc) {
        my ($beg, $end) = @$_;
        if($pos > $beg) {
            if($pos <= $end) {
                $len_crop += $pos-1 - $beg + 1;
                push @$locN, [$pos, $end];
            } else {
                $len_crop += $end - $beg + 1;
            }
        } else {
            push @$locN, [$beg, $end];
        }
    }
    return ($locN, $len_crop);
}
sub pipe_augustus {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    my $d02 = "$dir/02_raw";
    run_augustus_batch($fi, $d02);
    my $f11 = "$dir/11.gtb";
    sum_gff_batch($fi, $d02, $f11);
}

sub run_augustus_batch_simple {
    my ($fi, $dirO) = @_;
    make_path($dirO) unless -d $dirO;
    remove_tree($dirO, {keep_root => 1});
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("running Augustus(de novo mode)");

    my $f_bin = $ENV{"Augustus"}."/bin/augustus";
    die "$f_bin not there\n" unless -s $f_bin;
    my $d_cfg = $ENV{"Augustus"}."/config";
    my $t = readTable(-in=>$fi, -header=>1);
    my $f_fas = $ENV{"TMP_DIR"}."/augustus_simple_".int(rand(1000)).".fa";
    
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
        my $fo = "$dirO/$id";
#    next if $id != 30;
        writeFile($f_fas, ">tmp", $seq);
        my $cmd = "$f_bin --AUGUSTUS_CONFIG_PATH=$d_cfg --species=arabidopsis --gff3=on --strand=forward --noInFrameStop=true $f_fas > $fo";
        runCmd($cmd, 0);
        printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    system("rm $f_fas");
}
sub pipe_augustus_simple {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    my $d02 = "$dir/02_raw";
    run_augustus_batch_simple($fi, $d02);
    my $f11 = "$dir/11.gtb";
    sum_gff_batch($fi, $d02, $f11);
}

sub run_genemark_batch {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    remove_tree($dir, {keep_root=>1});
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("running GeneMark");
    
    my $f_bin = $ENV{"GeneMark"}."/gmhmme3";
    $log->error_die("$f_bin not there") unless -s $f_bin;

    my $f_mod;
    if($dir =~ /Athaliana/) {
        $f_mod = "a_thaliana.mod";
    } elsif($dir =~ /Mtruncatula/) {
        $f_mod = "m_truncatula.mod";
    } elsif($dir =~ /Osativa/) {
        $f_mod = "o_sativa.mod";
    } else {
        $log->warn("No proper MOD file to use -> using a_thaliana.mod instead");
        $f_mod = "a_thaliana.mod";
    }
    $f_mod = $ENV{"GeneMark"}."/".$f_mod;

    my $t = readTable(-in=>$fi, -header=>1);
    my $f_fas = $ENV{"TMP_DIR"}."/genemark_".int(rand(1000)).".fa";
    
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#        next if $id != 1;
        writeFile($f_fas, ">tmp", $seq);
        my $fo = "$dir/$id";
        runCmd("$f_bin -m $f_mod -f gff3 -o $fo $f_fas", 0);
        printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    system("rm $f_fas"); 
}
sub pipe_genemark {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    my $d02 = "$dir/02_raw";
    run_genemark_batch($fi, $d02);
    my $f11 = "$dir/11.gtb";
    sum_gff_batch($fi, $d02, $f11);
}

sub run_glimmerhmm_batch {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    remove_tree($dir, {keep_root=>1});
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("running GlimmerHMM");
    
    my $f_bin = $ENV{"GlimmerHMM"}."/bin/glimmerhmm";
    $log->error_die("$f_bin not there") unless -s $f_bin;

    my $d_train;
    if($dir =~ /Athaliana/) {
        $d_train = "arabidopsis";
    } elsif($dir =~ /Osativa/) {
        $d_train = "rice";
    } else {
        $log->warn("no training directory available -> using arabidopsis instead");
        $d_train = "arabidopsis";
    }
    $d_train = $ENV{"GlimmerHMM"}."/trained_dir/$d_train";

    my $t = readTable(-in=>$fi, -header=>1);
    my $f_fas = $ENV{"TMP_DIR"}."/glimmerhmm_".int(rand(1000)).".fa";
    
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 1;
        writeFile($f_fas, ">tmp", $seq);
        my $fo = "$dir/$id";
        runCmd("$f_bin $f_fas $d_train -g > $fo", 0);
        printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    system("rm $f_fas"); 
}
sub pipe_glimmerhmm {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    my $d02 = "$dir/02_raw";
    run_glimmerhmm_batch($fi, $d02);
    my $f11 = "$dir/11.gtb";
    sum_gff_batch($fi, $d02, $f11);
}

sub run_geneid_batch {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    remove_tree($dir, {keep_root=>1});
    my $log = Log::Log4perl->get_logger("Software");
    $log->info("running geneid");

    my $f_bin = $ENV{"GeneID"}."/bin/geneid";
    $log->error_die("$f_bin not there") unless -s $f_bin;

    my $t = readTable(-in=>$fi, -header=>1);
    my $f_fas = $ENV{"TMP_DIR"}."/geneid_".int(rand(1000)).".fa";
  
    my $f_param;
    if($dir =~ /Athaliana/) {
        $f_param = "arabidopsis";
    } elsif($dir =~ /Osativa/) {
        $f_param = "rice";
    } else {
        $log->warn("no param file available -> using arabidopsis instead");
        $f_param = "arabidopsis";
    }
    $f_param = $ENV{"GeneID"}."/param/$f_param.param";
  
    for my $i (0..$t->nofRow-1) {
        my ($id, $fam, $chr, $beg, $end, $srd, $locS, $begr, $endr, $begL, $endL, $locLS, $e, $seq_pro, $seq) = $t->row($i);
#    next if $id != 1;
        writeFile($f_fas, ">tmp", $seq);
        my $fo = "$dir/$id";
        runCmd("$f_bin -3 -W -P $f_param $f_fas > $fo", 0);
        printf "  %5d / %5d done...\n", $i+1, $t->nofRow if ($i+1) % 1000 == 0;
    }
    system("rm $f_fas"); 
}
sub pipe_geneid {
    my ($fi, $dir) = @_;
    make_path($dir) unless -d $dir;
    my $d02 = "$dir/02_raw";
    run_geneid_batch($fi, $d02);
    my $f11 = "$dir/11.gtb";
    sum_gff_batch($fi, $d02, $f11);
}


1;
__END__
