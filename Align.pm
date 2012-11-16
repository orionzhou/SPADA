package Align;
use strict;
use Bio::AlignIO;
use Bio::Matrix::IO;
use Common;
use Seq;
use Data::Dumper; 
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/read_aln_ids aln_fmt_convert
    run_clustalw run_tcoffee run_clustalo run_pal2nal
    pwAln run_dotplot
    aln_score_pw aln_score_vector aln_score_group aln_score_matrix/;
@EXPORT_OK = qw//;

sub read_aln_ids {
    my ($fi, $format) = @_;
    $format ||= "clustalw";
    my $ai = Bio::AlignIO->new(-file=>"<$fi", -format=>$format);
    my @ids;
    while(my $aln = $ai->next_aln()) {
        for my $seq ($aln->each_seq()) {
            my $id = $seq->id;
            push @ids, $id;
        }
    }
    return @ids;
}
sub aln_fmt_convert {
    my ($fi, $fo, $if, $of) = @_;
    my $ai = Bio::AlignIO->new(-file=>"<$fi", -format=>$if);
    my $ao = Bio::AlignIO->new(-file=>">$fo", -format=>$of);
    while(my $aln = $ai->next_aln()) {
        $ao->write_aln($aln);
    }
}
sub run_clustalw {
    my ($seqAry, $fOut, $type) = rearrange(['seqs', 'out', 'type'], @_);
    $type = uc($type);
    my $f_bin = $ENV{"ClustalW"}."/bin/clustalw2";
    my $f1 = $ENV{"TMP_DIR"}."/clustalw_".int(rand(1000)).".fa";
    writeSeq($seqAry, $f1);
    my $cmd = qq/$f_bin -INFILE=$f1 -ALIGN -OUTFILE=$fOut -OUTORDER=INPUT -TYPE=$type/;
    my $matrixP ||= "BLOSUM";
    my $matrixD ||= "CLUSTALW";
    $cmd .= " -PWMATRIX=$matrixP -MATRIX=$matrixP" if $type eq "PROTEIN";
    $cmd .= " -PWDNAMATRIX=$matrixD -DNAMATRIX=$matrixD" if $type eq "DNA";
    runCmd($cmd, 0);
    my $fTree = $f1;
    $fTree =~ s/\.\w+$/\.dnd/i;
    system("rm $fTree $f1"); 
}
sub run_clustalo {
    my ($seqs, $fi, $fo, $f_dis, $f_hmm) = rearrange(['seqs', 'in', 'out', 'dis', 'f_hmm'], @_);
    
    my $mode = defined($fi) ? 1 : 2;
    if($mode == 2) {
        if( defined $ENV{"TMP_DIR"} ) {
            $fi = $ENV{"TMP_DIR"}."/clustalo_".int(rand(1000)).".fa";
        } else {
            $fi = "$fo.tmp";
        }
        writeSeq($seqs, $fi);
    }
    my $f_bin = $ENV{"ClustalO"} ? $ENV{"ClustalO"}."/bin/clustalo" : "clustalo";

    my $cmd = qq/$f_bin -i $fi -o $fo --outfmt=clu --force/;
    $cmd .= " --hmm-in=$f_hmm" if $f_hmm;
    $cmd .= " --full --distmat-out=$f_dis" if $f_dis;
    runCmd($cmd, -1);
    die "error running clustalo\n" unless -s $fo;
    system("rm $fi") if $mode == 2; 
}
sub run_pal2nal {
    my ($fIn1, $fIn2, $fOut) = rearrange(['in1', 'in2', 'out'], @_);
    die "$fIn1 is not there\n" unless -s $fIn1;
    die "$fIn2 is not there\n" unless -s $fIn2;
    my $cmd = qq/pal2nal.pl $fIn2 $fIn1 > $fOut/;
    runCmd($cmd, 0);
}
sub run_tcoffee {
    my ($seqAry, $fo) = rearrange(['seqs', 'out'], @_);
    my $f1 = $ENV{"TMP_DIR"}."/tcoffee_".int(rand(1000)).".fa";
    writeSeq(-seqs=>$seqAry, -out=>$f1);
    my $cmd = qq/t_coffee $f1 -outfile $fo -gapopen -100 -gapext -5/;
    runCmd($cmd, 0);
    system("rm $f1 $fo.html"); 
    system("rm $fo.dnd") if -f "$fo.dnd";
    system("rm tcoffee.dnd") if -f "tcoffee.dnd";
}

sub pwAln {
    my ($seqAry, $fo, $format) = rearrange(['seqs', 'out', 'format'], @_);
    my $f_bin = "water";
    $format ||= "pair";
    my $fi1 = $ENV{"TMP_DIR"}."/water_1_".int(rand(1000)).".fa";
    my $fi2 = $ENV{"TMP_DIR"}."/water_2_".int(rand(1000)).".fa";
    writeSeq(-seqs=>$seqAry->[0], -out=>$fi1);
    writeSeq(-seqs=>$seqAry->[1], -out=>$fi2);
    my ($sid1, $sid2) = map {$_->id} @$seqAry;
    my $cmd = qq/$f_bin $fi1 $fi2 -gapopen 20 -gapextend 2 -aformat $format -outfile $fo/;
    runCmd($cmd, 0);
    system("rm $fi1 $fi2");
}
sub run_dotplot {
    my ($seqAry, $fo, $ws) = rearrange(['seqs', 'out', "ws"], @_);
    $ws ||= 10;
    my $fi1 = $ENV{"TMP_DIR"}."/dottup_1_".int(rand(1000)).".fa";
    my $fi2 = $ENV{"TMP_DIR"}."/dottup_2_".int(rand(1000)).".fa";
    writeSeq(-seqs=>$seqAry->[0], -out=>$fi1);
    writeSeq(-seqs=>$seqAry->[1], -out=>$fi2);
    my ($sid1, $sid2) = map {$_->id} @$seqAry;
    my $cmd = qq/dottup $fi1 $fi2 -graph png -word $ws/;
    runCmd($cmd, 0);
    system("rm $fi1 $fi2");
    if( -s "dottup.1.png" ) {
        system("mv dottup.1.png $fo");
    } else {
        die "dotplot failed\n";
    }
}

sub aln_score_pw {
    my ($seq1, $seq2) = @_;
    die "seq length not equal\n$seq1\n$seq2\n" unless length($seq1) == length($seq2);
    ($seq1, $seq2) = (uc($seq1), uc($seq2));
    $seq1 =~ s/\./\-/g;
    $seq2 =~ s/\./\-/g;
    
    my ($gap_open, $gap_ext) = (-8, -1);

    my $parser = Bio::Matrix::IO->new(-format=>'scoring', -file=>$ENV{"SPADA_SOURCE"}.'/BLOSUM62');
    my $matrix = $parser->next_matrix();

    my $score = 0;
    my $tag_gap_open = 1;
    for my $i (0..length($seq1)-1) {
        my $aa1 = substr($seq1, $i, 1);
        my $aa2 = substr($seq2, $i, 1);
        if($aa1 eq "-" && $aa2 eq "-") {
        } elsif($aa1 ne "-" && $aa2 ne "-") {
            $score += $matrix->get_entry($aa1, $aa2);
            $tag_gap_open = 1;
        } else {
            if($tag_gap_open == 1) {
                $score += $gap_open;
                $tag_gap_open = 0;
            } else {
                $score += $gap_ext;
            }
        }
    }
    return $score;
}
sub aln_score_vector {
    my ($fi, $id) = @_;
    my $h_seq = readSeq($fi, 2);
    die "$id not in $fi\n" unless exists $h_seq->{$id};
    my $seq1 = $h_seq->{$id};
    my @scores;
    for my $id2 (keys(%$h_seq)) {
        next if $id2 eq $id;
        my $seq2 = $h_seq->{$id2};
        push @scores, aln_score_pw($seq1, $seq2);
    }
    my $score = sprintf "%d", sum(@scores) / @scores;
    return $score;
}
sub aln_score_group {
    my ($fi, $idsQ, $idsT) = @_;
    my $h_seq = readSeq($fi, 2);

    my @idsQ_n = grep { !exists($h_seq->{$_}) } @$idsQ;
    die "idsQ not all found\n".Dumper($idsQ).Dumper($h_seq) if @idsQ_n;
    my $h = { map {$_ => 0} @$idsQ };
    $idsT ||= [ grep { !exists($h->{$_}) } keys(%$h_seq) ];

    for my $idQ (@$idsQ) {
        my @scores;
        for my $idT (@$idsT) {
            my $score = aln_score_pw($h_seq->{$idQ}, $h_seq->{$idT});
            push @scores, $score;
        }
        $h->{$idQ} = sprintf "%.01f", sum(@scores) / @scores;
    }
    return $h;
}
sub aln_score_matrix {
    my ($fi, $ids, $hi) = @_;

    my $h_seq = readSeq($fi, 2);
    my @seqs = map {$h_seq->{$_}} @$ids;
    my $nseq = @seqs;
    
    my ($hs, $m) = ({}, []);
    for my $i (0..$nseq-1) {
        $m->[$i] = [];
        my @js;
        if(defined($hi)) {
            @js = @{$hi->[$i]};
        } else {
            push @js, (0..$i-1) if $i >= 1;
            push @js, ($i+1, $nseq-1) if $i <= $nseq-2;
        }
        for my $j (@js) {
            if(exists($hs->{"$i.$j"})) {
                push @{$m->[$i]}, $hs->{"$i.$j"};
            } else {
                my $score = aln_score_pw($seqs[$i], $seqs[$j]);
                push @{$m->[$i]}, $score;
                $hs->{"$j.$i"} = $score;
            }
        }
    }
    return $m;
}


sub align {
# opt( 1 - dna, 2 - protein, 3 - dna2protein, 4 - cds )
    my ($seqAry, $fOut, $opt, $prog) = rearrange(['seqs', 'out', "opt", "program"], @_);
    my @uniqIds = uniq( map { $_->id } @$seqAry );
    my @idxs;
    for my $i (0..$#uniqIds) {
        my $idx = first_index {$seqAry->[$_]->id eq $uniqIds[$i]} (0..@$seqAry-1);
        push @idxs, $idx;
    }
    $seqAry = [ @$seqAry[@idxs] ];
    my $idH;
    for my $i (0..@$seqAry-1) {
        my $seq = $seqAry->[$i];
        my $id = $seq->id;
        if($id =~ /[\:]/ || length($id) > 13) {
            my $idN = $id;
            $idN =~ s/[\:]/_/g;
            $idN =~ s/Medtr//g;
            $idN =~ s/contig\_/ctg/ig;
            $idN =~ s/\.trim//;
            $idN = substr($idN, 0, 13) if length($idN) > 13;
            $idN =~ s/\W+$//;
            $seqAry->[$i] = Bio::Seq->new(-id=>$idN, -seq=>$seq->seq);
            $idH->{$idN} = $id;
        }
    }
    my $ref = undef;
    if(@$seqAry == 1) {
        print "only 1 sequence -> alignment not made\n";
        return undef;
    }
    if($opt == 1 || $opt == 2) {
        if($prog =~ /^(water)|(needle)$/) {
            die "not 2 seqs\n" unless @$seqAry == 2;
            pwAln(-seqs=>$seqAry, -out=>$fOut, -format=>'pair', -program=>$prog);
            $ref = getAlnDesc($fOut, 'emboss', $idH);
        } else {
            run_tcoffee(-seqs=>$seqAry, -out=>$fOut);
#      run_clustalw(-seqs=>$seqAry, -out=>$fOut, -type=> $opt==1 ? "DNA" : "PROTEIN");
            $ref = getAlnDesc($fOut, 'clustalw', $idH);
        }
    } else {
        my $seqAry2 = [];
        for (@$seqAry) {
            if($_->alphabet eq "dna") {
                my $seqP = $_->translate(-terminator=>'X')->seq;
                $seqP =~ s/X$//i;
                push @$seqAry2, Bio::Seq->new(-id=>$_->id, -seq=>$seqP);
            } elsif($_->alphabet eq "protein") {
                push @$seqAry2, $_;
            }
        }
        if($opt == 3) {
            if($prog && $prog =~ /^(water)|(needle)$/) {
                die "not 2 seqs\n" unless @$seqAry == 2;
                pwAln(-seqs=>$seqAry2, -out=>$fOut, -format=>'pair', -program=>$prog);
                $ref = getAlnDesc($fOut, 'emboss', $idH);
            } else {
                run_tcoffee(-seqs=>$seqAry2, -out=>$fOut);
#        run_clustalw(-seqs=>$seqAry2, -out=>$fOut,, -type=>"PROTEIN");
                $ref = getAlnDesc($fOut, 'clustalw', $idH);
            }
        } elsif($opt == 4) {
            my $f1 = "/tmp/cds.fa";
            writeSeq(-seqs=>$seqAry, -out=>$f1);
            my $f2 = "/tmp/protein.aln";
            if($prog =~ /^(water)|(needle)$/) {
                die "not 2 seqs\n" unless @$seqAry == 2;
                pwAln(-seqs=>$seqAry2, -out=>$f2, -format=>'clustal', -program=>$prog);
            } else {
                run_tcoffee(-seqs=>$seqAry2, -out=>$f2);
#        run_clustalw(-seqs=>$seqAry2, -out=>$f2,, -type=>"PROTEIN");
            }
            $ref = getAlnDesc($f2, 'clustalw', $idH);
            run_pal2nal(-in1=>$f1, -in2=>$f2, -out=>$fOut);
            system("rm $f1 $f2");
        }
    }
    return $ref;
}


1;
__END__
