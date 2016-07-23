package Common;
use strict; 
use Bio::Seq;
use Data::Table;
use File::Basename;
use Data::Dumper;
use IPC::Open3; use Symbol qw/gensym/;
use Clone qw/clone/;
use Graph;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/pretty tonum isnumber getDigits prettyStr 
  runCmd parse_gff_tags
  rearrange group bsearch readTable
  mergeArray aryCmp
  is_revsrd get_revsrd
  getIdxRange getPhase sample_serial
  rmRedPairs backOneLine scaleNumber writeFile/;
@EXPORT_OK = qw//;
sub rearrange {
    my( $order, @param ) = @_;
    return unless @param;
    my %param;
    if (ref $param[0] eq 'HASH') {
        %param = %{$param[0]};
    } else {
        return @param unless (defined($param[0]) && substr($param[0],0,1) eq '-');
        my $i;
        for ($i=0;$i<@param;$i+=2) {
            $param[$i]=~s/^\-//;  # get rid of initial - if present
            $param[$i]=~tr/a-z/A-Z/; # parameters are upper case
        }
        %param = @param; # convert into associative array
    }
    my (@return_array);
    local($^W) = 0;
    my ($key)='';
    foreach $key (@$order) {
        my ($value);
        if (ref($key) eq 'ARRAY') {
            foreach (@$key) {
            last if defined($value);
            $value = $param{uc($_)};
            delete $param{$_};
        }
        } else {
            $value = $param{uc($key)};
            delete $param{$key};
        }
        push(@return_array,$value);
    }
    push (@return_array,\%param) if %param;
    return @return_array;
}
sub pretty {
    my ($num, $sep) = @_;
    $sep = defined $sep ? $sep: ",";
    die "not a number: $num\n" if $num !~ /^\d+$/;
    $num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1$sep/g;
    return $num;
}
sub prettyStr {
  my ($str, $len, $sep) = @_;
  $len ||= 5;
  $sep ||= " ";
  my @ary;
  for(my $i = 0; $i*$len < length($str); $i++) {
    push @ary, substr($str, $i*$len, $len);
  }
  return join($sep, @ary);
}
sub tonum {
  my ($num) = @_;
  $num =~ s/[^\d]//g;
  return $num;
}
sub getDigits {
  my ($num) = @_;
  my $digit = 1;
  while(int($num/10) >= 1) {
    $num = int($num/10);
    $digit ++;
  }
  return $digit;
}
sub isnumber {
  my ($a) = @_;
  if( $a =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
    return 1;
  } else {
    return 0;
  }
}
sub runCmd {
  my ($cmd, $opt) = @_; # opt=1 [print to STDOUT]; opt=2 [return output]
  $opt = defined($opt) ? $opt : 1;

  print $cmd."\n" if $opt == 1;
  open( PS, $cmd." |" ) or die "Failed: $!\n";
  my @output;
  while( <PS> ) {
    if($opt == 1) {
      print STDOUT $_;
    } elsif($opt == 2) {
      chomp;
      push @output, $_;
    }
  }
  return \@output if $opt == 2;
}
sub parse_gff_tags {
  my ($str) = @_;
  my @tagAry = split(";", $str);
  my $h;
  for (@tagAry) {
    my ($tag, $value) = split "=";
    die $str if $value eq "";
    $value =~ s/\=/\:/g;
    $value =~ s/\;/\|/g;
    $h->{$tag} = $value;
  }
  return $h;
}
sub mergeArray {
  my @refAry = @_;
  my @rst;
  for my $ref (@refAry) {
    if(ref($ref) eq "ARRAY") {
      push @rst, mergeArray(@$ref);
    } elsif(!ref($ref)) {
      push @rst, $ref;
    } else{
      die("unknown type $ref\n");
    }
  }
  return @rst;
}
sub aryCmp {
  my ($ary1, $ary2) = @_;
  my @ary1 = uniq(@$ary1);
  my @ary2 = uniq(@$ary2);
  my @ary_merged = uniq(@ary1, @ary2);
  my %hash1 = map {$_=>1} @ary1;
  my %hash2 = map {$_=>1} @ary2;
  my (@ary_share, @ary_1, @ary_2);
  for my $ele (@ary_merged) {
    if(exists $hash1{$ele} && exists $hash2{$ele}) {
      push @ary_share, $ele;
    } elsif(exists $hash1{$ele} && !exists $hash2{$ele}) {
      push @ary_1, $ele;
    } elsif(!exists $hash1{$ele} && exists $hash2{$ele}) {
      push @ary_2, $ele;
    } else {
      die "$ele not in ary1 nor ary2\n";
    }
  }
  return (\@ary_share, \@ary_1, \@ary_2);
}
sub readTable {
  my ($fi, $fh, $header, $skip) = rearrange(['in', 'inh', 'header', 'skip'], @_);
  die "no file or handler passed" if !$fi && !$fh;
  $skip ||= 0;
  $header ||= 0;
  if(!$fh) {
    -e $fi || die "$fi is not there\n";
    open($fh, "<$fi") || die "cannot read $fi\n";
  }

  my $t;
  if($header ne 1 && $header ne 0) {
    $t = Data::Table::fromTSV($fh, 0, $header, {skip_lines=>$skip, skip_pattern=>'(^\s*#)|(^\s*$)'});
  } else {
    $t = Data::Table::fromTSV($fh, $header, {skip_lines=>$skip, skip_pattern=>'(^\s*#)|(^\s*$)'});
  }
#  printf "\t%s: cols[%d] rows[%d]\n", basename($fi), $t->nofCol, $t->nofRow;
  return $t;
}
sub scaleNumber {
  my ($ary, $down, $up) = rearrange(["value", "down", "up"], @_);
  my $rst;
  ($down, $up) = $down<=$up ? ($down, $up) : ($up, $down);
  my $minX = min(@$ary);
  my $maxX = max(@$ary);
  $down ||= 0;
  $up ||= 100;
  my $factor = ($up-$down) / ($maxX-$minX);
  for my $x (@$ary) {
    my $y = $down + ($x-$minX) * $factor;
    push @$rst, $y;
  }
  return $rst;
}
sub backOneLine {
  my ($fH) = @_;
  my $char;
  my $flag_nonblank = 0;
  while(seek($fH, -1, 1)) {
    read($fH, $char, 1);
    last if ($char eq "\n" && $flag_nonblank == 1);
    $flag_nonblank = 1 if $char ne "\n";
    seek($fH, -1, 1);
  }
  return tell($fH);
}
sub get_revsrd {
  my ($srdI) = @_;
  return "-" if $srdI =~ /^[\+1]$/;
  return "+" if $srdI =~ /^\-1?$/;
  die "unknonw strand: $srdI\n";
}
sub is_revsrd {
  my ($srd1, $srd2) = @_;
  $srd1 = 1 if $srd1 eq "+";
  $srd1 = -1 if $srd1 eq "-";
  $srd2 = 1 if $srd2 eq "+";
  $srd2 = -1 if $srd2 eq "-";
  return 1 if $srd1 * $srd2 == -1;
  return 0 if $srd1 * $srd2 == 1;
  die "unknown strands: $srd1  $srd2\n";
}

sub rmRedPairs {
    my ($edgeRef) = @_;
    my $g = Graph::Undirected->new();
    $g->add_edges(map {@$_} @$edgeRef);
    my @pairs = ();
    my @es;
    my @vs = $g->vertices();
    my $cnt = 0;
    while( @pairs < @vs/2 ) {
        die Dumper($edgeRef) if ++$cnt > 10;
        @es = $g->edges();
#    print join("\t", map {join("-", @$_)} @es)."\n";
        @pairs = ();
        for my $e (@es) {
            my ($v1, $v2) = @$e;
            push @pairs, $e if $g->degree($v1) == 1 || $g->degree($v2) == 1;
        }
        if(@pairs < @vs/2) {
            my $v2 = first_value {$g->degree($_) == 2} @vs;
            $g->delete_edge(@{[$g->edges_at($v2)]->[0]});
        }
    }
    print @pairs." pairs obtained for @vs vertices\n" unless @vs/2 == @pairs;
    return \@pairs;
}
sub group {
  my ($ary) = @_;
  my $ref = {};
  for my $i (0..@$ary-1) {
    my $ele = $ary->[$i];
    if(exists $ref->{$ele}) {
      $ref->{$ele}->[1] ++;
    } else {
      $ref->{$ele} = [$i, 1];
    }
  }
  return $ref;
}
sub bsearch {
    my ($ary, $word) = @_;
    my $low = 0;
    my $high = @$ary - 1;
    my $try;
    while( $low <= $high ) {
        $try = int( ($low + $high) / 2 );
        die $word."\n".join("\t", @$ary)."\n" if $try < 0 || $try > @$ary - 1;
        $low  = $try + 1, next if $ary->[$try] < $word;
        $high = $try - 1, next if $ary->[$try] > $word;
        return $try;
    }
    return $try;
}
sub getIdxRange {
    my @a = @_;
    my @idx;
    my $str = join("", @a);
    my $i = 0;
    while($str =~ /0(1*)/g) {
        my $len = length($1);
        push @idx, [$i, $i+$len];
        $i += $len + 1;
    }
    return @idx;
}
sub getPhase {
    my ($loc, $srd) = @_;
    $srd ||= "+";
    $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
    $loc = [ reverse @$loc ] if $srd =~ /^\-1?$/;

    my @phases;
    my $len = 0;
    for (@$loc) {
        my ($beg, $end) = @$_;
        push @phases, (3-$len%3) % 3;
        $len += $end - $beg + 1;
    }
    return \@phases;
}
sub sample_serial {
    my ($n, $m) = @_;
    my $inc = $n / $m;
    return map {1+$inc*$_-1} (0..$m-1);
}
sub parse_old_loc_str {
  my ($locS) = @_;
  my $srd;
  my $loc = [];
  if($locS =~ /^complement\(([\w\.\,]+)\)$/) {
    $srd = "-";
    $locS = $1;
  } else {
    $srd = "+";
  }
  while($locS =~ /(\d+)\.\.(\d+)/g) {
    push @$loc, [$1, $2];
  }
  $loc = [ sort {$a->[0] <=> $b->[0]} @$loc ];
  return ($loc, $srd);
}
sub writeFile {
    my ($fo, @strs) = @_;
    open(FH, ">$fo") or die "cannot open $fo for writing\n";
    print FH join("\n", @strs)."\n";
    close FH;
}



1;
__END__
