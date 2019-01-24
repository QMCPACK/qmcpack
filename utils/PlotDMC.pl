#!/usr/bin/env perl
use POSIX qw/floor/;
use strict;
use FileHandle;
use Getopt::Long;

my %config = do "/gprojects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";
my $gnuplot = $config{gnuplot};
                                                                                
my $epsfile;                                                                    
GetOptions('eps=s' => \$epsfile);       

my $file = shift;
my $start = shift;
(defined $start) || ($start = 0);
my $end = shift @ARGV;
if (!(defined $end)) {
    $end = -1;
}

open(FILE, "$file");
my @filedata;
while (defined (my $line = <FILE>)) {
   my @tmp = split('\s+', $line); 
   push(@filedata, [@tmp]); 
}
close(FILE);


my $min;
my $max;
($min, $max) = getMinMax(\@filedata, 2, $start, $end);
my $enTicString = getTicString($min,$max);

my $popmin;
my $popmax;
($popmin, $popmax) = getMinMax(\@filedata, 5, $start, $end);
my $popTicString = getTicString($popmin,$popmax);

my $plotstring = "set title \"$file\"\n";
if ($epsfile) {
    $plotstring .= "set term post color enhanced 20\n set output \"$epsfile\"\n";
}
$plotstring .= "set multiplot\n set ylabel \"Energy (Hartrees)\"\n ";
$plotstring .= "set origin 0,0.5\n set size 0.82,0.5\n";
$plotstring .= "set format x \"\"\n set ytics $enTicString \n unset y2tics\n set lmargin 12\n";
if ($end < 0) {
  $plotstring .= "plot [$start:][$min:$max] \"$file\" u 1:2 w l notitle, \"\" u 1:5 lt 3 notitle w l \n";
} else {
  $plotstring .= "plot [$start:$end][$min:$max] \"$file\" u 1:2 w l notitle, \"\" u 1:5 lt 3 notitle w l \n";
}
$plotstring .= "set format x \n unset title\n set xlabel \"DMC step\"\n set ylabel \"Population (Walkers)\"\n";
$plotstring .= "set xtics\n set ytics $popTicString\n set bmargin 4\n";
if ($end < 0) {
  $plotstring .= "set origin 0,0\n plot [$start:][$popmin:$popmax] \"\" u 1:5 w l notitle\n";
} else {
  $plotstring .= "set origin 0,0\n plot [$start:$end][$popmin:$popmax] \"\" u 1:5 w l notitle\n";
}
$plotstring .= "set format x \"\"\n";
$plotstring .= "set origin 0.8,0.5\n set size 0.2, 0.5\n set title \"histogram\"\n";
$plotstring .= "set lmargin 0\n unset xlabel\n unset ylabel \n set bmargin 1\n";
$plotstring .= "unset ytics \n set y2tics $enTicString \n set format y2 \"\"\n ";
$plotstring .= getGPLHistString(\@filedata, 2, $start, $end);
$plotstring .= "set origin 0.8,0.0\n set size 0.2, 0.5\n unset title \n set y2tics $popTicString \n";
$plotstring .= "set format y2 \"\"\n set xlabel \"counts\"\n set bmargin 4\n";
$plotstring .= getGPLHistString(\@filedata, 5, $start, $end);
$plotstring .= "unset multiplot\n";
unless($epsfile) {
    $plotstring .= "pause -1\n";
}

open(GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $plotstring;
unless($epsfile) {
    my $blah = <>; # Hack to leave graph up until user hits a key
}
close(GPL);


#############################################################################

sub getTicString {
    use POSIX qw/floor ceil/;
    my $low = shift;
    my $high = shift;
#    print "    $low   $high\n";
    my $range = $high - $low;
    my $scale = 10 ** floor(log($range)/log(10));
    my $ntics = int($range/$scale);
    if ($ntics == 1) {
	$scale *= 0.2;
	$ntics = int($range/$scale); 
    } elsif ($ntics == 2) {
	$scale *= 0.4;
	$ntics = int($range/$scale);
    } elsif ($ntics == 3) {
	$scale *= 0.4;
	$ntics = int($range/$scale);
    } elsif ($ntics == 4) {
	$scale *= 0.5;
	$ntics = int($range/$scale);
    } elsif ($ntics == 5) {
	$scale *= 0.8;
	$ntics = int($range/$scale);
    }    
    my $lowval = ceil($low/$scale)*$scale;
    my $highval = floor($high/$scale)*$scale;
    my $ticstring = "$lowval,$scale,$highval";
    return $ticstring;
}


sub getMinMax {
    my $ArrayRef = shift;
    my $column = shift;
    my $start = shift;
    my $stop = shift;

    if ($stop < 0) {
	$stop = $#{$ArrayRef};
    }
    
    my $min = 1000000000;
    my $max = -1000000000;
    my $numsamples = 0;
    for (my $i = $start+1; $i <= $stop; $i++) {
	my $data = ${$ArrayRef}[$i][$column];
	if ($data < $min) {
	    $min = $data;
	}
	if ($data > $max) {
	    $max = $data;
	}
	$numsamples++;
    }
    return ($min, $max, $numsamples);
}

sub getGPLHistString {
    my $ArrayRef = shift;
    my $column = shift;
    my $start = shift;
    my $stop = shift;

    my $min;
    my $max;
    my $numsamples;
    ($min, $max, $numsamples) = getMinMax($ArrayRef,$column,$start,$stop);    
    my @arr = getColumn($ArrayRef,$column,$start,$stop);
    
    my $numbins = floor($numsamples / 50);
    
    my @histdata;
    for (my $i = 0; $i < $numbins; $i++) {
	push @histdata, 0;
    }
    
    my $binwidth = ($max-$min)/($numbins+0.0);
    my $maxcount = 0;
    foreach my $i (@arr) {
	my $binnum = floor(($i-$min)/$binwidth);
	if ($binnum == $numbins || $binnum < 0) {
	    $binnum--;
	}
	$histdata[$binnum]++;
	if ($histdata[$binnum] > $maxcount) {
	    $maxcount = $histdata[$binnum];
	}
    }
    my $plotmax = $maxcount+2;
    
#    my $plotstring = "unset ytics; set y2tics\n";
    my $plotstring .= "plot [$plotmax:0][$min:$max] \"-\" u 1:2 w l notitle\n "; 
#    my $plotstring .= "plot [10:0][$min:$max] \"-\" u 1:2 w l notitle\n "; 
    for (my $i = 0; $i < $numbins; $i++) {
	my $boxxmin = $min+$i*$binwidth;
	my $boxxmax = $min+($i+1)*$binwidth;
	my $boxymax = $histdata[$i];
	my $tag = $i+1;
	$plotstring .= "0 $boxxmin\n $boxymax $boxxmin\n $boxymax $boxxmax\n 0 $boxxmax\n";
    }
    $plotstring .= "end \n";
    return $plotstring;
}

sub getColumn {
    my $matrix = shift;
    my $colnum = shift;
    my $startnum = shift;
    my $end = shift;

    if ($end < 0) {
	$end = $#{$matrix};
    }

    $startnum || ($startnum = 0);
    my @outmat;
    for (my $i = $startnum+1; $i <= $end; $i++) {
	push @outmat, ${$matrix}[$i][$colnum];
    }
    return @outmat;
}



