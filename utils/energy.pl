#!/usr/bin/env perl
use POSIX qw/floor/;
use Getopt::Long;
use strict;

## This is a perl script to replace Energy.py.
## This is created solely because numpy can be hard to come by on new machines
## and is a real pain to build by hand


my $inEv = 0;
GetOptions('ev' => \$inEv);

my $filename = shift @ARGV || die "First argument must be name of file to analyze!\n";
my $start = shift @ARGV;
((defined $start) && ($start >= 0)) || die "Second argument must be row for start of data to be analyzed!\n";
my $factor = shift @ARGV;
if (!(defined $factor)) {
    $factor = 1;
}
my $end = shift @ARGV;
if (!(defined $end)) {
    $end = -1;
}

if ($inEv) {
    $factor /= 27.2113966413;
}


open (FILE, "$filename") || die "Cannot open file $filename\n";

my @filedata;
while (defined (my $line = <FILE>)) {
   my @tmp = split('\s+', $line); 
   push(@filedata, [@tmp]); 
}
close(FILE);


isRectangular(\@filedata) || die "$filename is not rectangular\n";


my $MPC = 0;
my $Locen = 0;
my $kecorr = 0;
my $ewald = 0;
my $enerr = 0;

my $blockcpu=0;
my $blockweight=0;

for (my $i = 2; $i <= $#{$filedata[0]}; $i++) {    
    my $colname =  $filedata[0][$i];
    my @arr = getColumn(\@filedata,$i,$start,$end);

    my $avg;
    my $error;
    if ($colname eq 'AcceptRatio' || $colname eq 'BlockCPU' || $colname eq 'BlockWeight' || $colname eq 'Efficiency') {
	($avg, $error) = simplestats(\@arr);
    } else {
	($avg, $error) = stats(\@arr, $factor);
    }
    if ($colname eq 'LocalEnergy_sq') {
	$colname = 'Variance';
	my @arr2 = getColumn(\@filedata, $i-1, $start, $end);
	my $avg2;
	my $err2;
	($avg2, $err2) = stats(\@arr2, $factor);
	$avg /= $factor;
	$error /= $factor;
	$avg = $avg - $avg2*$avg2;
    }
    if ($colname eq 'LocalEnergy') {
	$Locen = $avg;
	$enerr = $error;
    } elsif ($colname eq 'ElecElec') {
	$ewald = $avg;
    } elsif ($colname eq 'MPC') {
	$MPC = $avg;
    } elsif ($colname eq 'KEcorr') {
	$kecorr = $avg;
    }
    my $formatStr = formatString($error);
    my $str = sprintf("%-21s =  $formatStr +/- $formatStr\n", $colname, $avg, $error);
    if($colname eq 'BlockCPU')
    {
      $blockcpu=$avg;
    }
    if($colname eq 'BlockWeight')
    {
      $blockweight=$avg;
    }
    print $str;
}

if($blockcpu>0)
{
    my $formatStr = formatString(0.01);
    my $str = sprintf("%-21s =  $formatStr +/- $formatStr\n", "Efficiency", $blockweight/$blockcpu, 0);
    print $str;
}

my $correction = $kecorr;
if ($ewald != 0 && $MPC != 0) {
    $correction += $MPC - $ewald;
}
if (abs($correction) > 1.0e-12) {
    print "--------------------------------------------------------------\n";
    my $formatStr = formatString($enerr);
    printf("%-21s =  $formatStr +/- $formatStr\n", "Corrected Energy", $Locen+$correction, $enerr);
}

sub simplestats {
    my $arr = shift;
    my $avg;
    my $sqavg;
    ($avg, $sqavg) = getAvgAndSqAv($arr);
    my $var = $sqavg-$avg*$avg;
    return ($avg, sqrt($var/$#{$arr}));
}


sub stats {
    my $arr = shift;
    my $locfac = shift;
    my $avg;
    my $sqav;
    ($avg, $sqav) = getAvgAndSqAv($arr);
    $avg /= $locfac;
    $sqav /= $locfac*$locfac;
    my $var = $sqav-$avg*$avg;
    my $c=0;
    my $N = $#{$arr}+1;
    my $tempC = 0.5;
    my $kappa = 0.0;
    while($tempC > 0 and $c < ($N-1)) {
	$kappa+=2*$tempC;
	$c += 1;
	$tempC=corr($c, $arr, $avg*$locfac, $var*$locfac*$locfac);
    }
    if ($kappa == 0.0) {
	$kappa = 1.0;
    }
    my $Neff = $N/$kappa;
    if ($Neff == 0.0) {
	$Neff = 1.0;
    }
    if ($var < 0.0) {
	$var = 0.0;
    }
    my $error = sqrt($var/$Neff);
    return ($avg, $error);
}

sub isRectangular {
    my $isRect = 1;
    my $matrix = shift;
    my $numcols = $#{${$matrix}[0]};
    for (my $i = 1; $i <= $#{$matrix}; $i++) {
	if ($numcols != $#{${$matrix}[$i]}) {
	    $isRect = 0;
	}
    }
    $isRect;
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

sub getSum {
    my $arr = shift;
    my $outval;
    for (my $i = 0; $i <= $#{$arr}; $i++) {
	$outval += ${$arr}[$i];
    }
    return $outval;
}


sub getAvgAndSqAv {
    my $arr = shift;
    my $sum;
    my $sumsq;
    for (my $i = 0; $i <= $#{$arr}; $i++) {
	$sum += ${$arr}[$i];
	$sumsq += ${$arr}[$i]*${$arr}[$i];
    }
    return ($sum/($#{$arr}+1), $sumsq/($#{$arr}+1));
}



sub corr {
    my $i = shift;
    my $arr = shift;
    my $mean = shift;
    my $var = shift;
    
    my $N = $#{$arr}+1;
    if ($var < 1.0e-10) {
	return 1;
    }
    my @newarr;
    for (my $ind = 0; $ind < $N-$i; $ind++) {
	my $val1 = ${$arr}[$ind]-$mean;
	my $val2 = ${$arr}[$i+$ind]-$mean;
	push @newarr, $val1*$val2;
    }
    return getSum(\@newarr)/$var/($N-$i);
}

sub formatString {
    my $error = shift;
    my $rightDigits=8;
    if ($error != 0.0) {
	$rightDigits = -floor(log($error)/log(10))+1;
    }
    if ($rightDigits < 0) {
	$rightDigits = 0;
    }
    my $formatstr = '%16.'.$rightDigits."f";
    return $formatstr;
}


