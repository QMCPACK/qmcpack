#!/usr/bin/env perl
use POSIX qw/floor/;
use strict;

# To use this to analyze a series of twists, the files should all follow the pattern
# basename-tw#.s???.scalar.dat.  All of these twists will be assumed to have equal
# weight.  Give the code the basename of the lowest numbered twist and then the second
# argument specifies how many sequentially numbered twists after that one should
# be included.


my $template = shift @ARGV || die "First argument must be pattern name for one of the data files!\n";
my $numfiles = shift @ARGV || die "Second argument must be number of data files (sequentially from the pattern one) to be analyzed\n";
my $start = shift @ARGV;
((defined $start) && ($start >= 0)) || die "Third argument must be row for start of data to be analyzed!\n";
my $factor = shift @ARGV;
if (!(defined $factor)) {
    $factor = 1;
}
my $end = shift @ARGV;
if (!(defined $end)) {
    $end = -1;
}

#print "template = $template\n";
#print "numfiles = $numfiles\n";

$template =~ /(.*-tw)(\d+)(\..*)?(\.s...\.scalar\.dat)/;
my $prefix = $1;
my $tnum = $2;
my $suffix = $4;

#print "prefix = $prefix,  tnum = $tnum,  suffix = $suffix\n";

opendir DIR, "."; 
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR; 
closedir DIR; 


my @twistens;
my @twisterrs;
my @rawfiles;
if ($numfiles > 1) {

    for (my $i = $tnum; $i < $tnum+$numfiles; $i++) {
	my @match = grep /$prefix$i(\..*)?$suffix/, @files;
	
	my $curfile;
	if ($match[0]) {
	    $curfile = shift @match;
	    push @rawfiles, $curfile;
	    push @twistens, 0.0;
	    push @twisterrs, 0.0;
	} else {
	    die "Cannot find enough raw data files to satsify this request\n";
	}
#    print "Looking for a file like $prefix" . $i . ".????$suffix:   found $curfile\n";
    }
} else {
    $tnum = 0;
    push @rawfiles, $template;
    push @twistens, 0.0;
    push @twisterrs, 0.0;
}

open (FILE, "$template") || die "Cannot open file $template\n";

my @filedata;
while (defined (my $line = <FILE>)) {
   my @tmp = split('\s+', $line); 
   push(@filedata, [@tmp]); 
}
close(FILE);

isRectangular(\@filedata) || die "$template is not rectangular\n";
my $numcomponents = $#{$filedata[0]}-1;
my @components;
my @componenterr;
my @componentnames;
for (my $i = 0; $i < $numcomponents; $i++) {
    push @components, 0.0;
    push @componenterr, 0.0;
    my $colname = $filedata[0][$i+2];
    if ($colname eq 'LocalEnergy_sq') {
	push @componentnames, "Variance";
    } else {
	push @componentnames, $colname;
    }
}


for (my $fnum = 0; $fnum < $numfiles; $fnum++) {
    
    open (FILE, "$rawfiles[$fnum]");
    # load in data
    my @newfiledata;
    while (defined (my $line = <FILE>)) {
	my @tmp = split('\s+',$line);
	push (@newfiledata, [@tmp]);
    }
    close(FILE);
    # loop over columns and do statistics
    for (my $i = 2; $i <= $#{$newfiledata[0]}; $i++) {    
	my $colname =  $newfiledata[0][$i];
	my @arr = getColumn(\@newfiledata,$i,$start,$end);
	for (my $j = 0; $j <= $#arr; $j++) {
	    if ($colname eq 'AcceptRatio' || $colname eq 'BlockCPU' || $colname eq 'BlockWeight') {
		$arr[$j] *= 1;
	    } else {
		$arr[$j] /= $factor;
	    }
	}
	my $avg;
	my $error;
	($avg, $error) = stats(\@arr);
	if ($colname eq 'LocalEnergy_sq') {
	    my @arr2 = getColumn(\@filedata, $i-1, $start, $end);
	    my $avg2;
	    my $err2;
	    ($avg2, $err2) = stats(\@arr2);
	    $avg = $avg - $avg2*$avg2/$factor;
	}
	if ($colname eq 'LocalEnergy') {
	    $twistens[$fnum] = $avg;
	    $twisterrs[$fnum] = $error;
	}
	$components[$i-2] += $avg * (1.0/$numfiles);
	$componenterr[$i-2] += $error*$error*(1.0/$numfiles/$numfiles);
    }
}

my $MPC = 0;
my $Locen = 0;
my $kecorr = 0;
my $ewald = 0;
my $enerr = 0;

for (my $i = 0; $i < $numcomponents; $i++) {
    my $colname = $componentnames[$i];
    my $avg = $components[$i];
    my $error = sqrt($componenterr[$i]);
    my $formatStr = formatString($error);
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
    my $str = sprintf("%-21s =  $formatStr +/- $formatStr\n", $colname, $avg, $error);
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


print "--------------------------------------------------------------\n";
print "Raw Energies for All Twists\n";
for (my $i = 0; $i < $numfiles; $i++) {
    my $avg = $twistens[$i];
    my $error = $twisterrs[$i];
    my $twistnum = $i + $tnum;
    my $colname = "Twist \#$twistnum";
    my $formatStr = formatString($error);
    my $str = sprintf("%-21s =  $formatStr +/- $formatStr\n", $colname, $avg, $error);
    print $str;
}    



#########################################################################################
sub stats {
    my $arr = shift;
    my $avg = getAvg($arr);
    my $sqav = getSqAvg($arr);
    my $var = $sqav-$avg*$avg;
    my $c=0;
    my $N = $#{$arr}+1;
    my $tempC = 0.5;
    my $kappa = 0.0;
    while($tempC > 0 and $c < ($N-1)) {
	$kappa+=2*$tempC;
	$c += 1;
	$tempC=corr($c, $arr, $avg, $var);
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

sub getAvg {
    my $arr = shift;
    my $size = $#{$arr};
    return getSum($arr)/($size+1);
}

sub getSqAvg {
    my $arr = shift;
    my @newarr;
    for (my $i = 0; $i <= $#{$arr}; $i++) {
	my $val = ${$arr}[$i];
	push @newarr, $val*$val;
    }
    return getAvg(\@newarr);
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
    } elsif ($rightDigits > 12) {
	$rightDigits = 12;
    }

    my $formatstr = '%16.'.$rightDigits."f";
    return $formatstr;
}


