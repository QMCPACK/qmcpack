#!/usr/bin/env perl 
use strict;
use POSIX qw/floor/;
use Getopt::Long;
use FileHandle;


# This tool uses energy.pl to plot the energy of a series of VMC optimizations.
# Simply give it the filename of one of the .scalar.dat files and it will plot
# the energy as a function of sequence.  The locations of gnuplot and energy.pl
# may need to be configured on the next lines

my $exactke;
my $epsfile;
GetOptions('ke=f' => \$exactke,
           'eps=s' => \$epsfile);


my %config = do "/gprojects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";

my $gnuplot = $config{gnuplot};
my $energytool = $config{energytool};

my $template = shift || die "Must give a template file name as the first argument here\n";
my $start = 50;
$start = shift;
(defined $start) || ($start = 0);
my $plotstart = shift @ARGV;
if (!(defined $plotstart)) {
    $plotstart = 1;
}


$template =~ /(.*f)(\d\.\d\d)(.*\.s\d\d\d\.scalar\.dat)/;
my $prefix = $1;
my $fact = $2;
my $suffix = $3;
#print "prefix = $prefix, factor = $fact, suffix = $suffix\n";


opendir DIR, "."; 
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR; 
closedir DIR; 

my @rawfiles;
foreach my $str (@files) {
    if ($str =~  /(.*f)(\d\.\d\d)(.*\.s\d\d\d\.scalar\.dat)/) {
	push @rawfiles, $str;
    }
}


my @factors;
my @fulldata;
my $varmin = 10000000000;
my $varmax = -1;
my $kemin = 1000000000;
my $kemax = -100000000;
my $minen = 1000000000;
my $maxen = -1000000000;
foreach my $file (sort byFactor @rawfiles) {
    $file =~ /(.*f)(\d\.\d\d)(.*\.s\d\d\d\.scalar\.dat)/;
    my $factor = $2;
    push @factors, $factor;

    my $str = `$energytool $file $start`;
#    print "String = $str\n";
    my @data = split(/\n/, $str);
#    print "data[0] = $data[0]\n";
    my $en;
    my $ke;
    my $var;
    foreach my $line (@data) {
	if ($line =~ /LocalEnergy/) {
#	    print "LocEnLine = $line\n";
	    my @linedata = split(/\s+/,$line);
	    $en = "$linedata[2]  $linedata[4]";
#	    print "en = $en\n";
	    if ($linedata[2] < $minen) {
		$minen = $linedata[2];
	    }
	    if ($linedata[2] > $maxen) {
		$maxen = $linedata[2];
	    }
	} elsif ($line =~ /Variance/) {
	    my @linedata = split(/\s+/,$line);
#	    print "VarLine = $line\n";
	    $var = "$linedata[2]  $linedata[4]";
	    if ($linedata[2] < $varmin) {
		$varmin = $linedata[2];
	    }
	    if ($linedata[2] > $varmax) {
		$varmax = $linedata[2];
	    }
	} elsif ($line =~ /Kinetic/) {
	    my @linedata = split(/\s+/,$line);
	    $ke = "$linedata[2]  $linedata[4]";
	    if ($linedata[2] < $kemin) {
		$kemin = $linedata[2];
	    }
	    if ($linedata[2] > $kemax) {
		$kemax = $linedata[2];
	    }
#	    print "KELine = $line\n";
	}
    }
    $str = "$factor   $en  $var  $ke";
#    print "data = $str\n";
    push(@fulldata, $str);

#    my $str = `$energytool $file $start | head -1`;
#    my @data = split(/\s+/,$str);
#    my $energy = "$data[2]  $data[4]";   
#    push @energies, $energy;

#    $str = `$energytool $file $start | grep Kinetic`;
#    @data = split(/\s+/, $str);
#    my $kenergy = "$data[2]   $data[4]";
#    push @kenergies, $kenergy;
}

my $plotstring;
if ($epsfile) {
    $plotstring .= "set term post color enhanced 16\n set output \"$epsfile\" \n";
}
$plotstring .= "set title \"B-spline factor vs energies for VMC\" \n";
$plotstring .= "set multiplot \n set origin 0,0.6\n set size 1,0.4\n";
$plotstring .= "set ylabel \"total energy (Ha)\"\n unset xtics\n";
$plotstring .= "set lmargin 15\n";
my $st = $factors[0]-0.05;
my $ed = $factors[$#factors]+0.05;
my $xtics = "set xtics (";
for (my $i = 0; $i <= $#factors; $i++) {
    $xtics .= "\"$factors[$i]\" $factors[$i]";
    if ($i < $#factors) {
	$xtics .= ", " ;
    }
}
$plotstring .= "$xtics )\n"; 
my $yincr = ($maxen-$minen)/4.0;
$plotstring .= "set ytics $yincr\n";
$plotstring .= "plot [$st:$ed] \"-\" u 1:2:3 notitle w e \n";
for (my $i = 0; $i <= $#factors; $i++) {
#    $plotstring .= "$factors[$i]  $energies[$i]\n";
    $plotstring .= "$fulldata[$i]\n"; 
   
}
$plotstring .= "end \n";

# Do plot of variance
$plotstring .= "set origin 0,0.33\n set size 1,0.3\n";
$plotstring .= "unset title\n";
$plotstring .= "set ylabel \"LocEn Variance\"\n set logscale y\n";
$yincr = ($varmax-$varmin)/4.0;
$plotstring .= "set ytics (";
for (my $yticval = $varmin; $yticval <= $varmax; $yticval += $yincr) {
    $plotstring .= "\"$yticval\" $yticval";
    if ($yticval + $yincr <= $varmax) {
	$plotstring .= ", ";
    }
}
$plotstring .= ")\n";


$plotstring .= "plot [$st:$ed][0.95*$varmin:1.05*$varmax] \"-\" u 1:4:5 notitle w e \n";

for (my $i = 0; $i <= $#factors; $i++) {
#    $plotstring .= "$factors[$i]   $kenergies[$i]\n";
    $plotstring .= "$fulldata[$i]\n";
}
$plotstring .= "end \n";
$plotstring .= "unset logscale y\n";

# Do plot of kinetic energy
$plotstring .= "set origin 0,0\n set size 1,0.35\n";
$plotstring .= "unset title\n";
$plotstring .= "set xlabel \"b-spline factor\"\n";
$plotstring .= "set ylabel \"kinetic energy (Ha)\"\n";
$yincr = ($kemax-$kemin)/4.0;
$plotstring .= "set ytics $yincr\n";
if ($exactke) {
    $plotstring .= "f(x) = $exactke\n";
    $plotstring .= "plot [$st:$ed] f(x) lw 3 lt 1 notitle, \"-\" u 1:6:7 notitle w e \n";
} else {
    $plotstring .= "plot [$st:$ed] \"-\" u 1:6:7 notitle w e \n";
}
for (my $i = 0; $i <= $#factors; $i++) {
#    $plotstring .= "$factors[$i]   $kenergies[$i]\n";
    $plotstring .= "$fulldata[$i]\n";
}
$plotstring .= "end \n";
$plotstring .= "unset multiplot\n";
unless($epsfile) {
    $plotstring .= "pause -1\n";
}

open(GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $plotstring;
unless($epsfile) {
    my $redundantString = <>; # Hack to leave graph up until user hits a key
}
close(GPL);


#############################################################################

sub byFactor {
    my $left = $a;
    my $right = $b;
    $left =~ /(.*f)(\d\.\d\d)(.*\.s\d\d\d\.scalar\.dat)/;
    my $leftseq = $2;
    $right =~ /(.*f)(\d\.\d\d)(.*\.s\d\d\d\.scalar\.dat)/;
    my $rightseq = $2;
    $leftseq <=> $rightseq;
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




