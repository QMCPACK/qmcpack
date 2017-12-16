#!/usr/bin/env perl
use strict;
use FileHandle;
use Getopt::Long;

my %config = do "/home/abenali/Work/src/qmcpack/utils/setup-qmc-conf.pl";
my $gnuplot = $config{gnuplot};
my $energytool = $config{energytool};

my $epsfile;
my $supercellsize = 1;
my $maxstepsize = 10;
my $quadratic;
GetOptions('eps=s' => \$epsfile,
           'supercellsize=i' => \$supercellsize,
           'maxstepsize=f' => \$maxstepsize,
           'quadratic' => \$quadratic);

my $fname = shift;
my $start = shift;
$start > 0 || die "second command line argument must be the starting value at which to take data from the DMC output\n";


open (XMLFILE, "$fname") || die "Cannot open template xml file: $fname\n";
my @temdata = <XMLFILE>;
close (XMLFILE);




my $projid;
my $seriesStart;
my $counter = -1;
my $startdmcsec = 0;
my $stopdmcsec = 0;
my $largesttstep = -1.0;
my %tsteps;
my %fnames;
foreach my $line (@temdata) {
  if ($line =~ /project/ && $line =~ /id/) {
      $line =~ /id\s*?=\s*?"(.*?)"/;
      $projid = $1;
      $line =~ /series\s*?=\s*?"(.*)"/;
      $seriesStart = $1;
      $counter = $seriesStart;
  }

  if ($line =~ /qmc/ && $line =~ /method/) {
      if ($line =~ /dmc/) {
	my $prettyseries;
	if ($counter < 10) {
	  $prettyseries = sprintf("00%d", $counter);
        } elsif ($counter < 100) {
          $prettyseries = sprintf("0%2d", $counter);
        } else {
          $prettyseries = $counter;
        }
        $fnames{$counter} = "$projid.s$prettyseries.scalar.dat";
#        print "series $counter is dmc, scalars in: $projid.s$prettyseries.scalar.dat\n";
        $startdmcsec = 1;
      } else {
#        print "series $counter is not dmc\n";
      }
  }
  if ($startdmcsec) {
      if ($line =~ /timestep/) {
        $line =~ />\s*(.*?)\s*</;
	my $tstep = $1;
	if ($tstep > $largesttstep) {
	    $largesttstep = $tstep;
	}
        $tsteps{$counter} = $tstep;
#        print "  For series $counter the timestep is $tstep\n";
      }
      if ($line =~ /<\/qmc>/) {
	$startdmcsec = 0;
      }
  }
  if ($line =~ /<\/qmc>/) {
     $counter++;
  }

}

if ($largesttstep > $maxstepsize) {
    $largesttstep = $maxstepsize;
}

#print "Project id = $projid\n";
#print "Starting series = $seriesStart\n";


my $gplstring = "set title \"timestep convergence of DMC\" \n";
if ($epsfile) {
    $gplstring .= "set term post color enhanced 20\n set output \"$epsfile\"\n";
}
$gplstring .= "set xlabel \"timestep\"; set ylabel \"energy (Ha)\"\n";
##$gplstring .= "set fit quiet\n";
$gplstring .= "set print \"tsteps\"\n";
if ($quadratic) {
    $gplstring .= "f(x) = a+b*x+c*x**2\n";
    $gplstring .= "fit [:$maxstepsize] f(x) \"-\" u 1:2:3 via a,b,c\n";
} else {
    $gplstring .= "f(x) = a+b*x\n";
    $gplstring .= "fit [:$maxstepsize] f(x) \"-\" u 1:2:3 via a,b\n";
}
foreach my $key ( sort by_key keys %fnames ) {
   my $str = `$energytool $fnames{$key} $start | head -1`;
   my @data = split(/\s+/,$str);
   if ($tsteps{$key} <= $maxstepsize) {
       my $energy = "$tsteps{$key}   " .  $data[2]/$supercellsize . "   " .  $data[4]/$supercellsize . "\n";
       $gplstring .= "$energy";  
   }
}
$gplstring .= "end\n";

foreach my $key ( sort by_key keys %fnames ) {
   my $str = `$energytool $fnames{$key} $start | head -1`;
   my @data = split(/\s+/,$str);
   my $energy = "$tsteps{$key}   " .  $data[2]/$supercellsize . "   " .  $data[4]/$supercellsize;
   $gplstring .= "print \"$energy\"\n ";  
}


my $plotmax = $largesttstep*1.05;
$gplstring .= "plot [0:$plotmax] \"-\" u 1:2:3 lw 3 ps 2 pt 7 notitle w e, f(x) w l lw 3 notitle \n";

#my $xtics = "set xtics (";
foreach my $key ( sort by_key keys %fnames ) {
   my $str = `$energytool $fnames{$key} $start | head -1`;
   my @data = split(/\s+/,$str);
   my $energy = "$tsteps{$key}   " .  $data[2]/$supercellsize . "   " .  $data[4]/$supercellsize . "\n";
   $gplstring .= "$energy";  
}

$gplstring .= "end \n";
unless($epsfile) {
  $gplstring .= "pause -1\n";
}
open(GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $gplstring;
#print $gplstring;
unless($epsfile) {
  my $redundantString = <>; # Hack to leave graph up until user hits a key
}
close(GPL);




sub by_key {
  return $a <=> $b;
}


