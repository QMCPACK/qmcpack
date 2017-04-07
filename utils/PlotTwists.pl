#!/usr/bin/env perl
use POSIX qw/floor/;
use strict;
use FileHandle;

# This is a utility to track the progress of several running twists with the following
# pattern: basename-tw#.g???.s???.scalar.dat.  Simply give the name of one of the files
# to the script and it will plot all matching files it can find.  Before using the
# script, the location of gnuplot may need to be set on the previous line

my %config = do "/gprojects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";
my $gnuplot = $config{gnuplot};

my $template = shift || die "Must give a template file name as the first argument here\n";
my $start = shift;
(defined $start) || ($start = 0);
my $end = shift @ARGV;
if (!(defined $end)) {
    $end = -1;
}


$template =~ /(.*-tw)(\d+)\..*(\.s...\.scalar\.dat)/;
my $prefix = $1;
my $tnum = $2;
my $suffix = $3;


opendir DIR, "."; 
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR; 
closedir DIR; 

my @rawfiles;
foreach my $str (@files) {
    if ($str =~  /$prefix(\d+)\..*$suffix/) {
	push @rawfiles, $str;
    }
}

my $plotstring = "plot [$start:] ";
my $twists = $#rawfiles;
my $blah = $twists+1;
print "Number of twists = $blah\n";
my $counter = 0;
foreach my $file (sort byTwist @rawfiles) {
    $file =~ /(.*-tw)(\d+)\..*(\.s...\.scalar\.dat)/;
    my $twist = $2;
    $plotstring .= " \"$file\" u 1:2 w l ti \"Twist $twist\"";
#    print "Plotting file: $file\n";
    $counter++;
    if ($counter <= $twists) {
#	print "Got here\n";
	$plotstring .=", ";
    }
}
$plotstring .= "\n";
$plotstring .= "pause -1\n";



open(GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $plotstring;
my $redundantString = <>; # Hack to leave graph up until user hits a key
close(GPL);


#############################################################################

sub byTwist {
    my $left = $a;
    my $right = $b;
    $left =~ /(.*-tw)(\d+)\..*(\.s...\.scalar\.dat)/;
    my $lefttwist = $2;
    $right =~ /(.*-tw)(\d+)\..*(\.s...\.scalar\.dat)/;
    my $righttwist = $2;
    $lefttwist <=> $righttwist;
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



