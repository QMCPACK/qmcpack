#!/usr/bin/env perl
use strict;
use Getopt::Long;
use POSIX qw/floor/;

my %config = do "/gprojects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";
my $getSupercell = $config{supercell};

print "Supercell command is at: $getSupercell\n";



my @ptvs;
my $mindet = -1;
my $maxdet = -1;
my $maxentry = 7;

GetOptions('ptvs=f{9}' => \@ptvs,
	   'min=i' => \$mindet,
	   'max=i' => \$maxdet, 
	   'maxentry=i' => \$maxentry);

if ($mindet < 0 || $maxdet < 0 || $maxdet < $mindet) {
    die "Must specify maximium and minum size for tilemat using --max and --min where max > min\n";
}



my $radius;
my $numbins = 32;
my $binsize = 0.1;
my $binmin = 0.45;

my @bins;
my @radii;
my @volratios;

for (my $i = 0; $i < $numbins; $i++) {
    $bins[$i] = 0.0;
    $radii[$i] = 0.0;
    $volratios[$i] = 1000.0;
}
my $primcelldet = getDet(\@ptvs);

for (my $i = $mindet; $i <= $maxdet; $i++) {
#    print "i = $i\n";
    my $out = `$getSupercell --ptvs @ptvs --target $i --maxentry $maxentry`;
#    print "output of supercell command is: $out\n";
    my @parsed = split(/\s+/, $out);
    $radius = $parsed[0];
    print "radius for i = $i is $radius\n";
    my $volratio = abs($primcelldet)*$i/4.18879020478639098458/$radius**3;
    my $binnum = floor(($radius-$binmin)/$binsize);
    if ($volratio < $volratios[$binnum]) {
	$bins[$binnum] = $i;
	$radii[$binnum] = $radius;
	$volratios[$binnum] = $volratio;
    }
}

for (my $i = 0; $i < $numbins; $i++) {
    if ($bins[$i] > 0) {
	print "Good value for tilemat with det = $bins[$i], giving radius $radii[$i], and volratio $volratios[$i]\n";
    }
}
			




######################################################################


	
sub getDet {
    my $matref = shift;
    my $val = $$matref[0]*($$matref[4]*$$matref[8] - $$matref[7]*$$matref[5])
	- $$matref[1]*($$matref[3]*$$matref[8] - $$matref[5]*$$matref[6])
	+ $$matref[2]*($$matref[3]*$$matref[7]-$$matref[4]*$$matref[6]);
    return $val;
}


