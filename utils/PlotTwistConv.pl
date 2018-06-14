#!/usr/bin/env perl
use strict;
use FileHandle;
use Getopt::Long;


my %config = do "/gprojects/qmcpack/redsky/bin/setup-qmc-conf.pl";

my $supercellsize;
my $start;
my $epsfile;
my $smallestTwist;
my $largestTwist;
my $withoutDFT = 0;
GetOptions('supercellsize=i' => \$supercellsize,
	   'start=i' => \$start,
	   'eps=s' => \$epsfile,
	   'smallesttwist=i' => \$smallestTwist,
	   'largesttwist=i' => \$largestTwist,
	   'withoutdft' => \$withoutDFT);

$#ARGV == 0 || die "Must give a pwscf input file as the argument to this script\n";
my $inputFile = $ARGV[0];

my $baseName = $inputFile;
$baseName =~ s/\.in//g;
$baseName =~ s/-scf//g;


unless ($supercellsize) {
    $supercellsize = 1;
}
unless ($start) {
    $start = 50;
}
unless ($smallestTwist) {
    $smallestTwist = 0;
}
unless ($largestTwist) {
    $largestTwist = 10000000000000000;
}

opendir DIR, ".";
my @files = grep { $_ ne '.' && $_ ne '..'} readdir DIR;
closedir DIR;

my %numToDir;
getDirs(\@files, $supercellsize, \%numToDir);

$inputFile =~ s/\.in/\.out/;
my $scfenout = `grep \"\\!\" $inputFile\n`;
my @arr = split('\s+', $scfenout);
my $scfen = $arr[4];


my %numToEn;
my %numToErr;
my %numToCorrection;
my $allWithCorrections = 1;

# Loop over different calcs and get twist averaged energy
foreach my $keyval ( sort {$a <=> $b} keys %numToDir ) {
    if ($keyval > $smallestTwist && $keyval < $largestTwist) {
	my $pscfname = $baseName;
	if ($keyval > 1) {
	    $pscfname .= "-${keyval}twists";
	}
	if ($numToDir{$keyval} =~ /-S(\d+)/) {
	    $pscfname .= "-S$1";
	}
	$pscfname .= "-pscf.out";
	my $pscfen;
	if (-e $pscfname) {
	    my $pscfenout = `grep \"total energy\" $pscfname`;
	    my @arr2 = split('\s+', $pscfenout);
	    $pscfen = $arr2[4];
	} else {
	    $allWithCorrections = 0;
	}
	my $dftcorr = ($scfen-$pscfen)/2.0;
	$numToCorrection{$keyval} = $dftcorr;
	
	
	my $templatefile = getTwistTemplateFile($keyval, $numToDir{$keyval}, $baseName);
        my $enln = `cd $numToDir{$keyval};   $config{twenergytool} $templatefile $keyval $start | head -1`;
	`cd ..`;
	chomp($enln);
	my @locdata = split('\s+', $enln);
	$numToEn{$keyval} = $locdata[2];
	$numToErr{$keyval} = $locdata[4];
    }
}

# Now write the input for gnuplot
my $gplstring;
if ($epsfile) {
    $gplstring .= "set term post color enhanced 20\n set output \"$epsfile\"\n";
}
$gplstring .= "set title \"Convergence of DMC energy with Twist Averaging\"\n";
$gplstring .= "set xlabel \"Number of Twists\"\n";
$gplstring .= "set ylabel \"Energy (Ha)\"\n";

$gplstring .= "set xtics (";
my $counter = 0;
foreach my $keyval (sort {$a <=> $b} keys %numToEn) {
    $gplstring .= "\"$keyval\" $counter,";
    $counter++;
}
$gplstring = substr($gplstring,0,-1);
$gplstring .= ")\n";

if ($withoutDFT) {
    $allWithCorrections = 0;
}

if ($allWithCorrections) {
    $gplstring .= "plot [-0.5:$counter-0.5] \"-\" u 0:2:3 w e ti \"DMC Energy with DFT correction\"\n";
} else {
    $gplstring .= "plot [-0.5:$counter-0.5] \"-\" u 0:2:3 w e ti \"DMC Energy\"\n";
}

foreach my $keyval (sort {$a <=> $b} keys %numToEn) {
    my $en = $numToEn{$keyval};
    if ($allWithCorrections) {
	$en += $numToCorrection{$keyval};
    }
    $gplstring .= "$keyval   $en  $numToErr{$keyval}\n";
} 
$gplstring .= "end \n";
unless ($epsfile) {
    $gplstring .= "pause -1\n";
}

my $gnuplot = $config{gnuplot};
open (GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $gplstring;

unless($epsfile) {
    <STDIN>;
} 
close(GPL);

#print $gplstring;







##########################################################################################
# subroutines
##########################################################################################

sub getTwistTemplateFile {
    my $numTwists = shift;
    my $dir = shift;
    my $baseName = shift;

    if ($dir =~ /-S(\d+)/) {
	$baseName .= "-S$1";
    }


    my $outfname;
    my $xmlbase;
    
    my $isProd = 0;
    my $isDMC = 0;
    

    if ($dir =~ /production/) {
	$isProd = 1;
    } else {
	$isDMC = 1;
    }

    ## Now we know which input file created the output
    if ($numTwists == 1) {
	if ($isDMC) {
	    $xmlbase = "$baseName-dmc";
	} else {
	    $xmlbase = "$baseName-production";
	}
    } else {
	if ($isDMC) {
	    $xmlbase = "$baseName-dmc-tw0";
	} else {
	    $xmlbase = "$baseName-production-tw0";
	}
    }

    ## Now find which output sequence from this was the dmc sequence
    open (XMLFILE, "$dir/$xmlbase.xml");
    my @temdata = <XMLFILE>;
    close (XMLFILE);

    my $projid;
    my $seriesStart;
    my $counter;
    my $prettyseries;
  
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
		if ($counter < 10) {
		    $prettyseries = sprintf("00%d", $counter);
		} elsif ($counter < 100) {
		    $prettyseries = sprintf("0%2d", $counter);
		} else {
		    $prettyseries = $counter;
		}
	    } 	    
	}
	if ($line =~ /<\/qmc>/) {
	    $counter++;
	}
    }

    ## now we look through the files in the directory for one whose name
    ## contains $projid and s$prettyseries.scalar.dat
    
    opendir DIR, "$dir";
    my @locfiles = grep { $_ ne '.' && $_ ne '..'} readdir DIR;
    closedir DIR;


    foreach my $file (@locfiles) {
	if ($file =~ /$projid.*s$prettyseries\.scalar\.dat/) {
	    $outfname = $file;
	}
    }

    return $outfname;
}




sub getDirs {
    my $fileref = shift;
    my $ssize = shift;
    my $ndref = shift;

    foreach my $str (@{$fileref}) {
	if ($ssize > 1) { 
	    if ($str =~ /^dmc-S$ssize-(\d+)twists$/) {
		my $numtwist = $1; 
		unless ($$ndref{$numtwist}) { 
		    #print "directory: $str, number of twists: $numtwist\n";
		    $$ndref{$numtwist} = $str;
		} 
		if ($$ndref{$numtwist} =~ /production/) {
		    $$ndref{$numtwist} = $str;
		}
	    } elsif ($str =~ /^dmc-S$ssize$/) {
		my $numtwist = 1;
		unless ($$ndref{$numtwist}) {
		    #print "directory: $str, number of twists: $numtwist\n";
		    $$ndref{$numtwist} = $str;
		}
		if ($$ndref{$numtwist} =~ /production/) {
		    $$ndref{$numtwist} = $str;
		}
	    } elsif ($str =~ /^production-S$ssize-(\d+)twists$/) {
		my $numtwist = $1;
		unless ($$ndref{$numtwist}) {
		    $$ndref{$numtwist} = $str;
		}
	    } elsif ($str =~/^production-S$ssize*/) {
		my $numtwist = 1;
		unless ($$ndref{$numtwist}) {
		    $$ndref{$numtwist} = $str;
		}
	    }   
	} else {
	    if ($str =~ /^dmc-(\d+)twists$/) {
		my $numtwist = $1;
		unless ($$ndref{$numtwist}) {
		    #print "directory: $str, number of twists: $numtwist\n";
		    $$ndref{$numtwist} = $str;
		}
		if ($$ndref{$numtwist} =~ /production/) {
		    $$ndref{$numtwist} = $str;
		}
	    }
	    if ($str =~ /^dmc$/) {
		my $numtwist = 1;
		unless ($$ndref{$numtwist}) {
		    #print "directory: $str, number of twists: $numtwist\n";
		    $$ndref{$numtwist} = $str;
		}
		if ($$ndref{$numtwist} =~ /production/) {
		    $$ndref{$numtwist} = $str;
		}
	    } elsif ($str =~ /^production-(\d+)twists$/) {
		my $numtwist = $1;
		unless ($$ndref{$numtwist}) {
		    $$ndref{$numtwist} = $str;
		}
	    } elsif ($str =~/^production$/) {
		my $numtwist = 1;
		unless ($$ndref{$numtwist}) {
		    $$ndref{$numtwist} = $str;
		}
	    }   
	    
	}
    }
}



