#!/usr/bin/env perl 
use POSIX qw/floor/;
use strict;
use FileHandle;
use Getopt::Long;

my %config = do "/gprojects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";

my $gnuplot = $config{gnuplot};
my $energytool = $config{energytool};

my $epsfile;
GetOptions('eps=s' => \$epsfile);

my $template = shift || die "Must give a template file name as the first argument here\n";
my $start = shift;
(defined $start) || ($start = 0);
my $plotstart = shift @ARGV;
if (!(defined $plotstart)) {
    $plotstart = 1;
}



$template =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
my $prefix = $1;
my $tnum = $2;
my $suffix = $3;
#print "prefix = $prefix, sequence number = $tnum, suffix = $suffix\n";


opendir DIR, "."; 
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR; 
closedir DIR; 

my @rawfiles;
foreach my $str (@files) {
    if ($str =~  /$prefix\d\d\d$suffix/) {
	push @rawfiles, $str;
    }
}


my @sequences;
my @energies;
my @correnergies;
foreach my $file (sort bySequence @rawfiles) {
    $file =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
    my $sequence = $2;
    push @sequences, $sequence;
    my $str = `$energytool $file $start | head -1`;
    my @data = split(/\s+/,$str);
    my $energy = "$data[2]  $data[4]";   
    print "$sequence  $energy\n";
    push @energies, $energy;
#    $str = `$energytool $file $start | tail -1`;
#    @data = split(/\s+/,$str);
#    my $correnergy = "$data[3]  $data[5]";   
#    push @correnergies, $correnergy;
}

my $plotstring = "set title \"Sequence vs Energy for VMC optimization\" \n";

if ($epsfile) {
   $plotstring .= "set term post color enhanced 20\n set output \"$epsfile\"\n";
}
$plotstring .= "set xlabel \"sequence\"; set ylabel \"energy (Ha)\"\n";
my $st = $plotstart-1.2;
my $ed = $#sequences+0.2;
my $xtics = "set xtics (";
$plotstring .= "plot [$st:$ed] \"-\" u 1:2:3 notitle w e \n";
for (my $i = 0; $i <= $#sequences; $i++) {
#    $plotstring .= "$i  $correnergies[$i]\n";
    $plotstring .= "$i  $energies[$i]\n"; 
   
    $xtics .= "\"$sequences[$i]\" $i";
    if ($i < $#sequences) {
	$xtics .= ", " ;
    }
}
$xtics .= ")\n";
$plotstring .= "end \n";
unless ($epsfile) {
  $plotstring .= "pause -1\n";
}
$plotstring = "$xtics $plotstring";

open(GPL, "|$gnuplot");
GPL->autoflush(1);
print GPL $plotstring;
unless($epsfile) {
  my $redundantString = <>; # Hack to leave graph up until user hits a key
}
close(GPL);


#############################################################################

sub bySequence {
    my $left = $a;
    my $right = $b;
    $left =~ /\.s(\d\d\d)\.scalar\.dat/;
    my $leftseq = $1;
    $right =~ /\.s(\d\d\d)\.scalar\.dat/;
    my $rightseq = $1;
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


