#!/usr/bin/env perl
use strict;
use Getopt::Long;
use POSIX qw/floor fmod ceil pow/;


# Input pwscf file must follow some rather strict conventions.  It must specify the pseudo
# directory directly, use ibrav=1 or ibrav=0, Decleare explicitly ntyp, nbnd and nspin.
# finally, the atomic positions should be given in terms of the PTV of the cell 
# (not real space coordinates although it will probably be worthwhile to relax this later)

# will allow only the use of ibrav=4 (hcp), ibrav=3 (bcc), ibrav=2 (fcc), ibrav=1 (simple cubic) and ibrav=0 (arbitrary)
#print "$0 @ARGV\n";

my %config = do "/gpfs/mira-home/projects/qmcpack/qmcpack/utils/setup-qmc-conf.pl";

#print "The location of ppconvert is: $config{ppconvert}\n";


###############################################################################
# Parse through the various options which can be given to the script
###############################################################################
my $testwvfcn;
my $getTilemat;
my $analyzeSym;
my $genNSCF;
my $genFSDFT;
my $convBspline;
my $optwvfcn;
my $convDMCTstep;
my $dmcCalc;
my $help;

## Supercell description keywords
my (@inSuperCellTwist, @inSuperCellKGrid, @inSuperCellKShift);
my @toptilingmatrix;
my $targetSsize;
my $isAe = 0;
my @bconds;
my $pbc;
my $obc;

## Wavefunction keywords
my $wvfcnfile;
my $splfactor;
my $sglprecsplines;
my $truncsplines;
my ($minfactor, $maxfactor, $factorinc);
my $twistnum = 0;
my $withjas;
my $toprcutTwoBody = 0;
my $toprcutOneBody = 0;
my $halfboxrcut;
my $jastrowPerIon;

## General MC keywords
my $walkers; # automatically sets both vmcwalkers and dmcwalkers
my $targetpop;

## VMC keywords
my ($vmcblocks, $vmcwalkers, $vmcwarmupsteps, $vmctimestep);
my ($vmcsteps, $vmcSubsteps, $vmcequiltime, $vmcdecorrtime);
my $vmcdrift;
my $numsamples;
$vmcwarmupsteps = 0;
$vmcsteps = 1;
$vmcSubsteps = 1;

## Optimization keywords
my ($optSamples, $optLoops);
my ($oneBodySplinePts, $twoBodySplinePts);

## DMC keywords
my ($dmcblocks, $dmcwalkers, $dmcwarmupsteps, $dmctimestep);
my ($dmcsteps, $dmctstep, $dmcUseTmoves, $mindmctstep, $maxdmctstep, $dmctstepinc);
my ($dmcequiltime, $dmcruntime, $dmcblocktime);
$dmcUseTmoves=1; 

## General execution keywords
my $useGPU = 0;

## Optimize for massive parallelism
my $massivelyParallel = 0;


GetOptions('testwvfcn' => \$testwvfcn,
   'gettilemat' => \$getTilemat,
   'analyzesym' => \$analyzeSym,
   'genwfn' => \$genNSCF,
   'genfsdft' => \$genFSDFT,
           'splconv' => \$convBspline,
   'optwvfcn' => \$optwvfcn,
   'convdmctstep' => \$convDMCTstep,
   'dmc' => \$dmcCalc,
   'wvfcnfile=s' => \$wvfcnfile,
   'tilemat=i{9}' => \@toptilingmatrix,
   'bconds=i{3}' => \@bconds,
   'pbc' => \$pbc,
   'obc' => \$obc,
   'walkers=i' => \$walkers,
   'targetpop=i' => \$targetpop,
   'vmcblocks=i' => \$vmcblocks,
   'withjas:1' => \$withjas,
   'vmcwalkers=i' => \$vmcwalkers,
   'vmcwarmupsteps=i' => \$vmcwarmupsteps,
   'vmctimestep=f' => \$vmctimestep,
   'vmcequiltime=f' => \$vmcequiltime,
   'vmcdecorrtime=f' => \$vmcdecorrtime,
   'numsamples=i' => \$numsamples,
   'vmcsteps=i' => \$vmcsteps,
   'vmcSubsteps=f' => \$vmcSubsteps,
   'vmcdrift' => \$vmcdrift,
   'dmcblocks=i' => \$dmcblocks,
   'dmcwalkers=i' => \$dmcwalkers,
   'dmcwarmupsteps=i' => \$dmcwarmupsteps,
   'dmctimestep=f' => \$dmctimestep,
   'dmcsteps=i' => \$dmcsteps,
   'dmctstep=f' => \$dmctstep,
   'dmcequiltime=f' => \$dmcequiltime,
   'dmcruntime=f' => \$dmcruntime,
   'dmcblocktime=f' => \$dmcblocktime,
   'dmcUseTmoves' => \$dmcUseTmoves,
   'mindmctstep=f' => \$mindmctstep,
   'maxdmctstep=f' => \$maxdmctstep,
   'dmctstepinc=f' => \$dmctstepinc,   
   'minfactor=f' => \$minfactor,
   'maxfactor=f' => \$maxfactor,
   'factorinc=f' => \$factorinc,
   'splfactor=f' => \$splfactor,
           'splprec' => \$sglprecsplines,
   'spltrunc' => \$truncsplines,
   'onebodysplinepts=i' => \$oneBodySplinePts,
   'twobodysplinepts=i' => \$twoBodySplinePts,
   'onebodyrcut=f' => \$toprcutOneBody,
   'twobodyrcut=f' => \$toprcutTwoBody,
   'halfboxrcut' => \$halfboxrcut,
   'jasperion' => \$jastrowPerIon,
   'optsamples=i' => \$optSamples,
   'optloops=i' => \$optLoops,
   'isAe=i' => \$isAe,
   'useGPU' => \$useGPU,
   'massivelyparallel' => \$massivelyParallel,
   'twistnum=i' => \$twistnum,
   'kpoint=f{3}' => \@inSuperCellTwist,
   'kgrid=i{3}' => \@inSuperCellKGrid,
   'kshift=f{3}' => \@inSuperCellKShift,
   'supercellsize=i' => \$targetSsize,
           'help' => \$help);


print "supercell size is: $targetSsize\n";

if (!$targetSsize) {
  $targetSsize = 1;
}

if (! @toptilingmatrix ) {
  @toptilingmatrix = (1, 0, 0, 0, 1, 0, 0, 0, 1);
}

if (! @inSuperCellKGrid) {
    @inSuperCellKGrid = (1, 1, 1);
}
if (! @inSuperCellKShift) {
    @inSuperCellKShift = (0, 0, 0);
}
if (! @bconds ) {
# in future should check to see which type of dft code we are using before doing this
# for instance, if it is a quantum chemistry type of code, the default boundary conditions
# should be open
    @bconds = (1, 1, 1);
}
if ($pbc) {
  @bconds = (1, 1, 1);
} elsif ($obc) {
  @bconds = (0, 0, 0);
}
my $hasOpenBC = 0;
for (my $i = 0; $i < 3; $i++) {
  if ($bconds[$i] == 0) {
    $hasOpenBC = 1;
  }
}

unless($testwvfcn || $genNSCF || $genFSDFT || $convBspline || $optwvfcn || $dmcCalc || $convDMCTstep || $getTilemat || $analyzeSym) {
    globalUsage();
}

for (my $i = 0; $i < 3; $i++) {
    if ($inSuperCellKShift[$i] != 0 && $inSuperCellKShift[$i] != 1) {
die "SuperCell K Shifts can only be either 0 or 1 for each direction\n This corresponds to the grid being centerd on 0 or on the center of the cell\n";
    }
}

($#ARGV == 0 || $help) || die "Must give a pwscf input file as the argument to this script\n";
my $inputFile = $ARGV[0];

$useGPU = 1;

if ($walkers) {
    $vmcwalkers = $walkers;
    $dmcwalkers = $walkers;
}

################################################################################
#
#
# Branch of code to test that given a supercell and supercell twists that the provided
# hdf wavefunction file contains the appropriate kpoints for the requested simulations
#
#
#  LNS TODO: Needs to be adapted to work with symmetry
#
################################################################################
if ($testwvfcn) {
###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);
#####################################################################################
    
    unless((@toptilingmatrix) && $wvfcnfile && ((@inSuperCellTwist) || (@inSuperCellKGrid))) {
die "Must give an input tilematrix, wavefunction file and either a specific\nkpoint or a kgrid for testwvfcn to make sense\n";
    }
    
    my $extractKvecsCommand = $config{kptlister};
    print "$extractKvecsCommand\n";
    my $kvecoutput = `$extractKvecsCommand $wvfcnfile`;
    my @kvecs = split('\s+', $kvecoutput);

    my $sstwists;
    if ((@inSuperCellTwist)) {
$sstwists = 1;
    } else {
$sstwists = $inSuperCellKGrid[0]*$inSuperCellKGrid[1]*$inSuperCellKGrid[2];
    }
    print "Looking for $sstwists supercell twists\n";


    analyzeTwists(\@toptilingmatrix, \@kvecs);
}




################################################################################
#
# Branch of code to get a good tile matrix for a given supercell size
#
################################################################################
if ($getTilemat) {
    unless($targetSsize) {
die "Must give a supercellsize for the gettilemat option to make sense!\n";
    }
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);
    print "@toptilingmatrix\n";
}


################################################################################
#
# Branch of code to generate dummy pwscf files for the supercell so that
# we can let pwscf do the heavy lifting of figuring out which twists are
# equivalent by symmetry and what their multiplicity is
#
# The resulting pwscf files should be run and other parts of this code
# such as genNSCF or genFSDFT will use the results to generate minimal 
# lists of k-points for the wavefunction generation step 
#
# LNS TODO: Also need to update dmcCalc section so that it only generates
#           input files for the required supercell twists
#
################################################################################
if ($analyzeSym) {
    unless($targetSsize || @toptilingmatrix ) {
die "Must give a supercellsize or specify tiling matrix for the" . 
            "supercell for the analyzesym option to make sense!\n";
    }

    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);

    my $numSuperCellTwists = 1;
    for (my $i = 0; $i < 3; $i++) {
$numSuperCellTwists *= $inSuperCellKGrid[$i];
    }    
    my $kptstring = "K_POINTS automatic\n";
    $kptstring .= "$inSuperCellKGrid[0]  $inSuperCellKGrid[1]  $inSuperCellKGrid[2]  $inSuperCellKShift[0]   $inSuperCellKShift[1]   $inSuperCellKShift[2]\n";

    open(IF, $inputFile) || die "cannot open input file $inputFile given on the command line\n";
    my @fdata = <IF>;
    close(IF);

    genPwscfSymm($inputFile, \@fdata, $kptstring);

    my $symmetryFileName = $inputFile;
    $symmetryFileName =~ s/\.in//g;
    $symmetryFileName =~ s/-scf//g;
    $symmetryFileName .= "-SYM";
    $symmetryFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $symmetryFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$symmetryFileName .= "-S$targetSsize";
    }
    $symmetryFileName .= ".in";

    open(SYMMETRY, ">$symmetryFileName");
    print SYMMETRY @fdata;
    close(SYMMETRY);
}


################################################################################
#
# Branch of code to generate NSCF input files
#
################################################################################
if ($genNSCF) {
###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);
#####################################################################################

    if ($help) {
NSCFUsage();
    }

    print "Generating files for nscf creation of wavefunctions\n";

    if (!$#inSuperCellTwist && (@inSuperCellKGrid || @inSuperCellKShift)) {
die "Must specify either kpoint or both kgrid and kshift, but not another combination\n";
    }
    if (@inSuperCellKGrid && !@inSuperCellKShift) {
die "Must specify both kgrid and kshift\n";
    }

    my @superCellTwists;
    getSuperCellTwists($inputFile, $targetSsize, \@inSuperCellKGrid, \@inSuperCellKShift, \@superCellTwists);
#    print "Before findSupercell call, toptilingmatrix = ";
#    for (my $i = 0; $i < 9; $i++) {
#      print "$toptilingmatrix[$i]  ";
#    }
#    print "\n";
#    print "targetSsize = $targetSsize\n";
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);
#    print "After findSupercell call, toptilingmatrix = ";
#    for (my $i = 0; $i < 9; $i++) {
#      print "$toptilingmatrix[$i]  ";
#    }
#    print "\n";    

    #####################################################################################
    # Parse the pwscf input file for most of the information I will need
    #####################################################################################
    open(IF, $inputFile) || die "cannot open input file $inputFile given on the command line\n";

    my @fdata = <IF>;
    close(IF);

    parsePwscfInput(\@fdata, \$calcPrefix, \$pseudoDir, \$outdir, \$celldim, \$numSpecies, \$numAts,
      \$numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials, \@ionIds, \%atNameToPP, 
                    \@posArray);
    #####################################################################################

    my @primCellKpts;
    print "About to get Primitive Cell K-Points\n";
    my $numsskpts = ($#superCellTwists+1)/4;
    print "Number of supercell k-points: $numsskpts\n"; 
    print "Outside getPrimCellKpts call, toptilingmatrix = ";
    for (my $i = 0; $i < 9; $i++) {
      print "$toptilingmatrix[$i]  ";
    }
    print "\n";
    getPrimCellKpts(\@cell_ptv, \@superCellTwists, \@toptilingmatrix, \@primCellKpts, \@inSuperCellKGrid, \@inSuperCellKShift);

#    my @superCellKpts;
#    getKVectors(\@superCellKpts, \@cell_ptv, \@toptilingmatrix, \@inSuperCellKGrid, \@inSuperCellKShift);
    

    my $kPointsSection = "K_POINTS {crystal}\n";
    my $numkpts = ($#primCellKpts+1)/4;
    $kPointsSection .= "$numkpts\n";
    for (my $i = 0; $i <= $#primCellKpts; $i+=4) {
$kPointsSection .= sprintf ("  %16.14f  %16.14f  %16.14f  %16.14f\n", $primCellKpts[$i],$primCellKpts[$i+1],$primCellKpts[$i+2],$primCellKpts[$i+3]);
    }

    my @nscfData = @fdata;
    genNscfPrep(\@nscfData, $kPointsSection);
    my $nscfFileName = $inputFile;
    $nscfFileName =~ s/\.in//g;
    $nscfFileName =~ s/-scf//g;
    $nscfFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $nscfFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$nscfFileName .= "-S$targetSsize";
    }
    $nscfFileName .= "-nscf.in";

    open(NSCF, ">$nscfFileName");
    print NSCF @nscfData;
    close(NSCF);
   
    my @pw2xData;
    genPw2qmcpack(\@pw2xData, $outdir, $calcPrefix);
    my $pw2xFileName = $inputFile;
    $pw2xFileName =~ s/\.in//g;
    $pw2xFileName =~ s/-scf//g;
    $pw2xFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $pw2xFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$pw2xFileName .= "-S$targetSsize";
    }
    $pw2xFileName .= "-pw2x.in";
    open(PW2X, ">$pw2xFileName");
    print PW2X @pw2xData;
    close(PW2X);
}


################################################################################
#
# Branch of code to generate PSCF and KZK input files, also gives pw2casino input
# so that Finite kinetic energy of the specific kpoints can be found.  Note that
# KZK is probably not using the right cell-size in choosing the DFT functional
#
################################################################################
if ($genFSDFT) {
###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);
#####################################################################################

    if ($help) {
FSDFTUsage();
    }

    print "Generating files for pscf calculation of energy and KZK correction\n";

    if (!$#inSuperCellTwist && (@inSuperCellKGrid || @inSuperCellKShift)) {
die "Must specify either kpoint or both kgrid and kshift, but not another combination\n";
    }
    if (@inSuperCellKGrid && !@inSuperCellKShift) {
die "Must specify both kgrid and kshift\n";
    }

    my @superCellTwists;
    getSuperCellTwists($inputFile, $targetSsize, \@inSuperCellKGrid, \@inSuperCellKShift, \@superCellTwists);
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);

    #####################################################################################
    # Parse the pwscf input file for most of the information I will need
    #####################################################################################
    open(IF, $inputFile) || die "cannot open input file $inputFile given on the command line\n";

    my @fdata = <IF>;
    close(IF);

    parsePwscfInput(\@fdata, \$calcPrefix, \$pseudoDir, \$outdir, \$celldim, \$numSpecies, \$numAts,
      \$numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials, \@ionIds, \%atNameToPP, 
                    \@posArray);
    #####################################################################################

    my @primCellKpts;
    print "About to get Primitive Cell K-Points\n";
    my $numsskpts = ($#superCellTwists+1)/4;
    print "Number of supercell k-points: $numsskpts\n"; 
    getPrimCellKpts(\@cell_ptv, \@superCellTwists, \@toptilingmatrix, \@primCellKpts, \@inSuperCellKGrid, \@inSuperCellKShift);

#    my @superCellKpts;
#    getKVectors(\@superCellKpts, \@cell_ptv, \@toptilingmatrix, \@inSuperCellKGrid, \@inSuperCellKShift);
    

    my $kPointsSection = "K_POINTS {crystal}\n";
    my $numkpts = ($#primCellKpts+1)/4;
    $kPointsSection .= "$numkpts\n";
    for (my $i = 0; $i <= $#primCellKpts; $i+=4) {
$kPointsSection .= sprintf ("  %16.14f  %16.14f  %16.14f  %16.14f\n", $primCellKpts[$i],$primCellKpts[$i+1],$primCellKpts[$i+2],$primCellKpts[$i+3]);
    }

    my @pscfData = @fdata;
    genPseudoScfPrep(\@pscfData, $kPointsSection);

    my $pscfFileName = $inputFile;
    $pscfFileName =~ s/\.in//g;
    $pscfFileName =~ s/-scf//g;
    $pscfFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $pscfFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$pscfFileName .= "-S$targetSsize";
    }
    $pscfFileName .= "-pscf.in";

    open(PSCF, ">$pscfFileName");
    print PSCF @pscfData;
    close(PSCF);
   
    my @kzkData = @fdata;
    genKzkPrep(\@kzkData, $kPointsSection);

    my $kzkFileName = $inputFile;
    $kzkFileName =~ s/\.in//g;
    $kzkFileName =~ s/-scf//g;
    $kzkFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $kzkFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$kzkFileName .= "-S$targetSsize";
    }
    $kzkFileName .= "-kzk.in";

    open(KZK, ">$kzkFileName");
    print KZK @kzkData;
    close(KZK);


    my @pw2casino;
    genPw2casino(\@pw2casino, $outdir, $calcPrefix);
    my $pw2casinoFileName = $inputFile;
    $pw2casinoFileName =~ s/\.in//g;
    $pw2casinoFileName =~ s/-scf//g;
    $pw2casinoFileName .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $pw2casinoFileName .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$pw2casinoFileName .= "-S$targetSsize";
    }
    $pw2casinoFileName .= "-pw2casino.in";

    open(PW2CASINO, ">$pw2casinoFileName");
    print PW2CASINO @pw2casino;
    close(PW2CASINO);
}


################################################################################
#
# Branch of code to generate qmcpack input files to test convergence of the
# spacing of the spline mesh.  These will be VMC calculations with no jastrow
# OR, if --withjas is given then we will start from an optimized jastrow as
# found in the optimization 
#
################################################################################
if ($convBspline) {
    my $calcType = 0;
    my $usage = "Must specify all of the keywords:\n";
    $usage .=   "   wvfcnfile, vmcequiltime, vmcdecorrtime, vmctimestep, numsamples, minfactor,\n";
    $usage .=   "   maxfactor and factorinc\n";   
    $usage .=   "when doing a spline convergence test\n\n";
    $usage .=   "Can also give withjas keyword to use the contents of the optimization directory to specify\n";
    $usage .=   "an initial jastrow\n";
        
    if(($wvfcnfile && $vmcequiltime && $vmcdecorrtime && $vmctimestep && $numsamples && $minfactor && $maxfactor && $factorinc)) {
$calcType = 2;
    } else {
die $usage;
    }

###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);

    my ($baseName, @fdata, $topSpinDependentWvfcn);
    my (@atomsCharge, @atomsValence, @atomsAtomicNumbers);
    my ($numUpElec, $numDownElec, $numelec);
#####################################################################################
# In an effort to refactor and simplify, put everything until the parts
# where we are getting qmc input sections into a single function
#####################################################################################

    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);

    getSystemInformation($inputFile, $calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies,
 $numAts, $numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials,
 \@ionIds, \%atNameToPP, \@posArray, $baseName, \@fdata, $topSpinDependentWvfcn,
 \@toptilingmatrix,, $targetSsize, \@atomsCharge, \@atomsValence, \@atomsAtomicNumbers,
 $numelec, $numUpElec, $numDownElec);
    


#########################################################################################
# get the hamiltonian string
#########################################################################################
    my $hamiltonianString = getHamiltonianString(0, \%atNameToPP, 0, 0, 0, $isAe, \@atomsCharge, $hasOpenBC);
#########################################################################################

#########################################################################################
# get the VMC section
#########################################################################################
    my $vmcUseDrift = 0;
    if ($vmcdrift) {
$vmcUseDrift = 1;
    }

    unless ($vmcwalkers) {
$vmcwalkers = 1;
    }
    unless ($vmcblocks) {
$vmcblocks = 32;
    }

    my $vmcSection;
    if ($calcType == 2) {
$vmcSection = getVMCSection($useGPU, $vmctimestep, $vmcequiltime, $vmcdecorrtime,
    $vmcblocks, $vmcwalkers, $numsamples, $vmcUseDrift, 0);
    }

#########################################################################################


#########################################################################################
# get the qmc footer
#########################################################################################
    my $qmcFooterString = getQmcpackFooter();
#########################################################################################

    my $splDirName;
    my $optBaseName = "opt-$baseName";
    if ($targetSsize) {
$splDirName = "bsplineConv-S$targetSsize";
$baseName = $baseName . "-S$targetSsize";
    } else {
$splDirName = "bsplineConv";
    }
    mkdir $splDirName;

    my $optDirSpecies;
    if ($targetSsize) {
$optDirSpecies = "optimization-S$targetSsize";
$optBaseName .= "-S$targetSsize.s001.scalar.dat";
    } else {
$optDirSpecies = "optimization";
$optBaseName .= ".s001.scalar.dat";
    }
    my $lowestJasPerSpeciesSeq = -1;
    my $lowestJasPerSpeciesEn = 100000;
    if (-d $optDirSpecies) {
($lowestJasPerSpeciesSeq, $lowestJasPerSpeciesEn) = findLowestEnergyWvfcn($optDirSpecies, $optBaseName);
    }
    my $lowestJasPerIonSeq = -1;
    my $lowestJasPerIonEn = 100000;

    my $optDirIon;
    if ($targetSsize) {
$optDirIon = "optimization-JastrowPerIon-S$targetSsize";
    } else {
$optDirIon = "optimization-JastrowPerIon";
    }
    if (-d $optDirIon) {
($lowestJasPerIonSeq, $lowestJasPerIonEn) = findLowestEnergyWvfcn($optDirIon, $optBaseName);
    }

    my $lowestseq;
    my $optDir;
    if ($lowestJasPerIonEn < $lowestJasPerSpeciesEn) {
$optDir = $optDirIon;
$lowestseq = $lowestJasPerIonSeq;
$jastrowPerIon = 1;
    } else {
$optDir = $optDirSpecies;
$lowestseq = $lowestJasPerSpeciesSeq;
$jastrowPerIon = 0;
    }

# Loop over mesh factors and write the appropriate qmcpack input files
    if ($massivelyParallel) {
open (FILELIST, ">$splDirName/files");
    }
    for (my $splineFactor = $minfactor; $splineFactor < $maxfactor+0.0001; $splineFactor+= $factorinc) {
$splineFactor = sprintf("%3.2f", $splineFactor);
my $jobid = "$baseName-f" . $splineFactor;
my $fname = "$splDirName/$baseName-f" . $splineFactor . ".xml";

my $jasOptBaseName;
my $wvfcnString;
if ($withjas) {
#    $wvfcnString = getWvfcnStringFromOptDir($optDir, $optBaseName, $twistnum, $splineFactor, $useGPU, $wvfcnfile, $lowestseq, \@toptilingmatrix);
    $wvfcnString = getWvfcnStringFromOptDir2($optDir, $optBaseName, $lowestseq, 0, !$jastrowPerIon, \@atoms_name, $numAts, \@ionIds, "../" . $wvfcnfile, 
     \@toptilingmatrix, $twistnum, $useGPU, 1, $splineFactor, $sglprecsplines, $truncsplines,
     $numUpElec, $numDownElec, $topSpinDependentWvfcn, $isAe, \@atomsCharge);
} else {
    my @dummyArr;
    $wvfcnString = getWavefunctionString(0, 0, 0, 0, $toprcutTwoBody, 
 \@dummyArr, \@dummyArr, \@atoms_name,
 $numAts, \@ionIds, 0, $toprcutOneBody, \@dummyArr, "../" . $wvfcnfile, 
 \@toptilingmatrix, $twistnum, $useGPU, 1,
 $splineFactor, $sglprecsplines, $truncsplines, $numUpElec, 
 $numDownElec, $topSpinDependentWvfcn, $isAe, \@atomsCharge);
}
my $qmcHeaderString = getQmcpackHeader($jobid, 1, "Mesh Size Convergence Test for $baseName with factor = $splineFactor", 49154);

########################################################################################
# Get the particleset string.  
# TODO : Make this so that an auxiliary particleset is not generated when generateComplexWavefunction == 0
########################################################################################
    my $ptclsetString = getPtclSetString(0, !$jastrowPerIon, \@cell_ptv, $numUpElec,
 $numDownElec, $numAts, $#atoms_name, \@atoms_name, \@atomsCharge,
 \@atomsValence, \@atomsAtomicNumbers, \@posArray,
 \@ionIds, \@bconds);
########################################################################################

my $qmcFile = $qmcHeaderString . $ptclsetString . $wvfcnString . $hamiltonianString . "\n\n". $vmcSection . $qmcFooterString;

open(CONVFILE, ">$fname");
print CONVFILE $qmcFile;
close(CONVFILE);
if ($massivelyParallel) {
    my $locfname = "$baseName-f" . $splineFactor . ".xml";
    print FILELIST "$locfname\n";
}
    }
    if ($massivelyParallel) {
close(FILELIST);
    }
}


####################################################################################
#
# Branch of code to generate dmc input file using an already completed optimization
#
# LNS TODO: Need to update this to take advantage of the symmetry information
#
####################################################################################

if ($dmcCalc) {
    my $calcType = 0;
    if (($wvfcnfile && $vmcequiltime && $vmcdecorrtime && $vmctimestep && $dmctstep && $dmcequiltime &&
$dmcruntime && $dmcblocktime && $targetpop)) {
$calcType = 2;
    } else {
die "Must specify all of the keywords: wvfcnfile, vmctimestep, vmcequiltime, vmcdecorrtime, dmctstsp, dmcequiltime, dmcruntime\ndmcblocktime and targetpop when doing a DMC calculation\n";
    }


###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);

    my ($baseName, @fdata, $topSpinDependentWvfcn);
    my (@atomsCharge, @atomsValence, @atomsAtomicNumbers);
    my ($numUpElec, $numDownElec, $numelec, $numSupercellTwists);
#####################################################################################
# In an effort to refactor and simplify, put everything until the parts
# where we are getting qmc input sections into a single function
#####################################################################################
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);
    
    print "Target supercell size is $targetSsize\n";

    getSystemInformation($inputFile, $calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies,
 $numAts, $numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials,
 \@ionIds, \%atNameToPP, \@posArray, $baseName, \@fdata, $topSpinDependentWvfcn,
 \@toptilingmatrix, $targetSsize, \@atomsCharge, \@atomsValence, \@atomsAtomicNumbers,
 $numelec, $numUpElec, $numDownElec);

    print "baseName = $baseName\n";
    
    my @superCellTwists;
    getSuperCellTwists($inputFile, $targetSsize, \@inSuperCellKGrid, \@inSuperCellKShift, \@superCellTwists);
    $numSupercellTwists = ($#superCellTwists+1)/4;
#########################################################################################
# Grab optimized wavefunction by inspecting the output of optimization runs
#########################################################################################

    my $optDirSpecies;
    my $optBaseName = "opt-$baseName";
    if ($targetSsize) {
$optDirSpecies = "optimization-S$targetSsize";
$optBaseName .= "-S$targetSsize.s001.scalar.dat";
    } else {
$optDirSpecies = "optimization";
$optBaseName .= ".s001.scalar.dat";
    }
    my $lowestJasPerSpeciesSeq = -1;
    my $lowestJasPerSpeciesEn = 100000;
    if (-d $optDirSpecies) {
($lowestJasPerSpeciesSeq, $lowestJasPerSpeciesEn) = findLowestEnergyWvfcn($optDirSpecies, $optBaseName);
    }
    my $lowestJasPerIonSeq = -1;
    my $lowestJasPerIonEn = 100000;

    my $optDirIon;
    if ($targetSsize) {
$optDirIon = "optimization-JastrowPerIon-S$targetSsize";
    } else {
$optDirIon = "optimization-JastrowPerIon";
    }
    if (-d $optDirIon) {
($lowestJasPerIonSeq, $lowestJasPerIonEn) = findLowestEnergyWvfcn($optDirIon, $optBaseName);
    }

    my $lowestseq;
    my $optDir;
    if ($lowestJasPerIonEn < $lowestJasPerSpeciesEn) {
$optDir = $optDirIon;
$lowestseq = $lowestJasPerIonSeq;
$jastrowPerIon = 1;
    } else {
$optDir = $optDirSpecies;
$lowestseq = $lowestJasPerSpeciesSeq;
$jastrowPerIon = 0;
    }

    
    my @wvfcnStrings;
    
    for (my $i = 0; $i < $numSupercellTwists; $i++) {
      my @kpoint = ($superCellTwists[$i*4], $superCellTwists[$i*4+1], $superCellTwists[$i*4+2]);
      #my $wvfcnString = getWvfcnStringFromOptDir($optDir, $optBaseName, $i, $splfactor, $useGPU, $wvfcnfile, $lowestseq, \@toptilingmatrix);
      my $wvfcnString = getWvfcnStringFromOptDir2($optDir, $optBaseName, $lowestseq, 0, !$jastrowPerIon, \@atoms_name, $numAts, \@ionIds, "../" . $wvfcnfile, 
  \@toptilingmatrix, -1, $useGPU, 1, $splfactor, $sglprecsplines, $truncsplines,
  $numUpElec, $numDownElec, $topSpinDependentWvfcn, $isAe, \@atomsCharge, \@kpoint);
      push(@wvfcnStrings, $wvfcnString);
    }
#########################################################################################
# get the hamiltonian string
#########################################################################################
    my $hamiltonianString = getHamiltonianString(0, \%atNameToPP, !$hasOpenBC, 0, 1, $isAe, \@atomsCharge);
#########################################################################################

#########################################################################################
# Get the particleset string.  
# TODO : Make this so that an auxiliary particleset is not generated when generateComplexWavefunction == 0
# TODO : Test that this works with a non-identity tiling matrix generated above
#########################################################################################
    my $ptclsetString = getPtclSetString(0, !$jastrowPerIon, \@cell_ptv, $numUpElec,
 $numDownElec, $numAts, $#atoms_name, \@atoms_name, \@atomsCharge,
 \@atomsValence, \@atomsAtomicNumbers, \@posArray,
 \@ionIds, \@bconds);
#########################################################################################

#########################################################################################
# get the VMC section
#########################################################################################
    my $vmcUseDrift = 0;
    if ($vmcdrift) {
$vmcUseDrift = 1;
    }
    unless ($vmcwalkers) {
$vmcwalkers = 1;
    }

    my $vmcSection;
    if ($calcType == 2) {
$vmcSection = getVMCSection($useGPU, $vmctimestep, $vmcequiltime, $vmcdecorrtime,
    1, $vmcwalkers, $vmcwalkers, $vmcUseDrift, 1);
    }
#########################################################################################


#########################################################################################
# get the qmc header and footer
#########################################################################################
    my @qmcHeaderStrings;
    if ($targetSsize) {
$baseName = $baseName . "-S$targetSsize";
    }

    for (my $i = 0; $i < $numSupercellTwists; $i++) {
my $jobid = "$baseName-dmc";
if ($numSupercellTwists > 1) {
    $jobid .= "-tw$i";
}
my $qmcHeaderString = getQmcpackHeader($jobid, 1, "DMC for $baseName-tw$i", 49154);
push(@qmcHeaderStrings, $qmcHeaderString);
    }
    my $qmcFooterString = getQmcpackFooter();
#########################################################################################

    my $dmcdirname = "dmc";
    $dmcdirname .= "-supertwist${inSuperCellKGrid[0]}${inSuperCellKGrid[1]}${inSuperCellKGrid[2]}";
    $dmcdirname .= "-supershift${inSuperCellKShift[0]}${inSuperCellKShift[1]}${inSuperCellKShift[2]}";
    if ($targetSsize) {
$dmcdirname .= "-S$targetSsize";
    } 
    mkdir "$dmcdirname";
    

    my $dmcSection;
    ## Do three warmup sections where we run for $dmcequiltime using stochastic reconfiguration
    ## larget timesteps and one walker per thread

    my $firsttstep = 0.1;
    my $secondtstep  = 0.05;
    my $thirdtstep = $dmctstep*2;
    
    $dmcSection .= getDMCSectionNew($useGPU, $firsttstep, $dmcequiltime, $dmcequiltime,
    $firsttstep, $dmcUseTmoves, 0, 0, 1); 
    $dmcSection .= getDMCSectionNew($useGPU, $secondtstep, $dmcequiltime, $dmcequiltime,
    $secondtstep, $dmcUseTmoves, 0, 0, 1); 
    $dmcSection .= getDMCSectionNew($useGPU, $thirdtstep, $dmcequiltime, $dmcequiltime,
    $thirdtstep, $dmcUseTmoves, 0, 0, 1); 

    ## Now do a final run where we have increased the number of walkers and use the normal
    ## population control algorithm

    if ($massivelyParallel) {
# for massively parallel simulation, use the one walker per thread provided.
# note this can be dangerous, if population size bias is not zero then different
# runs with different numbers of processors will produce different answers!!
$dmcSection .= getDMCSectionNew($useGPU, $dmctstep, $dmcequiltime, $dmcruntime,
$dmcblocktime, $dmcUseTmoves, 0, $targetpop, 0);
    } else {
$dmcSection .= getDMCSectionNew($useGPU, $dmctstep, $dmcequiltime, $dmcruntime,
$dmcblocktime, $dmcUseTmoves, $targetpop, 0, 0);
    }

    if ($massivelyParallel) {
open (FILELIST, ">$dmcdirname/files");
    }
    for (my $i = 0; $i < $numSupercellTwists; $i++) {
my $qmcFile = $qmcHeaderStrings[$i] . $ptclsetString . $wvfcnStrings[$i] . $hamiltonianString;
$qmcFile .= "\n\n" . $vmcSection . $dmcSection . $qmcFooterString;

my $fname;
if ($numSupercellTwists > 1) {
    $fname = "$dmcdirname/$baseName-dmc-tw$i.xml";
} else {
    $fname = "$dmcdirname/$baseName-dmc.xml";
}
open(DMC, ">$fname");
print DMC $qmcFile;
close(DMC);
if ($massivelyParallel) {
    my $locfname = "$baseName-dmc-tw$i.xml";
    print FILELIST "$locfname\n";
}
    }
    if ($massivelyParallel) {
close(FILELIST);
    }
}

    

################################################################################
#
#
# Branch of code to generate qmcpack input files to optimize jastrow factors
# currently blindly using Jeremy's recommended block for Al.
# TODO: Add facility to use rescaling for optimization rather than quartic.
#       This will require decoupling optsamples from the number of vmc steps
#       so that H and S matrices can be relatively well converged without having
#       too many samples.  This will be extermely useful on machines which have
#       small amounts of memory.
#
#
################################################################################
if ($optwvfcn) {
    my $usage = "Must specify all of the keywords:\n";
    $usage   .= "   wvfcnfile, splfactor, vmctimestep, vmcequiltime, vmcdecorrtime\n";
    $usage   .= "   optsamples, optloops, onebodysplinepts and twobodysplinepts\n";
    $usage   .= "when optimizing a wavefunction\n";

    my $calcType = 0;
    if (($wvfcnfile && $splfactor && $vmctimestep, $vmcequiltime, $vmcdecorrtime
 && $optSamples && $optLoops && $oneBodySplinePts && $twoBodySplinePts )) {
      $calcType = 2;
    } else {
      die $usage;
    }
    
    if ($hasOpenBC && (!$halfboxrcut && (!$toprcutTwoBody || !$toprcutOneBody))) {
      die "If one of the boundary conditions is open (n instead of p)\n Then you must specify either halfboxrcut or both onebodyrcut and twobodyrcut\n";
    }

###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);

    my ($baseName, @fdata, $topSpinDependentWvfcn);
    my (@atomsCharge, @atomsValence, @atomsAtomicNumbers);
    my ($numUpElec, $numDownElec, $numelec);
#####################################################################################
# In an effort to refactor and simplify, put everything until the parts
# where we are getting qmc input sections into a single function
#####################################################################################
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);
    print "targetSsize = $targetSsize\n";

    getSystemInformation($inputFile, $calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies,
 $numAts, $numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials,
 \@ionIds, \%atNameToPP, \@posArray, $baseName, \@fdata, $topSpinDependentWvfcn,
 \@toptilingmatrix, $targetSsize, \@atomsCharge, \@atomsValence, \@atomsAtomicNumbers,
 $numelec, $numUpElec, $numDownElec);

    print "got system information\n";
########################################################################################
# Get the particleset string.  
# TODO : Make this so that an auxiliary particleset is not generated when generateComplexWavefunction == 0?
# TODO : Test that this works with a non-identity tiling matrix generated above
# TODO : Add support so that this will generate a complex wavefunction
########################################################################################
    my $ptclsetString = getPtclSetString(0, !$jastrowPerIon, \@cell_ptv, $numUpElec,
 $numDownElec, $numAts, $#atoms_name, \@atoms_name, \@atomsCharge,
 \@atomsValence, \@atomsAtomicNumbers, \@posArray,
 \@ionIds, \@bconds);
########################################################################################

#########################################################################################
# get the hamiltonian string
#########################################################################################
    my $hamiltonianString = getHamiltonianString(0, \%atNameToPP, !$hasOpenBC, 0, 1, $isAe, \@atomsCharge);
#########################################################################################

#########################################################################################
# get the wavefunction string
# TODO: Add support for wavefunctions where there is an independent jastrow
#       per atom, not per species (complexWavefunction in the subroutine's nomenclature)
#########################################################################################
    my @topUpUpCoefs;
    my @topUpDownCoefs;
    
    
    for (my $i = 0; $i < $twoBodySplinePts; $i++) {
push(@topUpUpCoefs, 0.0);
push(@topUpDownCoefs, 0.0);
    }
    my $twoBodyRcut = 0.0;
    my $numDens = 0.0;
    getCellProperties($numelec, \@cell_ptv, $twoBodyRcut, $numDens);
    if ($halfboxrcut) {
      $toprcutTwoBody = $twoBodyRcut;
      $toprcutOneBody = $twoBodyRcut;
    }
    getTwoBodyRPAJastrow($numDens, $twoBodyRcut, $twoBodySplinePts, \@topUpUpCoefs, \@topUpDownCoefs);
    
#    for (my $i = 0; $i < $twoBodySplinePts; $i++) {
#print "uu[$i] = $topUpUpCoefs[$i], ud[$i] = $topUpDownCoefs[$i]\n";
#    }

    my @topJastrowStarts;
    for (my $i = 0; $i <= $#atoms_name; $i++) {
for (my $j = 0; $j < $oneBodySplinePts; $j++) {
    push (@topJastrowStarts, 0.0);
}
    }

    my $wvfcnString;
    if ($jastrowPerIon) {
# check to see if optimization has already been done with jastrow per species, if so, start from there
my $optDir;
my $optBaseName;
if ($targetSsize) {
    $optDir = "optimization-S$targetSsize";
    $optBaseName = "opt-$baseName-S$targetSsize.s001.scalar.dat";
} else {
    $optDir = "optimization";
    $optBaseName = "opt-$baseName.s001.scalar.dat";
}
if (-d $optDir) {
    $wvfcnString = getWvfcnStringFromOptDir2($optDir, $optBaseName, 0, 0, !$jastrowPerIon, \@atoms_name, $numAts, \@ionIds, "../" . $wvfcnfile, 
     \@toptilingmatrix, $twistnum, $useGPU, 1, $splfactor, $sglprecsplines, $truncsplines,
     $numUpElec, $numDownElec, $topSpinDependentWvfcn, $isAe, \@atomsCharge);
} else {
    $wvfcnString = getWavefunctionString(0, 1, !$jastrowPerIon, $twoBodySplinePts, $toprcutTwoBody, \@topUpUpCoefs,
 \@topUpDownCoefs, \@atoms_name, $numAts, \@ionIds,
 $oneBodySplinePts, $toprcutOneBody, \@topJastrowStarts,
 "../" . $wvfcnfile, \@toptilingmatrix, $twistnum, $useGPU,
 1, $splfactor, $sglprecsplines, $truncsplines, $numUpElec, $numDownElec, 
 $topSpinDependentWvfcn, $isAe, \@atomsCharge);
}
    } else {
$wvfcnString = getWavefunctionString(0, 1, !$jastrowPerIon, $twoBodySplinePts, $toprcutTwoBody, \@topUpUpCoefs,
     \@topUpDownCoefs, \@atoms_name, $numAts, \@ionIds,
     $oneBodySplinePts, $toprcutOneBody, \@topJastrowStarts,
     "../" . $wvfcnfile, \@toptilingmatrix, $twistnum, $useGPU,
     1, $splfactor, $sglprecsplines, $truncsplines, $numUpElec, $numDownElec, 
     $topSpinDependentWvfcn, $isAe, \@atomsCharge);
    }

#########################################################################################
# get the qmc header and footer
#########################################################################################
    my $qmcFooterString = getQmcpackFooter();
    if ($targetSsize) {
$baseName = "$baseName-S$targetSsize";
    }
    my $jobid = "opt-$baseName";
    my $qmcHeaderString = getQmcpackHeader($jobid, 1, "Optimization of jastrows for $baseName", 49154 );
#########################################################################################   


#########################################################################################
# Get the optimization section for this job
#########################################################################################
    my $vmcdrift = 0;
    if ($vmcdrift) {
$vmcdrift = 1;
    }
    
    unless($vmcwalkers) {
$vmcwalkers = 1;
    }
    unless($vmcblocks) {
$vmcblocks = 500;
    }

    my $optimizationString = getOptSection($useGPU, $vmcwalkers, $vmctimestep, $vmcequiltime,
   $vmcdecorrtime, $vmcblocks, $optSamples, $optLoops, $vmcdrift);

#########################################################################################
# Put it all together and write to a file
#########################################################################################
    my $dirname;
    if ($jastrowPerIon) {
if ($targetSsize) {
    $dirname = "optimization-JastrowPerIon-S$targetSsize";
} else {
    $dirname = "optimization-JastrowPerIon";
}
    } else {
if ($targetSsize) {
    $dirname = "optimization-S$targetSsize";
} else {
    $dirname = "optimization";
}
    }
    mkdir $dirname;
    my $optFileName = "$dirname/opt-$baseName.xml";
    my $optFileString = $qmcHeaderString . $ptclsetString . $wvfcnString;
    $optFileString .= $hamiltonianString . "\n\n" . $optimizationString . $qmcFooterString;
    
    open (OPTFILE, ">$optFileName") || die "Cannot open file $optFileName\n";
    print OPTFILE $optFileString;
    close(OPTFILE);
}



################################################################################
#
# Branch of code to generate qmcpack input files to test convergence of the
# DMC timestep.  
#
################################################################################
if ($convDMCTstep) {
    my $calcType = 0;
    if (($wvfcnfile && $vmcequiltime && $vmcdecorrtime && $vmctimestep && $dmcequiltime &&
 $dmcruntime && $dmcblocktime && $targetpop && $mindmctstep && $maxdmctstep && $dmctstepinc)) {
      $calcType = 2;
    } else {
die "Must specify all of the keywords: \nwvfcnfile, vmcequiltime, vmcdecorrtime, vmctimestep, dmcequiltime, dmcruntime\ndmcblocktime, targetpop, maxdmctstep, dmctstepinc\nwhen doing a DMC timestep convergence test\n";
    } 
    
###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);

    my ($baseName, @fdata, $topSpinDependentWvfcn);
    my (@atomsCharge, @atomsValence, @atomsAtomicNumbers);
    my ($numUpElec, $numDownElec, $numelec);
#####################################################################################
# In an effort to refactor and simplify, put everything until the parts
# where we are getting qmc input sections into a single function
#####################################################################################
    findSupercell($inputFile, $targetSsize, \@toptilingmatrix);

    getSystemInformation($inputFile, $calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies,
 $numAts, $numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials,
 \@ionIds, \%atNameToPP, \@posArray, $baseName, \@fdata, $topSpinDependentWvfcn,
 \@toptilingmatrix, $targetSsize, \@atomsCharge, \@atomsValence, \@atomsAtomicNumbers,
 $numelec, $numUpElec, $numDownElec);
    

#########################################################################################
# Grab optimized wavefunction by inspecting the output of optimization runs
#########################################################################################

    my $optDirSpecies;
    my $optBaseName = "opt-$baseName";
    if ($targetSsize) {
$optDirSpecies = "optimization-S$targetSsize";
$optBaseName .= "-S$targetSsize.s001.scalar.dat";
    } else {
$optDirSpecies = "optimization";
$optBaseName .= ".s001.scalar.dat";
    }
    my $lowestJasPerSpeciesSeq = -1;
    my $lowestJasPerSpeciesEn = 100000;
    if (-d $optDirSpecies) {
($lowestJasPerSpeciesSeq, $lowestJasPerSpeciesEn) = findLowestEnergyWvfcn($optDirSpecies, $optBaseName);
    }
    my $lowestJasPerIonSeq = -1;
    my $lowestJasPerIonEn = 100000;

    my $optDirIon;
    if ($targetSsize) {
$optDirIon = "optimization-JastrowPerIon-S$targetSsize";
    } else {
$optDirIon = "optimization-JastrowPerIon";
    }
    if (-d $optDirIon) {
($lowestJasPerIonSeq, $lowestJasPerIonEn) = findLowestEnergyWvfcn($optDirIon, $optBaseName);
    }

    my $lowestseq;
    my $optDir;
    if ($lowestJasPerIonEn < $lowestJasPerSpeciesEn) {
$optDir = $optDirIon;
$lowestseq = $lowestJasPerIonSeq;
$jastrowPerIon = 1;
    } else {
$optDir = $optDirSpecies;
$lowestseq = $lowestJasPerSpeciesSeq;
$jastrowPerIon = 0;
    }

    print "Using optDir = $optDir, and sequence = $lowestseq\n";

    #my $wvfcnString = getWvfcnStringFromOptDir($optDir, $optBaseName, $twistnum, $splfactor, $useGPU, $wvfcnfile, $lowestseq, \@toptilingmatrix);
    my $wvfcnString = getWvfcnStringFromOptDir2($optDir, $optBaseName, $lowestseq, 0, !$jastrowPerIon, \@atoms_name, $numAts, \@ionIds, "../". $wvfcnfile, 
\@toptilingmatrix, $twistnum, $useGPU, 1, $splfactor, $sglprecsplines, $truncsplines,
$numUpElec, $numDownElec, $topSpinDependentWvfcn, $isAe, \@atomsCharge);

#########################################################################################
# Get the particleset string.  
# TODO : Make this so that an auxiliary particleset is not generated when generateComplexWavefunction == 0
# TODO : Test that this works with a non-identity tiling matrix generated above
#########################################################################################
    my $ptclsetString = getPtclSetString(0, !$jastrowPerIon, \@cell_ptv, $numUpElec,
 $numDownElec, $numAts, $#atoms_name, \@atoms_name, \@atomsCharge,
 \@atomsValence, \@atomsAtomicNumbers, \@posArray,
 \@ionIds, \@bconds);
#########################################################################################

#########################################################################################
# get the hamiltonian string
#########################################################################################
    my $hamiltonianString = getHamiltonianString(0, \%atNameToPP, !$hasOpenBC, 0, 1, $isAe, \@atomsCharge);
#########################################################################################

#########################################################################################
# get the VMC section
#########################################################################################
    my $vmcUseDrift = 0;
    if ($vmcdrift) {
$vmcUseDrift = 1;
    }
    unless ($vmcwalkers) {
$vmcwalkers = 1;
    }
    unless ($vmcblocks) {
$vmcblocks = 32;
    }

    my $vmcSection;

    if ($massivelyParallel) {
$vmcSection = getVMCSection($useGPU, $vmctimestep, $vmcequiltime, $vmcdecorrtime,
    $vmcblocks, $vmcwalkers, 0, $vmcUseDrift, 1);
    } else {
$vmcSection = getVMCSection($useGPU, $vmctimestep, $vmcequiltime, $vmcdecorrtime,
    $vmcblocks, $vmcwalkers, $targetpop, $vmcUseDrift, 0);
    }


#########################################################################################


#########################################################################################
# get the qmc header and footer
#########################################################################################
    my $dmcdirname;
    if ($targetSsize) {
$baseName .= "-S$targetSsize";
$dmcdirname = "dmcTstepConv-S$targetSsize";
    } else {
$dmcdirname = "dmcTstepConv";
    }

#########################################################################################

    
    mkdir "$dmcdirname";

   

    if ($massivelyParallel) {
open (FILELIST, ">$dmcdirname/files");
for (my $locdmctstep = $maxdmctstep; $locdmctstep > $mindmctstep-0.000001; $locdmctstep -= $dmctstepinc) {
    my $jobid = "$baseName-dmcTsteps-tstep$locdmctstep";
    my $qmcHeaderString = getQmcpackHeader($jobid, 1, "DMC Timestep convergence test for $baseName, using timestep $locdmctstep", 49154);
    my $qmcFooterString = getQmcpackFooter();
    
    my $qmcFile = $qmcHeaderString . $ptclsetString . $wvfcnString. $hamiltonianString . "\n\n" . $vmcSection;    
    my $dmcSection = getDMCSectionNew($useGPU, $locdmctstep, $dmcequiltime, $dmcruntime,
       $dmcblocktime, $dmcUseTmoves, 0, $targetpop, 0);
    $qmcFile .= $dmcSection . $qmcFooterString;
    my $locfname = "$baseName-dmcTsteps-tstep$locdmctstep.xml";
    my $fname = "$dmcdirname/$locfname";
    open(DMCTSTEP, ">$fname");
    print DMCTSTEP $qmcFile;
    close(DMCTSTEP);
    print FILELIST "$locfname\n";
}
close(FILELIST);
    } else {
my $jobid = "$baseName-dmcTsteps";
my $qmcHeaderString = getQmcpackHeader($jobid, 1, "DMC Timestep convergence test for $baseName", 49154);
my $qmcFooterString = getQmcpackFooter();

my $qmcFile = $qmcHeaderString . $ptclsetString . $wvfcnString. $hamiltonianString . "\n\n" . $vmcSection;
for (my $locdmctstep = $maxdmctstep; $locdmctstep > $mindmctstep-0.000001; $locdmctstep -= $dmctstepinc) {
    my $dmcSection;
    if ($calcType == 2) {
$dmcSection = getDMCSectionNew($useGPU, $locdmctstep, $dmcequiltime, $dmcruntime,
       $dmcblocktime, $dmcUseTmoves, $targetpop, 0, 0);
    }
    $qmcFile .= $dmcSection;
}
$qmcFile .= $qmcFooterString;

my $fname = "$dmcdirname/$baseName-dmcTsteps.xml";
open(DMCTSTEP, ">$fname");
print DMCTSTEP $qmcFile;
close(DMCTSTEP);
    }
}


##########################################################################################
##########################################################################################
##########################################################################################
# End of main program, start of helper routines
##########################################################################################
##########################################################################################
##########################################################################################

#####################################################################################
# Subroutine to get the appropriate tiling matrix for the requested supercell size
#####################################################################################
sub findSupercell {
  my $ifile = shift;
  my $targetSsize = shift;
  my $tilmat = shift;
  
  `touch SYMMETRY`;
  my $line = `sed -n '/TILEMAT $targetSsize :/=' SYMMETRY`;
  # Check to see if symmetry file already contains the tiling matrix for the supercell size
  if ($line) {
    chomp($line);
    my $line2 = `sed -n ${line}p SYMMETRY`;
    my @data = split(/\s+/, $line2);
    for (my $i = 3; $i < 12; $i++) {
      $$tilmat[$i-3] = $data[$i];
    } 
    # otherwise, use the supercell tool to look for a tiling matrix
  } else {
    open(IF, $inputFile) || die "cannot open input file $inputFile given on the command line\n";
    
    my @fdata = <IF>;
    close(IF);

    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies, $numAts);
    my ($numSpins, @cell_ptv, @atoms_name, @pseudoPotentials, @ionIds);
    my (%atNameToPP, @posArray);

    parsePwscfInput(\@fdata, \$calcPrefix, \$pseudoDir, \$outdir, \$celldim, \$numSpecies, \$numAts,
    \$numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials, \@ionIds, \%atNameToPP, 
    \@posArray);


    my $getSupercell = $config{supercell};
    my $out = `$getSupercell --ptvs @cell_ptv --target $targetSsize --maxentry 7`;
    my @data = split(/\s+/, $out);
    for (my $i = 1; $i < 10; $i++) {
      $$tilmat[$i-1] = $data[$i];
    }

    open (file, '>>SYMMETRY');
    print file "TILEMAT $targetSsize : @{$tilmat}\n";
    close(file);
  }
}

#####################################################################################
# Subroutine to get all system information
#####################################################################################
sub getSystemInformation {
    my $ifile = \shift;
    my $calcprefix = \shift;
    my $psdir = \shift;
    my $odir = \shift;
    my $cdim = \shift;
    my $nspecies = \shift;
    my $nats = \shift;
    my $nspins = \shift;
    my $cptvref = shift;
    my $anameref = shift;
    my $ppsref = shift;
    my $iidref  = shift;
    my $atnmtoppref = shift;
    my $posarrref = shift;
    my $bname = \shift;
    my $fdataref = shift;
    my $tspdepwvfcn = \shift;
    my $ttilmatref = shift;
    my $tss = shift;
    my $atchgref = shift;
    my $atvalref = shift;
    my $atatnumref = shift;
    my $nelec = \shift;
    my $nupelec = \shift;
    my $ndnelec = \shift;

#####################################################################################
# Parse the pwscf input file for most of the information I will need
#####################################################################################
    $$bname = $$ifile;
    $$bname =~ s/\.in//g;
    $$bname =~ s/-scf//g;


    open(IF, $$ifile) || die "cannot open input file $inputFile given on the command line\n";
    @{$fdataref} = <IF>;
    close(IF);
    
    parsePwscfInput($fdataref, $calcprefix, $psdir, $odir, $cdim, $nspecies, $nats,
      $nspins, $cptvref, $anameref, $ppsref, $iidref, $atnmtoppref, 
                    $posarrref);


    if ($$nspins > 1) {
$$tspdepwvfcn = 1;
    } else {
$$tspdepwvfcn = 0;
    }

#####################################################################################
#  Figure out best tiling matrix for given supercell size or take tiling matrix
#  from input line.  If using a supercell, update relevant variables from above
#####################################################################################
    if ($tss) {
      
      `touch SYMMETRY`;
      my $line = `sed -n '/TILEMAT $targetSsize :/=' SYMMETRY`;
      if ((@{$ttilmatref})) {
      } elsif ($line) {
chomp($line);
my $line2 = `sed -n ${line}p SYMMETRY`;
my @data = split(/\s+/, $line2);
for (my $i = 3; $i < 12; $i++) {
  $$ttilmatref[$i-3] = $data[$i];
} 
      } else {
# otherwise, use the supercell tool to look for a tiling matrix
my $getSupercell = $config{supercell};
my $out = `$getSupercell --ptvs @{$cptvref} --target $tss --maxentry 7`;
my @data = split(/\s+/, $out);
## Set tilematrix from output of getSupercell
for (my $i = 1; $i < 10; $i++) {
  $$ttilmatref[$i-1] = $data[$i];
}
      }
      
## Get new arrays for iidref and posarrref
my @newIonIds;
my @newIonPositions;
for (my $i = 0; $i < $$nats; $i++) {
    my @tempIonPositions;

    ## Repeat the name of each ion in ionids tss times
    for (my $j = 0; $j < $tss; $j++) {
push(@newIonIds, $$iidref[$i]);
    }
    
    ## Get the positions of all copies of this particular ion in the supercell
    ## and push them onto the newIonPositions array
    my @locbasis;
    for (my $j = 0; $j < 3; $j++) {
$locbasis[$j] = $$posarrref[$i*3+$j];
    }
    my @positions;
    getSupercellPos(\@tempIonPositions, \@locbasis, $cptvref, $ttilmatref);
    push(@newIonPositions, @tempIonPositions);
}

@{$iidref} = @newIonIds;
@{$posarrref} = @newIonPositions;

## Update the number of ions in the supercell
$$nats *= $tss;
## Set ptv to be the supercell's lattice vectors from the output of getSupercell
#for (my $i = 10; $i < 19; $i++) {
#    $$cptvref[$i-10] = $data[$i];
#}
my @tmpssptv;
getSuperCell(\@tmpssptv, $cptvref, $ttilmatref);
@{$cptvref} = @tmpssptv;

    }



    if (!(@{$ttilmatref})) {
@{$ttilmatref} = (1, 0, 0, 0, 1, 0, 0, 0, 1);
    }


###########################################################################################
# Need to scrape out of the pseudopotentials the charge of each ion, the number of valence
# electrons, and the atomic numbers, do this for fhi and for ncpp pseudopotentials
# Also generate xml input files and update atNameToPP to point to these files.
###########################################################################################
    my $curdir = `pwd`;
    chomp($curdir);
    foreach my $at (@{$iidref}) {
my $ppname = $$atnmtoppref{$at};
my $atValenceCharge;
my $atAtomicNumber;

getPPInfo($ppname, $atValenceCharge, $atAtomicNumber);

$$nelec += $atValenceCharge;
#push(@{$atchgref}, $atValenceCharge);
#push(@{$atvalref}, $atValenceCharge);
#push(@{$atatnumref}, $atAtomicNumber);
    }   

    foreach my $at (@{$anameref}) {
my $ppname = $$atnmtoppref{$at};
my $atValenceCharge;
my $atAtomicNumber;

getPPInfo($ppname, $atValenceCharge, $atAtomicNumber);
push(@{$atchgref}, $atValenceCharge);
push(@{$atvalref}, $atValenceCharge);
push(@{$atatnumref}, $atAtomicNumber);

if ($ppname =~ /ncpp/i) {
    $ppname =~ s/ncpp/xml/;
    my $ppbasename = $ppname;
    $ppbasename =~ /.*\/(.*\.xml)/g;
    $ppbasename = $1;
    my $fname = $curdir . "/$ppbasename";
    unless (-e $fname) {
print "wfconvert cannot convert from ncpp potentials to xml, make sure that\npotentials of the same basename but with xml format are in this directory.\n";
    }
    $$atnmtoppref{$at} = $fname;
} elsif ($ppname =~ /upf/i) {
    # convert upf potential to xml
    print "converting upf potential to xml with ppconvert\n";
    my $ppbasename = $ppname;    
    $ppbasename =~ s/upf/xml/;
    $ppbasename =~ /.*\/(.*\.xml)/g;
    $ppbasename = $1;
    `$config{ppconvert} --upf_pot $ppname --xml $ppbasename`;
    $$atnmtoppref{$at} = $curdir . "/$ppbasename";
}
    }

###########################################################################################

###########################################################################################
# Figure out how many of each type of electron to use.  Eventually could have functionality
# so that after the scf calculation has run the magnetic state is scraped from the dft run
# for now with no input, assume no net magnetization.  Also allow nelUp - nelDown to be input
# on the command line
# TODO : Add support for nelUp-nelDown from command line.  
# TODO : Add functionality to get this from dft run when scf is run ahead of time
# TODO : Make this work properly when tiling matrix is not unity
###########################################################################################
    my $remainder = $$nelec % 2;

    if ($remainder) {
die "There are an odd number of electrons in the cell.  I don't know what to do!\n";
    } 
    $$nupelec = $$nelec/2;
    $$ndnelec = $$nelec/2;
###########################################################################################
}



##########################################################################################
# Subroutine to scrape the system information from the pwscf scf input file
##########################################################################################
sub parsePwscfInput {
    my $filedataRef = shift;

    my $cp = shift;
    my $psdir = shift;
    my $odir = shift;
    my $cdim = shift;
    my $numSp = shift;
    my $numAt = shift;
    my $numSpin = shift;
    my $cptvRef = shift;
    my $atNamRef = shift;
    my $ppsRef = shift;
    my $ionsRef = shift;
    my $atnmtoppRef = shift;
    my $posarrRef = shift;
    my $coadim;
    
    ${$cp} = getPwscfToken("prefix", $filedataRef);
    ${$psdir} = getPwscfToken("pseudo_dir", $filedataRef);
    ${$odir} = getPwscfToken("outdir", $filedataRef);
    ${$cdim} = getPwscfToken('celldm(1)', $filedataRef);
    if (${$cdim} < 0) {
${$cdim} = 1.0;
    }
    $coadim = getPwscfToken('celldm(3)', $filedataRef);
    ${$numSp} = getPwscfToken("ntyp", $filedataRef);
    ${$numAt} = getPwscfToken("nat", $filedataRef);
    ${$numSpin} = getPwscfToken("nspin", $filedataRef);

    my $ibrav = getPwscfToken("ibrav", $filedataRef);
    if ($ibrav == 0) {
my $str = getPwscfCard("CELL_PARAMETERS", $filedataRef);
@{$cptvRef} = split(/\s+/, removeWhiteSpaceAtBeginning($str));
for (my $i = 0; $i < 9; $i++) {
    $cptvRef->[$i] *= ${$cdim};
}
    } elsif ($ibrav == 1) {
@{$cptvRef} = (${$cdim}, 0, 0, 0, ${$cdim}, 0, 0, 0, ${$cdim});
    } elsif ($ibrav == 2) {
my $val = 0.5*${$cdim};
@{$cptvRef} = (-$val, 0, $val, 0, $val, $val, -$val, $val, 0);
    } elsif ($ibrav == 3) {
my $val = 0.5*${$cdim};
@{$cptvRef} = ($val, $val, $val, -$val, $val, $val, -$val, -$val, $val);
    } elsif ($ibrav == 4) {
my $val = ${$cdim};
@{$cptvRef} = ($val, 0, 0, -$val*0.5, sqrt(3)*0.5*$val, 0, 0, 0, $val*$coadim);
    } else {
die "Ibrav $ibrav is not recognized, please recreate with ibrav = 0\n";
    }

#    print "Primitive Translation Vectors:\n";
#    print "${$cptvRef}[0] ${$cptvRef}[1] ${$cptvRef}[2]\n";
#    print "${$cptvRef}[3] ${$cptvRef}[4] ${$cptvRef}[5]\n";
#    print "${$cptvRef}[6] ${$cptvRef}[7] ${$cptvRef}[8]\n";


    # now try to read in the atomic species
    my $speciesString = getPwscfCard("ATOMIC_SPECIES", $filedataRef);
     my @allSpecies = split(/\s+/, removeWhiteSpaceAtBeginning($speciesString));
     $#allSpecies+1 == $$numSp * 3 || die "Number of species specified in pwscf file does not match with number of species specified in system section\n";
    for (my $i = 0; $i < $$numSp; $i++) {
$atnmtoppRef->{$allSpecies[$i*3]} = $$psdir . $allSpecies[$i*3+2];
$atNamRef->[$i] = $allSpecies[$i*3];
$ppsRef->[$i] = $$psdir . $allSpecies[$i*3+2];
    }

    # now read in the atomic positions
    my $posString = getPwscfCard("ATOMIC_POSITIONS", $filedataRef);
    my @allPos = split(/\s+/, removeWhiteSpaceAtBeginning($posString));
    $#allPos+1 == $$numAt*4 || die "Number of positions in ATOMIC_POSITIONS does not match the nat in the system section\n";
    
    for (my $i = 0; $i < $$numAt; $i++) {
push(@{$ionsRef}, $allPos[4*$i]);
push(@{$posarrRef}, $allPos[4*$i+1]);
push(@{$posarrRef}, $allPos[4*$i+2]);
push(@{$posarrRef}, $allPos[4*$i+3]);
    }
}
#################################################################################################################


#################################################################################################################
# Subroutine to write a qmc simulation header
#################################################################################################################
sub getQmcpackHeader {
    my $projectId = shift;
    my $seriesId = shift;
    my $commentString = shift;
    my $seed = shift;

    
    my $outputString = '<?xml version="1.0"?>' . "\n";
    $outputString .= '<simulation xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.mcc.uiuc.edu/qmc/schema/molecu.xsd">' . "\n";
    $outputString .= "  <project id=\"$projectId\" series=\"$seriesId\">\n";
    $outputString .= '    <application name="qmcapp" role="molecu" class="serial" version="0.2">' . "\n";
    $outputString .= "      $commentString\n";
    $outputString .= '    </application>' . "\n";
    $outputString .= "  </project>\n";
    $outputString .= "\n  <random seed=\"$seed\"/>\n";
}
#################################################################################################################

#################################################################################################################
# Subroutine to write a qmc simulation footer
#################################################################################################################
sub getQmcpackFooter {
    my $outputString = "</simulation>\n";
}
#################################################################################################################


#################################################################################################################
# Subroutine to look through the optimization directory and figure out which wavefunction
# had the lowest energy
#################################################################################################################
sub findLowestEnergyWvfcn {
    my $optDir = shift;
    my $optTemplateFile = shift;

    my $lowesten=100000000000000000000000000.0;
    my $lowestseqnum = -1;

    if (-d $optDir) {

my $start = 16;

my $energytool = $config{energytool};

$optTemplateFile =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
my $prefix = $1;
my $tnum = $2;
my $suffix = $3;

# get list of files in the optimization directory
opendir DIR, "$optDir"; 
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR; 
closedir DIR; 


my $ofname = $prefix;
chop($ofname);
chop($ofname);
$ofname .= ".out";
$ofname = "$optDir/$ofname";
print "outfile = $ofname\n";
my $rcutline = `grep rcut $ofname | head -n 1`;
my @arr = split(/\s+/, $rcutline);
my $rcutval = $arr[3];
print "rcut = $rcutval\n";

my @rawfiles;
foreach my $str (@files) {
    if ($str =~  /$prefix\d\d\d$suffix/) {
push @rawfiles, $str;
    }
}

# Loop over optimization data files and figure out which one has the lowest average energy
foreach my $file (sort bySequence @rawfiles) {
    $file =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
    $file = $optDir . "/" . $file;
    my $sequence = $2;
    my $str = `$energytool $file $start | head -1`;
    my @data = split(/\s+/,$str);
    if ($data[2] < $lowesten) {
$lowesten = $data[2];
$lowestseqnum = $sequence;
    }
    print "$file:  $data[2]\n";

}

# Get properly formatted sequence number (always 3 digits) for the best .opt.xml file
my $optfileprettyseq = sprintf("%3d", $lowestseqnum-1);
if ($optfileprettyseq < 10) {
    $optfileprettyseq = sprintf("00%d", $optfileprettyseq);
} elsif ($optfileprettyseq < 100) {
    $optfileprettyseq = sprintf("0%2d", $optfileprettyseq);
}

# Now state which wavefunction is best
#my $bestfile = "$optDir/$prefix$optfileprettyseq.opt.xml";
#print "The file with the best wavefunction is: $bestfile\n";
    }
    return ($lowestseqnum, $lowesten);
}


#################################################################################################################
# Subroutine to grab an optimized wavefunction from a directory containing an 
# optimization run.  Will take the wvfcn from the .opt.xml file that corresponds
# to the wavefunction with the lowest energy
#################################################################################################################

sub getWvfcnStringFromOptDir2 {
    my $outputString;
    my $optDir = shift;
    my $optTemplateFile = shift;
    my $seqnum = shift;

    my $standAlone = shift;
    my $isSimple = shift;
#    print "isSimple = $isSimple\n";

    my $atomsNameRef = shift;
    my $totalIons = shift;
    my $ionIdsRef = shift;
    $#{$ionIdsRef}+1 == $totalIons || die "Must give a list of ion names (ionIds) that is the same length as totalIons ($totalIons)\n";

    my $locHdfFileName = shift;
    my $tileMatrixRef = shift;
    my $locTwistNum = shift;
    my $locUseGPU = shift;
#    print "In wavefunction, use gpu = $locUseGPU\n";
    my $locUseMeshFactor = shift;
    my $locMeshFactor = shift;
    my $locSglprecsplines = shift;
    my $locTruncSplines = shift;
    my $locUpElecs = shift;
    my $locDownElecs = shift;
    my $locSpinDependentWvfcn = shift;

    my $locIsAe = shift;
    my $atomsChargeRef = shift;
    my $kpointref = shift;

    $optTemplateFile =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
    my $prefix = $1;
    my $tnum = $2;
    my $suffix = $3;

    # grab rcut from the optimization run's output file
    my $ofname = $prefix;
    chop($ofname);
    chop($ofname);
    $ofname .= ".out";
    $ofname = "$optDir/$ofname";


    ## THIS IS ASSUMING THAT ALL VALUES OF RCUT ARE THE SAME!!!
    ## GETTING THIS FROM THE OUTFILE IS GOOD BUT SHOULD BE ADJUSTED SO THAT IT WORKS FOR EACH JASTROW INDIVIDUALLY
    my $rcutlinenos = `sed -n '/rcut/=' $ofname`;
    my @rcutlinenoarrs = split (/\n/, $rcutlinenos);
    my %rcuthash;
    foreach my $linenum (@rcutlinenoarrs) {
my $rcutval = `sed -n ${linenum}p $ofname | awk '{print \$3}'`;
chomp($rcutval);
my $nextline = $linenum+2;
my $token = `sed -n ${nextline}p $ofname | awk '{print \$1}'`;
#print "token before substr = $token\n";
my $undindex = rindex($token, '_');
#print "undindex = $undindex\n";
$token = substr($token,0,$undindex);
#print "token = $token\n";
$rcuthash{$token} = $rcutval;
    }
    
    my $lowestseqnum = $seqnum;
    my $lowesten;
    unless($seqnum) {
($lowestseqnum, $lowesten) = findLowestEnergyWvfcn($optDir, $optTemplateFile);
    }

    # Get properly formatted sequence number (always 3 digits) for the best .opt.xml file
    my $optfileprettyseq = sprintf("%3d", $lowestseqnum-1);
    if ($optfileprettyseq < 10) {
$optfileprettyseq = sprintf("00%d", $optfileprettyseq);
    } elsif ($optfileprettyseq < 100) {
$optfileprettyseq = sprintf("0%2d", $optfileprettyseq);
    }

    # Now state which wavefunction is best
    my $bestfile = "$optDir/$prefix$optfileprettyseq.opt.xml";
    ##print "The file with the best wavefunction is: $bestfile\n";

    my $detsetline = `grep \"<determinantset\" $bestfile`;
    unless($locMeshFactor) {
$detsetline =~ m/meshfactor=\"(.*?)\"/;
$locMeshFactor = $1;
    }
    
    $detsetline =~ m/tilematrix=\"(.*?)\"/;
    my $tilemat = $1;
    $detsetline =~ m/precision=\"(.*?)\"/;
    my $precision = $1;
    
     
    my $coeflines = `sed -n '/coefficients/=' $bestfile`;
    my @coeflinearr = split (/\n/, $coeflines);
    my %coefhash;
    my %numcoefhash;
    foreach my $linenum (@coeflinearr) {
my $coefline = `sed -n ${linenum}p $bestfile`;
chomp($coefline);
#print "coefline = $coefline\n";
$coefline =~ m/id=\"(.*?)\"/;
my $id = $1;
$coefline =~ m/\>(.*?)\</;
my $coefs = $1;
$linenum--;
my $corrline = `sed -n ${linenum}p $bestfile`;
$corrline =~ m/size=\"(.*?)\"/;
my $coefsize = $1;
$coefhash{$id} = $coefs;
$numcoefhash{$id} = $coefsize;
    }

    if ($standAlone) {
$outputString .= '<?xml version="1.0"?>';
$outputString .= "\n<qmcsystem>\n";
    }
    $outputString .= "  <wavefunction name=\"psi0\" target=\"e\">\n";
    
    if ($tilemat) {
$outputString .= "    <determinantset type=\"einspline\" source=\"i\" href=\"$locHdfFileName\" tilematrix=\"$tilemat\"";
    } else {
$outputString .= "    <determinantset type=\"einspline\" source=\"i\" href=\"$locHdfFileName\" tilematrix=\"@{$tileMatrixRef}\"";
    }
    $outputString .= " twistnum=\"$locTwistNum\"";
    
    if ($locUseGPU) { $outputString .= " gpu=\"yes\""; }
    if ($locUseMeshFactor) { $outputString .= " meshfactor=\"$locMeshFactor\""; }
    if ($locSglprecsplines) { $outputString .= " precision=\"single\""; }
    if ($locTruncSplines) { $outputString .= " truncate=\"yes\""; }
    my $twistx = -1.0* $$kpointref[0];
    my $twisty = -1.0* $$kpointref[1];
    my $twistz = -1.0* $$kpointref[2];
    if ($locTwistNum == -1) { $outputString .= " twist=\"$twistx  $twisty  $twistz\""; }
    $outputString .= ">\n";
    $outputString .= "      <basisset/>\n";
    $outputString .= "      <slaterdeterminant>\n";
    $outputString .= "        <determinant id=\"updet\" size=\"$locUpElecs\" ref=\"updet\">\n";
    $outputString .= "          <occupation mode=\"ground\" spindataset=\"0\">\n";
    $outputString .= "          </occupation>\n";
    $outputString .= "        </determinant>\n";
    $outputString .= "        <determinant id=\"downdet\" size=\"$locDownElecs\" ref=\"downdet\">\n";
    if ($locSpinDependentWvfcn) {
        $outputString .= "          <occupation mode=\"ground\" spindataset=\"1\">\n";
    } else {
        $outputString .= "          <occupation mode=\"ground\" spindataset=\"0\">\n";
    }
    $outputString .= "          </occupation>\n";
    $outputString .= "        </determinant>\n";
    $outputString .= "      </slaterdeterminant>\n";
    $outputString .= "    </determinantset>\n";

    $outputString .= "    <jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\">\n";
    if ($rcuthash{"uu"} > 0) {
my $coefsize = $numcoefhash{"uu"};
my $rcutsize = $rcuthash{"uu"};
$outputString .= "      <correlation speciesA=\"u\" speciesB=\"u\" size=\"$coefsize\" rcut=\"$rcutsize\">\n";
    } else {
my $coefsize = $numcoefhash{"uu"};
$outputString .= "      <correlation speciesA=\"u\" speciesB=\"u\" size=\"$coefsize\">\n";
    }
    my $curcoefs = $coefhash{"uu"};
    $outputString .= "        <coefficients id=\"uu\" type=\"Array\"> $curcoefs </coefficients>\n";
    $outputString .= "      </correlation>\n";
    if ($rcuthash{ud} > 0) {
my $coefsize = $numcoefhash{"ud"};
my $rcutsize = $rcuthash{"ud"};
$outputString .= "      <correlation speciesA=\"u\" speciesB=\"d\" size=\"$coefsize\" rcut=\"$rcutsize\">\n";
    } else {
my $coefsize = $numcoefhash{"ud"};
$outputString .= "      <correlation speciesA=\"u\" speciesB=\"d\" size=\"$coefsize\">\n";
    }
    my $curcoefs = $coefhash{"ud"};
    $outputString .= "        <coefficients id=\"ud\" type=\"Array\"> $curcoefs </coefficients>\n";
    $outputString .= "      </correlation>\n";
    $outputString .= "    </jastrow>\n";

    if ($isSimple) {
# In this case we have only one jastrow for each species
$outputString .= "    <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"i\">\n";
for (my $i = 0; $i <= $#{$atomsNameRef}; $i++) {
    my $id = ${$atomsNameRef}[$i];
    if ($rcuthash{$id} > 0) {
$outputString .= "      <correlation elementType=\"$id\" cusp=\"0.0\" size=\"$numcoefhash{$id}\" rcut=\"$rcuthash{$id}\">\n";
    } else {
$outputString .= "      <correlation elementType=\"$id\" cusp=\"0.0\" size=\"$numcoefhash{$id}\">\n";
    }
    $outputString .= "        <coefficients id=\"${$atomsNameRef}[$i]\" type=\"Array\"> $coefhash{$id} </coefficients>\n";
    $outputString .= "      </correlation>\n";
}
$outputString .= "    </jastrow>\n";
if ($locIsAe) {
    $outputString .= "    <jastrow name=\"J1S\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"i\">\n";
    for (my $j = 0; $j <= $#{$atomsNameRef}; $j++) {
my $id = "${$atomsNameRef}[$j]-sr";
if ($numcoefhash{$id}) {
    $outputString .= "      <correlation elementType=\"${$atomsNameRef}[$j]\" cusp=\"${$atomsChargeRef}[$j]\" size=\"$numcoefhash{$id}\" rcut=\"$rcuthash{$id}\">\n";
    $outputString .= "         <coefficients id=\"$id\" type=\"Array\"> $coefhash{$id} </coefficients>\n";
    $outputString .= "      </correlation>\n";
} else {
    $outputString .= "      <correlation elementType=\"${$atomsNameRef}[$j]\" cusp=\"${$atomsChargeRef}[$j]\" size=\"3\" rcut=\"0.5\">\n";
    $outputString .= "         <coefficients id=\"${$atomsNameRef}[$j]-sr\" type=\"Array\"> -0.1 -0.05 -0.01 </coefficients>\n";
    $outputString .= "      </correlation>\n";
}
    }
    $outputString .= "    </jastrow>\n";
}
    } else {
# In this case we have an independent jastrow for every ion
$outputString .= "    <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"centers\">\n";
for (my $i = 0; $i <= $#{$ionIdsRef}; $i++) {
    my $id = "center$i";
    if (!$coefhash{$id}) {
$id = ${$ionIdsRef}[$i];
    }
#    print "id = $id\n";
    if ($id) {
if ($rcuthash{$id} > 0) {
    $outputString .= "      <correlation elementType=\"center$i\" cusp=\"0.0\" size=\"$numcoefhash{$id}\" rcut=\"$rcuthash{$id}\">\n";
} else {
    $outputString .= "      <correlation elementType=\"center$i\" cusp=\"0.0\" size=\"$numcoefhash{$id}\">\n";
}
$outputString .= "        <coefficients id=\"center$i\" type=\"Array\"> $coefhash{$id} </coefficients>\n ";
$outputString .= "      </correlation>\n";

    } else {
$outputString .= "      <correlation elementType=\"center$i\" cusp=\"0.0\" size=\"4\">\n";
$outputString .= "        <coefficients id=\"center$i\" type=\"Array\"> 0 0 0 0 </coefficients>\n";
$outputString .= "      </correlation>\n";
    }
}
$outputString .= "    </jastrow>\n";
    }

    $outputString .= "  </wavefunction>\n";
    if ($standAlone) {
$outputString .= "</qmcsystem>\n";
    }
    $outputString;
}

#################################################################################################################
# Subroutine to grab an optimized wavefunction from a directory containing an 
# optimization run.  Will take the wvfcn from the .opt.xml file that corresponds
# to the wavefunction with the lowest energy
#################################################################################################################
sub getWvfcnStringFromOptDir {
    my $optDir = shift;
    my $optTemplateFile = shift;
    my $twistNum = shift;
    my $splFac = shift;
    my $uGPU = shift;
    my $wfile = shift;
    my $seqnum = shift;
    my $tileMatrixRef = shift;

    $optTemplateFile =~ /(.*\.s)(\d\d\d)(\.scalar\.dat)/;
    my $prefix = $1;
    my $tnum = $2;
    my $suffix = $3;

    # grab rcut from the optimization run's output file
    my $ofname = $prefix;
    chop($ofname);
    chop($ofname);
    $ofname .= ".out";
    $ofname = "$optDir/$ofname";
    my $rcutline = `grep rcut $ofname | head -n 1`;
    my @arr = split(/\s+/, $rcutline);
    my $rcutval = $arr[3];
    print "rcut = $rcutval\n";

    my $lowestseqnum = $seqnum;
    unless($seqnum) {
$lowestseqnum = findLowestEnergyWvfcn($optDir, $optTemplateFile);
    }

    # Get properly formatted sequence number (always 3 digits) for the best .opt.xml file
    my $optfileprettyseq = sprintf("%3d", $lowestseqnum-1);
    if ($optfileprettyseq < 10) {
$optfileprettyseq = sprintf("00%d", $optfileprettyseq);
    } elsif ($optfileprettyseq < 100) {
$optfileprettyseq = sprintf("0%2d", $optfileprettyseq);
    }

    # Now state which wavefunction is best
    my $bestfile = "$optDir/$prefix$optfileprettyseq.opt.xml";
    ##print "The file with the best wavefunction is: $bestfile\n";


    # open the best file and start parsing it
    my $wvfcnString;
    open(BESTFILE, "$bestfile") || die "Cannot open file $bestfile\n";
    my @bfdata = <BESTFILE>;
    close(BESTFILE);

    my $start = 0;
    my $stop = 0;

### tilematrix=\"@{$tileMatrixRef}\"


    foreach my $line (@bfdata) {
if ($line =~ /\/wavefunction/) {
    $stop = 1;
    $wvfcnString .= $line;
}
if ($line =~ /wavefunction/) {
    $start = 1;
}

if ($start && !$stop) {
    if ($line =~ /<correlation/) {
unless ($line =~ /rcut/) {
    $line =~ s/<correlation/<correlation rcut=\"$rcutval\"/;
}
    }
    if ($line =~ /<determinantset/) {
if ($line =~ /href/) {
    $line =~ s/href\s*=\s*\".*?\"/href=\"..\/$wfile\"/;
} else {
    $line =~ s/>/ href=\"..\/$wfile\">/;
}
    }
    if ($line =~ /<determinantset/) {
if ($line =~ /twistnum/) {
    $line =~ s/twistnum\s*=\s*\"\d+\"/twistnum=\"$twistNum\"/;
} else {
    $line =~ s/>/ twistnum=\"$twistNum\">/;
}
    }
    if ($line =~ /<determinantset/) {
if ($line =~ /tilematrix/) {
    $line =~ s/tilematrix\s*=\s*\"\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*-?\d\s*\"/tilematrix=\"@{$tileMatrixRef}\"/;
} else {
    $line =~ s/>/ tilematrix=\"@{$tileMatrixRef}\">/;
}
    }
    if ($splFac) {
if ($line =~ /<determinantset/) {
    if ($line =~ /meshfactor/) {
$line =~ s/meshfactor\s*=\s*\"(.*?)\"/meshfactor=\"$splFac\"/;
    } else {
$line =~ s/>/ meshfactor=\"$splFac\">/;
    }
}   
    }
    if ($useGPU) {
if ($line =~ /<determinantset/) {
    if ($line =~ /gpu/) {
$line =~ s/gpu\s*=\s*\"(.*?)\"/gpu=\"yes\"/;
    } else {
$line =~ s/>/ gpu=\"yes\">/;
    }
}
    }
    
    $wvfcnString .= $line;
}
    }
    $wvfcnString;



}
#################################################################################################################

#################################################################################################################
# Subroutine to write a qmcpack wavefunction section
#################################################################################################################
sub getWavefunctionString {
    my $outputString;

    my $standAlone = shift;
    my $hasJastrow = shift;
    my $isSimple = shift;
#    print "isSimple = $isSimple\n";

    my $sizeTwoBody = shift;
    my $rcutTwoBody = shift;
    my $upUpCoefsRef = shift;
    my $upDownCoefsRef = shift;

    my $atomsNameRef = shift;
    my $totalIons = shift;
    my $ionIdsRef = shift;
    $#{$ionIdsRef}+1 == $totalIons || die "Must give a list of ion names (ionIds) that is the same length as totalIons ($totalIons)\n";

    my $sizeOneBody = shift;
    my $rcutOneBody = shift;
    my $jastrowStartsRef = shift;
    if ($hasJastrow) {
$#{$jastrowStartsRef}+1 == $sizeOneBody*($#{$atomsNameRef}+1) || die
    "jastrowStarts must contain inital values for all types of one body jastrow.  As such it has\n to be $sizeOneBody (sizeOneBody) * ($#{$atomsNameRef} + 1) (Number of species in atomsName array) elements long\n";
    }

    # Now need to create a data structure that maps from the values in ionIds to the
    # index in jastrowStarts where we will grab the starting coefficients for the wfns file
    my %atomIdToIndexHash;
    if ($hasJastrow) {
for (my $i = 0; $i <= $#{$atomsNameRef}; $i++) {
    $atomIdToIndexHash{${$atomsNameRef}[$i]} = $i*$sizeOneBody;
}
    }
    # Now the starting index in jastrowStarts for each element in ionIds should be:
    # $atomIdToIndexHash{$atomsNameRef[$i]}

    my $locHdfFileName = shift;
    my $tileMatrixRef = shift;
    my $locTwistNum = shift;
    my $locUseGPU = shift;
#    print "In wavefunction, use gpu = $locUseGPU\n";
    my $locUseMeshFactor = shift;
    my $locMeshFactor = shift;
    my $locSglprecsplines = shift;
    my $locTruncSplines = shift;
    my $locUpElecs = shift;
    my $locDownElecs = shift;
    my $locSpinDependentWvfcn = shift;

    my $locIsAe = shift;
    my $atomsChargeRef = shift;

    if ($standAlone) {
$outputString .= '<?xml version="1.0"?>';
$outputString .= "\n<qmcsystem>\n";
    }
    $outputString .= "  <wavefunction name=\"psi0\" target=\"e\">\n";

    $outputString .= "    <determinantset type=\"einspline\" href=\"$locHdfFileName\" source=\"i\" tilematrix=\"@{$tileMatrixRef}\" twistnum=\"$locTwistNum\"";
    if ($locUseGPU) { $outputString .= " gpu=\"yes\""; }
    if ($locUseMeshFactor) { $outputString .= " meshfactor=\"$locMeshFactor\""; }
    if ($locSglprecsplines) { $outputString .= " precision=\"single\""; }
    if ($locTruncSplines) { $outputString .= " truncate=\"yes\""; }
    $outputString .= ">\n";
    $outputString .= "      <basisset/>\n";
    $outputString .= "      <slaterdeterminant>\n";
    $outputString .= "        <determinant id=\"updet\" size=\"$locUpElecs\" ref=\"updet\">\n";
    $outputString .= "          <occupation mode=\"ground\" spindataset=\"0\">\n";
    $outputString .= "          </occupation>\n";
    $outputString .= "        </determinant>\n";
    $outputString .= "        <determinant id=\"downdet\" size=\"$locDownElecs\" ref=\"downdet\">\n";
    if ($locSpinDependentWvfcn) {
        $outputString .= "          <occupation mode=\"ground\" spindataset=\"1\">\n";
    } else {
        $outputString .= "          <occupation mode=\"ground\" spindataset=\"0\">\n";
    }
    $outputString .= "          </occupation>\n";
    $outputString .= "        </determinant>\n";
    $outputString .= "      </slaterdeterminant>\n";
    $outputString .= "    </determinantset>\n";

    if ($hasJastrow) {
$outputString .= "    <jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\">\n";
if ($rcutTwoBody > 0) {
    $outputString .= "      <correlation speciesA=\"u\" speciesB=\"u\" size=\"$sizeTwoBody\" rcut=\"$rcutTwoBody\">\n";
} else {
    $outputString .= "      <correlation speciesA=\"u\" speciesB=\"u\" size=\"$sizeTwoBody\">\n";
}
$outputString .= "        <coefficients id=\"uu\" type=\"Array\">  ";
$#{$upUpCoefsRef} == $sizeTwoBody-1 || die "In getWavefunctionString, sizeTwoBody and dimensions of upUpCoefs do not agree\n";
for (my $i = 0; $i < $#{$upUpCoefsRef}+1; $i++) {
    $outputString .= "${$upUpCoefsRef}[$i]  ";
}
$outputString .= "</coefficients>\n";
$outputString .= "      </correlation>\n";
if ($rcutTwoBody > 0) {
    $outputString .= "      <correlation speciesA=\"u\" speciesB=\"d\" size=\"$sizeTwoBody\" rcut=\"$rcutTwoBody\">\n";
} else {
    $outputString .= "      <correlation speciesA=\"u\" speciesB=\"d\" size=\"$sizeTwoBody\">\n";
}
$outputString .= "        <coefficients id=\"ud\" type=\"Array\">  ";
$#{$upDownCoefsRef} == $sizeTwoBody-1 || die "In getWavefunctionString, sizeTwoBody and dimensions of upDownCoefs do not agree\n";
for (my $i = 0; $i < $#{$upDownCoefsRef}+1; $i++) {
    $outputString .= "${$upDownCoefsRef}[$i]  ";
}
$outputString .= "</coefficients>\n";
$outputString .= "      </correlation>\n";
$outputString .= "    </jastrow>\n";
    }

    if ($hasJastrow) {
if ($isSimple) {
    # In this case we have only one jastrow for each species
    $outputString .= "    <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"i\">\n";
    for (my $i = 0; $i <= $#{$atomsNameRef}; $i++) {
if ($rcutOneBody > 0) {
    $outputString .= "      <correlation elementType=\"${$atomsNameRef}[$i]\" cusp=\"0.0\" size=\"$sizeOneBody\" rcut=\"$rcutOneBody\">\n";
} else {
    $outputString .= "      <correlation elementType=\"${$atomsNameRef}[$i]\" cusp=\"0.0\" size=\"$sizeOneBody\">\n";
}
$outputString .= "        <coefficients id=\"${$atomsNameRef}[$i]\" type=\"Array\"> ";
for (my $j = 0; $j < $sizeOneBody; $j++) {
    $outputString .= ${$jastrowStartsRef}[$j+$atomIdToIndexHash{${$atomsNameRef}[$i]}] . "  ";
}
$outputString .= "</coefficients>\n";
$outputString .= "      </correlation>\n";
    }
    $outputString .= "    </jastrow>\n";
    if ($locIsAe) {
$outputString .= "    <jastrow name=\"J1S\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"i\">\n";
for (my $j = 0; $j <= $#{$atomsNameRef}; $j++) {
    $outputString .= "      <correlation elementType=\"${$atomsNameRef}[$j]\" cusp=\"${$atomsChargeRef}[$j]\" size=\"$sizeOneBody\" rcut=\"0.5\">\n";
    $outputString .= "         <coefficients id=\"${$atomsNameRef}[$j]-sr\" type=\"Array\"> -0.1 -0.05 -0.01 </coefficients>\n";
    $outputString .= "      </correlation>\n";
}
$outputString .= "    </jastrow>\n";
    }
} else {
    # In this case we have an independent jastrow for every ion
    $outputString .= "    <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"centers\">\n";
    for (my $i = 0; $i <= $#{$ionIdsRef}; $i++) {
if ($rcutOneBody > 0) {
    $outputString .= "      <correlation elementType=\"center$i\" cusp=\"0.0\" size=\"$sizeOneBody\" rcut=\"$rcutOneBody\">\n";
} else {
    $outputString .= "      <correlation elementType=\"center$i\" cusp=\"0.0\" size=\"$sizeOneBody\">\n";
}
$outputString .= "        <coefficients id=\"center$i\" type=\"Array\"> ";
for (my $j = 0; $j < $sizeOneBody; $j++) {
    $outputString .= ${$jastrowStartsRef}[$j+$atomIdToIndexHash{${$ionIdsRef}[$i]}] . "  ";
}
$outputString .= "</coefficients>\n";
$outputString .= "      </correlation>\n";
    }
    $outputString .= "    </jastrow>\n";
}
    }


    $outputString .= "  </wavefunction>\n";
    if ($standAlone) {
$outputString .= "</qmcsystem>\n";
    }
    $outputString;
}
#################################################################################################################


#################################################################################################################
# Subroutine to write a qmcpack particleset
# first argument is a flag which is 1 if this is a simple particleset.
# if it is 0 an auxiliary particleset is created called centers which can
# be useful for attaching independent jastrows to each ion
#################################################################################################################
sub getPtclSetString {
    my $outputString;

    my $standAlone = shift;
    my $isSimple = shift;
    my $lattRef = shift;

    my $numUpElectrons = shift;
    my $numDownElectrons = shift;

    my $numIons = shift;
    my $numIonSpecies = shift;
    my $atomsNameRef = shift;
    my $atomsChargeRef = shift;
    my $atomsValenceRef = shift;
    my $atomsAtomicNumbersRef = shift;

    my $ionPosRef = shift;
    $#{$ionPosRef}+1 == $numIons*3 || die "Ion Positions array handed to writePtclset supposed to have\n 3 * $numIons (3 * \$numIons) entries, but it has $#{$ionPosRef}+1\n";
    my $ionIdsRef = shift;
    $#{$ionIdsRef}+1 == $numIons || die "Ion ids array handed to writePtclset supposed to have \n$numIons entries, but it has $#{$ionIdsRef}+1\n";

    my $bcondsRef = shift;    

    if ($standAlone) {
$outputString .= "<?xml version=\"1.0\"?>\n";
    }

    my $bcondsString;
    for (my $i = 0; $i < 3; $i++) {
      if ($$bcondsRef[$i] == 1) {
$bcondsString .= "p ";
      } else {
$bcondsString .= "n ";
      }
    }

    $outputString .= "  <qmcsystem>\n";    
    $outputString .= "  <simulationcell>\n";
    $outputString .= "    <parameter name=\"lattice\">\n";
    $outputString .= "      ${$lattRef}[0]   ${$lattRef}[1]   ${$lattRef}[2]\n";
    $outputString .= "      ${$lattRef}[3]   ${$lattRef}[4]   ${$lattRef}[5]\n";
    $outputString .= "      ${$lattRef}[6]   ${$lattRef}[7]   ${$lattRef}[8]\n";
    $outputString .= "      </parameter>\n";
    $outputString .= "      <parameter name=\"bconds\">$bcondsString</parameter>\n";  
    if ($$bcondsRef[0] == 1 || $$bcondsRef[1] == 1 || $$bcondsRef[2] == 1) {
      $outputString .= "      <parameter name=\"LR_dim_cutoff\">15</parameter>\n";
    }
    $outputString .= "  </simulationcell>\n";
    $outputString .= "  </qmcsystem>\n";

    $outputString .= "  <particleset name=\"i\" size=\"$numIons\">\n";
    for (my $i = 0; $i <= $numIonSpecies; $i++) {
        $outputString .= "    <group name=\"${$atomsNameRef}[$i]\">\n";
        $outputString .= "      <parameter name=\"charge\">${$atomsChargeRef}[$i]</parameter>\n";
        $outputString .= "      <parameter name=\"valence\">${$atomsValenceRef}[$i]</parameter>\n";
        $outputString .= "      <parameter name=\"atomicnumber\">${$atomsAtomicNumbersRef}[$i]</parameter>\n";
        $outputString .= "    </group>\n";
    }
    $outputString .= "    <attrib name=\"position\" datatype=\"posArray\" condition=\"1\">\n";
    for (my $i = 0; $i < $numIons; $i++) {
        my @locarr = (${$ionPosRef}[3*$i], ${$ionPosRef}[3*$i+1], ${$ionPosRef}[3*$i+2]);
        $outputString .= "    $locarr[0]   $locarr[1]   $locarr[2]\n";
    }
    $outputString .= "    </attrib>\n";

    
    $outputString .= "    <attrib name=\"ionid\" datatype=\"stringArray\">\n";
    
    for (my $i = 0; $i*8+8 <= $numIons; $i++) {
        $outputString .= "      ";
        for (my $j = 0; $j < 8; $j++) {
            my $id = $i*8+$j;
            $outputString .= "${$ionIdsRef}[$id]  ";
        }
        $outputString .= "\n";
    }    
    if ($numIons % 8 > 0) {
        $outputString .= "      ";
        for (my $i = 0; $i < ($numIons % 8); $i++) {
            my $id = 8*floor($numIons/8)+$i;
            $outputString .= "${$ionIdsRef}[$id]  ";
        }
        $outputString .= "\n";
    }
    $outputString .= "    </attrib>\n";
    $outputString .= "  </particleset>\n";

    $outputString .= "  <particleset name=\"e\" random=\"yes\" randomsrc=\"i\">\n";
    $outputString .= "    <group name=\"u\" size=\"$numUpElectrons\">\n";
    $outputString .= "      <parameter name=\"charge\">-1</parameter>\n";
    $outputString .= "    </group>\n";
    $outputString .= "    <group name=\"d\" size=\"$numDownElectrons\">\n";
    $outputString .= "      <parameter name=\"charge\">-1</parameter>\n";
    $outputString .= "    </group>\n";
    $outputString .= "  </particleset>\n";
    if (!$isSimple) {
        $outputString .= "  <particleset name=\"centers\" size=\"$numIons\">\n";
        for (my $i = 0; $i < $numIons; $i++) {
            $outputString .= "    <group name=\"center$i\"></group>\n";
        }
        $outputString .= "    <attrib name=\"position\" datatype=\"posArray\" condition=\"1\">\n";
        for (my $i = 0; $i < $numIons; $i++) {
            my @locarr = (${$ionPosRef}[3*$i], ${$ionPosRef}[3*$i+1], ${$ionPosRef}[3*$i+2]);
            $outputString .= "    $locarr[0]   $locarr[1]   $locarr[2]\n";
        }
        $outputString .= "    </attrib>\n";
        $outputString .= "    <attrib name=\"ionid\" datatype=\"stringArray\">\n";
        for (my $i = 0; $i*8+8 <= $numIons; $i++) {
            $outputString .= "      ";
            for (my $j = 0; $j < 8; $j++) {
                my $id = $i*8+$j;
                $outputString .= "center$id  ";
            }
            $outputString .= "\n";
        }
        if ($numIons % 8 > 0) {
            $outputString .= "      ";
            for (my $i = 0; $i < ($numIons % 8); $i++) {
                my $id = 8*floor($numIons/8)+$i;
                $outputString .= "center$id  ";
            }
            $outputString .= "\n";
        }


        $outputString .= "    </attrib>\n";
        $outputString .= "  </particleset>\n\n";
    }

    if ($standAlone) {
    }
    $outputString;
}
#################################################################################################################




#################################################################################################################
# Subroutine to write qmcpack Hamiltonian section
#################################################################################################################
sub getHamiltonianString {
    my $standAlone = shift;
    my $atNameToPPRef = shift;
    my $locUseMPC = shift;
    my $locMPCIsPhysical = shift;
    my $locUseKECorr = shift;
    my $locIsAe = shift;
    my $atomsChargeRef = shift;
    my $ppName;
    my $ppFile;

    my $outputString;

    if ($standAlone) {
$outputString .= '<?xml version="1.0"?>' . "\n";
$outputString .= '<casing> <!-- for whatever reason all included sections must be in some kind of arbitrary tag -->' . "\n";
    }
    $outputString .= '  <hamiltonian name="h0" type="generic" target="e">' . "\n";
    $outputString .= '    <pairpot type="pseudo" name="PseudoPot" source="i" wavefunction="psi0" format="xml">' ."\n";
    if (!$locIsAe) {
while (($ppName, $ppFile) = each(%{$atNameToPPRef})) {
    $outputString .= "      <pseudo elementType=\"$ppName\" href=\"$ppFile\"/>\n";
}
    } else {
my $i = 0;
while (($ppName, $ppFile) = each(%{$atNameToPPRef})) {
#    my $chg = ${$atomsChargeRef}[$i];
    my $chg = AtNameToNumber($ppName);
    $i++;
    $outputString .= "      <pseudo elementType=\"$ppName\">\n";
            $outputString .= "        <header symbol=\"$ppName\" atomic-number=\"$chg\" zval=\"$chg\" />\n";
            $outputString .= "        <local>\n";
            $outputString .= "          <grid type=\"linear\" ri=\"0.0\" rf=\"4.0\" npts=\"201\" />\n";
            $outputString .= "        </local>\n";
            $outputString .= "      </pseudo>\n";
}
    }
    $outputString .= '    </pairpot>' ."\n";
    $outputString .= '    <constant name="IonIon" type="coulomb" source="i" target="i"/>' . "\n";
    if ($locUseMPC) {
if ($locMPCIsPhysical) {
    $outputString .= '    <pairpot name="MPC" type="MPC" source="e" target="e" ecut="60.0" physical="true"/>' . "\n";
    $outputString .= '    <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="false"/>' . "\n";
} else {
    $outputString .= '    <pairpot name="MPC" type="MPC" source="e" target="e" ecut="60.0" physical="false"/>' . "\n";
    $outputString .= '    <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="true"/>' . "\n";
}
    } else {
$outputString .= '    <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="true"/>' . "\n";
    }
    if ($locUseKECorr) {
$outputString .= '    <estimator name="KEcorr" type="chiesa" source="e" psi="psi0"/>' . "\n";
    }
    $outputString .= '  </hamiltonian>' . "\n";
    if ($standAlone) {
$outputString .= '</casing>' . "\n";
    }
    $outputString;
}
#################################################################################################################


#################################################################################################################
# Subroutine to write a vmc section
#################################################################################################################
sub getVMCSection {
    my $useGPU_ = shift;
    my $vmctimestep_ = shift;
    my $vmcequiltime_ = shift;
    my $vmcdecorrtime_ = shift;
    my $vmcblocks_ = shift;
    my $vmcwalkers_ = shift;
    my $targetpop_ = shift;
    my $useDrift_ = shift;
    my $samplesPerThread = shift;
    
    my $warmupSteps_ = floor($vmcequiltime_/$vmctimestep_);
    my $StBS_ = floor($vmcdecorrtime_/$vmctimestep_);

    my $outputString;
    if ($useGPU_) {
$outputString .= "  <qmc method=\"vmc\" move=\"pbyp\" gpu=\"yes\">\n";
    } else {
$outputString .= "  <qmc method=\"vmc\" move=\"pbyp\">\n";
    }
    $outputString .= "    <estimator name=\"LocalEnergy\" hdf5=\"no\"/>\n";
    $outputString .= "    <parameter name=\"walkers\">    $vmcwalkers_ </parameter>\n";
    if ($samplesPerThread) {
      $outputString .= "    <parameter name=\"samplesperthread\">    $targetpop_ </parameter>\n";
    } else {
      $outputString .= "    <parameter name=\"samples\">    $targetpop_ </parameter>\n";    
    }
    $outputString .= "    <parameter name=\"stepsbetweensamples\">    1 </parameter>\n";
    $outputString .= "    <parameter name=\"substeps\">  $StBS_ </parameter>\n";
    $outputString .= "    <parameter name=\"warmupSteps\">  $warmupSteps_ </parameter>\n";
    $outputString .= "    <parameter name=\"blocks\">  $vmcblocks_ </parameter>\n";
    $outputString .= "    <parameter name=\"timestep\">  $vmctimestep_ </parameter>\n";
    if ($useDrift_) {
$outputString .= "    <parameter name=\"usedrift\">  yes </parameter>\n";
    } else {
    $outputString .= "    <parameter name=\"usedrift\">   no </parameter>\n";
    }
    $outputString .= "  </qmc>\n";
}

#################################################################################################################
# New Subroutine to write a dmc section
#################################################################################################################
sub getDMCSectionNew {
    my $useGPU_ = shift;
    my $dmctstep_ = shift;
    my $dmcequiltime_ = shift;
    my $dmcruntime_ = shift;
    my $dmcblocktime_ = shift;
    my $dmcUseTmoves_ = shift;
    my $targetPop_ = shift;
    my $minimumTargetPop_ = shift;
    my $warmupFlag_ = shift;
    

    my $dmcWarmupSteps_ = $dmcequiltime_ / $dmctstep_;
    my $totDMCSteps_ = ($dmcruntime+$dmcequiltime_) / $dmctstep_;
    my $dmcStepsPerBlock_ = floor($dmcblocktime_ / $dmctstep_);
    my $dmcBlocks_ = floor($totDMCSteps_ / $dmcStepsPerBlock_);


    my $outputString;
    if ($useGPU_) {
$outputString .= "  <qmc method=\"dmc\" move=\"pbyp\" checkpoint=\"-1\" gpu=\"yes\">\n";
    } else {
$outputString .= "  <qmc method=\"dmc\" move=\"pbyp\" checkpoint=\"-1\">\n";
    }
    $outputString .= "    <estimator name=\"LocalEnergy\" hdf5=\"no\"/>\n";
    if ($warmupFlag_) {
      $outputString .= "    <parameter name=\"reconfiguration\">   yes </parameter>\n";
    }
    if (!$warmupFlag_) {
      if ($targetPop_) {
  $outputString .= "    <parameter name=\"targetwalkers\">  $targetPop_ </parameter>\n";
      } elsif ($minimumTargetPop_) {
  $outputString .= "    <parameter name=\"minimumtargetwalkers\">  $minimumTargetPop_ </parameter>\n";
      }

      $outputString .= "    <parameter name=\"reconfiguration\">   no </parameter>\n";
      $outputString .= "    <parameter name=\"warmupSteps\">  $dmcWarmupSteps_ </parameter>\n";
    }
    $outputString .= "    <parameter name=\"timestep\">  $dmctstep_ </parameter>\n";
    $outputString .= "    <parameter name=\"steps\">   $dmcStepsPerBlock_ </parameter>\n";
    $outputString .= "    <parameter name=\"blocks\">  $dmcBlocks_ </parameter>\n";
    if ($dmcUseTmoves_) {
$outputString .= "    <parameter name=\"nonlocalmoves\">  yes </parameter>\n";
    } else {
    $outputString .= "    <parameter name=\"nonlocalmoves\">   no </parameter>\n";
    }
    $outputString .= "  </qmc>\n";
}


################################################################################################################
# Subroutine to write an optimization section
# Using a routine from Jeremy that worked for Al
# TODO: Add support for pure VMC variance minimization
################################################################################################################
sub getOptSection {
    my $useGPU = shift;
    my $walkers = shift;
    my $vmctimestep = shift;
    my $vmcequiltime = shift;
    my $vmcdecorrtime = shift;
    my $vmcblocks = shift;
    my $numOptSamples = shift;
    my $numOptLoops = shift;
    my $useDrift = shift;

    my $warmupSteps = floor($vmcequiltime/$vmctimestep/10);
    my $sbs = floor($vmcdecorrtime/$vmctimestep);
    print " Warmup Steps = (VMC Equil time / VMC Time Step) / 10\n";
    print " Sub steps = VMC Decorrelation Time / VMC Time Step\n";
    print "$useDrift\n";
    my $outputString;
   
    
    if ($numOptLoops) {
$outputString .= "<loop max=\"$numOptLoops\">\n";
    }
    if ($useGPU) {
$outputString .= "  <qmc method=\"linear\" move=\"pbyp\" checkpoint=\"-1\" gpu=\"yes\">\n";
    } else {
$outputString .= "  <qmc method=\"linear\" move=\"pbyp\" checkpoint=\"-1\" gpu=\"no\">\n";
    }
    $outputString .= "    <parameter name=\"blocks\">   100 </parameter>\n";
    $outputString .= "    <parameter name=\"warmupsteps\"> 10 </parameter>\n";
    $outputString .= "    <parameter name=\"stepsbetweensamples\">  1 </parameter>\n";
    $outputString .= "    <parameter name=\"substeps\">$sbs     </parameter>\n";
    $outputString .= "    <parameter name=\"timestep\">  $vmctimestep  </parameter>\n";
    $outputString .= "    <parameter name=\"walkers\">  $walkers </parameter>\n";
    $outputString .= "    <parameter name=\"samples\">  $numOptSamples  </parameter>\n";
    $outputString .= "    <parameter name=\"minwalkers\">  0.5 </parameter>\n";
    $outputString .= "    <parameter name=\"maxWeight\">    1e9 </parameter>\n";
    if ($useDrift) {
$outputString .= "    <parameter name=\"useDrift\">  yes </parameter>\n";
    } else {
$outputString .= "    <parameter name=\"useDrift\">   no </parameter>\n";
    }
    $outputString .= "    <estimator name=\"LocalEnergy\" hdf5=\"no\"/>\n";
    $outputString .= "    <cost name=\"energy\">                   0.0 </cost>\n";
    $outputString .= "    <cost name=\"unreweightedvariance\">     1.0 </cost>\n";
    $outputString .= "    <cost name=\"reweightedvariance\">       0.0 </cost>\n";
    $outputString .= "    <parameter name=\"MinMethod\">rescale</parameter>\n";
    $outputString .= "    <parameter name=\"GEVMethod\">mixed</parameter>\n";
    $outputString .= "    <parameter name=\"beta\">  0.0 </parameter>\n";
    $outputString .= "    <parameter name=\"exp0\"> -16 </parameter>\n";
    $outputString .= "    <parameter name=\"nonlocalpp\">no</parameter>\n";
    $outputString .= "    <parameter name=\"useBuffer\">no</parameter>\n";
    $outputString .= "    <parameter name=\"bigchange\">9.0</parameter>\n";
    $outputString .= "    <parameter name=\"alloweddifference\"> 1.0e-3 </parameter>\n";
    $outputString .= "    <parameter name=\"stepsize\">4.0e-1</parameter>\n";
    $outputString .= "    <parameter name=\"stabilizerscale\">  1.0 </parameter>\n";
    $outputString .= "    <parameter name=\"nstabilizers\"> 3 </parameter>\n";
    $outputString .= "    <parameter name=\"max_its\"> 1 </parameter>\n";
    $outputString .= "  </qmc>\n";
    if ($numOptLoops) {
$outputString .= "</loop>\n";
    }
    $outputString .= "<loop max=\"6\">\n";
    if ($useGPU) {
$outputString .= "  <qmc method=\"linear\" move=\"pbyp\" checkpoint=\"-1\" gpu=\"yes\">\n";
    } else {
$outputString .= "  <qmc method=\"linear\" move=\"pbyp\" checkpoint=\"-1\" gpu=\"no\">\n";
    }
    $outputString .= "    <parameter name=\"blocks\">   100 </parameter>\n";
    $outputString .= "    <parameter name=\"warmupsteps\"> 10 </parameter>\n";
    $outputString .= "    <parameter name=\"stepsbetweensamples\">   1 </parameter>\n";
    $outputString .= "    <parameter name=\"substeps\">$sbs     </parameter>\n";
    $outputString .= "    <parameter name=\"timestep\">  $vmctimestep  </parameter>\n";
    $outputString .= "    <parameter name=\"walkers\">  $walkers </parameter>\n";
    $outputString .= "    <parameter name=\"samples\">  $numOptSamples  </parameter>\n";
    $outputString .= "    <parameter name=\"minwalkers\">  0.5 </parameter>\n";
    $outputString .= "    <parameter name=\"maxWeight\">    1e9 </parameter>\n";
    if ($useDrift) {
$outputString .= "    <parameter name=\"useDrift\">  yes </parameter>\n";
    } else {
$outputString .= "    <parameter name=\"useDrift\">   no </parameter>\n";
    }
    $outputString .= "    <estimator name=\"LocalEnergy\" hdf5=\"no\"/>\n";
    $outputString .= "    <cost name=\"energy\">                   0.3</cost>\n";
    $outputString .= "    <cost name=\"unreweightedvariance\">     0.0 </cost>\n";
    $outputString .= "    <cost name=\"reweightedvariance\">       0.7</cost>\n";
    $outputString .= "    <parameter name=\"MinMethod\">quartic</parameter>\n";
    $outputString .= "    <parameter name=\"GEVMethod\">mixed</parameter>\n";
    $outputString .= "    <parameter name=\"beta\">  0.0 </parameter>\n";
    $outputString .= "    <parameter name=\"exp0\"> -16 </parameter>\n";
    $outputString .= "    <parameter name=\"nonlocalpp\">no</parameter>\n";
    $outputString .= "    <parameter name=\"useBuffer\">no</parameter>\n";
    $outputString .= "    <parameter name=\"bigchange\">5.0</parameter>\n";
    $outputString .= "    <parameter name=\"alloweddifference\"> 1.0e-4 </parameter>\n";
    $outputString .= "    <parameter name=\"stepsize\">4.0e-1</parameter>\n";
    $outputString .= "    <parameter name=\"stabilizerscale\">  1.0 </parameter>\n";
    $outputString .= "    <parameter name=\"nstabilizers\"> 3 </parameter>\n";
    $outputString .= "    <parameter name=\"max_its\"> 1 </parameter>\n";
    $outputString .= "  </qmc>\n";
    if ($numOptLoops) {
$outputString .= "</loop>\n";
    }
    $outputString;
}


#################################################################################################################
# Subroutine to get the value of a token from a pwscf input file
#################################################################################################################
sub getPwscfToken {
    my $tokenName = shift;
    my $pwscfGetTokenFileDataRef = shift;
   
    my $outval = -999999999999999;
    for (my $i = 0; $i <= $#{$pwscfGetTokenFileDataRef}; $i++) {
my $line = $pwscfGetTokenFileDataRef->[$i];
if (index($line, $tokenName) >= 0) {
    # need to handle two cases here.  One where the token we want is in 
    # quotes and one where it is not
    # first remove all whitespace from the line
     my $cleanedLine = removeQuotes(removeWhiteSpace($line));
    my $startDataCharindex = index($cleanedLine,$tokenName)+length($tokenName)+1; 
    my $dataToParse = substr($cleanedLine, $startDataCharindex);
    my $nextCommaLoc = index($dataToParse, ',');
    
    if ($nextCommaLoc > 0) {
$outval = substr($dataToParse, 0, $nextCommaLoc);
    } else {
$outval = $dataToParse;
    }
}
    }
    $outval;
}
#################################################################################################################


#################################################################################################################
# Subroutine to change the value of a token from a pwscf input file
#################################################################################################################
sub changePwscfToken {
    my $tokenName = shift;
    my $newTokenValue = shift;
    my $pwscfGetTokenFileDataRef = shift;
   
    my $retval = 0;
    for (my $i = 0; $i <= $#{$pwscfGetTokenFileDataRef}; $i++) {
my $line = $pwscfGetTokenFileDataRef->[$i];
if ($line =~ /$tokenName/) {
    $line =~ s/($tokenName(\s*)?)=([^,\n]*)/$tokenName=$newTokenValue/;
    $retval = 1;
    $pwscfGetTokenFileDataRef->[$i] = $line;
}
    }
    $retval;
}

sub removePwscfToken {
    my $tokenName = shift;
    my $pwscfGetTokenFileDataRef = shift;
   
    my $retval = 0;
    for (my $i = 0; $i <= $#{$pwscfGetTokenFileDataRef}; $i++) {
my $line = $pwscfGetTokenFileDataRef->[$i];
if ($line =~ /$tokenName/) {
    $line =~ s/($tokenName(\s*)?)=([^,]*,)//;
    $retval = 1;
    $pwscfGetTokenFileDataRef->[$i] = $line;
}
    }
    $retval;
}

#################################################################################################################


#################################################################################################################
# Helper routine to remove all of the whitespace at the beginning of a string
#################################################################################################################
sub removeWhiteSpaceAtBeginning {
    my $inputString = shift;
    my $outputString;
    my $strlen = length($inputString);
    my $copy = 0;
    for (my $i = 0; $i <= $strlen; $i++) {
my $curchar = substr($inputString, $i, 1);
if ($curchar =~ /\s/) {
} else {
    $copy = 1;
}
if ($copy) {
    $outputString .= $curchar;
}
    }
    $outputString;
}
#################################################################################################################


#################################################################################################################
# Subroutine to remove all whitespace from a string
#################################################################################################################
sub removeWhiteSpace {
    my $inputString = shift;
    my $outputString;

    my $strlen = length($inputString);
    for (my $i = 0; $i <= $strlen; $i++) {
my $curchar = substr($inputString, $i, 1);
if ($curchar =~ /\s/) {
} else {
    $outputString .= $curchar;
}
    }
    $outputString;
}
#################################################################################################################


#################################################################################################################
# Helper routine to remove all quotes from a string
#################################################################################################################
sub removeQuotes {
    my $inputString = shift;
    my $outputString;
    
    my $strlen = length($inputString);
    for (my $i = 0; $i <= $strlen; $i++) {
my $curchar = substr($inputString, $i, 1);
if ($curchar =~ /\'/) {
} elsif ($curchar =~ /\"/) {
} else {
    $outputString .= $curchar;
}
    }
    $outputString;
}
#################################################################################################################


#################################################################################################################
# Subroutine to replace a card of a given name in a pwscf input file
#################################################################################################################
sub replacePwscfCard {
    my $dataref = shift;
    my $secname = shift;
    my $replaceString = shift;
    
    my @replaceArray;
    if (ref($replaceString) eq 'ARRAY') {
@replaceArray = @{$replaceString};
    } elsif (ref($replaceString) eq 'SCALAR') {
@replaceArray = split(/(\n)/, removeWhiteSpaceAtBeginning(${$replaceString}));
    }

#    print "in replacePwscfCard, first line of \@replaceArray is: $replaceArray[0]\n";

    my $start = -1;
    my $stop = -1;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
        if ($line =~ $secname) {
            $start = $i;
} elsif ($start != -1 && $set == 0 && ((!($line =~ /[a-z]/) && !($line =~ /[0-9]/) && !($line =~ /\//)) || $line =~ /\{/)) {
            $stop = $i-1;
            $set = 1;
        }
        $i++
    }
    if ($start != -1 && $stop == -1) {
        $stop = $i+1;
    }
             
    if ($start == -1 && $stop == -1) {
#print "got here\n";
push(@{$dataref}, @replaceArray);
    } else {
splice @{$dataref}, $start, ($stop-$start), @replaceArray;
    }
    return $dataref; 
}
#################################################################################################################


#################################################################################################################
# Subroutine to replace a section of a pwscf input file
#################################################################################################################
sub replacePwscfSection {
    my $dataref = shift;
    my $secname = shift;
    my $replaceString = shift;
    
    my $start = -1;
    my $stop;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
        if ($line =~ $secname) {
            $start = $i;
        } elsif ($start != -1 && $set == 0 && $line =~ "\/" && !($line =~ "\=")) {
            $stop = $i;
            $set = 1;
        }
        $i++;
    }
    splice @{$dataref}, $start, ($stop-$start+1), @{$replaceString};
    return $dataref;
}
#################################################################################################################


#################################################################################################################
# Subroutine to get a section of a pwscf input file
#################################################################################################################
sub getPwscfSection {
    my $secname = shift;
    my $dataref = shift;
    my $outputString = "";
    my $start = -1;
    my $stop = -1;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
if ($line =~ /$secname/i) {
    $start = $i-1;
} elsif ($start != -1 && $set == 0 && $line =~ /\//) {
    $stop = $i;
    $set = 1;
}
$i++;
    }
    if ($start != -1 && $stop == -1) {
        $stop = $i;
    }
    for (my $i = $start+1; $i <= $stop; $i++) {
$outputString .= $dataref->[$i];
    }

    return $outputString; 
}
#################################################################################################################


#################################################################################################################
# Subroutine to get a card from a pwscf input file (does not include card name)
#################################################################################################################
sub getPwscfCard {
    my $secname = shift;
    my $dataref = shift;
    my $outputString = "";
    
    my $start = -1;
    my $stop = -1;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
if ($line =~ /$secname/i) {
            $start = $i;
#} elsif ($start != -1 && $set == 0 && ( (!($line =~ /[a-z]/) && !($line =~ /[0-9]/) && !($line =~ /\//) ) || $line =~ /\{/)) {
} elsif ($start != -1 && $set == 0 && ( ( !($line =~ /[0-9]/) && !($line =~ /\//) ) || $line =~ /\{/)) {
            $stop = $i;
            $set = 1;
        }
        $i++;
    }
    if ($start != -1 && $stop == -1) {
        $stop = $i;
    }
    for (my $i = $start+1; $i < $stop; $i++) {
$outputString .= $dataref->[$i];
    }

    return $outputString; 
}
#################################################################################################################


#################################################################################################################
# Subroutine to get a card from a pwscf input file (includes card name)
#################################################################################################################
sub getPwscfWholeCard {
    my $secname = shift;
    my $dataref = shift;
    my $outputString = "";
    
    my $start = -1;
    my $stop = -1;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
if ($line =~ /$secname/i) {
            $start = $i;
} elsif ($start != -1 && $set == 0 && ((!($line =~ /[a-z]/) && !($line =~ /[0-9]/) && !($line =~ /\//)) || $line =~ /\{/)) {
            $stop = $i;
            $set = 1;
        }
        $i++;
    }
    if ($start != -1 && $stop == -1) {
        $stop = $i;
    }
    for (my $i = $start; $i < $stop; $i++) {
$outputString .= $dataref->[$i];
    }

    return $outputString; 
}
#################################################################################################################


#################################################################################################################
# Subroutine to add lines to the end of a section in a pwscf input file
#################################################################################################################
sub addToPwscfSection {
    my $dataref = shift;
    my $secname = shift;
    my $replaceString = shift;
    
    my $start = -1;
    my $stop;
    my $i = 0;
    my $set = 0;
    foreach my $line (@{$dataref}) {
        if ($line =~ $secname) {
            $start = $i;
        } elsif ($start != -1 && $set == 0 && $line =~ "\/" && !($line =~ "\=")) {
            $stop = $i;
            $set = 1;
        }
        $i++;
    }
    splice @{$dataref}, $stop, 0, @{$replaceString};
    return $dataref;
}
#################################################################################################################



#################################################################################################################
# Subroutine to start with a pwscf input file and a list of kpoints and 
# write a pwscf input file to do a non self consistent calculation at
# those k-points (useful in preparing wavefunctions for qmc)
#################################################################################################################
sub genNscfPrep {
    my $dataref = shift;
    my $kptstring = shift;
    
    my $i;
    my $hasnosym = 0;
    for ($i = 0; $i < $#{$dataref}; $i++) {
        if (${$dataref}[$i] =~ /calculation/) {
            ${$dataref}[$i] = "    calculation = \'nscf\'\n";
        }
        if (${$dataref}[$i] =~ /nosym/) {
            $hasnosym=1;
        }
    }
    

    replacePwscfCard($dataref, "K_POINTS", \$kptstring);
    if (!$hasnosym) {
        my @addToSystem = ("    nosym=.true.\n    noinv=.true.\n");
        addToPwscfSection($dataref, "\&system", \@addToSystem);
    }
}
#################################################################################################################

#################################################################################################################
# Subroutine to start with a pwscf input file, new (super)cell translation vectors,
# new list of atom positions and a grid of k-points.  Use this to extract the independent
# supercell twists that we will use
#################################################################################################################
sub genPwscfSymm {
    my $ifile = shift;
    my $dataref = shift;
    my $kptstring = shift;
    
###################################################################################
# Variable Declarations to be used later
##################################################################################
    my ($calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies);
    my ($numAts, $numSpins, @cell_ptv, @atoms_name, @pseudoPotentials);
    my (@ionIds, %atNameToPP, @posArray);

    my ($baseName, @fdata, $topSpinDependentWvfcn);
    my (@atomsCharge, @atomsValence, @atomsAtomicNumbers);
    my ($numUpElec, $numDownElec, $numelec);
#####################################################################################

    # now use getSystemInformation to figure out the supercell atom positions
    # and primitive translation vectors
    getSystemInformation($ifile, $calcPrefix, $pseudoDir, $outdir, $celldim, $numSpecies,
 $numAts, $numSpins, \@cell_ptv, \@atoms_name, \@pseudoPotentials,
 \@ionIds, \%atNameToPP, \@posArray, $baseName, \@fdata, $topSpinDependentWvfcn,
 \@toptilingmatrix,, $targetSsize, \@atomsCharge, \@atomsValence, \@atomsAtomicNumbers,
 $numelec, $numUpElec, $numDownElec);

    changePwscfToken("ibrav", 0, $dataref);
    removePwscfToken('celldm\(1\)', $dataref);
    removePwscfToken('celldm\(2\)', $dataref);
    removePwscfToken('celldm\(3\)', $dataref);
    removePwscfToken('celldm\(4\)', $dataref);
    removePwscfToken('celldm\(5\)', $dataref);
    removePwscfToken('celldm\(6\)', $dataref);
    changePwscfToken("ecutwfc", 1, $dataref);
    changePwscfToken("ecutrho", 4, $dataref);
    changePwscfToken("nat", $numAts, $dataref);
    changePwscfToken("disk_io", "\'none\'", $dataref);


    my @addToSystem = ("    verbosity = \'high\'\n");
    addToPwscfSection($dataref, "\&control", \@addToSystem);

    replacePwscfCard($dataref, "K_POINTS", \$kptstring);

    my $newatomscard = "ATOMIC_POSITIONS crystal\n";
    for (my $i = 0; $i < $numAts; $i++) {
$newatomscard .= "$ionIds[$i]   $posArray[$i*3]  $posArray[$i*3+1]  $posArray[$i*3+2]\n";
    }  
    replacePwscfCard($dataref, "ATOMIC_POSITIONS", \$newatomscard);

    my $newcellparamscard = "CELL_PARAMETERS\n";
    $newcellparamscard .=  "$cell_ptv[0]  $cell_ptv[1]  $cell_ptv[2]\n";
    $newcellparamscard .=  "$cell_ptv[3]  $cell_ptv[4]  $cell_ptv[5]\n";
    $newcellparamscard .=  "$cell_ptv[6]  $cell_ptv[7]  $cell_ptv[8]\n";
    replacePwscfCard($dataref, "CELL_PARAMETERS", \$newcellparamscard);
}
#################################################################################################################



#################################################################################################################
# Subroutine to start with a pwscf input file and a list of kpoints and 
# write a pwscf input file to do a non self consistent calculation at
# those k-points.  Contrary to the above section, this one does a single
# converged step for the electronic iteration and prints the energy.
# this can be useful to get the DFT energy for a set of k-points to be
# used in a finite size correction
#################################################################################################################
sub genPseudoScfPrep {
    my $dataref = shift;
    my $kptstring = shift;

    my $i;
    my $hasnosym = 0;
    my $haswfcollect = 0;
    my $hasdiskio = 0;
    for ($i = 0; $i < $#{$dataref}; $i++) {
        if (${$dataref}[$i] =~ /nosym/) {
            $hasnosym=1;
        }
        if (${$dataref}[$i] =~ /disk_io/) {
            $hasdiskio=1;
            ${$dataref}[$i] = "    disk_io = \'none\'\n";
        }
        if (${$dataref}[$i] =~ /wf_collect/) {
            $haswfcollect=1;
            ${$dataref}[$i] = "    wf_collect =.false.\n";
        }
    }
    
    if (!$hasdiskio) {
        my @addToSystem = ("    disk_io = \'none\'\n");
        addToPwscfSection($dataref, "\&control", \@addToSystem);
    }
    if (!$haswfcollect) {
        my @addToSystem = ("    wf_collect =.false.\n");
        addToPwscfSection($dataref, "\&control", \@addToSystem);
    }
     
    my @pscfelec;
    push(@pscfelec, "&electrons \n");
    push(@pscfelec, "   mixing_mode = 'plain'\n");
    push(@pscfelec, "   mixing_beta = 0.0\n");
    push(@pscfelec, "   electron_maxstep=1\n");
    push(@pscfelec, "   diago_thr_init=1E-12\n");
    push(@pscfelec, "   startingpot='file'\n");
    push(@pscfelec, "   conv_thr = 1.0d-10\n");
    push(@pscfelec, "   mixing_fixed_ns = 5\n");
    push(@pscfelec, "/\n");

     
    if (!$hasnosym) {
        my @addToSystem = ("    nosym=.true.\n    noinv=.true.\n");
        addToPwscfSection($dataref, "\&system", \@addToSystem);
    }
    
    replacePwscfCard($dataref, "K_POINTS", \$kptstring);
    replacePwscfSection($dataref, "\&electrons", \@pscfelec);
}
#################################################################################################################

#################################################################################################################
# Subroutine to start with a pwscf input file and a list of kpoints and 
# write a pwscf input file to do a self consistent calculation with the 
# KZK functional.  Note that this will be confused to say the least by
# the multiple k-points (PWSCF will misidentify the cell volume).  Use
# with CAUTION!!!!
#################################################################################################################
sub genKzkPrep {
    my $dataref = shift;
    my $kptstring = shift;

    my $i;
    my $hasnosym = 0;
    my $haswfcollect = 0;
    my $hasdiskio = 0;
    for ($i = 0; $i < $#{$dataref}; $i++) {
        if (${$dataref}[$i] =~ /nosym/) {
            $hasnosym=1;
        }
        if (${$dataref}[$i] =~ /disk_io/) {
            $hasdiskio=1;
            ${$dataref}[$i] = "    disk_io = \'none\'\n";
        }
        if (${$dataref}[$i] =~ /wf_collect/) {
            $haswfcollect=1;
            ${$dataref}[$i] = "    wf_collect =.false.\n";
        }
    }
    
    if (!$hasdiskio) {
        my @addToSystem = ("    disk_io = \'none\'\n");
        addToPwscfSection($dataref, "\&control", \@addToSystem);
    }
    if (!$haswfcollect) {
        my @addToSystem = ("    wf_collect =.false.\n");
        addToPwscfSection($dataref, "\&control", \@addToSystem);
    }

    my @addToSystem = ("    input_dft='KZK'\n");
    addToPwscfSection($dataref, "\&system", \@addToSystem);

    if (!$hasnosym) {
        my @addToSystem = ("    nosym=.true.\n    noinv=.true.\n");
        addToPwscfSection($dataref, "\&system", \@addToSystem);
    }
    
    replacePwscfCard($dataref, "K_POINTS", \$kptstring);
}
#################################################################################################################



#################################################################################################################
# Subroutine to write an input file for pw2qmcpack
#################################################################################################################
sub genPw2qmcpack {
    my $dataref = shift;
    # output directory from pwscf calculation
    my $odir = shift;
    # calculation prefix from pwscf calculation
    my $cp = shift;
    
    $dataref->[0] = "&inputpp\n";
    $dataref->[1] = "   outdir=\'$odir\'\n";
    $dataref->[2] = "   prefix=\'$cp\'\n";
    $dataref->[3] = "   write_psir=.false.\n";
    $dataref->[4] = "/\n";
}
#################################################################################################################

#################################################################################################################
# Subroutine to write an input file for pw2casino
#################################################################################################################
sub genPw2casino {
    my $dataref = shift;
    # output directory from pwscf calculation
    my $odir = shift;
    # calculation prefix from pwscf calculation
    my $cp = shift;
    
    $dataref->[0] = "&inputpp\n";
    $dataref->[1] = "   outdir=\'$odir\'\n";
    $dataref->[2] = "   prefix=\'$cp\'\n";
    $dataref->[4] = "/\n";
}
#################################################################################################################

######################################################################
#
# Subroutines to get the supercell twists from either the SYMMETRY file,
# the output of the specialized pwscf calculations or the arguments to
# setup-qmc.pl as a last resort.
#
######################################################################

sub getSuperCellTwists {
    my $inFile = shift;
    my $targetSize = shift;
    my $superCellGridRef = shift;
    my $superCellShiftRef = shift;
    my $superCellTwists = shift;

    my $concatgrid = "$$superCellGridRef[0]-$$superCellGridRef[1]-$$superCellGridRef[2]";
    my $concatshift = "$$superCellShiftRef[0]-$$superCellShiftRef[1]-$$superCellShiftRef[2]";
    my $sstheader = "SCTWISTS $targetSsize $concatgrid $concatshift";

    my $pwscfsymfilename = $inputFile;
    $pwscfsymfilename =~ s/\.in//g;
    $pwscfsymfilename =~ s/-scf//g;
    $pwscfsymfilename .= "-SYM-supertwist$$superCellGridRef[0]$$superCellGridRef[1]$$superCellGridRef[2]";
    $pwscfsymfilename .= "-supershift$$superCellShiftRef[0]$$superCellShiftRef[1]$$superCellShiftRef[2]";
    if ($targetSize) {
$pwscfsymfilename .= "-S$targetSize";
    }
    $pwscfsymfilename .= ".out";


    # first check to see if there is an entry in the SYMMETRY file
    `touch SYMMETRY`;
    my $sstmarkerline = `sed -n '/$sstheader/=' SYMMETRY`;
    if ($sstmarkerline) {
print "Getting entry from SYMMETRY file\n";
chomp($sstmarkerline);
my $numtwistline = $sstmarkerline+1;
my $numtwists = `sed -n ${numtwistline}p SYMMETRY`;
# print "   Number of twists: $numtwists\n";
for (my $i = $numtwistline+1; $i < $numtwistline+1+$numtwists; $i++) {
  my $twistline = `sed -n ${i}p SYMMETRY`;
#  print "  Parsed line: $twistline\n";
  my @data = split(/\s+/, $twistline);
  push(@{$superCellTwists}, $data[0]);
  push(@{$superCellTwists}, $data[1]);
  push(@{$superCellTwists}, $data[2]);
  push(@{$superCellTwists}, $data[3]);
}
    } elsif (-e $pwscfsymfilename) {
# there is not an entry in the SYMMETRY file, see if there is a 
# PWSCF symmetry run that we can look at
print "Getting entry from pwscf output file\n";
my $line = `grep \"number of k points=\" $pwscfsymfilename`;
my @data = split(/\s+/, $line);
my $numkpts = $data[5];

my $kptline = `sed -n '/number of k points=/=' $pwscfsymfilename`;
chomp $kptline;
my $firstredkptlinenum = $kptline + 4 + $numkpts;
my $firstredkptline = `sed -n ${firstredkptlinenum}p $pwscfsymfilename`;

open (SYM, ">>SYMMETRY");
print SYM "$sstheader\n";
print SYM "$numkpts\n";
for (my $i = $firstredkptlinenum; $i < $firstredkptlinenum+$numkpts; $i++) {
    my $line = `sed -n ${i}p $pwscfsymfilename`;
    my @data = split(/\s+/, $line);
    push (@{$superCellTwists}, $data[5]);
    push (@{$superCellTwists}, $data[6]);
    my $processme = $data[7];
    chop($processme);
    chop($processme);
    push (@{$superCellTwists}, $processme);
    push (@{$superCellTwists}, $data[10]);
    print SYM "$data[5]  $data[6]  $processme   $data[10]\n";
}
close(SYM);
    } else {
print "No symmetry information found, doing the naive thing!\n";
    
## Loop over k-points in the supercell k-pt mesh
my @ns;
for ($ns[0] = 0; $ns[0] < $$superCellGridRef[0]; $ns[0]++) {
    for ($ns[1] = 0; $ns[1] < $$superCellGridRef[1]; $ns[1]++) {
for ($ns[2] = 0; $ns[2] < $$superCellGridRef[2]; $ns[2]++) {
    
    for (my $i = 0; $i < 3; $i++) {
if ($$superCellShiftRef[$i]) {
    $ns[$i] += 0.5;
}
    }

    my @shift;
    for (my $i = 0; $i < 3; $i++) {
$shift[$i] = $ns[$i] / $$superCellGridRef[$i];
    }
    push(@{$superCellTwists}, $shift[0]);
    push(@{$superCellTwists}, $shift[1]);
    push(@{$superCellTwists}, $shift[2]);
    push(@{$superCellTwists}, 1.0);
}
    }
}
    }
 
    # note the kpoints should alread be unique, just using this function to sort them really
    combineKpts($superCellTwists);

#    my $numtwists = ($#{$superCellTwists}+1)/4;
#    print "numtwists = $numtwists\n";
#    for (my $i = 0; $i < $numtwists; $i++) {
#print "k $i:   $$superCellTwists[$i*4]   $$superCellTwists[$i*4+1]   $$superCellTwists[$i*4+2],  weight = $$superCellTwists[$i*4+3]\n";
#    }
}


######################################################################
#
# Subroutines to get the k-points at which to run a given simulation to
# produce a given supercell with a given set of supercell k-points
#
######################################################################

sub getPrimCellKpts {
  my $ptvref = shift;
  my $superCellTws = shift;
  my $tilematref = shift;
  my $primCellTws = shift;
  my $superCellGridRef = shift;
  my $superCellShiftRef = shift;
  
  my @superrlvs;
  my @supercell;
  my @n;
  my $eps = 1e-11;

  `touch SYMMETRY`;
  my $concatgrid = "$$superCellGridRef[0]-$$superCellGridRef[1]-$$superCellGridRef[2]";
  my $concatshift = "$$superCellShiftRef[0]-$$superCellShiftRef[1]-$$superCellShiftRef[2]";  
  my $pctwistheader = "PCTWISTS $targetSsize $concatgrid $concatshift";
  my $pctwmarkerline = `sed -n '/$pctwistheader/=' SYMMETRY`;
  if ($pctwmarkerline) {
    print "Getting primitive cell k-points from SYMMETRY file\n";
    chomp ($pctwmarkerline);
    my $numkptline = $pctwmarkerline+1;
    my $numtwists = `sed -n ${numkptline}p SYMMETRY`;
    for (my $i = $numkptline+1; $i < $numkptline+1+$numtwists; $i++) {
      my $kptline = `sed -n ${i}p SYMMETRY`;
      my @data = split(/\s+/, $kptline);
      push (@{$primCellTws}, $data[0]);
      push (@{$primCellTws}, $data[1]);
      push (@{$primCellTws}, $data[2]);
      push (@{$primCellTws}, $data[3]);
    }
  } else {
    # There is not an entry in the SYMMETRY file, figure out the appropriate twists and weights locally
    print "Doing search to get primitive cell kpoints\n";
    print "tilemat = $$tilematref[0]  $$tilematref[1]  $$tilematref[2]  $$tilematref[3]  $$tilematref[4]  $$tilematref[5]  $$tilematref[6]  $$tilematref[7]  $$tilematref[8]\n";  
    print "PTVs = $$ptvref[0]  $$ptvref[1]  $$ptvref[2]  $$ptvref[3]  $$ptvref[4]  $$ptvref[5]  $$ptvref[6]  $$ptvref[7]  $$ptvref[8]\n";  
    getSuperCell(\@supercell, $ptvref, $tilematref);
    print "Got Supercell\n";
    print "Supercell = $supercell[0]  $supercell[1]  $supercell[2]  $supercell[3]  $supercell[4]  $supercell[5]  $supercell[6]  $supercell[7]  $supercell[8]\n";  
    getSuperRlvs(\@superrlvs, \@supercell);
    print "Finished getSuperRlvs\n";

    ## Loop over the requested supercell k-points
    for (my $kp = 0; $kp < $#{$superCellTws}; $kp += 4) {
      my @locprimCellTws;
#      print "Looking for primitive cell k-points corresponding to super cell kpoint: $$superCellTws[$kp]  $$superCellTws[$kp+1]  $$superCellTws[$kp+2]\n";
      my @shift;
      for (my $j = 0; $j < 3; $j++) {
$shift[$j] = $$superCellTws[$kp+$j];
      }
      
      ##Search through multiples of the supercell RLV's (with shifts)
      ## and find if they belong to the FBZ of the primitive cell
      my $nmax = 14;
      for ($n[0] = -$nmax; $n[0] <= $nmax; $n[0]++) { 
for ($n[1] = -$nmax; $n[1] <= $nmax; $n[1]++) { 
  for ($n[2] = -$nmax; $n[2] <= $nmax; $n[2]++) { 
    
    ## Get location of the G vector
    my @G = (0, 0, 0);
    for (my $i = 0; $i < 3; $i++) {
      for (my $j = 0; $j < 3; $j++) {
$G[$j] += ($n[$i] + $shift[$i]) * $superrlvs[$i*3+$j];
      }
    }
    
    ## Check if it is in the FBZ of the primitive lattice
    my $inFBZ = 1;
    for (my $i = 0; $i < 3; $i++) {
      my $dotval = 0;
      for (my $j = 0; $j < 3; $j++) {
$dotval += $$ptvref[$i*3+$j]*$G[$j];
      }
      if ($dotval < -1-$eps || $dotval > $eps) {
$inFBZ = 0;
      }
    }
    
    ## If it's already in the FBZ of the primitive lattice,
    ## check whether we have already included it in our list
    ## of kpoints.  If not, add it to primCellTws
    if ($inFBZ) {
      my @twist = MatVec3($ptvref, \@G);
      ## Get only the fractional part of the twist
      for (my $i = 0; $i < 3; $i++) {
$twist[$i] -= round($twist[$i]);
      }
      ## Turn all -0.5 components into 0.5's
      ## Also get turn small roundoff errors around 0 into 0
      for (my $i = 0; $i < 3; $i++) {
  if ($twist[$i] < -0.49999999999) {
      $twist[$i] += 1;
  }
  if (abs($twist[$i]) < $eps) {
      $twist[$i] = 0;
  }
      }
      
      for (my $i = 0; $i < 3; $i++) {
push(@locprimCellTws, $twist[$i]);
      }
      push (@locprimCellTws, $$superCellTws[$kp+3]);
    }
  }
}
      } # ends search for primitive cell k-vectors related to this supercell twist
      pruneKpts(\@locprimCellTws);
      for (my $i = 0; $i <= $#locprimCellTws; $i++) {
push (@{$primCellTws}, $locprimCellTws[$i]);
      }
      combineKpts($primCellTws);
    }
    
    # now write the results to the symmetry file so that they can be reused
    open (SYM, ">>SYMMETRY");
    print SYM "$pctwistheader\n";
    my $numkpts = ($#{$primCellTws}+1) / 4;
    print SYM "$numkpts\n";
    for (my $i = 0; $i < $#{$primCellTws}; $i += 4) {
      print SYM "${$primCellTws}[$i]  ${$primCellTws}[$i+1]  ${$primCellTws}[$i+2]  ${$primCellTws}[$i+3]\n";
    }
    close(SYM);
  }
}
    

######################################################################
#
# Subroutines to get the k-points at which to run a given simulation to
# produce a given supercell with a given supercell k-grid
#
######################################################################


sub getKVectors {
    my $kvecref = shift;
    my $ptvref = shift;
    my $tilematref = shift;
    my $supermeshref = shift;
    my $supershiftref = shift;
    my $eps = 1e-8;

    my @superrlvs;
    my @supercell;
    my @ns;
    my @n;
    my $eps = 0.00000000001;
    
    getSuperCell(\@supercell, $ptvref, $tilematref);
#    print "supercell is:\n";
#    print "$supercell[0]  $supercell[1]  $supercell[2]\n";
#    print "$supercell[3]  $supercell[4]  $supercell[5]\n";
#    print "$supercell[6]  $supercell[7]  $supercell[8]\n";
    getSuperRlvs(\@superrlvs, \@supercell);
#    print "supercell reciprocal lattice vectors:\n";
#    print "$superrlvs[0]  $superrlvs[1]  $superrlvs[2]\n";
#    print "$superrlvs[3]  $superrlvs[4]  $superrlvs[5]\n";
#    print "$superrlvs[6]  $superrlvs[7]  $superrlvs[8]\n";

    
    ## Loop over k-points in the supercell k-pt mesh
    for ($ns[0] = 0; $ns[0] < $$supermeshref[0]; $ns[0]++) {
for ($ns[1] = 0; $ns[1] < $$supermeshref[1]; $ns[1]++) {
    for ($ns[2] = 0; $ns[2] < $$supermeshref[2]; $ns[2]++) {

my @shift;
for (my $i = 0; $i < 3; $i++) {
    $shift[$i] = ($ns[$i] + $$supershiftref[$i]) / $$supermeshref[$i];
}

##Search through multiples of the supercell RLV's (with shifts)
## and find if they belong to the FBZ of the primitive cell
my $nmax = 14;
for ($n[0] = -$nmax; $n[0] <= $nmax; $n[0]++) { 
    for ($n[1] = -$nmax; $n[1] <= $nmax; $n[1]++) { 
for ($n[2] = -$nmax; $n[2] <= $nmax; $n[2]++) { 
    
    ## Get location of the G vector
    my @G = (0, 0, 0);
    for (my $i = 0; $i < 3; $i++) {
for (my $j = 0; $j < 3; $j++) {
    $G[$j] += ($n[$i] + $shift[$i]) * $superrlvs[$i*3+$j];
}
    }
    
    ## Check if it is in the FBZ of the primitive lattice
    my $inFBZ = 1;
    for (my $i = 0; $i < 3; $i++) {
my $dotval = 0;
for (my $j = 0; $j < 3; $j++) {
    $dotval += $$ptvref[$i*3+$j]*$G[$j];
}
if ($dotval < -1-$eps || $dotval > $eps) {
    $inFBZ = 0;
}
    }


    ## If it's already in the FBZ of the primitive lattice,
    ## check whether we have already included it in our list
    ## of kpoints.  If not, add it to kvecref
    if ($inFBZ) {
my $found = 0;
my @twist = MatVec3($ptvref, \@G);
#print "@G is in the first BZ\n";
#print "   Twist = @twist\n";

## Get only the fractional part of the twist
for (my $i = 0; $i < 3; $i++) {
    $twist[$i] -= round($twist[$i]);
}
#print "   Its fractional part is @twist\n";


## Check and see whether we have already found this one
for (my $j = 0; $j < ($#{$kvecref}+1)/3.0; $j++) {
    my @diff;
    my $diffsz = 0;
    for (my $i = 0; $i < 3; $i++) {
$diff[$i] = $twist[$i] - $$kvecref[$j*3+$i];
$diff[$i] -= round($diff[$i]);
$diffsz += $diff[$i]*$diff[$i];
    }
    $found = $found || $diffsz < $eps;
}
if (! $found) {
    for (my $i = 0; $i < 3; $i++) {
push(@$kvecref, $twist[$i]);
    }
}       
    }

} ## matches innermost loop over potential k-vectors 
    }
}
    } ## matches innermost loop over supercell twists
}
    }
}


##################################################################
# Subroutine to get atomic positions in the supercell given 
# the basis.  Note this is not differentiating between species,
# will probably have to call this multiple times with a "basis"
# for each speces.
##################################################################

sub getSupercellPos {
    my $posref = shift;
    my $basisref = shift;
    my $ptvref = shift;
    my $tilematref = shift;

    my @supercell;
    my @superrlvs;
    getSuperCell(\@supercell, $ptvref, $tilematref);
    getSuperRlvs(\@superrlvs, \@supercell);


#    print "supercell is:\n";
#    print "$supercell[0]  $supercell[1]  $supercell[2]\n";
#    print "$supercell[3]  $supercell[4]  $supercell[5]\n";
#    print "$supercell[6]  $supercell[7]  $supercell[8]\n";
#    print "supercell reciprocal lattice vectors:\n";
#    print "$superrlvs[0]  $superrlvs[1]  $superrlvs[2]\n";
#    print "$superrlvs[3]  $superrlvs[4]  $superrlvs[5]\n";
#    print "$superrlvs[6]  $superrlvs[7]  $superrlvs[8]\n\n";

    my $size = -1;
    while ($#{$posref} <= abs(getDet($tilematref))* ($#{$basisref}+1) && $size < 8) {
#print "size of posref = $#{$posref}\n";

$size++;
#print "in getSupercellPos, size = $size\n";
my @trial = (0, 0, 0);

for (my $basiselem = 0; $basiselem < ($#{$basisref}+1); $basiselem += 3) {
    my $added = 0;
    my @curpos = @{$basisref}[$basiselem .. $basiselem+2];
    
    my @n;
    for ($n[0] = -$size; $n[0] <= $size; $n[0]++) {
for ($n[1] = -$size; $n[1] <= $size; $n[1]++) {
    for ($n[2] = -$size; $n[2] <= $size; $n[2]++) {

my @pctrial;
## position of candidate in primitive cell reduced coordinates
for (my $i = 0; $i < 3; $i++) {
    $pctrial[$i] = $curpos[$i] + $n[$i];
}
## now in real space
@trial = (0, 0, 0);
for (my $i = 0; $i < 3; $i++) {
    for (my $j = 0; $j < 3; $j++) {
$trial[$i] += $$ptvref[$j*3+$i] * $pctrial[$j];
    }
}
## now get reduced coordinates in the supercell
my @ssRedCoord = MatVec3(\@superrlvs, \@trial);

my $eps = 1e-10;
if ($ssRedCoord[0] > -$eps && $ssRedCoord[0] < 1-$eps && $ssRedCoord[1] > -$eps && $ssRedCoord[1] < 1-$eps && $ssRedCoord[2] > -$eps && $ssRedCoord[2] < 1-$eps) {
    my $found = 0;
    
    for (my $j = 0; $j < ($#{$posref}+1)/3.0; $j++) {
my $diffsz;
for (my $i = 0; $i < 3; $i++) {
    $diffsz += abs($ssRedCoord[$i] - $$posref[$j*3+$i]);
}
if ($diffsz < $eps) {
    $found = 1;
}
    }
    
    if (! $found) {
$added++;
for (my $i = 0; $i < 3; $i++) {
    if (abs($ssRedCoord[$i]) < 1.0e-12) {
$ssRedCoord[$i] = 0.0;
    }
    push(@$posref, $ssRedCoord[$i]);
}
    }
}
    }
}
    }
}
    } 
}




##################################################################################
# Helper Subroutines for manipulating supercells and k-vectors
#################################################################################
sub getSuperCell {
    my $supercellref = shift;
    my $ptvref = shift;
    my $tilematref = shift;

    $$supercellref[0] = $$tilematref[0]*$$ptvref[0]+$$tilematref[1]*$$ptvref[3]+$$tilematref[2]*$$ptvref[6];
    $$supercellref[1] = $$tilematref[0]*$$ptvref[1]+$$tilematref[1]*$$ptvref[4]+$$tilematref[2]*$$ptvref[7];
    $$supercellref[2] = $$tilematref[0]*$$ptvref[2]+$$tilematref[1]*$$ptvref[5]+$$tilematref[2]*$$ptvref[8];
    $$supercellref[3] = $$tilematref[3]*$$ptvref[0]+$$tilematref[4]*$$ptvref[3]+$$tilematref[5]*$$ptvref[6];
    $$supercellref[4] = $$tilematref[3]*$$ptvref[1]+$$tilematref[4]*$$ptvref[4]+$$tilematref[5]*$$ptvref[7];
    $$supercellref[5] = $$tilematref[3]*$$ptvref[2]+$$tilematref[4]*$$ptvref[5]+$$tilematref[5]*$$ptvref[8];
    $$supercellref[6] = $$tilematref[6]*$$ptvref[0]+$$tilematref[7]*$$ptvref[3]+$$tilematref[8]*$$ptvref[6];
    $$supercellref[7] = $$tilematref[6]*$$ptvref[1]+$$tilematref[7]*$$ptvref[4]+$$tilematref[8]*$$ptvref[7];
    $$supercellref[8] = $$tilematref[6]*$$ptvref[2]+$$tilematref[7]*$$ptvref[5]+$$tilematref[8]*$$ptvref[8];
}


sub getSuperRlvs {
    my $superrlvref = shift;
    my $ptvref = shift;

    my $invdet = 1.0/getDet($ptvref);
    $$superrlvref[0] = $invdet * ($$ptvref[4]*$$ptvref[8] - $$ptvref[5]*$$ptvref[7]);
    $$superrlvref[1] = $invdet * ($$ptvref[5]*$$ptvref[6] - $$ptvref[3]*$$ptvref[8]);
    $$superrlvref[2] = $invdet * ($$ptvref[3]*$$ptvref[7] - $$ptvref[4]*$$ptvref[6]);
    $$superrlvref[3] = $invdet * ($$ptvref[2]*$$ptvref[7] - $$ptvref[1]*$$ptvref[8]);
    $$superrlvref[4] = $invdet * ($$ptvref[0]*$$ptvref[8] - $$ptvref[2]*$$ptvref[6]);
    $$superrlvref[5] = $invdet * ($$ptvref[1]*$$ptvref[6] - $$ptvref[0]*$$ptvref[7]);
    $$superrlvref[6] = $invdet * ($$ptvref[1]*$$ptvref[5] - $$ptvref[2]*$$ptvref[4]);
    $$superrlvref[7] = $invdet * ($$ptvref[2]*$$ptvref[3] - $$ptvref[0]*$$ptvref[5]);
    $$superrlvref[8] = $invdet * ($$ptvref[0]*$$ptvref[4] - $$ptvref[1]*$$ptvref[3]);
}


sub getDet {
    my $matref = shift;
    my $val = $$matref[0]*($$matref[4]*$$matref[8] - $$matref[7]*$$matref[5])
- $$matref[1]*($$matref[3]*$$matref[8] - $$matref[5]*$$matref[6])
+ $$matref[2]*($$matref[3]*$$matref[7]-$$matref[4]*$$matref[6]);
    return $val;
}


sub getCellProperties {
    my $locNumElec = shift;
    my $matref = shift;
    my $locTwoBodyRcut = \shift;
    my $locNumDens = \shift;

    my $loccellvol = abs(getDet($matref));
    $$locNumDens = $locNumElec / $loccellvol;
    
    my $rmin = 10000000000000000000000;
    for (my $i = -1; $i <= 1; $i++) {
for (my $j = -1; $j <= 1; $j++) {
    for (my $k = -1; $k <= 1; $k++) {
if ( ($i != 0) || ($j != 0) || ($k != 0) ) {
    my @d = (0.0, 0.0, 0.0);
    $d[0] = $i*$$matref[0] + $j*$$matref[3] + $k*$$matref[6];
    $d[1] = $i*$$matref[1] + $j*$$matref[4] + $k*$$matref[7];
    $d[2] = $i*$$matref[2] + $j*$$matref[5] + $k*$$matref[8];
    my $dist = 0.5 * sqrt($d[0]*$d[0] + $d[1]*$d[1] + $d[2]*$d[2]);
    if ($dist < $rmin) {
$rmin = $dist;
    }
}
    }
}
    }
    $$locTwoBodyRcut = $rmin;
}

sub getTwoBodyRPAJastrow {
    my $locNumDens = shift;
    my $locTwoBodyRcut = shift;
    my $locTwoBodySplinePts = shift;
    my $tuuc = shift;
    my $tudc = shift;

    my $wp = sqrt(4.0*3.14159265358979*$locNumDens);
    print "wp = $wp\n";
    print "numDens = $locNumDens\n";
    
    my $dr = $locTwoBodyRcut / ($locTwoBodySplinePts);
    my $i = 0;
    for (my $r = 0.02; $r <= $locTwoBodyRcut+0.000001; $r += $dr) {
$$tuuc[$i] = (0.5 / $wp / $r) * ( 1.0 - exp(-$r * sqrt($wp / 2.0)) ) * exp(-($r*2.0/$locTwoBodyRcut)**2);
$$tudc[$i] = (0.5 / $wp / $r) * ( 1.0 - exp(-$r * sqrt($wp)) ) * exp(-($r*2.0/$locTwoBodyRcut)**2);
$i++;
    }
}

sub FracPart {
    my $vecref = shift;
    my @intpart = IntPart($vecref);
    
    my @outarr;
    for (my $i = 0; $i <= $#{$vecref}; $i++) {
push(@outarr, $$vecref[$i] - $intpart[$i]);
    }
    return @outarr;
}

sub round {
    my $val = shift;
    
    if ($val > 0) { 
return floor($val + 0.5);
    } else {
return ceil($val - 0.5);
    }
}


## This only works with a 3x3 matrix and a 3 element vector
sub MatVec3 {
    my $matref = shift;
    my $vecref = shift;
    my @outarr = (0, 0, 0);

    for (my $i = 0; $i < 3; $i++) {
for (my $j = 0; $j < 3; $j++) {
    $outarr[$i] += $$matref[$i*3+$j] * $$vecref[$j];
}
    }
    return @outarr;
}

## This only works with 3 element vectors
sub dot {
    my $lvec = shift;
    my $rvec = shift;
    
    my $val = $$lvec[0] * $$rvec[0] + $$lvec[1] * $$rvec[1] + $$lvec[2] * $$rvec[2];
    return $val;
}

##################################################################################
# Helper Subroutine for analyzing a list of kpoints with respect to a tilematrix
#################################################################################

sub analyzeTwists {
    my $tilematref = shift;
    my $kvecref = shift;


    my @superIndex;
    my @superFracs;

    print "File contained " . ($#{$kvecref}+1)/3 . " primitive cell kvectors\n";
    print "Determinant of tilematrix = " . getDet($tilematref) . "\n";
    print "Hoping to find " . (($#{$kvecref}+1)/3 / getDet($tilematref)) . " supercell twists\n";

    for (my $i = 0; $i <= $#{$kvecref}; $i += 3) {
my @primTwist = ($$kvecref[$i],  $$kvecref[$i+1],  $$kvecref[$i+2]);
my @superTwist = MatVec3($tilematref, \@primTwist);
my @frac = FracPart(\@superTwist);

my $found = 0;
for (my $j = 0; $j < ($#superFracs+1); $j += 3) {
    my @diff;
    for (my $k = 0; $k < 3; $k++) {
push(@diff, $frac[$k] - $superFracs[$j+$k]);
    }
    my $diffsz = dot(\@diff, \@diff);
    if ($diffsz < 1.0e-6) {
$found = 1;
push(@superIndex, $j);
    } 
}
if (!$found) {
    push(@superIndex, $#superFracs+1);
    push(@superFracs, $frac[0]);
    push(@superFracs, $frac[1]);
    push(@superFracs, $frac[2]);
}
    }

    my $numSuperTwists = ($#superFracs+1)/3;
    print "Found $numSuperTwists distinct supercell twists\n";

    ## For each supercell twist, count how many primitive cell twists belong to it
    my %PrimTwistsPerSuperTwistIndex;
    for (my $i = 0; $i <= $#superIndex; $i ++) {
$PrimTwistsPerSuperTwistIndex{$superIndex[$i]}++;
    }
    my %FreqTwists;
    foreach (sort (keys %PrimTwistsPerSuperTwistIndex)) {
$FreqTwists{$PrimTwistsPerSuperTwistIndex{$_}}++
    }
    foreach (sort (keys %FreqTwists)) {
print "There are $FreqTwists{$_} SuperCell Twists that have $_ associated Primitive Cell Twists in the Wavefunction\n";
    } 
}




#################################################################################################################
# Subroutine to scrape valence charge and atom number from upf or ncpp pseudopotential files
#################################################################################################################
sub getPPInfo {
    my $ppfile = shift;
    my $valChg = \shift;
    my $atNum = \shift;
    
    my $pptype;
    if ($ppfile =~ /ncpp/i) {
$pptype = "NCPP";
    } elsif ($ppfile =~ /upf/i) {
$pptype = "UPF";
    }

    my @pplines;
    open(PPFILE, $ppfile) || die "Cannot open pseudopotential file $ppfile\n";;
    @pplines = <PPFILE>;
    close(PPFILE);

    if ($pptype eq "ncpp" || $pptype eq "NCPP") {
my @fields = split(/,/, $pplines[1]);
$$valChg = $fields[1];
        @fields = split(/\s+/, $pplines[2]);
$$atNum = $fields[0];
#print "   Got here with $$valChg valence electrons, and atomic Number $$atNum\n";
    } elsif ($pptype eq "upf" || $pptype eq "UPF") {
#print "this is a upf pseudopotential\n";
foreach my $line (@pplines) {
    if ($line =~ /z_valence/) {
my @fields = split(/\"/, removeWhiteSpaceAtBeginning($line));
$$valChg = $fields[1];
#print "In upf parser, valence is $$valChg\n";
    }
    if ($line =~ /Z valence/) {
my @fields = split(/\s+/, removeWhiteSpaceAtBeginning($line));
$$valChg = $fields[0];
#print "In upf parser, valence is $$valChg\n";
    }
    if ($line =~ /element/) {
my @fields = split(/\"/, removeWhiteSpaceAtBeginning($line));
my $elementName = $fields[1];
#print "In upf parser, element name is: $elementName\n";
$$atNum = AtNameToNumber($elementName);
    }
    if ($line =~ /Element/) {
my @fields = split(/\s+/, removeWhiteSpaceAtBeginning($line));
my $elementName = $fields[0];
#print "In upf parser, element name is: $elementName\n";
$$atNum = AtNameToNumber($elementName);
    }
}
    } else {
die ("Pseudopotential type $pptype is not recognized for pseudopotential $ppfile!\n");
    }
    
}
#################################################################################################################


#################################################################################################################
# Subroutine to go from element name given in a upf file to an atomic number.
# Note:  The list stops at #103, Lr
# Note:  This is only for the names as given in the UPF files.  Capitalization is rather specific
#################################################################################################################
sub AtNameToNumber {
    my $ename = shift;
    

    my $atNum = -1;
    my %hash = (
"H" => "1",     "He" => "2","Li" => "3","Be" => "4",
"B" => "5","C" => "6","N" => "7","O" => "8",
"F" => "9","Ne" => "10","Na" => "11","Mg" => "12",
"Al" => "13","Si" => "14","P" => "15","S" => "16",
"Cl" => "17","Ar" => "18","K" => "19","Ca" => "20",
"Sc" => "21","Ti" => "22","V" => "23","Cr" => "24",
"Mn" => "25","Fe" => "26","Co" => "27","Ni" => "28",
"Cu" => "29","Zn" => "30","Ga" => "31","Ge" => "32",
"As" => "33","Se" => "34","Br" => "35","Kr" => "36",
"Rb" => "37",   "Sr" => "38",   "Y" => "39",    "Zr" => "40",
"Nb" => "41",   "Mo" => "42",   "Tc" => "43",   "Ru" => "44",
"Rh" => "45",   "Pd" => "46",   "Ag" => "47",   "Cd" => "48",
"In" => "49",   "Sn" => "50",   "Sb" => "51",   "Te" => "52",
"I" => "53",    "Xe" => "54",   "Cs" => "55",   "Ba" => "56",
"La" => "57",   "Ce" => "58",   "Pr" => "59",   "Nd" => "60",
"Pm" => "61",   "Sm" => "62",   "Eu" => "63",   "Gd" => "64",
"Tb" => "65",   "Dy" => "66",   "Ho" => "67",   "Er" => "68",
"Tm" => "69",   "Yb" => "70",   "Lu" => "71",   "Hf" => "72",   
"Ta" => "73",   "W" => "74",    "Re" => "75",   "Os" => "76",
"Ir" => "77",   "Pt" => "78",   "Au" => "79",   "Hg" => "80",
"Tl" => "81",   "Pb" => "82",   "Bi" => "83",   "Po" => "84",
"At" => "85",   "Rn" => "86",   "Fr" => "87",   "Ra" => "88",
"Ac" => "89",   "Th" => "90",   "Pa" => "91",   "U" => "92",
"Np" => "93",   "Pu" => "94",   "Am" => "95",   "Cm" => "96",
"Bk" => "97",   "Cf" => "98",   "Es" => "99",   "Fm" => "100",
"Md" => "101",  "No" => "102",  "Lr" => "103");

    $atNum = $hash{$ename};
    if ($atNum < 0) {
die "Cannot find atomic number for $ename derived from upf pseudopotential.  Check AtNameToNumber routine\n";
    }
    $atNum;
}
#################################################################################################################

sub bySequence {
    my $left = $a;
    my $right = $b;
    $left =~ /\.s(\d\d\d)\.scalar\.dat/;
    my $leftseq = $1;
    $right =~ /\.s(\d\d\d)\.scalar\.dat/;
    my $rightseq = $1;
    $leftseq <=> $rightseq;
}

sub sortKpts {
# assumes it is getting a list where kpoints are stored as four consecutive entries 
# (first three are the kpt and fourth is the weight)
# This routine will sort these k-points numerically first on their first entry, next
# on the second and finally on their third;
  my $arrayref = shift;
  my @locarray;
  for (my $i = 0; $i <= $#{$arrayref}; $i+=4) {
    my $token = "$$arrayref[$i]   $$arrayref[$i+1]   $$arrayref[$i+2]   $$arrayref[$i+3]";
    push(@locarray, $token);
  }
  my @sortedarr = (sort byKpt @locarray);
  for (my $i = 0; $i <= $#sortedarr; $i++) {
    my $line = $sortedarr[$i];
    my @data = split(/\s+/,$line);
    $$arrayref[$i*4] = $data[0];
    $$arrayref[$i*4+1] = $data[1];
    $$arrayref[$i*4+2] = $data[2];
    $$arrayref[$i*4+3] = $data[3];
  }
}

sub byKpt {
  my $left = $a;
  my $right = $b;
  my @ldata = split(/\s+/,$left);
  my @rdata = split(/\s+/,$right);
  my $eps = 1e-7;
  if ($ldata[0] < ($rdata[0] - $eps)) {
    return -1;
  } elsif ($ldata[0] > ($rdata[0] + $eps)) {
    return 1;
  } elsif ($ldata[1] < ($rdata[1] - $eps)){
    return -1;
  } elsif ($ldata[1] > ($rdata[1] + $eps)) {
    return 1;
  } elsif ($ldata[2] < ($rdata[2] - $eps)) {
    return -1;
  } elsif ($ldata[2] > ($rdata[2] + $eps)) {
    return 1;
  } elsif ($ldata[3] < ($rdata[3] - $eps)) {
    return -1;
  } elsif ($ldata[3] > ($rdata[3] + $eps)) {
    return 1;
  }
  return 0;
}

sub combineKpts {
#returns a sorted list or kpts where any duplicates have been combined
  my $arrayref = shift;
  my $eps = 1e-7;
  sortKpts($arrayref);
  
  my @newarr;
  for (my $i = 0; $i < 4; $i++) {
    push (@newarr, $$arrayref[$i]);
  }
  for (my $i = 4; $i <= $#{$arrayref}; $i+=4) {
    if ((abs($$arrayref[$i] - $$arrayref[$i-4]) < $eps) && 
(abs($$arrayref[$i+1] - $$arrayref[$i-3]) < $eps) && 
(abs($$arrayref[$i+2] - $$arrayref[$i-2]) < $eps)) {
      # then this kpt is equal to the one before it
      # rather than add it to newarr, just add it's weight to
      # the last element of newarr
      $newarr[$#newarr] += $$arrayref[$i+3];
    } else {
      for (my $j = 0; $j < 4; $j++) {
push(@newarr, $$arrayref[$i+$j]);
      }
    }
  }
  @{$arrayref} = @newarr;
} 

sub pruneKpts {
#returns a sorted list or kpts where any duplicates have been pruned
  my $arrayref = shift;
  my $eps = 1e-7;
  sortKpts($arrayref);
  
  my @newarr;
  for (my $i = 0; $i < 4; $i++) {
    push (@newarr, $$arrayref[$i]);
  }
  for (my $i = 4; $i <= $#{$arrayref}; $i+=4) {
    if ((abs($$arrayref[$i] - $$arrayref[$i-4]) < $eps) && 
(abs($$arrayref[$i+1] - $$arrayref[$i-3]) < $eps) && 
(abs($$arrayref[$i+2] - $$arrayref[$i-2]) < $eps)) {
      # then this kpt is equal to the one before it
      # just skip it
    } else {
      for (my $j = 0; $j < 4; $j++) {
push(@newarr, $$arrayref[$i+$j]);
      }
    }
  }
  @{$arrayref} = @newarr;
} 

#################################################################################################################

sub globalUsage {
    my $usage;
    $usage .= "\nsetup-qmc.pl handles all aspects of generating the input files for various\n";
    $usage .= "aspects of qmc calculations of solids starting from pwscf input files\n";
    $usage .= "\n";
    $usage .= "setup-qmc.pl is invoked as:\n";
    $usage .= "  setup-qmc.pl --Function [--suboptions] pwscf-infile.in\n\n";
    $usage .= "Where --Function is one of the following:\n";
    $usage .= "   --testwvfcn (analyze eshdf wavefunction with regards to a particular supercell and twists)\n";
    $usage .= "   --gettilemat (get a tilematrix with optimizes the supercell shape with respect to simulation cell radius)\n";
    $usage .= "   --genwfn (create input files to generate suitable trial wavefunctions for qmcpack)\n";
    $usage .= "   --genfsdft (create input files to get non-self consistent energy, pw2casino and kzk)\n";
    $usage .= "   --splconv (test convergence of spline spacing)\n";
    $usage .= "   --optwvfcn (Optimize wavefunction (just jastrow for now))\n";
    $usage .= "   --convdmctstep (test convergence of the dmc timestep) \n";
    $usage .= "   --dmc (dmc calculation)\n\n";
    die ($usage);
}

sub NSCFUsage {
    my $usage;
    $usage .= "\n--genwfn  generates input files for pw.x, and pw2qmcpack.x so that a\n";
    $usage .= "          non self consistent wavefunction can be generated after a\n";
    $usage .= "          fully converged dft calculation has been done\n\n";
    $usage .= "To specify the size of the supercell (if there is to be one) give one \n";
    $usage .= "of these options:\n";
    $usage .= "   --supercellsize n (performs search for optimal supercell with n copies\n";
    $usage .= "                      of the primitive cell)\n";
    $usage .= "   --tilemat i i i i i i i i i (uses the specified tile matrix to get supercell\n";
    $usage .= "                                where the 9 i's are integers)\n\n";
    $usage .= "Can also optionally run at a series of twists in the supercell specified with\n";
    $usage .= "the following keywords:\n";
    $usage .= "   --kpoint f f f (Do calculation at twist given by 3 floats (in reduced coordinates))\n";
    $usage .= "   --kshift f f f (Shift the grid of supercell k-points by f f f)\n";
    $usage .= "   --kmesh i i i  (Generate wfn for ixixi mesh of k-points of the supercell)\n\n";
    die ($usage);
}

sub FSDFTUsage {
    my $usage;
    $usage .= "\n--genfsdft  generates input files for pw.x, and pw2qmcpack.x so that the\n";
    $usage .= "            energy of the non self-consistend wavefunction can found within dft\n";
    $usage .= "            Also -pw2casino.in and a pw.x input file to attempt to get the kzk\n";
    $usage .= "            energy for this size cell is written\n\n";
    $usage .= "To specify the size of the supercell (if there is to be one) give one \n";
    $usage .= "of these options:\n";
    $usage .= "   --supercellsize n (performs search for optimal supercell with n copies\n";
    $usage .= "                      of the primitive cell)\n";
    $usage .= "   --tilemat i i i i i i i i i (uses the specified tile matrix to get supercell\n";
    $usage .= "                                where the 9 i's are integers)\n\n";
    $usage .= "Can also optionally run at a series of twists in the supercell specified with\n";
    $usage .= "the following keywords:\n";
    $usage .= "   --kpoint f f f (Do calculation at twist given by 3 floats (in reduced coordinates))\n";
    $usage .= "   --kshift f f f (Shift the grid of supercell k-points by f f f)\n";
    $usage .= "   --kmesh i i i  (Generate wfn for ixixi mesh of k-points of the supercell)\n\n";
    die ($usage);
}

