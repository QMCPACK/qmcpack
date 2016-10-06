#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory 
#//////////////////////////////////////////////////////////////////////////////////////


SRCDIR = src
BUILDDIR = build
COMMONDIR = src/common

#Intel Compiler with MKL build
#CXX = icpc
#CXXFLAGS = -O3 -g -restrict -mkl
#BLASLIBS =

# g++ build using installed MKL
CXX = g++
CXXFLAGS = -O3 -g -Drestrict=__restrict__ 
BLASLIBS = -L/usr/local/apps/intel/Compiler/14.0/106/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -pthread




# dummy rule to build everything
all: $(BUILDDIR)/ppconvert

# useful macros
PPCONVERT_OBJS = $(BUILDDIR)/CubicSpline.o $(BUILDDIR)/NLPPClass.o $(BUILDDIR)/ParseCommand.o $(BUILDDIR)/ParserClass.o \
$(BUILDDIR)/XMLWriterClass2.o

PPCONVERT_COMMON_OBJS = $(BUILDDIR)/IO.o $(BUILDDIR)/IOASCII.o \
$(BUILDDIR)/Functionals.o $(BUILDDIR)/Fitting.o $(BUILDDIR)/MatrixOps.o $(BUILDDIR)/Grid.o $(BUILDDIR)/GKIntegration.o \
$(BUILDDIR)/CubicSplineCommon.o $(BUILDDIR)/RadialWF.o $(BUILDDIR)/DFTAtom.o $(BUILDDIR)/Potential.o $(BUILDDIR)/CoulombPot.o \
$(BUILDDIR)/SplinePot.o $(BUILDDIR)/ScreenedPot.o $(BUILDDIR)/PotentialBase.o $(BUILDDIR)/GeneralPot.o $(BUILDDIR)/NLPP.o

COMMON_IO_HEADERS = $(COMMONDIR)/IO.h $(COMMONDIR)/IOBase.h $(COMMONDIR)/IOVar.h $(COMMONDIR)/IOVarBase.h \
$(COMMONDIR)/IOVarASCII.h $(COMMONDIR)/IOASCII.h

# rules for each executable
$(BUILDDIR)/ppconvert: $(PPCONVERT_OBJS) $(PPCONVERT_COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(BLASLIBS) -o $@

# rules for objects in src/
$(BUILDDIR)/CubicSpline.o: $(SRCDIR)/CubicSpline.cc $(SRCDIR)/CubicSpline.h $(SRCDIR)/GeneralGrid.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(SRCDIR) -c $(SRCDIR)/CubicSpline.cc -o $@
$(BUILDDIR)/NLPPClass.o: $(SRCDIR)/NLPPClass.cc $(SRCDIR)/NLPPClass.h $(SRCDIR)/CubicSpline.h $(SRCDIR)/GeneralGrid.h $(SRCDIR)/XMLWriterClass2.h $(SRCDIR)/ParserClass.h $(SRCDIR)/ParseCommand.h $(COMMON_IO_HEADERS) $(COMMONDIR)/DFTAtom.h $(COMMONDIR)/AtomBase.h $(COMMONDIR)/RadialWF.h $(COMMONDIR)/Potential.h $(COMMONDIR)/CoulombPot.h $(COMMONDIR)/ScreenedPot.h $(COMMONDIR)/SplinePot.h $(COMMONDIR)/GeneralPot.h $(COMMONDIR)/NLPP.h $(COMMONDIR)/PotentialBase.h $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(SRCDIR) -c $(SRCDIR)/NLPPClass.cc -o $@
$(BUILDDIR)/ParseCommand.o: $(SRCDIR)/ParseCommand.cc $(SRCDIR)/ParseCommand.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(SRCDIR) -c $(SRCDIR)/ParseCommand.cc -o $@
$(BUILDDIR)/ParserClass.o: $(SRCDIR)/ParserClass.cc $(SRCDIR)/ParserClass.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(SRCDIR) -c $(SRCDIR)/ParserClass.cc -o $@
$(BUILDDIR)/XMLWriterClass2.o: $(SRCDIR)/XMLWriterClass2.cc $(SRCDIR)/XMLWriterClass2.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(SRCDIR) -c $(SRCDIR)/XMLWriterClass2.cc -o $@


# rules for objects in src/common (mostly derived from libcommon)
$(BUILDDIR)/CoulombPot.o: $(COMMONDIR)/CoulombPot.cc $(COMMONDIR)/CoulombPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/CoulombPot.cc -o $@
$(BUILDDIR)/CubicSplineCommon.o: $(COMMONDIR)/CubicSplineCommon.cc $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMON_IO_HEADERS) $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/CubicSplineCommon.cc -o $@
$(BUILDDIR)/DFTAtom.o: $(COMMONDIR)/DFTAtom.cc $(COMMONDIR)/DFTAtom.h $(COMMONDIR)/AtomBase.h $(COMMONDIR)/RadialWF.h $(COMMONDIR)/Potential.h $(COMMONDIR)/CoulombPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/ScreenedPot.h $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h $(COMMONDIR)/SplinePot.h $(COMMONDIR)/GeneralPot.h $(COMMONDIR)/NLPP.h $(COMMONDIR)/RungeKutta.h $(COMMONDIR)/Functionals.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/DFTAtom.cc -o $@
$(BUILDDIR)/Fitting.o: $(COMMONDIR)/Fitting.cc $(COMMONDIR)/Fitting.h $(COMMONDIR)/Blitz.h $(COMMONDIR)/MatrixOps.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/Fitting.cc -o $@
$(BUILDDIR)/Functionals.o: $(COMMONDIR)/Functionals.cc $(COMMONDIR)/Functionals.h $(COMMONDIR)/Blitz.h $(SRCDIR)/config.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/Functionals.cc -o $@
$(BUILDDIR)/GeneralPot.o: $(COMMONDIR)/GeneralPot.cc $(COMMONDIR)/GeneralPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/GeneralPot.cc -o $@
$(BUILDDIR)/Grid.o: $(COMMONDIR)/Grid.cc $(COMMONDIR)/Grid.h $(COMMON_IO_HEADERS) $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/Grid.cc -o $@
$(BUILDDIR)/GKIntegration.o: $(COMMONDIR)/GKIntegration.cc $(COMMONDIR)/GKIntegration.h $(COMMONDIR)/Standard.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/GKIntegration.cc -o $@
$(BUILDDIR)/IO.o: $(COMMONDIR)/IO.cc $(COMMON_IO_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/IO.cc -o $@
#$(BUILDDIR)/IOHDF5.o: $(COMMONDIR)/IOHDF5.cc $(COMMON_IO_HEADERS)
#	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(HDF_INCLUDEDIR) -c $(COMMONDIR)/IOHDF5.cc -o $@
#$(BUILDDIR)/IOVarHDF5.o: $(COMMONDIR)/IOVarHDF5.cc $(COMMON_IO_HEADERS)
#	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -I$(HDF_INCLUDEDIR) -c $(COMMONDIR)/IOVarHDF5.cc -o $@
$(BUILDDIR)/IOASCII.o: $(COMMONDIR)/IOASCII.cc $(COMMON_IO_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/IOASCII.cc -o $@
$(BUILDDIR)/MatrixOps.o: $(COMMONDIR)/MatrixOps.cc $(COMMONDIR)/MatrixOps.h $(COMMONDIR)/Blitz.h $(SRCDIR)/config.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/MatrixOps.cc -o $@
$(BUILDDIR)/NLPP.o: $(COMMONDIR)/NLPP.cc $(COMMONDIR)/NLPP.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h $(COMMONDIR)/Standard.h $(COMMONDIR)/MatrixOps.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/NLPP.cc -o $@
$(BUILDDIR)/Potential.o: $(COMMONDIR)/Potential.cc $(COMMONDIR)/Potential.h $(COMMONDIR)/CoulombPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/ScreenedPot.h $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h $(COMMONDIR)/SplinePot.h $(COMMONDIR)/GeneralPot.h $(COMMONDIR)/NLPP.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/Potential.cc -o $@
$(BUILDDIR)/PotentialBase.o: $(COMMONDIR)/PotentialBase.cc $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS)
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/PotentialBase.cc -o $@
$(BUILDDIR)/RadialWF.o: $(COMMONDIR)/RadialWF.cc $(COMMONDIR)/RadialWF.h $(COMMONDIR)/Potential.h $(COMMONDIR)/CoulombPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/ScreenedPot.h $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h $(COMMONDIR)/SplinePot.h $(COMMONDIR)/GeneralPot.h $(COMMONDIR)/NLPP.h $(COMMONDIR)/RungeKutta.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/RadialWF.cc -o $@
$(BUILDDIR)/ScreenedPot.o: $(COMMONDIR)/ScreenedPot.cc $(COMMONDIR)/ScreenedPot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/ScreenedPot.cc -o $@
$(BUILDDIR)/SplinePot.o: $(COMMONDIR)/SplinePot.cc $(COMMONDIR)/SplinePot.h $(COMMONDIR)/PotentialBase.h $(COMMON_IO_HEADERS) $(COMMONDIR)/CubicSplineCommon.h $(COMMONDIR)/Grid.h $(COMMONDIR)/Blitz.h
	$(CXX) $(CXXFLAGS) -I$(COMMONDIR) -c $(COMMONDIR)/SplinePot.cc -o $@

clean:
	rm -f $(BUILDDIR)/*
