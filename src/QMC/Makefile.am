

bin_PROGRAMS = qmc
qmc_SOURCES = qmc.cpp DOMProcessor.cpp \
     QMCDriver.cpp VMC.cpp DMC.cpp \
     EstimatorManager.cpp \
  ../Utilities/RandRandom.cpp ../Utilities/RandomGenerator.cpp \
  ../Utilities/OhmmsObject.cpp \
  ../Particle/ParticleBase.cpp \
  ../Particle/MCWalkerConfiguration.cpp \
  ../Particle/DistanceTable.cpp \
  ../QMCHamiltonians/QMCHamiltonianBase.cpp \
  ../QMCWaveFunctions/HFAtomicSTOSet.cpp \
  ../QMCWaveFunctions/HFAtomicSTOSet.IO.cpp \
  ../QMCWaveFunctions/HeSTOClementiRottie.cpp \
  ../QMCWaveFunctions/RadialGridFunctorBuilder.cpp \
  ../QMCWaveFunctions/MOSonRadialGrid.cpp \
  ../QMCWaveFunctions/MOSwithSTO.cpp \
  ../QMCWaveFunctions/MolecularOrbitalBuilder.cpp \
  ../QMCWaveFunctions/TrialWaveFunction.cpp \
  ../QMCWaveFunctions/TrialWaveFunctionBuilder.cpp

INCLUDES = -I$(top_srcdir)/src
LDADD = -llapack -lblas -lg2c -lxml2 -lhdf5 -lz

