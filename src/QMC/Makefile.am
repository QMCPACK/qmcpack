QMC_COMMON_SRCS =  QMCDriver.cpp\
VMC.cpp\
VMC_OPT.cpp\
MolecuDMC.cpp\
WaveFunctionTester.cpp\
../Estimators/ScalarEstimatorManager.cpp \
QMCApps.cpp

if USE_MOLECUAPPS
QMC_MOLECU =\
    ../QMCWaveFunctions/OrbitalBuilderBase.cpp \
    ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp \
    ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp\
    ../QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp\
    ../QMCWaveFunctions/MolecularOrbitals/RadialGridFunctorBuilder.cpp\
    ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp\
    ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp\
    ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp\
    MolecuApps.cpp
endif

#if USE_QDAPPS
#QMC_QD = \
#   ../QMCWaveFunctions/OrbitalBuilderBase.cpp \
#    ../Numerics/Spline3D/Grid1D.cpp \
#    ../Numerics/Spline3D/Grid3D.cpp \
#    ../Numerics/Spline3D/CubicSpline.cpp \
#    ../Numerics/Spline3D/TriCubicSpline.cpp \
#    ../Numerics/Spline3D/TriCubicSplineSet.cpp \
#    ../QMCWaveFunctions/TriCubicSplineBuilder.cpp \
#    ../QMCWaveFunctions/QDwfBuilder.cpp \
#    ../QMCWaveFunctions/QDwf.cpp \
#    QDApps.cpp \
#endif

bin_PROGRAMS = qmcapp
qmcapp_SOURCES = \
$(QMC_MOLECU) \
$(QMC_QD) \
$(QMC_COMMON_SRCS)\
main.cpp

INCLUDES = -I$(top_srcdir)/src
qmcapp_LDADD = @QMCBASE_LIBS@ 
qmcapp_LDFLAGS = @LDFLAGS@ 
qmcapp_DEPENDENCIES = ${top_builddir}/src/libqmcbase.a

