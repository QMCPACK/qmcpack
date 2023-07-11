#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2022 QMCPACK developers.
#//
#// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
#//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//////////////////////////////////////////////////////////////////////////////////////

add_subdirectory(detail)

set(WFBASE_SRCS
    OptimizableFunctorBase.cpp
    VariableSet.cpp
    WaveFunctionPool.cpp
    WaveFunctionComponent.cpp
    WaveFunctionComponentBuilder.cpp
    SPOSetBuilder.cpp
    SPOInfo.cpp
    SPOSetInfo.cpp
    SPOSetInputInfo.cpp
    SPOSet.cpp
    CompositeSPOSet.cpp
    HarmonicOscillator/SHOSet.cpp
    HarmonicOscillator/SHOSetBuilder.cpp
    ExampleHeBuilder.cpp
    ExampleHeComponent.cpp)

if(NOT QMC_COMPLEX)
  set(WFBASE_SRCS ${WFBASE_SRCS} RotatedSPOs.cpp)
endif(NOT QMC_COMPLEX)

if(QMC_COMPLEX)
  set(WFBASE_SRCS ${WFBASE_SRCS} SpinorSet.cpp)
endif(QMC_COMPLEX)
########################
# build jastrows
########################
#common jastrows
set(JASTROW_SRCS
    Jastrow/JastrowBuilder.cpp
    Jastrow/kSpaceJastrow.cpp
    Jastrow/kSpaceJastrowBuilder.cpp
    Jastrow/RadialJastrowBuilder.cpp
    Jastrow/CountingJastrowBuilder.cpp
    Jastrow/RPAJastrow.cpp
    Jastrow/J1OrbitalSoA.cpp
    LatticeGaussianProduct.cpp
    LatticeGaussianProductBuilder.cpp)
set(JASTROW_OMPTARGET_SRCS
    Jastrow/TwoBodyJastrow.cpp
    Jastrow/BsplineFunctor.cpp)
set(FERMION_SRCS ${FERMION_SRCS} ElectronGas/FreeOrbital.cpp ElectronGas/FreeOrbitalBuilder.cpp)

# wavefunctions only availbale to 3-dim problems
if(OHMMS_DIM MATCHES 3)

  set(JASTROW_SRCS ${JASTROW_SRCS} Jastrow/eeI_JastrowBuilder.cpp Jastrow/CountingJastrowBuilder.cpp)

  set(FERMION_SRCS ${FERMION_SRCS} LCAO/LCAOrbitalSet.cpp LCAO/LCAOrbitalBuilder.cpp LCAO/MultiQuinticSpline1D.cpp
                   LCAO/AOBasisBuilder.cpp LCAO/SoaLocalizedBasisSet.cpp)
  if(QMC_COMPLEX)
    set(FERMION_SRCS ${FERMION_SRCS} LCAO/LCAOSpinorBuilder.cpp)
  else(QMC_COMPLEX)
    #LCAO cusp correction is not ready for complex
    set(FERMION_SRCS ${FERMION_SRCS} LCAO/LCAOrbitalSetWithCorrection.cpp
                     LCAO/CuspCorrectionConstruction.cpp LCAO/SoaCuspCorrection.cpp)
  endif(QMC_COMPLEX)

  if(HAVE_EINSPLINE)
    if(NOT MIXED_PRECISION)
      set(FERMION_SRCS ${FERMION_SRCS} EinsplineSet.cpp)
    endif(NOT MIXED_PRECISION)
    set(FERMION_SRCS
        ${FERMION_SRCS}
        EinsplineSetBuilderCommon.cpp
        EinsplineSetBuilderOld.cpp
        AtomicOrbital.cpp
        EinsplineSetBuilderReadBands_ESHDF.cpp
        EinsplineSetBuilderESHDF.fft.cpp
        EinsplineSetBuilder_createSPOs.cpp
        BsplineFactory/createComplexDouble.cpp
        BsplineFactory/createComplexSingle.cpp
        BsplineFactory/HybridRepCenterOrbitals.cpp
        BandInfo.cpp
        BsplineFactory/BsplineReaderBase.cpp)
    set(FERMION_OMPTARGET_SRCS Fermion/DiracDeterminantBatched.cpp Fermion/MultiDiracDeterminant.2.cpp)
    if(QMC_COMPLEX)
      set(FERMION_SRCS ${FERMION_SRCS} EinsplineSpinorSetBuilder.cpp BsplineFactory/SplineC2C.cpp)
      set(FERMION_OMPTARGET_SRCS ${FERMION_OMPTARGET_SRCS} BsplineFactory/SplineC2COMPTarget.cpp)
    else(QMC_COMPLEX)
      set(FERMION_SRCS ${FERMION_SRCS} BsplineFactory/createRealSingle.cpp BsplineFactory/createRealDouble.cpp
        BsplineFactory/SplineC2R.cpp BsplineFactory/SplineR2R.cpp)
      set(FERMION_OMPTARGET_SRCS ${FERMION_OMPTARGET_SRCS} BsplineFactory/SplineC2ROMPTarget.cpp)
    endif(QMC_COMPLEX)

  endif(HAVE_EINSPLINE)

  # plane wave SPO
  set(FERMION_SRCS ${FERMION_SRCS} PlaneWave/PWBasis.cpp PlaneWave/PWParameterSet.cpp PlaneWave/PWOrbitalBuilder.cpp)
  if(QMC_COMPLEX)
    set(FERMION_SRCS ${FERMION_SRCS} PlaneWave/PWOrbitalSet.cpp)
  else()
    set(FERMION_SRCS ${FERMION_SRCS} PlaneWave/PWRealOrbitalSet.cpp)
  endif(QMC_COMPLEX)

  if(NOT QMC_COMPLEX)
    set(FERMION_SRCS ${FERMION_SRCS} AGPDeterminant.cpp AGPDeterminantBuilder.cpp)
  endif(NOT QMC_COMPLEX)

endif(OHMMS_DIM MATCHES 3)

set(FERMION_SRCS
    ${FERMION_SRCS}
    Fermion/DiracDeterminant.cpp
    Fermion/MultiDiracDeterminant.cpp
    Fermion/SlaterDet.cpp
    Fermion/SlaterDetBuilder.cpp
    Fermion/BackflowBuilder.cpp
    Fermion/BackflowTransformation.cpp
    Fermion/DiracDeterminantWithBackflow.cpp
    Fermion/SlaterDetWithBackflow.cpp
    SPOSetBuilderFactory.cpp
    TrialWaveFunction.cpp
    TWFdispatcher.cpp
    TWFFastDerivWrapper.cpp
    TWFGrads.cpp
    WaveFunctionFactory.cpp)
set(FERMION_OMPTARGET_SRCS ${FERMION_OMPTARGET_SRCS} Fermion/MultiSlaterDetTableMethod.cpp)

####################################
# create libqmcwfs
####################################
if(USE_OBJECT_TARGET)
  add_library(qmcwfs OBJECT ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
else(USE_OBJECT_TARGET)
  add_library(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
endif(USE_OBJECT_TARGET)

target_include_directories(qmcwfs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcwfs PUBLIC qmcutil qmcparticle platform_runtime platform_LA)
target_link_libraries(qmcwfs PRIVATE einspline Math::FFTW3)

if(USE_OBJECT_TARGET)
  add_library(qmcwfs_omptarget OBJECT ${JASTROW_OMPTARGET_SRCS} ${FERMION_OMPTARGET_SRCS})
else(USE_OBJECT_TARGET)
  add_library(qmcwfs_omptarget ${JASTROW_OMPTARGET_SRCS} ${FERMION_OMPTARGET_SRCS})
endif(USE_OBJECT_TARGET)

target_include_directories(qmcwfs_omptarget PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcwfs_omptarget PUBLIC qmcutil qmcparticle containers platform_LA)

target_compile_options(qmcwfs_omptarget PRIVATE "$<$<BOOL:${ENABLE_OFFLOAD_CLANG_DEBUG_O3}>:$<$<CONFIG:DEBUG>:-O3>>")

target_link_libraries(qmcwfs PUBLIC qmcwfs_omptarget)

if(ENABLE_CUDA)
  target_link_libraries(qmcwfs PRIVATE qmcwfs_cuda)
elseif(ENABLE_SYCL)
  target_link_libraries(qmcwfs PRIVATE qmcwfs_sycl)
endif()

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
