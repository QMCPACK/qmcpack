#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 Jeongnim Kim and QMCPACK developers.
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
    Jastrow/J2OrbitalSoA.cpp
    LatticeGaussianProduct.cpp
    LatticeGaussianProductBuilder.cpp)

if(ENABLE_OFFLOAD)
  set(JASTROW_SRCS ${JASTROW_SRCS} Jastrow/J2OMPTarget.cpp)
endif()

if(QMC_COMPLEX)
  set(FERMION_SRCS ${FERMION_SRCS} ElectronGas/ElectronGasComplexOrbitalBuilder.cpp)
else(QMC_COMPLEX)
  set(FERMION_SRCS ${FERMION_SRCS} ElectronGas/ElectronGasOrbitalBuilder.cpp)

endif(QMC_COMPLEX)

# wavefunctions only availbale to 3-dim problems
if(OHMMS_DIM MATCHES 3)

  set(JASTROW_SRCS ${JASTROW_SRCS} Jastrow/eeI_JastrowBuilder.cpp Jastrow/CountingJastrowBuilder.cpp)

  set(FERMION_SRCS
      ${FERMION_SRCS}
      LCAO/LCAOrbitalSet.cpp
      LCAO/LCAOrbitalBuilder.cpp
      LCAO/SoaSphericalTensor.cpp
      LCAO/MultiQuinticSpline1D.cpp
      LCAO/AOBasisBuilder.cpp
      LCAO/SoaLocalizedBasisSet.cpp)
  if(QMC_COMPLEX)
    set(FERMION_SRCS ${FERMION_SRCS} LCAO/LCAOSpinorBuilder.cpp)
  else(QMC_COMPLEX)
    #LCAO cusp correction is not ready for complex
    set(FERMION_SRCS ${FERMION_SRCS} LCAO/CuspCorrection.cpp LCAO/LCAOrbitalSetWithCorrection.cpp
                     LCAO/CuspCorrectionConstruction.cpp LCAO/SoaCuspCorrection.cpp)
  endif(QMC_COMPLEX)

  if(QMC_CUDA)
    set(WFSSRCS_CUDA detail/CUDA_legacy/determinant_update.cu detail/CUDA_legacy/delayed_update.cu
                     Fermion/DiracDeterminantCUDA.cpp Fermion/SlaterDetCUDA.cpp TrialWaveFunction_CUDA.cpp)
    set(WFSSRCS_CUDA
        ${WFSSRCS_CUDA} detail/CUDA_legacy/BsplineJastrowCuda.cu detail/CUDA_legacy/BsplineJastrowCudaPBC.cu
        Jastrow/OneBodyJastrowOrbitalBspline.cpp Jastrow/TwoBodyJastrowOrbitalBspline.cpp)
    if(HAVE_EINSPLINE)
      set(WFSSRCS_CUDA ${WFSSRCS_CUDA} EinsplineSetCuda.cpp detail/CUDA_legacy/PhaseFactors.cu
                       detail/CUDA_legacy/AtomicOrbitalCuda.cu)
    endif(HAVE_EINSPLINE)
  endif(QMC_CUDA)

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
    if(QMC_COMPLEX)
      set(FERMION_SRCS ${FERMION_SRCS} EinsplineSpinorSetBuilder.cpp BsplineFactory/SplineC2C.cpp)
      if(ENABLE_OFFLOAD)
        set(FERMION_SRCS ${FERMION_SRCS} BsplineFactory/SplineC2COMPTarget.cpp)
      endif(ENABLE_OFFLOAD)
    else(QMC_COMPLEX)
      set(FERMION_SRCS ${FERMION_SRCS} BsplineFactory/createRealSingle.cpp BsplineFactory/createRealDouble.cpp
                       BsplineFactory/SplineC2R.cpp BsplineFactory/SplineR2R.cpp)
      if(ENABLE_OFFLOAD)
        set(FERMION_SRCS ${FERMION_SRCS} BsplineFactory/SplineC2ROMPTarget.cpp)
      endif(ENABLE_OFFLOAD)
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
    Fermion/DiracDeterminantBatched.cpp
    Fermion/SlaterDet.cpp
    Fermion/SlaterDetBuilder.cpp
    Fermion/MultiSlaterDetTableMethod.cpp
    Fermion/MultiDiracDeterminant.cpp
    Fermion/MultiDiracDeterminant.2.cpp
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

if(ENABLE_CUDA)
  set(WFSSRCS_CUDA ${WFSSRCS_CUDA} detail/CUDA/cuBLAS_LU.cu detail/CUDA/delayed_update_helper.cu
                   detail/CUDA/matrix_update_helper.cu)
endif(ENABLE_CUDA)

####################################
# create libqmcwfs
####################################
if(USE_OBJECT_TARGET)
  add_library(qmcwfs OBJECT ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
else(USE_OBJECT_TARGET)
  add_library(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
endif(USE_OBJECT_TARGET)

if(QMC_CUDA OR ENABLE_CUDA)
  if(NOT QMC_CUDA2HIP)
    add_library(qmcwfs_cuda ${WFSSRCS_CUDA})
  else()
    hip_add_library(qmcwfs_cuda ${WFSSRCS_CUDA})
  endif(NOT QMC_CUDA2HIP)
  if(QMC_CUDA)
    target_link_libraries(qmcwfs_cuda PRIVATE qmcparticle)
  endif()
  target_link_libraries(qmcwfs_cuda PRIVATE platform_LA qmc_external_thrust)
  target_include_directories(qmcwfs_cuda PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(qmcwfs PRIVATE qmcwfs_cuda)
endif(QMC_CUDA OR ENABLE_CUDA)

target_include_directories(qmcwfs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcwfs PUBLIC qmcutil qmcparticle platform_runtime)
target_link_libraries(qmcwfs PRIVATE einspline platform_LA Math::FFTW3)

# do not merge with qmcwfs target unless all qmcwfs unit tests pass santitizers.
if(ENABLE_CUDA)
  set(DIRECT_INVERSION_SRCS detail/CUDA/cuBLAS_LU.cu)
  if(NOT QMC_CUDA2HIP)
    add_library(qmcwfs_direct_inversion_cuda "${DIRECT_INVERSION_SRCS}")
  else()
    hip_add_library(qmcwfs_direct_inversion_cuda "${DIRECT_INVERSION_SRCS}")
  endif()
  target_link_libraries(qmcwfs_direct_inversion_cuda PUBLIC platform_LA)
  # CUDA_ADD_LIBRARY is from deprecated FindCuda and ignores TARGET_INCLUDE_DIRECTORIES so
  # headers in this target contain full relative path from project root.
  target_link_libraries(qmcwfs PUBLIC qmcwfs_direct_inversion_cuda)
endif(ENABLE_CUDA)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
