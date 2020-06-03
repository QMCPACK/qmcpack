#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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
#//
#// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//////////////////////////////////////////////////////////////////////////////////////
    
    
SET(WFBASE_SRCS
  WaveFunctionPool.cpp
  WaveFunctionComponent.cpp
  DiffWaveFunctionComponent.cpp
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
  ExampleHeComponent.cpp
  )

IF(NOT QMC_COMPLEX)
  SET(WFBASE_SRCS ${WFBASE_SRCS}
    RotatedSPOs.cpp
    FDLRWfn.cpp
    )
ENDIF(NOT QMC_COMPLEX)

IF(QMC_COMPLEX)
  SET(WFBASE_SRCS ${WFBASE_SRCS}
    SpinorSet.cpp
  )
ENDIF(QMC_COMPLEX)
########################
# build jastrows
########################
#common jastrows
SET(JASTROW_SRCS
  Jastrow/JastrowBuilder.cpp
  Jastrow/kSpaceJastrow.cpp
  Jastrow/kSpaceJastrowBuilder.cpp
  Jastrow/RadialJastrowBuilder.cpp
  Jastrow/CountingJastrowBuilder.cpp
  Jastrow/RPAJastrow.cpp
  LatticeGaussianProduct.cpp
  LatticeGaussianProductBuilder.cpp
  Fermion/SPOSetProxy.cpp
  Fermion/SPOSetProxyForMSD.cpp
  )

IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasComplexOrbitalBuilder.cpp
    )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasOrbitalBuilder.cpp
    )

ENDIF(QMC_COMPLEX)


# wavefunctions only availbale to 3-dim problems
IF(OHMMS_DIM MATCHES 3)

  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/eeI_JastrowBuilder.cpp
    Jastrow/CountingJastrowBuilder.cpp
    )

  SET(FERMION_SRCS ${FERMION_SRCS} lcao/LCAOrbitalSet.cpp lcao/LCAOrbitalBuilder.cpp)
  IF(NOT QMC_COMPLEX)
    #lcao cusp correction is not ready for complex
    SET(FERMION_SRCS ${FERMION_SRCS}
      lcao/CuspCorrection.cpp
      lcao/LCAOrbitalSetWithCorrection.cpp
      lcao/CuspCorrectionConstruction.cpp
      )
  ENDIF(NOT QMC_COMPLEX)

  IF(QMC_CUDA)
    SET(WFSSRCS_CUDA
      Fermion/determinant_update.cu
      Fermion/delayed_update.cu
      Fermion/DiracDeterminantCUDA.cpp
      Fermion/SlaterDetCUDA.cpp
      TrialWaveFunction_CUDA.cpp
      )
    SET(WFSSRCS_CUDA ${WFSSRCS_CUDA}
      Jastrow/BsplineJastrowCuda.cu
      Jastrow/BsplineJastrowCudaPBC.cu
      Jastrow/OneBodyJastrowOrbitalBspline.cpp
      Jastrow/TwoBodyJastrowOrbitalBspline.cpp
      )
    IF(HAVE_EINSPLINE)
      SET(WFSSRCS_CUDA ${WFSSRCS_CUDA}
        EinsplineSetCuda.cpp
        PhaseFactors.cu
        AtomicOrbitalCuda.cu
        )
    ENDIF(HAVE_EINSPLINE)
  ENDIF(QMC_CUDA)

  IF(HAVE_EINSPLINE)
    IF(NOT MIXED_PRECISION)
      SET(FERMION_SRCS ${FERMION_SRCS} EinsplineSet.cpp)
    ENDIF(NOT MIXED_PRECISION)
    SET(FERMION_SRCS ${FERMION_SRCS}
      EinsplineSetBuilderCommon.cpp
      EinsplineSetBuilderOld.cpp
      MuffinTin.cpp
      AtomicOrbital.cpp
      EinsplineSetBuilderReadBands_ESHDF.cpp
      EinsplineSetBuilderESHDF.fft.cpp
      EinsplineSetBuilder_createSPOs.cpp
      BsplineFactory/createComplexDouble.cpp
      BsplineFactory/createComplexSingle.cpp
      BandInfo.cpp
      BsplineFactory/BsplineReaderBase.cpp
      )
    IF(QMC_COMPLEX)
      SET(FERMION_SRCS ${FERMION_SRCS}
        EinsplineSpinorSetBuilder.cpp
        BsplineFactory/SplineC2C.cpp
      )
    ENDIF(QMC_COMPLEX)
    IF(NOT QMC_COMPLEX)
      SET(FERMION_SRCS ${FERMION_SRCS}
        BsplineFactory/createRealSingle.cpp
        BsplineFactory/createRealDouble.cpp
        BsplineFactory/SplineC2R.cpp
        BsplineFactory/SplineR2R.cpp
        )
      IF(ENABLE_OFFLOAD)
        SET(FERMION_SRCS ${FERMION_SRCS} BsplineFactory/SplineC2ROMP.cpp)
      ENDIF(ENABLE_OFFLOAD)
    ENDIF(NOT QMC_COMPLEX)

  ENDIF(HAVE_EINSPLINE)

  # plane wave SPO
  SET(FERMION_SRCS ${FERMION_SRCS}
    PlaneWave/PWBasis.cpp
    PlaneWave/PWParameterSet.cpp
    PlaneWave/PWOrbitalBuilder.cpp
    )
  IF(QMC_COMPLEX)
    SET(FERMION_SRCS ${FERMION_SRCS}
      PlaneWave/PWOrbitalSet.cpp
    )
  ELSE()
    SET(FERMION_SRCS ${FERMION_SRCS}
      PlaneWave/PWRealOrbitalSet.cpp
    )
  ENDIF(QMC_COMPLEX)

  IF(NOT QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    AGPDeterminant.cpp AGPDeterminantBuilder.cpp
    )
  ENDIF(NOT QMC_COMPLEX)

ENDIF(OHMMS_DIM MATCHES 3)

SET(FERMION_SRCS ${FERMION_SRCS}
  Fermion/DiracDeterminant.cpp
  Fermion/DiracDeterminantBatched.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  Fermion/MultiSlaterDeterminant.cpp
  Fermion/MultiSlaterDeterminantFast.cpp
  Fermion/MultiDiracDeterminant.cpp
  Fermion/MultiDiracDeterminant.2.cpp
  Fermion/BackflowBuilder.cpp 
  Fermion/DiracDeterminantWithBackflow.cpp
  Fermion/SlaterDetWithBackflow.cpp
  Fermion/MultiSlaterDeterminantWithBackflow.cpp
  SPOSetBuilderFactory.cpp
  TrialWaveFunction.cpp
  WaveFunctionFactory.cpp
  )

IF(ENABLE_CUDA)
  SET(WFSSRCS_CUDA ${WFSSRCS_CUDA}
      Fermion/delayed_update_helper.cu
      Fermion/matrix_update_helper.cu
      )
ENDIF(ENABLE_CUDA)

####################################
# create libqmcwfs
####################################
IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcwfs OBJECT ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ELSE(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ENDIF(USE_OBJECT_TARGET)

IF(QMC_CUDA OR ENABLE_CUDA)
  CUDA_ADD_LIBRARY(qmcwfs_cuda ${WFSSRCS_CUDA})
  CUDA_ADD_CUBLAS_TO_TARGET(qmcwfs_cuda)
  TARGET_LINK_LIBRARIES(qmcwfs_cuda PUBLIC qmcparticle)
  TARGET_LINK_LIBRARIES(qmcwfs PRIVATE qmcwfs_cuda ${CUDA_cusolver_LIBRARY})
ENDIF(QMC_CUDA OR ENABLE_CUDA)

TARGET_LINK_LIBRARIES(qmcwfs PUBLIC qmcparticle platform_device)
TARGET_LINK_LIBRARIES(qmcwfs PRIVATE einspline Math::BLAS_LAPACK Math::FFTW3 Math::scalar_vector_functions)

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
