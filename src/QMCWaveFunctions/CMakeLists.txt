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
#//
#// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//////////////////////////////////////////////////////////////////////////////////////
    
    
SET(WFBASE_SRCS
  OrbitalBase.cpp
  DiffOrbitalBase.cpp
  OrbitalBuilderBase.cpp
  SPOSetBuilder.cpp
  SPOInfo.cpp
  SPOSetInfo.cpp
  SPOSetInputInfo.cpp
  SPOSetBase.cpp
  OptimizableSPOSet.cpp
  CompositeSPOSet.cpp
  HarmonicOscillator/SHOSet.cpp
  HarmonicOscillator/SHOSetBuilder.cpp
  )

IF(NOT QMC_COMPLEX)
  SET(WFBASE_SRCS ${WFBASE_SRCS}
    FDLRWfn.cpp
    )
ENDIF(NOT QMC_COMPLEX)


########################
# build jastrows
########################
#common jastrows
SET(JASTROW_SRCS
  Jastrow/LRTwoBodyJastrow.cpp
  Jastrow/PadeJastrowBuilder.cpp
  Jastrow/JastrowBuilder.cpp
  Jastrow/BsplineJastrowBuilder.cpp
  Jastrow/kSpaceJastrow.cpp
  Jastrow/kSpaceJastrowBuilder.cpp
  Jastrow/RPAJastrow.cpp
  Jastrow/singleRPAJastrowBuilder.cpp
  Jastrow/JABBuilder.cpp
  IonOrbital.cpp
  IonOrbitalBuilder.cpp
  OptimizableSPOBuilder.cpp
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
    )

  IF(NOT QMC_COMPLEX AND NOT ENABLE_SOA)
    SET(FERMION_SRCS ${FERMION_SRCS}
      MolecularOrbitals/STOBuilder.cpp
      MolecularOrbitals/GTOBuilder.cpp
      MolecularOrbitals/NGOBuilder.cpp
      )
  ENDIF(NOT QMC_COMPLEX AND NOT ENABLE_SOA)

  IF(NOT QMC_COMPLEX AND ENABLE_SOA)
    #lcao is deigned for localized orbitals in bulk and complex
    SET(FERMION_SRCS ${FERMION_SRCS}  lcao/LCAOrbitalSet.cpp lcao/LCAOrbitalBuilder.cpp)
  ENDIF(NOT QMC_COMPLEX AND ENABLE_SOA)

  IF(QMC_CUDA)
    SET(FERMION_SRCS ${FERMION_SRCS}
      Fermion/determinant_update.cu 
      Fermion/DiracDeterminantCUDA.cpp
      Fermion/SlaterDetCUDA.cpp
      TrialWaveFunction_CUDA.cpp
      )
    IF(ENABLE_SOA)
      SET(JASTROW_SRCS ${JASTROW_SRCS}
        Jastrow/BsplineJastrowCuda.cu
        Jastrow/BsplineJastrowCudaPBC.cu
        Jastrow/OneBodyJastrowOrbitalBspline.cpp
        Jastrow/TwoBodyJastrowOrbitalBspline.cpp
        )
    ELSE(ENABLE_SOA)
      SET(JASTROW_SRCS ${JASTROW_SRCS}
        Jastrow/BsplineJastrowCuda.cu
        Jastrow/BsplineJastrowCudaPBC.cu
        Jastrow/OneBodyJastrowOrbitalBsplineAoS.cpp
        Jastrow/TwoBodyJastrowOrbitalBsplineAoS.cpp
        )
    ENDIF(ENABLE_SOA)
    IF(HAVE_EINSPLINE)
      SET(FERMION_SRCS ${FERMION_SRCS}
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

    IF(NOT QMC_COMPLEX)
      SET(FERMION_SRCS ${FERMION_SRCS}
        BsplineFactory/createRealSingle.cpp
        BsplineFactory/createRealDouble.cpp
        )
    ENDIF()

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

  #only experimental version
  IF(QMC_BUILD_LEVEL GREATER 2)
    IF(NOT QMC_COMPLEX)
    SET(FERMION_SRCS ${FERMION_SRCS}
      AGPDeterminant.cpp AGPDeterminantBuilder.cpp
      )
    ENDIF(NOT QMC_COMPLEX)
  ENDIF(QMC_BUILD_LEVEL GREATER 2)

ENDIF(OHMMS_DIM MATCHES 3)

SET(FERMION_SRCS ${FERMION_SRCS}
  Fermion/DiracDeterminantBase.cpp
  Fermion/DiracDeterminantSoA.cpp
  Fermion/DiracDeterminantOpt.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  Fermion/MultiSlaterDeterminant.cpp
  Fermion/MultiSlaterDeterminantFast.cpp
  Fermion/MultiDiracDeterminantBase.cpp
  Fermion/MultiDiracDeterminantBase.2.cpp
  Fermion/BackflowBuilder.cpp 
  Fermion/DiracDeterminantWithBackflow.cpp
  Fermion/SlaterDetWithBackflow.cpp
  Fermion/MultiSlaterDeterminantWithBackflow.cpp
  SPOSetBuilderFactory.cpp
  TrialWaveFunction.cpp
  WaveFunctionFactory.cpp
  )

IF(NOT QMC_COMPLEX AND NOT ENABLE_SOA)
  SET(FERMION_SRCS ${FERMION_SRCS}
    Fermion/SlaterDetOpt.cpp
    )
ENDIF(NOT QMC_COMPLEX AND NOT ENABLE_SOA)

####################################
# create libqmcwfs
####################################
IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ELSE(QMC_CUDA)
  ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ENDIF(QMC_CUDA)
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
