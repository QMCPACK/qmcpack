#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
#//                    D.C. Yang, University of Illinois at Urbana-Champaign
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(HAMSRCS
    HamiltonianPool.cpp
    ACForce.cpp
    OperatorBase.cpp
    QMCHamiltonian.cpp
    Hdispatcher.cpp
    BareKineticEnergy.cpp
    CoulombPBCAA.cpp
    CoulombPBCAB.cpp
    EwaldRef.cpp
    NonLocalTOperator.cpp
    ForwardWalking.cpp
    PairCorrEstimator.cpp
    DensityEstimator.cpp
    SkPot.cpp
    SkEstimator.cpp
    SkAllEstimator.cpp
    MomentumEstimator.cpp
    ForceBase.cpp
    BareForce.cpp
    HamiltonianFactory.cpp
    CoulombPotentialFactory.cpp
    HarmonicExternalPotential.cpp
    GridExternalPotential.cpp
    StaticStructureFactor.cpp
    OrbitalImages.cpp
    SpinDensity.cpp
    SpeciesKineticEnergy.cpp
    LatticeDeviationEstimator.cpp
    SpaceWarpTransformation.cpp
    ObservableHelper.cpp)

if(OHMMS_DIM MATCHES 3)
  set(HAMSRCS
      ${HAMSRCS}
      ChiesaCorrection.cpp
      ForceCeperley.cpp
      ForceChiesaPBCAA.cpp
      StressPBC.cpp
      LocalECPotential.cpp
      NonLocalECPComponent.cpp
      NonLocalECPotential.deriv.cpp
      NonLocalECPotential.cpp
      L2Potential.cpp
      SOECPComponent.cpp
      SOECPotential.cpp
      ECPotentialBuilder.cpp
      ECPComponentBuilder.cpp
      ECPComponentBuilder.1.cpp
      ECPComponentBuilder.2.cpp
      ECPComponentBuilder_L2.cpp)

  if(NOT REMOVE_TRACEMANAGER)
    set(HAMSRCS ${HAMSRCS} ../Estimators/TraceManager.cpp ReferencePoints.cpp SpaceGrid.cpp EnergyDensityEstimator.cpp
                DensityMatrices1B.cpp)
  endif()

  if(HAVE_LIBFFTW)
    set(HAMSRCS ${HAMSRCS} MPC.cpp)
  endif(HAVE_LIBFFTW)

endif(OHMMS_DIM MATCHES 3)

if(USE_OBJECT_TARGET)
  add_library(qmcham OBJECT ${HAMSRCS})
else()
  add_library(qmcham ${HAMSRCS})
endif()

target_include_directories(qmcham PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(qmcham PUBLIC qmcwfs)

target_link_libraries(qmcham PRIVATE einspline platform_LA Math::FFTW3)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
