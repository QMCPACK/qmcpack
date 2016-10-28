aos branch
----------


Main changes in DistanceTableData instantiation
  Particle/DistanceTable.cpp
  Particle/DistanceTable.h
  Particle/DistanceTableData.h
  Particle/ParticleSet.cpp
  Particle/ParticleSet.h
  ParticleBase/ParticleBase.cpp
  ParticleBase/ParticleBase.h

Important to note that DistanceTableData supports both AoS and SoA data
structures. Once the major classes in QMCWaveFunctions and QMCHamiltonians
adopt the SoA structure and there is no conflict between the two data types,
DistanceTableData has to be modified to eliminate inefficient code sections.

I'm not aware of any use of these in the performance critical sections and
advise everyone to use them only during the development stage.

\code
  /// return the distance |R[iadj(i,nj)]-R[i]|
  inline RealType distance(int i, int nj) const
  {
    return (DTType)? r_m2(i,nj): r_m[M[i]+nj];
  }

  /// return the displacement R[iadj(i,nj)]-R[i]
  inline PosType displacement(int i, int nj) const
  {
    return (DTType)? dr_m2(i,nj): dr_m[M[i]+nj];
  }

  //!< Returns a number of neighbors of the i-th ptcl.
  inline IndexType nadj(int i) const
  {
    return (DTType)? M[i]:M[i+1]-M[i];
  }

  //!< Returns the id of nj-th neighbor for i-th ptcl
  inline IndexType iadj(int i, int nj) const
  {
    return (DTType)? J2(i,nj):J[M[i] +nj];
  }
\endcode

These files are modifed to use DT_AOS. These classes need to be modified to use DT_SOA later.

  Particle/tests/test_distance_table.cpp
  Particle/tests/test_particle.cpp
  QMCApp/InitMolecularSystem.cpp
  QMCDrivers/WaveFunctionTester.cpp
  QMCDrivers/tests/test_vmc.cpp
  QMCTools/MolecularOrbitalBasis.h

  QMCHamiltonians/GaussianPot.h
  QMCHamiltonians/HFDHE2Potential.h
  QMCHamiltonians/HFDHE2Potential_tail.h
  QMCHamiltonians/HFDHE2_Moroni1995.cpp
  QMCHamiltonians/HardSphere.h
  QMCHamiltonians/HeEPotential.h
  QMCHamiltonians/HePressure.h
  QMCHamiltonians/HeSAPT_smoothed.cpp
  QMCHamiltonians/HusePot.h
  QMCHamiltonians/JelliumPotential.h
  QMCHamiltonians/LennardJones_smoothed.cpp
  QMCHamiltonians/LocalCorePolPotential.cpp
  QMCHamiltonians/LocalECPotential.cpp
  QMCHamiltonians/LocalMomentEstimator.cpp
  QMCHamiltonians/MPC.cpp
  QMCHamiltonians/ModPosTelPot.h
  QMCHamiltonians/NonLocalECPotential.cpp
  QMCHamiltonians/NumberFluctuations.h
  QMCHamiltonians/NumericalRadialPotential.cpp
  QMCHamiltonians/OscillatoryPot.h
  QMCHamiltonians/PulayForce.cpp
  QMCHamiltonians/StressPBC.cpp
  QMCHamiltonians/ZeroVarianceForce.cpp
  QMCHamiltonians/tests/test_bare_kinetic.cpp
  QMCHamiltonians/tests/test_coulomb_pbcAB.cpp

  QMCWaveFunctions/EinsplineSetBuilder_createSPOs.cpp
  QMCWaveFunctions/Experimental/CorrectingOrbitalBasisSet.h
  QMCWaveFunctions/Fermion/BackflowTransformation.h
  QMCWaveFunctions/Fermion/Backflow_eI.h
  QMCWaveFunctions/Fermion/Backflow_eI_spin.h
  QMCWaveFunctions/Fermion/Backflow_ee.h
  QMCWaveFunctions/IonOrbital.cpp
  QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h
  QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h
  QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h
  QMCWaveFunctions/Jastrow/OneBodySpinJastrowOrbital.h
  QMCWaveFunctions/Jastrow/ThreeBodyBlockSparse.cpp
  QMCWaveFunctions/Jastrow/ThreeBodyGeminal.cpp
  QMCWaveFunctions/Jastrow/eeI_JastrowOrbital.h
  QMCWaveFunctions/LocalizedBasisSet.h
  QMCWaveFunctions/SparseLocalizedBasisSet.h
  QMCWaveFunctions/tests/test_bspline_jastrow.cpp
  QMCWaveFunctions/tests/test_einset.cpp
  QMCWaveFunctions/tests/test_pw.cpp
  QMCWaveFunctions/tests/test_wf.cpp
 
