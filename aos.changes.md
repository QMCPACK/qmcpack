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

These files are modifed to use DT_AOS. These classes need to be modified to use DT_SOA later.

  Particle/tests/test_distance_table.cpp
  Particle/tests/test_particle.cpp
  QMCApp/InitMolecularSystem.cpp
  QMCDrivers/WaveFunctionTester.cpp
  QMCDrivers/tests/test_vmc.cpp

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
