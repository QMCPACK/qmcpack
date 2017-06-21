//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DummyQMC.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

/// Constructor.
DummyQMC::DummyQMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h)
{
  RootName = "dummy";
  QMCType ="dummy";
}

bool DummyQMC::run()
{
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = sqrt(Tau);
  std::cout << "Lattice of ParticleSet " << std::endl;
  W.Lattice.print(std::cout);
  //property container to hold temporary properties, such as a local energy
  //MCWalkerConfiguration::PropertyContainer_t Properties;
  MCWalkerConfiguration::Walker_t& thisWalker(**W.begin());
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandom(deltaR);
  //new poosition
  W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
  //apply boundary condition: put everything back to a unit cell
  W.applyBC(W.R);
  //update the distance table associated with W
  W.update();
  //evaluate wave function
  //update the properties: note that we are getting \f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
  ValueType logpsi(Psi.evaluateLog(W));
  RealType eloc=H.evaluate(W);
  return true;
}

bool
DummyQMC::put(xmlNodePtr q)
{
  //nothing to do
  return true;
}

}

