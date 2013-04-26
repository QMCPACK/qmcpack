//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  cout << "Lattice of ParticleSet " << endl;
  W.Lattice.print(cout);
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

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
