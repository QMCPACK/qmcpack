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
    
    
#include "QMCDrivers/VMC.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

/// Constructor.
VMC::VMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h)
{
  RootName = "vmc";
  QMCType ="VMC";
}

/** Run the VMC algorithm.
 *
 * Advance the walkers nblocks*nsteps timesteps.
 * For each timestep:
 * <ul>
 * <li> Move all the particles of a walker.
 * <li> Calculate the properties for the new walker configuration.
 * <li> Accept/reject the new configuration.
 * <li> Accumulate the estimators.
 * </ul>
 * For each block:
 * <ul>
 * <li> Flush the estimators and print to file.
 * <li> (Optional) Print the ensemble of walker configurations.
 * </ul>
 *
 * Default mode: Print the ensemble of walker configurations
 * at the end of the run.
 */
bool VMC::run()
{
  Estimators->reportHeader(AppendRun);
  Estimators->reset();
  IndexType block = 0;
  double wh=0.0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  do
  {
    IndexType step = 0;
    nAccept = 0;
    nReject=0;
    //start a block
    Estimators->startBlock();
    do
    {
      advanceWalkerByWalker();
      //cerr << "HEY!!  I'M NOT MOVING ANYBODY; JUST FOR TESTING. -FROM VMC.CPP" << std::endl;
      //cerr << "going to evaluate Psi.  it is ";
      //RealType logpsi(Psi.evaluateLog(W));
      //cerr << logpsi << std::endl;
      //cerr << "this is for coords" << std::endl;
      //for (int w=0; w<W.R.size(); w++)
      //	cerr << w << ": " << W.R[w] << std::endl;
      step++;
      CurrentStep++;
      Estimators->accumulate(W);
    }
    while(step<nSteps);
    //stop a block, pass the acceptance rate of this block
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    block++;
    nAccept = 0;
    nReject = 0;
    recordBlock(block);
  }
  while(block<nBlocks);
  //finalize a qmc section
  return finalize(block);
}

/**  Advance all the walkers one timstep.
 *
 Propose a move for each walker from its old
 position \f${\bf R'}\f$ to a new position \f${\bf R}\f$
 \f[
 {\bf R'} + {\bf \chi} +
 \tau {\bf v_{drift}}({\bf R'}) =  {\bf R},
 \f]
 where \f$ {\bf \chi} \f$ is a 3N-diminsional Gaussian
 of mean zero and variance \f$ \tau \f$ and
 \f$ {\bf v_{drift}} \f$ is the drift velocity
 \f[
 {\bf v_{drift}}({\bf R'}) = {\bf \nabla}
 \ln |\Psi_T({\bf R'})| = \Psi_T({\bf R'})^{-1}
 {\bf \nabla} \Psi_T({\bf R'}).
 \f]
 Metropolis accept/reject with probability
 \f[
 P_{accept}(\mathbf{R'}\rightarrow\mathbf{R}) =
 \min\left[1,\frac{G(\mathbf{R}\rightarrow\mathbf{R'})
 \Psi_T(\mathbf{R})^2}{G(\mathbf{R'}\rightarrow\mathbf{R})
 \Psi_T(\mathbf{R'})^2}\right],
 \f]
 where \f$ G \f$ is the drift-diffusion Green's function
 \f[
 G(\mathbf{R'} \rightarrow
 \mathbf{R}) = (2\pi\tau)^{-3/2}\exp \left[ -
 (\mathbf{R}-\mathbf{R'}-\tau \mathbf{v_{drift}}
 (\mathbf{R'}))^2/2\tau \right].
 \f]
 If the move is accepted, update the walker configuration and
 properties.  For rejected moves, do not update except for the
 Age which needs to be incremented by one.
 */
void
VMC::advanceWalkerByWalker()
{
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = sqrt(Tau);
  //property container to hold temporary properties, such as a local energy
  //MCWalkerConfiguration::PropertyContainer_t Properties;
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  while(it != it_end)
  {
    MCWalkerConfiguration::Walker_t& thisWalker(**it);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(deltaR);
    W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
    //update the distance table associated with W
    //DistanceTable::update(W);
    W.update();
    //evaluate wave function
    //update the properties: note that we are getting \f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
    RealType logpsi(Psi.evaluateLog(W));
    //deltaR = W.R - (*it)->R - (*it)->Drift;
    //RealType forwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
    //RealType forwardGF = exp(-0.5*Dot(deltaR,deltaR));
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    //converting gradients W.G to drifts, D = tau*G (reuse G)
    //RealType scale = getDriftScale(Tau,W.G);
    //drift = scale*W.G;
    setScaledDrift(Tau,W.G,drift);
    //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
    deltaR = thisWalker.R - W.R - drift;
    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
    RealType g= exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if(Random() > g)
    {
      thisWalker.Age++;
      ++nReject;
    }
    else
    {
      //thisWalker.Age=0;
      //Properties(LOCALENERGY) = H.evaluate(W);
      //Properties(LOCALPOTENTIAL) = H.getLocalPotential();
      RealType eloc=H.evaluate(W);
      thisWalker.R = W.R;
      thisWalker.Drift = drift;
      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.saveProperty(thisWalker.getPropertyBase());
      ++nAccept;
    }
    ++it;
  }
}

bool
VMC::put(xmlNodePtr q)
{
  //nothing to do
  return true;
}

}

