//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "QMCDrivers/DriftOperators.h"
#include "LongRange/StructFact.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/SymmetryOperations.h"
#include "Numerics/Blasf.h"
#include <sstream>

namespace qmcplusplus
{


WaveFunctionTester::WaveFunctionTester(MCWalkerConfiguration& w,
                                       TrialWaveFunction& psi,
                                       QMCHamiltonian& h,
                                       ParticleSetPool &ptclPool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),checkRatio("no"),checkClone("no"), checkHamPbyP("no"),
  PtclPool(ptclPool), wftricks("no"),checkEloc("no"), checkBasic("yes"), checkRatioV("no"),
  deltaParam(0.0), toleranceParam(0.0), outputDeltaVsError(false), checkSlaterDet(true)
{
  m_param.add(checkRatio,"ratio","string");
  m_param.add(checkClone,"clone","string");
  m_param.add(checkHamPbyP,"hamiltonianpbyp","string");
  m_param.add(sourceName,"source","string");
  m_param.add(wftricks,"orbitalutility","string");
  m_param.add(checkEloc,"printEloc","string");
  m_param.add(checkBasic,"basic","string");
  m_param.add(checkRatioV,"virtual_move","string");
  m_param.add(deltaParam,"delta","none");
  m_param.add(toleranceParam,"tolerance","none");
  m_param.add(checkSlaterDetOption,"sd","string");

  deltaR.resize(w.getTotalNum());
  makeGaussRandom(deltaR);
}

WaveFunctionTester::~WaveFunctionTester()
{
}

/*!
 * \brief Test the evaluation of the wavefunction, gradient and laplacian
 by comparing to the numerical evaluation.
 *
 Use the finite difference formulas formulas
 \f[
 \nabla_i f({\bf R}) = \frac{f({\bf R+\Delta r_i}) - f({\bf R})}{2\Delta r_i}
 \f]
 and
 \f[
 \nabla_i^2 f({\bf R}) = \sum_{x,y,z} \frac{f({\bf R}+\Delta x_i)
 - 2 f({\bf R}) + f({\bf R}-\Delta x_i)}{2\Delta x_i^2},
 \f]
 where \f$ f = \ln \Psi \f$ and \f$ \Delta r_i \f$ is a
 small displacement for the ith particle.
*/

bool
WaveFunctionTester::run()
{
  //DistanceTable::create(1);
  char fname[16];
  sprintf(fname,"wftest.%03d",OHMMS::Controller->rank());
  fout.open(fname);
  fout.precision(15);

  app_log() << "Starting a Wavefunction tester.  Additional information in "  << fname << std::endl;

  put(qmcNode);
  if (checkSlaterDetOption=="no") checkSlaterDet = false;
  if (checkRatio == "yes")
  {
    //runRatioTest();
    runRatioTest2();
  }
  else if (checkClone == "yes")
    runCloneTest();
  else if(checkEloc != "no")
    printEloc();
  else if (sourceName.size() != 0)
  {
    runGradSourceTest();
    runZeroVarianceTest();
  }
  else if (checkRatio =="deriv")
  {
    makeGaussRandom(deltaR);
    deltaR *=0.2;
    runDerivTest();
    runDerivNLPPTest();
  }
  else if (checkRatio =="derivclone")
    runDerivCloneTest();
  else if (wftricks =="rotate")
    runwftricks();
  else if (wftricks =="plot")
    runNodePlot();
  else if (checkBasic == "yes")
    runBasicTest();
  else if (checkRatioV == "yes")
    runRatioV();
  else
    app_log() << "No wavefunction test specified" << std::endl;

  //RealType ene = H.evaluate(W);
  //app_log() << " Energy " << ene << std::endl;
  return true;
}

void WaveFunctionTester::runCloneTest()
{
  for (int iter=0; iter<4; ++iter)
  {
    app_log() << "Clone" << iter << std::endl;
    ParticleSet* w_clone = new MCWalkerConfiguration(W);
    TrialWaveFunction *psi_clone = Psi.makeClone(*w_clone);
    QMCHamiltonian *h_clone = H.makeClone(*w_clone,*psi_clone);
    h_clone->setPrimary(false);
    IndexType nskipped = 0;
    RealType sig2Enloc=0, sig2Drift=0;
    RealType delta = 0.00001;
    RealType delta2 = 2*delta;
    ValueType c1 = 1.0/delta/2.0;
    ValueType c2 = 1.0/delta/delta;
    int nat = W.getTotalNum();
    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
    //copy the properties of the working walker
    Properties = awalker->Properties;
    W.R = awalker->R;
    W.update();
    ValueType logpsi1 = Psi.evaluateLog(W);
    RealType eloc1  = H.evaluate(W);
    w_clone->R=awalker->R;
    w_clone->update();
    ValueType logpsi2 = psi_clone->evaluateLog(*w_clone);
    RealType eloc2  = h_clone->evaluate(*w_clone);
    app_log() << "Testing walker-by-walker functions " << std::endl;
    app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << std::endl;
    app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << std::endl;
    app_log() << "Testing pbyp functions " << std::endl;
    Walker_t::WFBuffer_t &wbuffer(awalker->DataSet);
    wbuffer.clear();
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << std::endl;
    Psi.registerData(W,wbuffer);
    wbuffer.allocate();
    Psi.copyFromBuffer(W,wbuffer);
    Psi.evaluateLog(W);
    logpsi1 = Psi.updateBuffer(W,wbuffer,false);
    eloc1= H.evaluate(W);
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << std::endl;
    wbuffer.clear();
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << std::endl;
    psi_clone->registerData(W,wbuffer);
    wbuffer.allocate();
    Psi.copyFromBuffer(W,wbuffer);
    Psi.evaluateLog(W);
    logpsi2 = Psi.updateBuffer(W,wbuffer,false);
    eloc2= H.evaluate(*w_clone);
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << std::endl;
    app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << std::endl;
    app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << std::endl;
    delete h_clone;
    delete psi_clone;
    delete w_clone;
  }
}

void WaveFunctionTester::printEloc()
{
  ParticleSetPool::PoolType::iterator p;
  for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
    app_log() << "ParticelSet = " << p->first << std::endl;
  // Find source ParticleSet
  ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
  if(pit == PtclPool.getPool().end())
     APP_ABORT("Unknown source \"" + sourceName + "\" in printEloc in WaveFunctionTester.");
  ParticleSet& source = *((*pit).second);
  app_log() << "Source = " <<sourceName <<"  " <<(*pit).first << std::endl;
  int nel = W.getTotalNum();
  int ncenter = source.getTotalNum();
//    int ncenter = 3;
//    std::cout <<"number of centers: " <<source.getTotalNum() << std::endl;
//    std::cout <<"0: " <<source.R[0] << std::endl;
//    std::cout <<"1: " <<source.R[1] << std::endl;
//    std::cout <<"2: " <<source.R[2] << std::endl;
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
  //copy the properties of the working walker
  Properties = awalker->Properties;
  W.R = awalker->R;
  W.update();
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);
  app_log() << "  Logpsi: " <<logpsi  << std::endl;
  app_log() << "  HamTest " << "  Total " <<  eloc << std::endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  //RealType psi = Psi.evaluateLog(W);
  //int iat=0;
  double maxR = 1000000.0;
  std::vector<int> closestElectron(ncenter);
  for(int iat=0; iat<ncenter; iat++)
  {
    maxR=10000000;
    for(int k=0; k<nel; k++)
    {
      double dx = std::sqrt( (W.R[k][0]-source.R[iat][0])*(W.R[k][0]-source.R[iat][0])
                             +(W.R[k][1]-source.R[iat][1])*(W.R[k][1]-source.R[iat][1])
                             +(W.R[k][2]-source.R[iat][2])*(W.R[k][2]-source.R[iat][2]));
      if(dx < maxR)
      {
        maxR = dx;
        closestElectron[iat]=k;
      }
    }
  }
//    closestElectron[iat]=1;
  std::ofstream out("eloc.dat");
  double x,dx=1.0/499.0;
  for (int k=0; k<500; k++)
  {
    x=-0.5+k*dx;
    out<<x <<"  ";
    for(int iat=0; iat<ncenter; iat++)
    {
      PosType tempR = W.R[closestElectron[iat]];
      W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;
      W.R[closestElectron[iat]][0] += x;
      W.update();
      ValueType logpsi_p = Psi.evaluateLog(W);
      ValueType ene = H.evaluate(W);
      out<<ene <<"  ";
      W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;
      W.R[closestElectron[iat]][1] += x;
      W.update();
      logpsi_p = Psi.evaluateLog(W);
      ene = H.evaluate(W);
      out<<ene <<"  ";
      W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;
      W.R[closestElectron[iat]][2] += x;
      W.update();
      logpsi_p = Psi.evaluateLog(W);
      ene = H.evaluate(W);
      out<<ene <<"  ";
      W.R[closestElectron[iat]] = tempR;
    }
    out<< std::endl;
  }
  out.close();
}


class FiniteDifference : public QMCTraits
{
public:
  enum FiniteDiffType {
    FiniteDiff_LowOrder,  // use simplest low-order formulas
    FiniteDiff_Richardson  // use Richardson extrapolation
  };
  FiniteDifference(FiniteDiffType fd_type=FiniteDiff_Richardson) : m_RichardsonSize(10), m_fd_type(fd_type) {}

  int m_RichardsonSize;


  FiniteDiffType m_fd_type;

  struct PositionChange {
    int index; // particle index
    PosType r;
  };
  typedef std::vector<PositionChange> PosChangeVector;
  typedef std::vector<ValueType> ValueVector;


  /** Generate points to evaluate */
  void finiteDifferencePoints(RealType delta, MCWalkerConfiguration& W,
                              PosChangeVector &positions);

  /** Compute finite difference after log psi is computed for each point */
  void computeFiniteDiff(RealType delta,
                         PosChangeVector &positions,
                         ValueVector &values,
                         ParticleSet::ParticleGradient_t &G_fd,
                         ParticleSet::ParticleLaplacian_t &L_fd);

  void computeFiniteDiffLowOrder(RealType delta,
                         PosChangeVector &positions,
                         ValueVector &values,
                         ParticleSet::ParticleGradient_t &G_fd,
                         ParticleSet::ParticleLaplacian_t &L_fd);

  void computeFiniteDiffRichardson(RealType delta,
                         PosChangeVector &positions,
                         ValueVector &values,
                         ParticleSet::ParticleGradient_t &G_fd,
                         ParticleSet::ParticleLaplacian_t &L_fd);
};


void FiniteDifference::finiteDifferencePoints(RealType delta, MCWalkerConfiguration& W,
                                              PosChangeVector &positions)
{
  // First position is the central point
  PositionChange p;
  p.index = 0;
  p.r = W.R[0];
  positions.push_back(p);

  int nat = W.getTotalNum();
  for (int iat=0; iat<nat; iat++)
  {
    PositionChange p;
    p.index = iat;
    PosType r0 = W.R[iat];

    for (int idim=0; idim<OHMMS_DIM; idim++)
    {
      p.r = r0;
      p.r[idim] = r0[idim] - delta;
      positions.push_back(p);

      p.r = r0;
      p.r[idim] = r0[idim] + delta;
      positions.push_back(p);

      if (m_fd_type == FiniteDiff_Richardson)
      {
        RealType dd = delta/2;
        for (int nr = 0; nr < m_RichardsonSize; nr++)
        {
          p.r = r0;
          p.r[idim] = r0[idim] - dd;
          positions.push_back(p);

          p.r = r0;
          p.r[idim] = r0[idim] + dd;
          positions.push_back(p);

          dd = dd/2;
        }
      }
    }
  }
}


void FiniteDifference::computeFiniteDiff(RealType delta,
                                   PosChangeVector &positions,
                                   ValueVector &values,
                                   ParticleSet::ParticleGradient_t &G_fd,
                                   ParticleSet::ParticleLaplacian_t &L_fd)
{
    assert(positions.size() == values.size());
    if (positions.size() == 0)
        return;

    ValueType logpsi = values[0];

    if (m_fd_type == FiniteDiff_LowOrder)
    {
      computeFiniteDiffLowOrder(delta, positions, values, G_fd, L_fd);
    }
    else if (m_fd_type == FiniteDiff_Richardson)
    {
      computeFiniteDiffRichardson(delta, positions, values, G_fd, L_fd);
    }
}

void FiniteDifference::computeFiniteDiffLowOrder(RealType delta,
                                   PosChangeVector &positions,
                                   ValueVector &values,
                                   ParticleSet::ParticleGradient_t &G_fd,
                                   ParticleSet::ParticleLaplacian_t &L_fd)
{
    ValueType logpsi = values[0];

    // lowest order derivative formula
    ValueType c1 = 1.0/delta/2.0;
    ValueType c2 = 1.0/delta/delta;

    const RealType twoD(2*OHMMS_DIM);
    const int pt_per_deriv = 2; // number of points per derivative
    for (int pt_i = 1; pt_i < values.size(); pt_i += pt_per_deriv*OHMMS_DIM)
    {
      GradType g0;
      ValueType lap0 = 0.0;
      for (int idim=0; idim<OHMMS_DIM; idim++)
      {
        int idx = pt_i + idim*pt_per_deriv;
        ValueType logpsi_m = values[idx];
        ValueType logpsi_p = values[idx+1];

        g0[idim] = logpsi_p - logpsi_m;
        lap0 += logpsi_p + logpsi_m;
      }

      int iat = positions[pt_i].index;
      GradType g = c1*g0;
      G_fd[iat] = g;

      ValueType lap = c2*(lap0 - twoD*logpsi);
      //ValueType lap = c2*(lap0 - 2.0*OHMMS_DIM*logpsi);
      L_fd[iat] = lap;
    }
}

// Use Richardson extrapolation to compute the derivatives.
// The 'delta' parameter should not be small, as with fixed
// order finite difference methods.  The algorithm  will zoom
// in on the right size of parameter.
void FiniteDifference::computeFiniteDiffRichardson(RealType delta,
                                           PosChangeVector &positions,
                                           ValueVector &values,
                                           ParticleSet::ParticleGradient_t &G_fd,
                                           ParticleSet::ParticleLaplacian_t &L_fd)
{
  RealType tol = 1e-7;
  ValueType logpsi = values[0];

  const int pt_per_deriv = 2*(m_RichardsonSize+1); // number of points per derivative
  for (int pt_i = 1; pt_i < values.size(); pt_i += pt_per_deriv*OHMMS_DIM)
  {
    GradType g0;
    GradType gmin;
    ValueType lmin;
    std::vector<GradType> g_base(m_RichardsonSize+1);
    std::vector<GradType> g_rich(m_RichardsonSize+1);
    std::vector<GradType> g_prev(m_RichardsonSize+1);

    std::vector<ValueType> l_base(m_RichardsonSize+1);
    std::vector<ValueType> l_rich(m_RichardsonSize+1);
    std::vector<ValueType> l_prev(m_RichardsonSize+1);

    // Initial gradients and Laplacians at different deltas.
    RealType dd = delta;
    const ValueType ctwo(2);
    for (int inr = 0; inr < m_RichardsonSize+1; inr++)
    {
      RealType twodd = 2*dd;
      RealType ddsq = dd*dd;
        l_base[inr] = 0.0;
        for (int idim=0; idim<OHMMS_DIM; idim++)
        {
           int idx = pt_i + idim*pt_per_deriv + 2*inr;
           ValueType logpsi_m = values[idx];
           ValueType logpsi_p = values[idx + 1];
           g_base[inr][idim] = (logpsi_p - logpsi_m)/twodd;
           l_base[inr] += (logpsi_p + logpsi_m - ctwo*logpsi)/ddsq;
           //g_base[inr][idim] = (logpsi_p - logpsi_m)/dd/2.0;
           //l_base[inr] += (logpsi_p + logpsi_m - 2.0*logpsi)/(dd*dd);
        }
        dd = dd/2;
    }

    // Gradient

    g_prev[0] = g_base[0];
    RealType fac = 1;
    bool found = false;
    for (int inr = 1; inr < m_RichardsonSize+1; inr++)
    {
        g_rich[0] = g_base[inr];

        fac *= 4;
        for (int j = 1; j < inr+1; j++) {
          g_rich[j] = g_rich[j-1] + ( g_rich[j-1] - g_prev[j-1])/(fac-1);
        }

        RealType err1 = 0.0;
        RealType norm = 0.0;
        for (int idim=0; idim<OHMMS_DIM; idim++)
        {
            err1 += std::abs(g_rich[inr][idim] - g_prev[inr-1][idim]);
            norm += std::abs(g_prev[inr-1][idim]);
        }

        RealType err_rel = err1/norm;

        // Not sure about the best stopping criteria
        if (err_rel < tol)
        {
          gmin = g_rich[inr];
          found = true;
          break;
        }
        g_prev = g_rich;
    }

    if (!found)
    {
      gmin = g_rich[m_RichardsonSize];
    }


    // Laplacian
    // TODO: eliminate the copied code between the gradient and Laplacian
    //       computations.

    l_prev[0] = l_base[0];

    fac = 1;
    found = false;
    for (int inr = 1; inr < m_RichardsonSize+1; inr++)
    {
        l_rich[0] = l_base[inr];

        fac *= 4;
        for (int j = 1; j < inr+1; j++) {
          l_rich[j] = l_rich[j-1] + ( l_rich[j-1] - l_prev[j-1])/(fac-1);
        }

        RealType err1 = std::abs(l_rich[inr] - l_prev[inr-1]);
        RealType err_rel = std::abs(err1/l_prev[inr-1]);

        if (err_rel < tol)
        {
          lmin = l_rich[inr];
          found = true;
          break;
        }
        l_prev = l_rich;
    }

    if (!found)
    {
      lmin = l_rich[m_RichardsonSize];
    }

    int iat = positions[pt_i].index;
    G_fd[iat] = gmin;
    L_fd[iat] = lmin;
  }
}

// Compute numerical gradient and Laplacian
void WaveFunctionTester::computeNumericalGrad(RealType delta,
                                              ParticleSet::ParticleGradient_t &G_fd, // finite difference
                                              ParticleSet::ParticleLaplacian_t &L_fd)
{
  FiniteDifference fd(FiniteDifference::FiniteDiff_LowOrder);
  //FiniteDifference fd(FiniteDifference::FiniteDiff_Richardson);
  FiniteDifference::PosChangeVector positions;

  fd.finiteDifferencePoints(delta, W, positions);

  FiniteDifference::ValueVector logpsi_vals;
  FiniteDifference::PosChangeVector::iterator it;

  for (it = positions.begin(); it != positions.end(); it++)
  {
    PosType r0 = W.R[it->index];
    W.R[it->index] = it->r;
    W.update();
    RealType logpsi0 = Psi.evaluateLog(W);
    RealType phase0 = Psi.getPhase();
#if defined(QMC_COMPLEX)
    ValueType logpsi = std::complex<OHMMS_PRECISION>(logpsi0,phase0);
#else
    ValueType logpsi = logpsi0;
#endif
    logpsi_vals.push_back(logpsi);

    W.R[it->index] = r0;
    W.update();
    Psi.evaluateLog(W);
  }

  fd.computeFiniteDiff(delta, positions, logpsi_vals, G_fd, L_fd);
}

// Usually
// lower_iat = 0
// upper_iat = nat
bool WaveFunctionTester::checkGradients(int lower_iat, int upper_iat,
                                           ParticleSet::ParticleGradient_t &G,
                                           ParticleSet::ParticleLaplacian_t &L,
                                           ParticleSet::ParticleGradient_t &G_fd,
                                           ParticleSet::ParticleLaplacian_t &L_fd,
                                           std::stringstream &log,
                                           int indent /* = 0 */)
{

  ParticleSet::Scalar_t rel_tol = 1e-3;
  ParticleSet::Scalar_t abs_tol = 1e-7;
  if (toleranceParam > 0.0)
  {
    rel_tol = toleranceParam;
  } 

  bool all_okay = true;
  std::string pad(4*indent, ' ');

  
  for (int iat=lower_iat; iat<upper_iat; iat++)
  {
    ParticleSet::Scalar_t L_err = std::abs(L[iat]-L_fd[iat]);
    ParticleSet::Scalar_t L_rel_denom = std::max( std::abs(L[iat]), std::abs(L_fd[iat]) );
    ParticleSet::Scalar_t L_err_rel = std::abs( L_err / L_rel_denom );
    
    if (L_err_rel > rel_tol && L_err > abs_tol)
    {
      if (L_err_rel > rel_tol)
      {
        log << pad << "Finite difference Laplacian exceeds relative tolerance (" << rel_tol << ") for particle " << iat << std::endl;
      } 
      else
      {
        log << pad << "Finite difference Laplacian exceeds absolute tolerance (" << abs_tol << ") for particle " << iat << std::endl;
      }
      log << pad << "  Analytic    = " << L[iat] << std::endl;
      log << pad << "  Finite diff = " << L_fd[iat] << std::endl;
      log << pad << "  Error       = " << L_err  << "  Relative Error = " << L_err_rel << std::endl;
      all_okay = false;
    }

    ParticleSet::Scalar_t G_err_rel[OHMMS_DIM];
    for (int idim=0; idim<OHMMS_DIM; idim++)
    {
      ParticleSet::Scalar_t G_err = std::abs(G[iat][idim]-G_fd[iat][idim]);
      ParticleSet::Scalar_t G_rel_denom = std::max( std::abs(G[iat][idim]), std::abs(G_fd[iat][idim]) );
      G_err_rel[idim] = std::abs( G_err / G[iat][idim] );

      if (G_err_rel[idim] > rel_tol && G_err > abs_tol)
      {
        if (G_err_rel[idim] > rel_tol)
        {
          log << pad << "Finite difference gradient exceeds relative tolerance (" << rel_tol << ") for particle " << iat;
        }
        else
        {
          log << pad << "Finite difference gradient exceeds absolute tolerance (" << abs_tol << ") for particle " << iat;
        }
        log << " component " << idim << std::endl;
        log << pad << "  Analytic    = " << G[iat][idim] << std::endl;
        log << pad << "  Finite diff = " << G_fd[iat][idim] << std::endl;
        log << pad << "  Error       = " << G_err <<  "  Relative Error = " << G_err_rel[idim] << std::endl;
        all_okay = false;
      }
    }

    fout << pad << "For particle #" << iat << " at " << W.R[iat] << std::endl;
    fout << pad << "Gradient      = " << std::setw(12) << G[iat] << std::endl;
    fout << pad << "  Finite diff = " << std::setw(12) << G_fd[iat] << std::endl;
    fout << pad << "  Error       = " << std::setw(12) << G[iat]-G_fd[iat] << std::endl;
    fout << pad << "  Relative Error = ";
    for (int idim = 0; idim<OHMMS_DIM; idim++)
    {
        fout << G_err_rel[idim] << " ";
    }
    fout << std::endl << std::endl;
    fout << pad << "Laplacian     = " << std::setw(12) << L[iat] << std::endl;
    fout << pad << "  Finite diff = " << std::setw(12) << L_fd[iat] << std::endl;
    fout << pad << "  Error       = " << std::setw(12) << L[iat]-L_fd[iat] << "  Relative Error = " << L_err_rel << std::endl << std::endl;
  }
  return all_okay;
}

bool WaveFunctionTester::checkGradientAtConfiguration(MCWalkerConfiguration::Walker_t *W1, std::stringstream &fail_log, bool &ignore)
{

  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);

  W.loadWalker(*W1, true);

  // compute analytic values
  Psi.evaluateLog(W);
  G = W.G;
  L = W.L;

  // Use a single delta with a fairly large tolerance
  //computeNumericalGrad(delta, G1, L1);

  RealType delta = 1.0e-4;
  if (deltaParam > 0.0)
  {
    delta = deltaParam;
  } 
  FiniteDifference fd(FiniteDifference::FiniteDiff_LowOrder);
  //RealType delta = 1.0;
  //FiniteDifference fd(FiniteDifference::FiniteDiff_Richardson);

  FiniteDifference::PosChangeVector positions;

  fd.finiteDifferencePoints(delta, W, positions);

  FiniteDifference::ValueVector logpsi_vals;
  FiniteDifference::PosChangeVector::iterator it;

  for (it = positions.begin(); it != positions.end(); it++)
  {
    PosType r0 = W.R[it->index];
    W.R[it->index] = it->r;
    W.update();
    RealType logpsi0 = Psi.evaluateLog(W);
    RealType phase0 = Psi.getPhase();
#if defined(QMC_COMPLEX)
    ValueType logpsi = std::complex<OHMMS_PRECISION>(logpsi0,phase0);
#else
    ValueType logpsi = logpsi0;
#endif
    logpsi_vals.push_back(logpsi);

    W.R[it->index] = r0;
    W.update();
    Psi.evaluateLog(W);
  }

  fd.computeFiniteDiff(delta, positions, logpsi_vals, G1, L1);

  fout << "delta = " << delta << std::endl;

  // TODO - better choice of tolerance
  // TODO - adjust delta and tolerance based on precision of wavefunction

  bool all_okay = checkGradients(0, nat, G, L, G1, L1, fail_log);
  RealType tol = 1e-3;
  if (toleranceParam> 0.0)
  {
    tol = toleranceParam;
  } 

  for (int iorb = 0; iorb < Psi.getOrbitals().size(); iorb++)
  {
    OrbitalBase *orb = Psi.getOrbitals()[iorb];

    ParticleSet::ParticleGradient_t G(nat), tmpG(nat), G1(nat);
    ParticleSet::ParticleLaplacian_t L(nat), tmpL(nat), L1(nat);


    RealType logpsi1 = orb->evaluateLog(W, G, L);

    fail_log << "Orbital " << iorb << " " << orb->OrbitalName << " log psi = " << logpsi1 << std::endl;

    FiniteDifference::ValueVector logpsi_vals;
    FiniteDifference::PosChangeVector::iterator it;
    for (it = positions.begin(); it != positions.end(); it++)
    {
      PosType r0 = W.R[it->index];
      W.R[it->index] = it->r;
      W.update();
      ParticleSet::SingleParticlePos_t zeroR;
      W.makeMove(it->index,zeroR);

      RealType logpsi0 = orb->evaluateLog(W, tmpG, tmpL);
      RealType phase0 = orb->PhaseValue;
#if defined(QMC_COMPLEX)
      ValueType logpsi = std::complex<OHMMS_PRECISION>(logpsi0,phase0);
#else
      ValueType logpsi = logpsi0;
#endif
      logpsi_vals.push_back(logpsi);
      W.rejectMove(it->index);

      W.R[it->index] = r0;
      W.update();
    }
    fd.computeFiniteDiff(delta, positions, logpsi_vals, G1, L1);

    fout << "  Orbital " << iorb << " " << orb->OrbitalName << std::endl;

    if (!checkGradients(0, nat, G, L, G1, L1, fail_log, 1))
    {
      all_okay = false;
    }

    if (!checkSlaterDet) continue; // skip SlaterDet check if <backflow> is present
    // DiracDeterminantWithBackflow::evaluateLog requires a call to BackflowTransformation::evaluate in its owning SlaterDetWithBackflow to work correctly.
    SlaterDet *sd = dynamic_cast<SlaterDet *>(orb);
    if (sd)
    {
      for (int isd = 0; isd < sd->Dets.size(); isd++)
      {
        ParticleSet::ParticleGradient_t G(nat), tmpG(nat), G1(nat);
        ParticleSet::ParticleLaplacian_t L(nat), tmpL(nat), L1(nat);
        DiracDeterminantBase *det = sd->Dets[isd];
        RealType logpsi2 = det->evaluateLog(W, G, L); // this won't work with backflow
        fail_log << "  Slater Determiant " << isd << " (for particles " << det->FirstIndex << " to " << det->LastIndex << ") log psi = " << logpsi2 << std::endl;
        // Should really check the condition number on the matrix determinant.
        // For now, just ignore values that too small.
        if (logpsi2 < -40.0)
        {
            ignore = true;
        }
        FiniteDifference::ValueVector logpsi_vals;
        FiniteDifference::PosChangeVector::iterator it;
        for (it = positions.begin(); it != positions.end(); it++)
        {
          PosType r0 = W.R[it->index];
          W.R[it->index] = it->r;
          W.update();

          RealType logpsi0 = det->evaluateLog(W, tmpG, tmpL);
          RealType phase0 = det->PhaseValue;
    #if defined(QMC_COMPLEX)
          ValueType logpsi = std::complex<OHMMS_PRECISION>(logpsi0,phase0);
    #else
          ValueType logpsi = logpsi0;
    #endif
          logpsi_vals.push_back(logpsi);

          W.R[it->index] = r0;
          W.update();
        }
        fd.computeFiniteDiff(delta, positions, logpsi_vals, G1, L1);

        if (!checkGradients(det->FirstIndex, det->LastIndex, G, L, G1, L1, fail_log, 2))
        {
          all_okay = false;
        }

#if 0
        // Testing single particle orbitals doesn't work yet - probably something
        // with setup after setting the position.
        std::map<std::string, SPOSetBasePtr>::iterator spo_it = sd->mySPOSet.begin();
        for (; spo_it != sd->mySPOSet.end(); spo_it++)
        {
          SPOSetBasePtr spo = spo_it->second;
          fail_log << "      SPO set = " << spo_it->first <<  " name = " << spo->className;
          fail_log << " orbital set size = " << spo->size();
          fail_log << " basis set size = " << spo->getBasisSetSize() << std::endl;

          ParticleSet::ParticleGradient_t G(nat), tmpG(nat), G1(nat);
          ParticleSet::ParticleLaplacian_t L(nat), tmpL(nat), L1(nat);
          RealType logpsi3 = det->evaluateLog(W, G, L);
          FiniteDifference::ValueVector logpsi_vals;
          FiniteDifference::PosChangeVector::iterator it;
          for (it = positions.begin(); it != positions.end(); it++)
          {
            PosType r0 = W.R[it->index];
            W.R[it->index] = it->r;
            W.update();
            ParticleSet::SingleParticlePos_t zeroR;
            W.makeMove(it->index,zeroR);

            SPOSetBase::ValueVector_t psi(spo->size());

            spo->evaluate(W, it->index, psi);
            ValueType logpsi = psi[0];
            logpsi_vals.push_back(logpsi);

            W.rejectMove(it->index);
            W.R[it->index] = r0;
            W.update();
          }
          fd.computeFiniteDiff(delta, positions, logpsi_vals, G1, L1);

          if (!checkGradients(det->FirstIndex, det->LastIndex, G, L, G1, L1, fail_log, 3))
          {
            all_okay = false;
          }
        }
#endif
      }
    }
  }
  return all_okay;
}

void WaveFunctionTester::runBasicTest()
{
  RealType sig2Enloc=0, sig2Drift=0;

  int nat = W.getTotalNum();
  fout << "Numerical gradient and Laplacian test" << std::endl;

  std::stringstream fail_log;
  bool all_okay = true;
  int fails = 0;
  int nconfig = 0;
  int nignore = 0;

  MCWalkerConfiguration::iterator Wit(W.begin());
  for (; Wit != W.end(); Wit++)
  {
    fout << "Walker # " << nconfig << std::endl;
    std::stringstream fail_log1;
    bool ignore = false;
    bool this_okay =  checkGradientAtConfiguration(*Wit, fail_log1, ignore);
    if (ignore)
    {
      nignore++;
    }
    if (!this_okay && !ignore)
    {
      fail_log << "Walker # " << nconfig << std::endl;
      fail_log << fail_log1.str();
      fail_log << std::endl;
      fails++;
      all_okay = false;
    }
    nconfig++;
  }

  app_log() << "Number of samples = " << nconfig << std::endl;
  app_log() << "Number ignored (bad positions) = " << nignore << std::endl << std::endl;
  app_log() << "Number of fails = " << fails << std::endl << std::endl;

  if (!all_okay)
  {
    std::string fail_name("wf_fails.dat");
    app_log() << "More detail on finite difference failures in " << fail_name << std::endl;
    std::ofstream eout(fail_name.c_str());
    eout << fail_log.str();
    eout.close();
  }

  app_log() << "Finite difference test: " << (all_okay?"PASS":"FAIL") << std::endl;

  // Compute approximation error vs. delta.
  if (outputDeltaVsError)
  {
    double delta = 1.0;
    bool inputOkay = true;

    std::ofstream dout(DeltaVsError.outputFile.c_str());

    int iat = DeltaVsError.particleIndex;
    int ig = DeltaVsError.gradientComponentIndex;

    if (iat < 0 || iat >= nat)
    {
      dout << "# Particle index (" << iat << ") is out of range (0 - " << nat-1 << ")" << std::endl;
      inputOkay = false;
    }

    if (ig < 0 || ig >= OHMMS_DIM)
    {
      dout << "# Gradient component index (" << ig << ") is out of range (0 - " << OHMMS_DIM-1 << ")" << std::endl;
      inputOkay = false;
    }

    if (inputOkay)
    {
      dout << "# Particle = " << iat << " Gradient component = " << ig << std::endl;
      dout << "#" << std::setw(11) <<  "delta" << std::setw(14) << "G_err_rel" << std::setw(14) << "L_err_rel" << std::endl;
      ParticleSet::ParticleGradient_t G(nat), G1(nat);
      ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
      for (int i = 0; i < 20; i++) {
        // compute analytic values
        G = W.G;
        L = W.L;
        Psi.evaluateLog(W);

        computeNumericalGrad(delta, G1, L1);
        ParticleSet::Scalar_t L_err = std::abs(L[iat]-L1[iat]);
        ParticleSet::Scalar_t L_err_rel = std::abs( L_err/L[iat] );
        ParticleSet::Scalar_t G_err = std::abs(G[iat][ig]-G1[iat][ig]);
        ParticleSet::Scalar_t G_err_rel = std::abs(G_err/G[iat][ig]);
        dout << std::setw(12) << delta;
        dout << std::setw(14) << std::abs(G_err_rel) << std::setw(14) << std::abs(L_err_rel);
        dout << std::endl;
        delta *= std::sqrt(0.1);
      }
    }
    dout.close();
  }

  fout << "Ratio test" << std::endl;

  RealType tol = 1e-3;
  RealType ratio_tol = 1e-9;
  bool any_ratio_fail = false;
  makeGaussRandom(deltaR);
  fout << "deltaR:" << std::endl;
  fout << deltaR << std::endl;
  fout << "Particle       Ratio of Ratios     Computed Ratio   Internal Ratio" << std::endl;
  //testing ratio alone
  for (int iat=0; iat<nat; iat++)
  {
    W.update();
    //ValueType psi_p = log(std::abs(Psi.evaluate(W)));
    RealType psi_p = Psi.evaluateLog(W);
    RealType phase_p=Psi.getPhase();
    W.makeMove(iat,deltaR[iat]);
    //W.update();
    RealType aratio = Psi.ratio(W,iat);
    RealType phaseDiff = Psi.getPhaseDiff();
    W.rejectMove(iat);
    Psi.rejectMove(iat);
    W.R[iat] += deltaR[iat];
    W.update();
    //ValueType psi_m = log(std::abs(Psi.evaluate(W)));
    RealType psi_m = Psi.evaluateLog(W);
    RealType phase_m=Psi.getPhase();
    
#if defined(QMC_COMPLEX)
    RealType ratioMag = std::exp(psi_m-psi_p);
    RealType dphase = phase_m-phase_p;
    if(phaseDiff < 0.0)
      phaseDiff += 2.0*M_PI;
    if(phaseDiff > 2.0*M_PI)
      phaseDiff -= 2.0*M_PI;
    if(dphase < 0.0)
      dphase += 2.0*M_PI;
    if(dphase > 2.0*M_PI)
      dphase -= 2.0*M_PI;
    ValueType ratDiff=std::complex<OHMMS_PRECISION>(ratioMag*std::cos(dphase),ratioMag*std::sin(dphase)) ;
    // TODO - test complex ratio against a tolerance
    fout << iat << " " << aratio*std::complex<OHMMS_PRECISION>(std::cos(phaseDiff),std::sin(phaseDiff))/ratDiff << " " << ratDiff << " " << aratio*std::complex<OHMMS_PRECISION>(std::cos(phaseDiff),std::sin(phaseDiff)) << std::endl;
    fout << "     ratioMag " << aratio/ratioMag << " " << ratioMag << std::endl;
    fout << "     PhaseDiff " << phaseDiff/dphase << " " << phaseDiff  <<" " << dphase << std::endl;
#else
    RealType ratDiff=std::exp(psi_m-psi_p)*std::cos(phase_m-phase_p);
    fout << iat << " " << aratio/ratDiff << " " << ratDiff << " " << aratio << std::endl;
    if (std::abs(aratio/ratDiff - 1.0) > ratio_tol)
    {
      app_log() << "Wavefunction ratio exceeds tolerance " << tol << ") for particle " << iat << std::endl;
      app_log() << "  Internally computed ratio = " << aratio << std::endl;
      app_log() << "  Separately computed ratio = " << ratDiff << std::endl;
      any_ratio_fail = true;
    }
#endif
  }
  app_log() << "Ratio test: " << (any_ratio_fail?"FAIL":"PASS") << std::endl;
}

void WaveFunctionTester::runRatioTest()
{
#if 0
  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
  ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);
  bool checkHam=(checkHamPbyP == "yes");
  Tau=0.025;
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while (it != it_end)
  {
    makeGaussRandom(deltaR);
    Walker_t::WFBuffer_t tbuffer;
    W.R = (**it).R+Tau*deltaR;
    (**it).R=W.R;
    W.update();
    RealType logpsi=Psi.registerData(W,tbuffer);
    RealType ene;
    if (checkHam)
      ene = H.registerData(W,tbuffer);
    else
      ene = H.evaluate(W);
    (*it)->DataSet=tbuffer;
    //RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
    H.saveProperty((*it)->getPropertyBase());
    ++it;
    app_log() << "  HamTest " << "  Total " <<  ene << std::endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  }
  fout << "  Update using drift " << std::endl;
  bool pbyp_mode=true;
  for (int iter=0; iter<4; ++iter)
  {
    int iw=0;
    it=W.begin();
    while (it != it_end)
    {
      fout << "\nStart Walker " << iw++ << std::endl;
      Walker_t& thisWalker(**it);
      W.loadWalker(thisWalker,pbyp_mode);
      Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
      Psi.copyFromBuffer(W,w_buffer);
      H.copyFromBuffer(W,w_buffer);
//             Psi.evaluateLog(W);
      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);
      makeGaussRandom(deltaR);
      //mave a move
      RealType ratio_accum(1.0);
      for (int iat=0; iat<nat; iat++)
      {
        PosType dr(Tau*deltaR[iat]);
        PosType newpos(W.makeMove(iat,dr));
        //RealType ratio=Psi.ratio(W,iat,dGp,dLp);
        W.dG=0;
        W.dL=0;
        RealType ratio=Psi.ratio(W,iat,W.dG,W.dL);
        Gp = W.G + W.dG;
        //RealType enew = H.evaluatePbyP(W,iat);
        if (checkHam)
          enew = H.evaluatePbyP(W,iat);
        if (ratio > Random())
        {
          fout << " Accepting a move for " << iat << std::endl;
          fout << " Energy after a move " << enew << std::endl;
          W.G += W.dG;
          W.L += W.dL;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          if (checkHam)
            H.acceptMove(iat);
          ratio_accum *= ratio;
        }
        else
        {
          fout << " Rejecting a move for " << iat << std::endl;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
          //H.rejectMove(iat);
        }
      }
      fout << " Energy after pbyp = " << H.getLocalEnergy() << std::endl;
      RealType newlogpsi_up = Psi.evaluateLog(W,w_buffer);
      W.saveWalker(thisWalker);
      RealType ene_up;
      if (checkHam)
        ene_up= H.evaluate(W,w_buffer);
      else
        ene_up = H.evaluate(W);
      Gp=W.G;
      Lp=W.L;
      W.R=thisWalker.R;
      W.update();
      RealType newlogpsi=Psi.updateBuffer(W,w_buffer,false);
      RealType ene = H.evaluate(W);
      thisWalker.resetProperty(newlogpsi,Psi.getPhase(),ene);
      //thisWalker.resetProperty(std::log(psi),Psi.getPhase(),ene);
      fout << iter << "  Energy by update = "<< ene_up << " " << ene << " "  << ene_up-ene << std::endl;
      fout << iter << " Ratio " << ratio_accum*ratio_accum
            << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " "
            << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << std::endl
            << " new log(psi) updated " << newlogpsi_up
            << " new log(psi) calculated " << newlogpsi
            << " old log(psi) " << logpsi << std::endl;
      fout << " Gradients " << std::endl;
      for (int iat=0; iat<nat; iat++)
        fout << W.G[iat]-Gp[iat] << W.G[iat] << std::endl; //W.G[iat] << G[iat] << std::endl;
      fout << " Laplacians " << std::endl;
      for (int iat=0; iat<nat; iat++)
        fout << W.L[iat]-Lp[iat] << " " << W.L[iat] << std::endl;
      ++it;
    }
  }
  fout << "  Update without drift : for VMC useDrift=\"no\"" << std::endl;
  for (int iter=0; iter<4; ++iter)
  {
    it=W.begin();
    int iw=0;
    while (it != it_end)
    {
      fout << "\nStart Walker " << iw++ << std::endl;
      Walker_t& thisWalker(**it);
      W.loadWalker(thisWalker,pbyp_mode);
      Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
      //Psi.updateBuffer(W,w_buffer,true);
      Psi.copyFromBuffer(W,w_buffer);
      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);
      //mave a move
      RealType ratio_accum(1.0);
      for (int substep=0; substep<3; ++substep)
      {
        makeGaussRandom(deltaR);
        for (int iat=0; iat<nat; iat++)
        {
          PosType dr(Tau*deltaR[iat]);
          PosType newpos(W.makeMove(iat,dr));
          RealType ratio=Psi.ratio(W,iat);
          RealType prob = ratio*ratio;
          if (prob > Random())
          {
            fout << " Accepting a move for " << iat << std::endl;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            ratio_accum *= ratio;
          }
          else
          {
            fout << " Rejecting a move for " << iat << std::endl;
            W.rejectMove(iat);
            Psi.rejectMove(iat);
          }
        }
        RealType logpsi_up = Psi.updateBuffer(W,w_buffer,false);
        W.saveWalker(thisWalker);
        RealType ene = H.evaluate(W);
        thisWalker.resetProperty(logpsi_up,Psi.getPhase(),ene);
      }
      Gp=W.G;
      Lp=W.L;
      W.update();
      RealType newlogpsi=Psi.evaluateLog(W);
      fout << iter << " Ratio " << ratio_accum*ratio_accum
            << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " "
            << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << std::endl
            << " new log(psi) " << newlogpsi
            << " old log(psi) " << logpsi << std::endl;
      fout << " Gradients " << std::endl;
      for (int iat=0; iat<nat; iat++)
      {
        fout << W.G[iat]-Gp[iat] << W.G[iat] << std::endl; //W.G[iat] << G[iat] << std::endl;
      }
      fout << " Laplacians " << std::endl;
      for (int iat=0; iat<nat; iat++)
      {
        fout << W.L[iat]-Lp[iat] << " " << W.L[iat] << std::endl;
      }
      ++it;
    }
  }
  //for(it=W.begin();it != it_end; ++it)
  //{
  //  Walker_t& thisWalker(**it);
  //  Walker_t::WFBuffer_t& w_buffer((*it)->DataSet);
  //  w_buffer.rewind();
  //  W.updateBuffer(**it,w_buffer);
  //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
  //}
 #endif
}

void WaveFunctionTester::runRatioTest2()
{
  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
  ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);
  Tau=0.025;
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  for (; it != it_end; ++it)
  {
    makeGaussRandom(deltaR);
    Walker_t::WFBuffer_t tbuffer;
    (**it).R  +=  Tau*deltaR;
    W.loadWalker(**it,true);
    Psi.registerData(W,tbuffer);
    tbuffer.allocate();
    Psi.copyFromBuffer(W,tbuffer);
    Psi.evaluateLog(W);
    RealType logpsi = Psi.updateBuffer(W,tbuffer,false);
    RealType ene = H.evaluate(W);
    (*it)->DataSet=tbuffer;
    //RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
    H.saveProperty((*it)->getPropertyBase());
    app_log() << "  HamTest " << "  Total " <<  ene << std::endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  }
  for (int iter=0; iter<20; ++iter)
  {
    int iw=0;
    it=W.begin();
    //while(it != it_end)
    for (; it != it_end; ++it)
    {
      fout << "\nStart Walker " << iw++ << std::endl;
      Walker_t& thisWalker(**it);
      W.loadWalker(thisWalker,true);
      Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
      Psi.copyFromBuffer(W,w_buffer);
      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);
      Psi.evaluateLog(W);
      ParticleSet::ParticleGradient_t realGrad(W.G);
      makeGaussRandom(deltaR);
      //mave a move
      RealType ratio_accum(1.0);
      for (int iat=0; iat<nat; iat++)
      {
        TinyVector<ParticleSet::ParticleValue_t,OHMMS_DIM> grad_now=Psi.evalGrad(W,iat);
        GradType grad_new;
        for(int sds=0; sds<3; sds++)
          fout<< realGrad[iat][sds]-grad_now[sds]<<" ";
        PosType dr(Tau*deltaR[iat]);
        PosType newpos(W.makeMove(iat,dr));
        RealType ratio2 = Psi.ratioGrad(W,iat,grad_new);
        W.rejectMove(iat);
        Psi.rejectMove(iat);
        newpos=W.makeMove(iat,dr);
        RealType ratio1 = Psi.ratio(W,iat);
        //Psi.rejectMove(iat);
        W.rejectMove(iat);
        fout << "  ratio1 = " << ratio1 << " ration2 = " << ratio2 << std::endl;
      }
    }
  }
  //for(it=W.begin();it != it_end; ++it)
  //{
  //  Walker_t& thisWalker(**it);
  //  Walker_t::WFBuffer_t& w_buffer((*it)->DataSet);
  //  w_buffer.rewind();
  //  W.updateBuffer(**it,w_buffer);
  //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
  //}
}

template<typename T, unsigned D>
inline void randomize(ParticleAttrib<TinyVector<T,D> >& displ, T fac)
{
  T* rv=&(displ[0][0]);
  assignUniformRand(rv, displ.size()*D, Random);
  for(int i=0; i<displ.size()*D; ++i) rv[i] =fac*(rv[i]-0.5);
}

void WaveFunctionTester::runRatioV()
{
#if 0
  app_log() << "WaveFunctionTester::runRatioV " << std::endl;
  int nat = W.getTotalNum();
  Tau=0.025;

  //create a VP with 8 virtual moves
  VirtualParticleSet vp(&W,8);
  W.enableVirtualMoves();

  //cheating
  const ParticleSet& ions=W.DistTables[1]->origin();
  DistanceTableData* dt_ie=W.DistTables[1];
  double Rmax=2.0;

  ParticleSet::ParticlePos_t sphere(8);
  std::vector<RealType> ratio_1(8), ratio_v(8);
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while (it != it_end)
  {
    makeGaussRandom(deltaR);
    Walker_t::WFBuffer_t tbuffer;
    W.R = (**it).R+Tau*deltaR;
    (**it).R=W.R;
    W.update();
    RealType logpsi=Psi.registerData(W,tbuffer);

    W.initVirtualMoves();

    for(int iat=0; iat<ions.getTotalNum(); ++iat)
    {
      for(int nn=dt_ie->M[iat],iel=0; nn<dt_ie->M[iat+1]; nn++,iel++)
      {
        register RealType r(dt_ie->r(nn));
        if(r>Rmax) continue;
        randomize(sphere,(RealType)0.5);
        
        for(int k=0; k<sphere.size(); ++k)
        {
          W.makeMoveOnSphere(iel,sphere[k]);
          ratio_1[k]=Psi.ratio(W,iel);
          W.rejectMove(iel);
          Psi.resetPhaseDiff();
        }

        vp.makeMoves(iel,sphere);

        Psi.evaluateRatios(vp,ratio_v);

        app_log() << "IAT = " << iat << " " << iel << std::endl;
        for(int k=0; k<sphere.size(); ++k)
        {
          app_log() << ratio_1[k]/ratio_v[k] << " " << ratio_1[k] << std::endl;
        }
        app_log() << std::endl;
      }
    }
    ++it;
  }
#endif
}

void WaveFunctionTester::runGradSourceTest()
{
  ParticleSetPool::PoolType::iterator p;
  for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
    app_log() << "ParticelSet = " << p->first << std::endl;
  // Find source ParticleSet
  ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
  app_log() << pit->first << std::endl;
  // if(pit == PtclPool.getPool().end())
  //   APP_ABORT("Unknown source \"" + sourceName + "\" WaveFunctionTester.");
  ParticleSet& source = *((*pit).second);
  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 0.00001;
  RealType delta2 = 2*delta;
  ValueType c1 = 1.0/delta/2.0;
  ValueType c2 = 1.0/delta/delta;
  int nat = W.getTotalNum();
  ParticleSet::ParticlePos_t deltaR(nat);
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
  //copy the properties of the working walker
  Properties = awalker->Properties;
  //sample a new walker configuration and copy to ParticleSet::R
  //makeGaussRandom(deltaR);
  W.R = awalker->R;
  //W.R += deltaR;
  W.update();
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);
  app_log() << "  HamTest " << "  Total " <<  eloc << std::endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;
  for (int isrc=0; isrc < 1/*source.getTotalNum()*/; isrc++)
  {
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad;
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad_FD;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_FD;
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      grad_grad[dim].resize(nat);
      lapl_grad[dim].resize(nat);
      grad_grad_FD[dim].resize(nat);
      lapl_grad_FD[dim].resize(nat);
    }
    Psi.evaluateLog(W);
    GradType grad_log = Psi.evalGradSource(W, source, isrc, grad_grad, lapl_grad);
    ValueType log = Psi.evaluateLog(W);
    //grad_log = Psi.evalGradSource (W, source, isrc);
    for (int iat=0; iat<nat; iat++)
    {
      PosType r0 = W.R[iat];
      GradType gFD[OHMMS_DIM];
      GradType lapFD = ValueType();
      for (int eldim=0; eldim<3; eldim++)
      {
        W.R[iat][eldim] = r0[eldim]+delta;
        W.update();
        ValueType log_p = Psi.evaluateLog(W);
        GradType gradlogpsi_p =  Psi.evalGradSource(W, source, isrc);
        W.R[iat][eldim] = r0[eldim]-delta;
        W.update();
        ValueType log_m = Psi.evaluateLog(W);
        GradType gradlogpsi_m = Psi.evalGradSource(W, source, isrc);
        lapFD    += gradlogpsi_m + gradlogpsi_p;
        gFD[eldim] = gradlogpsi_p - gradlogpsi_m;
        W.R[iat] = r0;
        W.update();
        //Psi.evaluateLog(W);
      }
      const ValueType six(6);
      for (int iondim=0; iondim<OHMMS_DIM; iondim++)
      {
        for (int eldim=0; eldim<OHMMS_DIM; eldim++)
          grad_grad_FD[iondim][iat][eldim] = c1*gFD[eldim][iondim];
        lapl_grad_FD[iondim][iat] = c2*(lapFD[iondim]-six*grad_log[iondim]);
      }
    }
    for (int dimsrc=0; dimsrc<OHMMS_DIM; dimsrc++)
    {
      for (int iat=0; iat<nat; iat++)
      {
        fout << "For particle #" << iat << " at " << W.R[iat] << std::endl;
        fout << "Gradient      = " << std::setw(12) << grad_grad[dimsrc][iat] << std::endl
             << "  Finite diff = " << std::setw(12) << grad_grad_FD[dimsrc][iat] << std::endl
             << "  Error       = " << std::setw(12)
             <<  grad_grad_FD[dimsrc][iat] - grad_grad[dimsrc][iat] << std::endl << std::endl;
        fout << "Laplacian     = " << std::setw(12) << lapl_grad[dimsrc][iat] << std::endl
             << "  Finite diff = " << std::setw(12) << lapl_grad_FD[dimsrc][iat] << std::endl
             << "  Error       = " << std::setw(12)
             << lapl_grad_FD[dimsrc][iat] - lapl_grad[dimsrc][iat] << std::endl << std::endl;
      }
    }
  }
}


void WaveFunctionTester::runZeroVarianceTest()
{
  ParticleSetPool::PoolType::iterator p;
  for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
    app_log() << "ParticelSet = " << p->first << std::endl;
  // Find source ParticleSet
  ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
  app_log() << pit->first << std::endl;
  // if(pit == PtclPool.getPool().end())
  //   APP_ABORT("Unknown source \"" + sourceName + "\" WaveFunctionTester.");
  ParticleSet& source = *((*pit).second);
  int nat = W.getTotalNum();
  ParticleSet::ParticlePos_t deltaR(nat);
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
  //copy the properties of the working walker
  Properties = awalker->Properties;
  //sample a new walker configuration and copy to ParticleSet::R
  //makeGaussRandom(deltaR);
  W.R = awalker->R;
  //W.R += deltaR;
  W.update();
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);
  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;
  PosType r1(5.0, 2.62, 2.55);
  W.R[1] = PosType(4.313, 5.989, 4.699);
  W.R[2] = PosType(5.813, 4.321, 4.893);
  W.R[3] = PosType(4.002, 5.502, 5.381);
  W.R[4] = PosType(5.901, 5.121, 5.311);
  W.R[5] = PosType(5.808, 4.083, 5.021);
  W.R[6] = PosType(4.750, 5.810, 4.732);
  W.R[7] = PosType(4.690, 5.901, 4.989);
  for (int i=1; i<8; i++)
    W.R[i] -= PosType(2.5, 2.5, 2.5);
  char fname[32];
  sprintf(fname,"ZVtest.%03d.dat",OHMMS::Controller->rank());
  FILE *fzout = fopen(fname, "w");
  TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad;
  TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad;
  TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad_FD;
  TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_FD;
  for (int dim=0; dim<OHMMS_DIM; dim++)
  {
    grad_grad[dim].resize(nat);
    lapl_grad[dim].resize(nat);
    grad_grad_FD[dim].resize(nat);
    lapl_grad_FD[dim].resize(nat);
  }
  for (r1[0]=0.0; r1[0]<5.0; r1[0]+=1.0e-4)
  {
    W.R[0] = r1;
    fprintf(fzout, "%1.8e %1.8e %1.8e ", r1[0], r1[1], r1[2]);
    RealType log = Psi.evaluateLog(W);
//        ValueType psi = std::cos(Psi.getPhase())*std::exp(log);//*W.PropertyList[SIGN];
#if defined(QMC_COMPLEX)
    RealType ratioMag = std::exp(log);
    ValueType psi=std::complex<OHMMS_PRECISION>(ratioMag*std::cos(Psi.getPhase()),ratioMag*std::sin(Psi.getPhase())) ;
#else
    ValueType psi = std::cos(Psi.getPhase())*std::exp(log);//*W.PropertyList[SIGN];
#endif
    double E = H.evaluate(W);
    //double KE = E - W.PropertyList[LOCALPOTENTIAL];
    double KE = -0.5*(Sum(W.L) + Dot(W.G,W.G));
#if defined(QMC_COMPLEX)
    fprintf(fzout, "%16.12e %16.12e %16.12e ", psi.real(), psi.imag(),KE);
#else
    fprintf(fzout, "%16.12e %16.12e ", psi, KE);
#endif
    for (int isrc=0; isrc < source.getTotalNum(); isrc++)
    {
      GradType grad_log = Psi.evalGradSource(W, source, isrc, grad_grad, lapl_grad);
      for (int dim=0; dim<OHMMS_DIM; dim++)
      {
        double ZV = 0.5*Sum(lapl_grad[dim]) + Dot(grad_grad[dim], W.G);
#if defined(QMC_COMPLEX)
        fprintf(fzout, "%16.12e %16.12e %16.12e ", ZV, grad_log[dim].real(), grad_log[dim].imag());
#else
        fprintf(fzout, "%16.12e %16.12e ", ZV, grad_log[dim]);
#endif
      }
    }
    fprintf(fzout, "\n");
  }
  fclose(fzout);
}



bool
WaveFunctionTester::put(xmlNodePtr q)
{
  myNode = q;
  xmlNodePtr tcur=q->children;
  while(tcur != NULL)
  {
    std::string cname((const char*)(tcur->name));
    if(cname == "delta_output")
    {
      outputDeltaVsError = true;
      DeltaVsError.put(tcur);
    }
    tcur = tcur->next;
  }

  bool success = putQMCInfo(q);

  return success;
}

void WaveFunctionTester::runDerivTest()
{
  app_log()<<" Testing derivatives"<< std::endl;
  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 1e-6;
  RealType delta2 = 2*delta;
  ValueType c1 = 1.0/delta/2.0;
  ValueType c2 = 1.0/delta/delta;
  int nat = W.getTotalNum();
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
  //copy the properties of the working walker
  Properties = awalker->Properties;
  //sample a new walker configuration and copy to ParticleSet::R
  W.R = awalker->R+deltaR;

  fout << "Position " << std::endl << W.R << std::endl;

  //W.R += deltaR;
  W.update();
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);
  app_log() << "  HamTest " << "  Total " <<  eloc << std::endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;
  fout<<"Gradients"<< std::endl;
  for (int iat=0; iat<W.R.size(); iat++)
  {
    for (int i=0; i<3 ; i++)
      fout<<W.G[iat][i]<<"  ";
    fout<< std::endl;
  }
  fout<<"Laplaians"<< std::endl;
  for (int iat=0; iat<W.R.size(); iat++)
  {
    fout<<W.L[iat]<<"  ";
    fout<< std::endl;
  }
  opt_variables_type wfVars,wfvar_prime;
//build optimizables from the wavefunction
  wfVars.clear();
  Psi.checkInVariables(wfVars);
  wfVars.resetIndex();
  Psi.checkOutVariables(wfVars);
  wfvar_prime= wfVars;
  wfVars.print(fout);
  int Nvars= wfVars.size();
  std::vector<RealType> Dsaved(Nvars);
  std::vector<RealType> HDsaved(Nvars);
  std::vector<RealType> PGradient(Nvars);
  std::vector<RealType> HGradient(Nvars);
  Psi.resetParameters(wfVars);
  logpsi = Psi.evaluateLog(W);

  //reuse the sphere
  H.setPrimary(false);

  eloc=H.evaluate(W);
  Psi.evaluateDerivatives(W, wfVars, Dsaved, HDsaved);
  RealType FiniteDiff = 1e-6;
  QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++)
      wfvar_prime[j]=wfVars[j];
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
//     Psi.checkOutVariables(wfvar_prime);
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiPlus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocPlus=H.getLocalEnergy()-H.getLocalPotential();
    wfvar_prime[i] = wfVars[i]- FiniteDiff;
//     Psi.checkOutVariables(wfvar_prime);
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiMinus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocMinus = H.getLocalEnergy()-H.getLocalPotential();
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }
  Psi.resetParameters(wfVars);
  fout<< std::endl<<"Deriv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i]) << std::endl;
  fout<< std::endl<<"Hderiv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i]) << std::endl;
}


void WaveFunctionTester::runDerivNLPPTest()
{
  char fname[16];
  sprintf(fname,"nlpp.%03d",OHMMS::Controller->rank());
  std::ofstream nlout(fname);
  nlout.precision(15);

  app_log()<<" Testing derivatives"<< std::endl;
  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 1e-6;
  RealType delta2 = 2*delta;
  ValueType c1 = 1.0/delta/2.0;
  ValueType c2 = 1.0/delta/delta;
  int nat = W.getTotalNum();
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
  //copy the properties of the working walker
  Properties = awalker->Properties;
  //sample a new walker configuration and copy to ParticleSet::R
  W.R = awalker->R+deltaR;

  //W.R += deltaR;
  W.update();
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);

  app_log() << "  HamTest " << "  Total " <<  eloc << std::endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;

  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;
  nlout<<"Gradients"<< std::endl;
  for (int iat=0; iat<W.R.size(); iat++)
  {
    for (int i=0; i<3 ; i++)
      nlout<<W.G[iat][i]<<"  ";
    nlout<< std::endl;
  }
  nlout<<"Laplaians"<< std::endl;
  for (int iat=0; iat<W.R.size(); iat++)
  {
    nlout<<W.L[iat]<<"  ";
    nlout<< std::endl;
  }
  opt_variables_type wfVars,wfvar_prime;
//build optimizables from the wavefunction
  wfVars.clear();
  Psi.checkInVariables(wfVars);
  wfVars.resetIndex();
  Psi.checkOutVariables(wfVars);
  wfvar_prime= wfVars;
  wfVars.print(nlout);
  int Nvars= wfVars.size();
  std::vector<RealType> Dsaved(Nvars);
  std::vector<RealType> HDsaved(Nvars);
  std::vector<RealType> PGradient(Nvars);
  std::vector<RealType> HGradient(Nvars);
  Psi.resetParameters(wfVars);

  logpsi = Psi.evaluateLog(W);

  //reuse the sphere for non-local pp
  H.setPrimary(false);

  std::vector<RealType> ene(4), ene_p(4), ene_m(4);
  Psi.evaluateDerivatives(W, wfVars, Dsaved, HDsaved);
  ene[0]=H.evaluateValueAndDerivatives(W,wfVars,Dsaved,HDsaved,true);
  app_log() << "Check the energy " << eloc << " " << H.getLocalEnergy() << " " << ene[0] << std::endl;

  RealType FiniteDiff = 1e-6;
  QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++)
      wfvar_prime[j]=wfVars[j];
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiPlus = Psi.evaluateLog(W);
    RealType elocPlus=H.evaluateVariableEnergy(W,true);

    //H.evaluate(W);
    //RealType elocPlus=H.getLocalEnergy()-H.getLocalPotential();

    wfvar_prime[i] = wfVars[i]- FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiMinus = Psi.evaluateLog(W);
    RealType elocMinus=H.evaluateVariableEnergy(W,true);

    //H.evaluate(W);
    //RealType elocMinus = H.getLocalEnergy()-H.getLocalPotential();
    
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }

  nlout<< std::endl<<"Deriv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    nlout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i]) << std::endl;
  nlout<< std::endl<<"Hderiv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    nlout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i]) << std::endl;
}


void WaveFunctionTester::runDerivCloneTest()
{
  app_log()<<" Testing derivatives clone"<< std::endl;
  RandomGenerator_t* Rng1= new RandomGenerator_t();
  RandomGenerator_t* Rng2= new RandomGenerator_t();
  (*Rng1) = (*Rng2);
  MCWalkerConfiguration* w_clone = new MCWalkerConfiguration(W);
  TrialWaveFunction *psi_clone = Psi.makeClone(*w_clone);
  QMCHamiltonian *h_clone = H.makeClone(*w_clone,*psi_clone);
  h_clone->setRandomGenerator(Rng2);
  H.setRandomGenerator(Rng1);
  h_clone->setPrimary(true);
  int nat = W.getTotalNum();
  ParticleSet::ParticlePos_t deltaR(nat);
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
//   MCWalkerConfiguration::Walker_t* bwalker = *(w_clone->begin());
//   bwalker->R = awalker->R;
  W.R = awalker->R;
  W.update();
  w_clone->R=awalker->R;
  w_clone->update();
  opt_variables_type wfVars;
  //build optimizables from the wavefunction
//   wfVars.clear();
  Psi.checkInVariables(wfVars);
  wfVars.resetIndex();
  Psi.checkOutVariables(wfVars);
  wfVars.print(fout);
  int Nvars= wfVars.size();
  opt_variables_type wfvar_prime;
//   wfvar_prime.insertFrom(wfVars);
//   wfvar_prime.clear();
  psi_clone->checkInVariables(wfvar_prime);
  wfvar_prime.resetIndex();
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  psi_clone->checkOutVariables(wfvar_prime);
  wfvar_prime.print(fout);
  psi_clone->resetParameters(wfvar_prime);
  Psi.resetParameters(wfVars);
  std::vector<RealType> Dsaved(Nvars,0), og_Dsaved(Nvars,0);
  std::vector<RealType> HDsaved(Nvars,0), og_HDsaved(Nvars,0);
  std::vector<RealType> PGradient(Nvars,0), og_PGradient(Nvars,0);
  std::vector<RealType> HGradient(Nvars,0), og_HGradient(Nvars,0);
  ValueType logpsi2 = psi_clone->evaluateLog(*w_clone);
  RealType eloc2  = h_clone->evaluate(*w_clone);
  psi_clone->evaluateDerivatives(*w_clone, wfvar_prime, Dsaved, HDsaved);
  ValueType logpsi1 = Psi.evaluateLog(W);
  RealType eloc1  = H.evaluate(W);
  Psi.evaluateDerivatives(W, wfVars, og_Dsaved, og_HDsaved);
  app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << std::endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << std::endl;
  app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << std::endl;
  for (int i=0; i<h_clone->sizeOfObservables(); i++)
    app_log() << "  HamTest " << h_clone->getObservableName(i) << " " << h_clone->getObservable(i) << std::endl;
//   app_log()<<" Saved quantities:"<< std::endl;
//   for(int i=0;i<Nvars;i++) app_log()<<Dsaved[i]<<" "<<og_Dsaved[i]<<"         "<<HDsaved[i]<<" "<<og_HDsaved[i]<< std::endl;
  RealType FiniteDiff = 1e-6;
  QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++)
      wfvar_prime[j]=wfVars[j];
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
    psi_clone->resetParameters(wfvar_prime);
    psi_clone->reset();
    w_clone->update();
    w_clone->G=0;
    w_clone->L=0;
    RealType logpsiPlus = psi_clone->evaluateLog(*w_clone);
    h_clone->evaluate(*w_clone);
    RealType elocPlus=h_clone->getLocalEnergy()-h_clone->getLocalPotential();
    wfvar_prime[i] = wfVars[i]- FiniteDiff;
    psi_clone->resetParameters(wfvar_prime);
    psi_clone->reset();
    w_clone->update();
    w_clone->G=0;
    w_clone->L=0;
    RealType logpsiMinus = psi_clone->evaluateLog(*w_clone);
    h_clone->evaluate(*w_clone);
    RealType elocMinus = h_clone->getLocalEnergy()-h_clone->getLocalPotential();
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }
  fout<<"CLONE"<< std::endl;
  fout<< std::endl<<"   Deriv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i])/PGradient[i] << std::endl;
  fout<< std::endl<<"   Hderiv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i])/HGradient[i] << std::endl;
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++)
      wfvar_prime[j]=wfVars[j];
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiPlus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocPlus=H.getLocalEnergy()-H.getLocalPotential();
    wfvar_prime[i] = wfVars[i]- FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiMinus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocMinus = H.getLocalEnergy()-H.getLocalPotential();
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }
  fout<<"ORIGINAL"<< std::endl;
  fout<< std::endl<<"   Deriv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i])/PGradient[i] << std::endl;
  fout<< std::endl<<"   Hderiv  Numeric Analytic"<< std::endl;
  for (int i=0; i<Nvars ; i++)
    fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i])/HGradient[i] << std::endl;
}
void WaveFunctionTester::runwftricks()
{
  std::vector<OrbitalBase*>& Orbitals=Psi.getOrbitals();
  app_log()<<" Total of "<<Orbitals.size()<<" orbitals."<< std::endl;
  int SDindex(0);
  for (int i=0; i<Orbitals.size(); i++)
    if ("SlaterDet"==Orbitals[i]->OrbitalName)
      SDindex=i;
  SPOSetBasePtr Phi= dynamic_cast<SlaterDet *>(Orbitals[SDindex])->getPhi();
  int NumOrbitals=Phi->getBasisSetSize();
  app_log()<<"Basis set size: "<<NumOrbitals<< std::endl;
  std::vector<int> SPONumbers(0,0);
  std::vector<int> irrepRotations(0,0);
  std::vector<int> Grid(0,0);
  xmlNodePtr kids=myNode->children;
  std::string doProj("yes");
  std::string doRotate("yes");
  std::string sClass("C2V");
  ParameterSet aAttrib;
  aAttrib.add(doProj,"projection","string");
  aAttrib.add(doRotate,"rotate","string");
  aAttrib.put(myNode);
  while(kids != NULL)
  {
    std::string cname((const char*)(kids->name));
    if(cname == "orbitals")
    {
      putContent(SPONumbers,kids);
    }
    else
      if(cname == "representations")
      {
        putContent(irrepRotations,kids);
      }
      else
        if(cname=="grid")
          putContent(Grid,kids);
    kids=kids->next;
  }
  ParticleSet::ParticlePos_t R_cart(1);
  R_cart.setUnit(PosUnit::CartesianUnit);
  ParticleSet::ParticlePos_t R_unit(1);
  R_unit.setUnit(PosUnit::LatticeUnit);
//       app_log()<<" My crystals basis set is:"<< std::endl;
//       std::vector<std::vector<RealType> > BasisMatrix(3, std::vector<RealType>(3,0.0));
//
//       for (int i=0;i<3;i++)
//       {
//         R_unit[0][0]=0;
//         R_unit[0][1]=0;
//         R_unit[0][2]=0;
//         R_unit[0][i]=1;
//         W.convert2Cart(R_unit,R_cart);
//         app_log()<<"basis_"<<i<<":  ("<<R_cart[0][0]<<", "<<R_cart[0][1]<<", "<<R_cart[0][2]<<")"<< std::endl;
//         for (int j=0;j<3;j++) BasisMatrix[j][i]=R_cart[0][j];
//       }
  int Nrotated(SPONumbers.size());
  app_log()<<" Projected orbitals: ";
  for(int i=0; i<Nrotated; i++)
    app_log()<< SPONumbers[i] <<" ";
  app_log()<< std::endl;
  //indexing trick
//       for(int i=0;i<Nrotated;i++) SPONumbers[i]-=1;
  SymmetryBuilder SO;
  SO.put(myNode);
  SymmetryGroup symOp(*SO.getSymmetryGroup());
//       SO.changeBasis(InverseBasisMatrix);
  OrbitalSetTraits<ValueType>::ValueVector_t values;
  values.resize(NumOrbitals);
  RealType overG0(1.0/Grid[0]);
  RealType overG1(1.0/Grid[1]);
  RealType overG2(1.0/Grid[2]);
  RealType overNpoints=  overG0*overG1*overG2;
  std::vector<RealType> NormPhi(Nrotated, 0.0);
  int totsymops = symOp.getSymmetriesSize();
  Matrix<RealType> SymmetryOrbitalValues;
  SymmetryOrbitalValues.resize(Nrotated,totsymops);
  int ctabledim = symOp.getClassesSize();
  Matrix<double> projs(Nrotated,ctabledim);
  Matrix<double> orthoProjs(Nrotated,Nrotated);
  std::vector<RealType> brokenSymmetryCharacter(totsymops);
  for(int k=0; k<Nrotated; k++)
    for(int l=0; l<totsymops; l++)
      brokenSymmetryCharacter[l] += irrepRotations[k]*symOp.getsymmetryCharacter(l,irrepRotations[k]-1);
//       app_log()<<"bsc: ";
//       for(int l=0;l<totsymops;l++) app_log()<<brokenSymmetryCharacter[l]<<" ";
//       app_log()<< std::endl;
//       for(int l=0;l<totsymops;l++) brokenSymmetryCharacter[l]+=0.5;
  if ((doProj=="yes")||(doRotate=="yes"))
  {
    OrbitalSetTraits<ValueType>::ValueVector_t identityValues(values.size());
    //Loop over grid
    for(int i=0; i<Grid[0]; i++)
      for(int j=0; j<Grid[1]; j++)
        for(int k=0; k<Grid[2]; k++)
        {
          //Loop over symmetry classes and small group operators
          for(int l=0; l<totsymops; l++)
          {
            R_unit[0][0]=overG0*RealType(i);// R_cart[0][0]=0;
            R_unit[0][1]=overG1*RealType(j);// R_cart[0][1]=0;
            R_unit[0][2]=overG2*RealType(k);// R_cart[0][2]=0;
//                 for(int a=0; a<3; a++) for(int b=0;b<3;b++) R_cart[0][a]+=BasisMatrix[a][b]*R_unit[0][b];
            W.convert2Cart(R_unit,R_cart);
            symOp.TransformSinglePosition(R_cart,l);
            W.R[0]=R_cart[0];
            values=0.0;
            //evaluate orbitals
//                 Phi->evaluate(W,0,values);
            Psi.evaluateLog(W);
            // YYYY: is the following two lines still maintained?
            //for(int n=0; n<NumOrbitals; n++)
            //  values[n] = Phi->t_logpsi(0,n);
            if (l==0)
            {
              identityValues=values;
#if defined(QMC_COMPLEX)
              for(int n=0; n<Nrotated; n++)
                NormPhi[n] += totsymops*real(values[SPONumbers[n]]*values[SPONumbers[n]]);
#else
              for(int n=0; n<Nrotated; n++)
                NormPhi[n] += totsymops*(values[SPONumbers[n]]*values[SPONumbers[n]]);
#endif
            }
            //now we have phi evaluated at the rotated/inverted/whichever coordinates
            for(int n=0; n<Nrotated; n++)
            {
              int N=SPONumbers[n];
#if defined(QMC_COMPLEX)
              RealType phi2 = real(values[N]*identityValues[N]);
#else
              RealType phi2 = (values[N]*identityValues[N]);
#endif
              SymmetryOrbitalValues(n,l) += phi2;
            }
            for(int n=0; n<Nrotated; n++)
              for(int p=0; p<Nrotated; p++)
              {
                int N=SPONumbers[n];
                int P=SPONumbers[p];
#if defined(QMC_COMPLEX)
                orthoProjs(n,p) += 0.5*real(identityValues[N]*values[P]+identityValues[P]*values[N])*brokenSymmetryCharacter[l];
#else
                orthoProjs(n,p) +=0.5*(identityValues[N]*values[P]+identityValues[P]*values[N])*brokenSymmetryCharacter[l];
#endif
              }
          }
        }
    for(int n=0; n<Nrotated; n++)
      for(int l=0; l<totsymops; l++)
        SymmetryOrbitalValues(n,l)/= NormPhi[n];
    for(int n=0; n<Nrotated; n++)
      for(int l=0; l<Nrotated; l++)
        orthoProjs(n,l) /= std::sqrt(NormPhi[n]*NormPhi[l]);
//       if (true){
//         app_log()<< std::endl;
//         for(int n=0;n<Nrotated;n++) {
//           for(int l=0;l<totsymops;l++) app_log()<<SymmetryOrbitalValues(n,l)<<" ";
//           app_log()<< std::endl;
//         }
//       app_log()<< std::endl;
//       }
    for(int n=0; n<Nrotated; n++)
    {
      if (false)
        app_log()<<" orbital #"<<SPONumbers[n]<< std::endl;
      for(int i=0; i<ctabledim; i++)
      {
        double proj(0);
        for(int j=0; j<totsymops; j++)
          proj+=symOp.getsymmetryCharacter(j,i)*SymmetryOrbitalValues(n,j);
        if (false)
          app_log()<<"  Rep "<<i<< ": "<<proj;
        projs(n,i)=proj<1e-4?0:proj;
      }
      if (false)
        app_log()<< std::endl;
    }
    if (true)
    {
      app_log()<<"Printing Projection Matrix"<< std::endl;
      for(int n=0; n<Nrotated; n++)
      {
        for(int l=0; l<ctabledim; l++)
          app_log()<<projs(n,l)<<" ";
        app_log()<< std::endl;
      }
      app_log()<< std::endl;
    }
    if (true)
    {
      app_log()<<"Printing Coefficient Matrix"<< std::endl;
      for(int n=0; n<Nrotated; n++)
      {
        for(int l=0; l<ctabledim; l++)
          app_log()<<std::sqrt(projs(n,l))<<" ";
        app_log()<< std::endl;
      }
      app_log()<< std::endl;
    }
    if (doRotate=="yes")
    {
//         app_log()<<"Printing Broken Symmetry Projection Matrix"<< std::endl;
//           for(int n=0;n<Nrotated;n++) {
//             for(int l=0;l<Nrotated;l++) app_log()<<orthoProjs(n,l)<<" ";
//             app_log()<< std::endl;
//           }
      char JOBU('A');
      char JOBVT('A');
      int vdim=Nrotated;
      Vector<double> Sigma(vdim);
      Matrix<double> U(vdim,vdim);
      Matrix<double> VT(vdim,vdim);
      int lwork=8*Nrotated;
      std::vector<double> work(lwork,0);
      int info(0);
      dgesvd(&JOBU, &JOBVT, &vdim, &vdim,
             orthoProjs.data(), &vdim, Sigma.data(), U.data(),
             &vdim, VT.data(), &vdim, &(work[0]),
             &lwork, &info);
      app_log()<<"Printing Rotation Matrix"<< std::endl;
      for(int n=0; n<vdim; n++)
      {
        for(int l=0; l<vdim; l++)
          app_log()<<VT(l,n)<<" ";
        app_log()<< std::endl;
      }
      app_log()<< std::endl<<"Printing Eigenvalues"<< std::endl;
      for(int n=0; n<vdim; n++)
        app_log()<<Sigma[n]<<" ";
      app_log()<< std::endl;
    }
  }
}

void  WaveFunctionTester::runNodePlot()
{
  xmlNodePtr kids=myNode->children;
  std::string doEnergy("no");
  ParameterSet aAttrib;
  aAttrib.add(doEnergy,"energy","string");
  aAttrib.put(myNode);
  std::vector<int> Grid;
  while(kids != NULL)
  {
    std::string cname((const char*)(kids->name));
    if(cname=="grid")
      putContent(Grid,kids);
    kids=kids->next;
  }
  ParticleSet::ParticlePos_t R_cart(1);
  R_cart.setUnit(PosUnit::CartesianUnit);
  ParticleSet::ParticlePos_t R_unit(1);
  R_unit.setUnit(PosUnit::LatticeUnit);
  Walker_t& thisWalker(**(W.begin()));
  W.loadWalker(thisWalker,true);
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  Psi.copyFromBuffer(W,w_buffer);
#if OHMMS_DIM==2
  assert(Grid.size()==2);
  char fname[16];
//       sprintf(fname,"loc.xy");
//       std::ofstream e_out(fname);
//       e_out.precision(6);
//       e_out<<"#e  x  y"<< std::endl;
  int nat = W.getTotalNum();
  int nup = W.getTotalNum()/2;//std::max(W.getSpeciesSet().findSpecies("u"),W.getSpeciesSet().findSpecies("d"));
//       for(int iat(0);iat<nat;iat++)
//         e_out<<iat<<" "<<W[0]->R[iat][0]<<" "<<W[0]->R[iat][1]<< std::endl;
  RealType overG0(1.0/Grid[0]);
  RealType overG1(1.0/Grid[1]);
  for(int iat(0); iat<nat; iat++)
  {
    W.update();
    std::stringstream fn;
    fn<<RootName.c_str()<<".ratios."<<iat<<".py";
    std::ofstream plot_out(fn.str().c_str());
    plot_out.precision(6);
//         plot_out<<"#e  x  y  ratio"<< std::endl;
    R_unit[0][0]=1.0;
    R_unit[0][1]=1.0;
    W.convert2Cart(R_unit,R_cart);
    RealType xmax=R_cart[0][0];
    RealType ymax=R_cart[0][1];
    plot_out<<"import matplotlib\n";
    plot_out<<"import numpy as np\n";
    plot_out<<"import matplotlib.cm as cm\n";
    plot_out<<"import matplotlib.mlab as mlab\n";
    plot_out<<"import matplotlib.pyplot as plt\n";
    plot_out<< std::endl;
    plot_out<<"matplotlib.rcParams['xtick.direction'] = 'out'\n";
    plot_out<<"matplotlib.rcParams['ytick.direction'] = 'out'\n";
    plot_out<< std::endl;
    plot_out<<"x = np.arange(0, "<<xmax<<", "<< xmax*overG0 <<")\n";
    plot_out<<"y = np.arange(0, "<<ymax<<", "<< ymax*overG1 <<")\n";
    plot_out<<"X, Y = np.meshgrid(x, y)\n";
    plot_out<<"Z = [";
    for(int i=0; i<Grid[0]; i++)
    {
      plot_out<<"[ ";
      for(int j=0; j<Grid[1]; j++)
      {
        R_unit[0][0]=overG0*RealType(i);
        R_unit[0][1]=overG1*RealType(j);
        W.convert2Cart(R_unit,R_cart);
        PosType dr(R_cart[0]-W.R[iat]);
        W.makeMove(iat,dr);
        RealType aratio = Psi.ratio(W,iat);
        W.rejectMove(iat);
        Psi.rejectMove(iat);
        plot_out<<aratio<<", ";
      }
      plot_out<<"], ";
    }
    plot_out<<"]\n";
    plot_out<<"up_y=[";
    for(int ix(0); ix<nup; ix++)
    {
      RealType yy(W[0]->R[ix][0]);
      while (yy>xmax)
        yy-=xmax;
      while (yy<0)
        yy+=xmax;
      plot_out<<yy<<", ";
    }
    plot_out<<"]\n";
    plot_out<<"up_x=[";
    for(int ix(0); ix<nup; ix++)
    {
      RealType yy(W[0]->R[ix][1]);
      while (yy>ymax)
        yy-=ymax;
      while (yy<0)
        yy+=ymax;
      plot_out<<yy<<", ";
    }
    plot_out<<"]\n";
    plot_out<<"dn_y=[";
    for(int ix(nup); ix<nat; ix++)
    {
      RealType yy(W[0]->R[ix][0]);
      while (yy>xmax)
        yy-=xmax;
      while (yy<0)
        yy+=xmax;
      plot_out<<yy<<", ";
    }
    plot_out<<"]\n";
    plot_out<<"dn_x=[";
    for(int ix(nup); ix<nat; ix++)
    {
      RealType yy(W[0]->R[ix][1]);
      while (yy>ymax)
        yy-=ymax;
      while (yy<0)
        yy+=ymax;
      plot_out<<yy<<", ";
    }
    plot_out<<"]\n";
    plot_out<<"matplotlib.rcParams['contour.negative_linestyle'] = 'solid'\n";
    plot_out<<"plt.figure()\n";
    plot_out<<"CS = plt.contourf(X, Y, Z, 5, cmap=cm.gray)\n";
    plot_out<<"CS2 = plt.contour(X, Y, Z, colors='k',levels=[0])\n";
    plot_out<<"PTu = plt.scatter(up_x,up_y, c='r', marker='o')\n";
    plot_out<<"PTd = plt.scatter(dn_x,dn_y, c='b', marker='d')\n";
    plot_out<<"plt.clabel(CS2, fontsize=9, inline=1)\n";
    plot_out<<"plt.title('2D Nodal Structure')\n";
    plot_out<<"plt.xlim(0,"<< ymax*(1.0-overG1) <<")\n";
    plot_out<<"plt.ylim(0,"<< xmax*(1.0-overG0) <<")\n";
    fn.str("");
    fn<<RootName.c_str()<<".ratios."<<iat<<".png";
    plot_out<<"plt.savefig('"<<fn.str().c_str()<<"', bbox_inches='tight', pad_inches=0.01 )\n";
  }
#elif OHMMS_DIM==3
//         assert(Grid.size()==3);
//
//         RealType overG0(1.0/Grid[0]);
//         RealType overG1(1.0/Grid[1]);
//         RealType overG2(1.0/Grid[2]);
//         int iat(0);
//         W.update();
//         plot_out<<"#e  x  y  z  ratio"<< std::endl;
//
//         for(int i=0;i<Grid[0];i++)
//           for(int j=0;j<Grid[1];j++)
//           for(int k=0;k<Grid[2];k++)
//           {
//             R_unit[iat][0]=overG0*RealType(i);
//             R_unit[iat][1]=overG1*RealType(j);
//             R_unit[iat][2]=overG2*RealType(k);
//             W.convert2Cart(R_unit,R_cart);
//             PosType dr(R_cart[iat]-W.R[iat]);
//
//             W.makeMove(iat,dr);
//             RealType aratio = Psi.ratio(W,iat);
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             plot_out<<iat<<" "<<R_cart[iat][0]<<" "<<R_cart[iat][1]<<" "<<R_cart[iat][2]<<" "<<aratio<<" "<< std::endl;
//           }
#endif
}

FiniteDiffErrData::FiniteDiffErrData() : particleIndex(0), gradientComponentIndex(0), outputFile("delta.dat") {}

bool
FiniteDiffErrData::put(xmlNodePtr q)
{
  ParameterSet param;
  param.add(outputFile,"file","string");
  param.add(particleIndex,"particle_index","none");
  param.add(gradientComponentIndex,"gradient_index","none");
  bool s=param.put(q);
  return s;
}



}

