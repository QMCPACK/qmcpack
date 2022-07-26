//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////


#include "MagDensityEstimator.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/LRCoulombSingleton.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/DistanceTable.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/ParticleBase/RandomSeqGeneratorGlobal.h"
#include "Utilities/ModernStringUtils.hpp"

namespace qmcplusplus
{
using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
using GridType       = LRCoulombSingleton::GridType;
using RadFunctorType = LRCoulombSingleton::RadFunctorType;

MagDensityEstimator::MagDensityEstimator(ParticleSet& elns, TrialWaveFunction& psi):
refPsi(psi),nSamples_(9),integrator_(MD_INT_SIMPSONS)
{
  update_mode_.set(COLLECTABLE, 1);
  Periodic = (elns.getLattice().SuperCellEnum != SUPERCELL_OPEN);
  for (int dim = 0; dim < OHMMS_DIM; ++dim)
  {
    density_max[dim] = elns.getLattice().Length[dim];
    ScaleFactor[dim] = 1.0 / elns.getLattice().Length[dim];
  }
}

void MagDensityEstimator::resetTargetParticleSet(ParticleSet& P) {}

ParticleSet::ParticleScalar MagDensityEstimator::generateUniformGrid(RealType start, RealType stop, int nGridPoints)
{
  ParticleSet::ParticleScalar sgrid;
  sgrid.resize(nGridPoints);
  RealType delta=(stop-start)/(nGridPoints-1.0);
  for(int i=0; i<nGridPoints; i++)
    sgrid[i]=start+i*delta;
  return sgrid;
}
    
ParticleSet::ParticleScalar MagDensityEstimator::generateRandomGrid(RealType start, RealType stop, int nSamples)
{
  ParticleSet::ParticleScalar sgrid;
  sgrid.resize(nSamples);
  makeUniformRandom(sgrid); //returns random number between [0,1)
  //scale and shift the random grid to be between [start,stop)
  sgrid=(stop-start)*sgrid;  
  sgrid=sgrid+start;
  return sgrid;
}

MagDensityEstimator::RealType MagDensityEstimator::integrateBySimpsonsRule(const std::vector<RealType>& fgrid, RealType gridDx)
{
  RealType sint(0.0);
  int gridsize=fgrid.size();
  for (int is = 1; is < gridsize - 1; is += 2)
  {
    sint += RealType(4. / 3.) * gridDx * fgrid[is];
    app_log()<<"A is="<<is<<" s="<<gridDx<<std::endl;
  }
  for (int is = 2; is < gridsize - 1; is += 2)
  {
    sint += RealType(2. / 3.) * gridDx * fgrid[is];
    app_log()<<"B is="<<is<<" s="<<gridDx<<std::endl;
  }
  sint += RealType(1. / 3.) * gridDx * fgrid[0];
  sint += RealType(1. / 3.) * gridDx * fgrid[gridsize-1];
  app_log()<<"Final indices are 0 and "<<gridsize-1<<" at 0 and "<<(gridsize-1)*gridDx<<std::endl;
  return sint;
}

MagDensityEstimator::RealType MagDensityEstimator::integrateByTrapzRule(const std::vector<RealType>& fgrid, RealType gridDx)
{
  RealType sint(0.0);
  int gridsize=fgrid.size();
  for (int is = 1; is < gridsize - 1; is++)
  {
    sint +=  gridDx * fgrid[is];
    app_log()<<"A is="<<is<<" s="<<gridDx<<std::endl;
  }
  sint += 0.5 * gridDx * fgrid[0];
  sint += 0.5 * gridDx * fgrid[gridsize-1];
  app_log()<<"Final indices are 0 and "<<gridsize-1<<" at 0 and "<<(gridsize-1)*gridDx<<std::endl;
  return sint;
}

MagDensityEstimator::RealType MagDensityEstimator::average(const std::vector<RealType>& fgrid)
{
  int gridsize=fgrid.size();
  RealType sum(0.0);
  for (int i = 0; i < gridsize; i++)
    sum+=fgrid[i];

  return sum/RealType(gridsize);
}

MagDensityEstimator::Return_t MagDensityEstimator::evaluate(ParticleSet& P)
{
  ParticleSet::ParticleScalar sgrid;
  sgrid.resize(nSamples_);
  std::vector<RealType> sxgrid;
  std::vector<RealType> sygrid;
  std::vector<RealType> szgrid;
  sxgrid.resize(nSamples_);
  sygrid.resize(nSamples_);
  szgrid.resize(nSamples_);
  const PosType dr(0.0); //Integration over spin variable doesn't change particle position.  
			 //This is the argument for the calcRatio calls, which require both.

  std::complex<RealType> eye(0,1.0);
  RealType wgt = t_walker_->Weight;

  app_log()<<"wgt = "<<wgt<<std::endl;
  app_log()<<" Periodic = "<<Periodic<<std::endl;
  if (Periodic)
  {
    for (int ig = 0; ig < P.groups(); ++ig)
    {
      for (int iat = P.first(ig); iat < P.last(ig); ++iat)  
      {
        PosType ru;
        ru    = P.getLattice().toUnit(P.R[iat]);
        int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
        int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
        int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
        app_log()<<"P.R[iat] = "<<P.R[iat]<<std::endl;  
        //This has to be complex type for spinors.  If not true, other part of code will
        //handle it.
        ValueType sx(0.0);
        ValueType sy(0.0);
        ValueType sz(0.0);
        
        //If Monte Carlo, we generate a random s grid between 0 and 2pi.,     
        if( integrator_ == MD_INT_MC )
          sgrid=generateRandomGrid(0,TWOPI,nSamples_);
        else //Otherwise evenly divide the interval.
          sgrid=generateUniformGrid(0,TWOPI,nSamples_);
        
        for(int samp=0; samp<nSamples_; samp++)
        {
          ValueType ratio;
          RealType ds=sgrid[samp]-P.spins[iat]; 
          P.makeMoveWithSpin(iat, 0.0, ds);
          ratio = refPsi.calcRatio(P, iat);
          P.rejectMove(iat); //reject the move
          refPsi.resetPhaseDiff();
          sxgrid[samp] = std::real(2.0*std::cos(sgrid[samp]+P.spins[iat])*ratio);
          sygrid[samp] = std::real(2.0*std::sin(sgrid[samp]+P.spins[iat])*ratio);
          szgrid[samp] = std::real(-2.0*eye*std::sin(ds)*ratio);
        }
        if(integrator_ == MD_INT_MC )
        {
          //Monte Carlo, so the 2pi normalization is already handled.
          sx=average(sxgrid);
          sy=average(sygrid);
          sz=average(szgrid);
        }
        else
        {
          RealType gridDs=TWOPI/(nSamples_ - 1.0);
        
          sx=integrateBySimpsonsRule(sxgrid,gridDs)/TWOPI;
          sy=integrateBySimpsonsRule(sygrid,gridDs)/TWOPI;
          sz=integrateBySimpsonsRule(szgrid,gridDs)/TWOPI;
        }
        app_log()<<" sx="<<sx<<" sy="<<sy<<" sz="<<sz<<std::endl;
        P.Collectables[getMagGridIndex(i, j, k, 0)] += wgt*std::real(sx); //1.0;
        P.Collectables[getMagGridIndex(i, j, k, 1)] += wgt*std::real(sy); //1.0;
        P.Collectables[getMagGridIndex(i, j, k, 2)] += wgt*std::real(sz); //1.0;
     
 
      }
    }
  }
  else
  {
//    for (int iat = 0; iat < P.getTotalNum(); ++iat)
//    {
//      PosType ru;
//      for (int dim = 0; dim < OHMMS_DIM; dim++)
//      {
//        ru[dim] = (P.R[iat][dim] - density_min[dim]) * ScaleFactor[dim];
//      }
//      if (ru[0] > 0.0 && ru[1] > 0.0 && ru[2] > 0.0 && ru[0] < 1.0 && ru[1] < 1.0 && ru[2] < 1.0)
//      {
//        int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
//        int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
//        int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
//        for( int idim = 0; idim < OHMMS_DIM; idim++)
//          P.Collectables[getMagGridIndex(i, j, k, idim)] += wgt; //1.0;
//      }
//    }
  }
  return 0.0;
}

void MagDensityEstimator::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  int nw = W.WalkerList.size();
  int N  = W.getTotalNum();
  if (Periodic)
  {
    for (int iw = 0; iw < nw; iw++)
    {
      Walker_t& w     = *W.WalkerList[iw];
      RealType weight = w.Weight / nw;
      for (int iat = 0; iat < N; iat++)
      {
        PosType ru;
        ru = W.getLattice().toUnit(w.R[iat]);
        // for (int dim=0; dim<OHMMS_DIM; dim++)
        // {
        //   ru[dim]=(w.R[iat][dim]-density_min[dim])*ScaleFactor[dim];
        // }
        int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
        int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
        int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
        W.Collectables[getGridIndex(i, j, k)] += weight;
      }
    }
  }
}

void MagDensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  //current index
  my_index_ = collectables.current();
  app_log()<<"addObservables = "<<NumGrids[OHMMS_DIM]<<std::endl;
  std::vector<RealType> tmp(OHMMS_DIM*NumGrids[OHMMS_DIM]);
  collectables.add(tmp.begin(), tmp.end());
  //potentialIndex=collectables.current();
  //vector<RealType> tmp2(NumGrids[OHMMS_DIM]);
  //collectables.add(tmp2.begin(),tmp2.end());
}

void MagDensityEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  int loc = h5desc.size();
  std::vector<int> ng(OHMMS_DIM+1);
  for (int i = 0; i < OHMMS_DIM; ++i)
    ng[i] = NumGrids[i];
  ng[OHMMS_DIM]=OHMMS_DIM;
  h5desc.emplace_back(name_);
  auto& h5o = h5desc.back();
  h5o.set_dimensions(ng, my_index_);
  h5o.open(gid);
}

void MagDensityEstimator::setObservables(PropertySetType& plist)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex);
}

void MagDensityEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex+offset);
}

/** check xml elements
 *
 * <estimator name="density" debug="no" delta="0.1 0.1 0.1"/>
 */
bool MagDensityEstimator::put(xmlNodePtr cur)
{
  Delta = 0.1;
  std::vector<double> delta;
  std::string debug("no");
  std::string potential("no");
  std::string integrator_string("simpson");
  OhmmsAttributeSet attrib;
  attrib.add(debug, "debug");
  attrib.add(potential, "potential");
  attrib.add(density_min[0], "x_min");
  attrib.add(density_min[1], "y_min");
  attrib.add(density_min[2], "z_min");
  attrib.add(density_max[0], "x_max");
  attrib.add(density_max[1], "y_max");
  attrib.add(density_max[2], "z_max");
  attrib.add(nSamples_, "nsamples");
  attrib.add(integrator_string,"spin_integral");
  attrib.add(Delta, "delta");
  attrib.put(cur);

  std::string tmp = lowerCase(integrator_string);
  size_t found = tmp.find("simp");
  if (tmp.find("simp") != std::string::npos)
  {
    app_log()<<"Setting up Simpson's 3/8 spin integration rule.\n";
    if (nSamples_%2 == 0)
    {
      nSamples_++;
      app_log()<<"Simpson's rule requires odd number of integration points.  Incrementing nSamples up by one.\n";
      integrator_=MD_INT_SIMPSONS; 
    }
  }
  else if(tmp.find("trap") != std::string::npos)
  {
    app_log()<<"Setting up Trapezoidal spin integration rule.\n";
    integrator_=MD_INT_TRAP;
  }
  else if (tmp.find("mc") != std::string::npos)
  {
    app_log()<<"Setting up Monte Carlo integration rule.\n";
    integrator_=MD_INT_MC;
  }
  else
  {
    APP_ABORT("DIE");
  }
   
  //trap
  //mc
  if (!Periodic)
  {
    for (int dim = 0; dim < OHMMS_DIM; ++dim)
      ScaleFactor[dim] = 1.0 / (density_max[dim] - density_min[dim]);
  }
  resize();
  return true;
}

bool MagDensityEstimator::get(std::ostream& os) const
{
  os << name_ << " bin =" << Delta << " bohrs " << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> MagDensityEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //default constructor is sufficient
  std::unique_ptr<MagDensityEstimator> myClone = std::make_unique<MagDensityEstimator>(qp,psi);
  myClone->NumGrids = NumGrids;
  return myClone;
}

void MagDensityEstimator::resize()
{
  for (int i = 0; i < OHMMS_DIM; ++i)
  {
    DeltaInv[i] = 1.0 / Delta[i];
    NumGrids[i] = static_cast<int>(DeltaInv[i]);
    if (NumGrids[i] < 2)
    {
      APP_ABORT("MagDensityEstimator::resize invalid bin size");
    }
  }
  app_log() << " MagDensityEstimator bin_size= " << NumGrids << " delta = " << Delta << std::endl;
  NumGrids[OHMMS_DIM] = NumGrids[0] * NumGrids[1] * NumGrids[2];
}

} // namespace qmcplusplus
