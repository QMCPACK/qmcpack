//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////
#include "Estimators/MagnetizationDensity.h"
#include "Estimators/MagnetizationDensityInput.h"

namespace qmcplusplus
{

MagnetizationDensity::MagnetizationDensity(MagnetizationDensityInput&& minput,
					   const Lattice & lat):
OperatorEstBase(DataLocality::crowd), input_(minput), lattice_(lat) 
{
  //Pull consistent corner, grids, etc., from already inititalized input.  
  MagnetizationDensityInput::DerivedParameters derived = minput.calculateDerivedParameters(lat);
  npoints_ = derived.npoints;
  grid_= derived.grid;
  gdims_ = derived.gdims;

  rcorner_ = derived.corner;
  center_  = rcorner_ + lattice_.Center;

  nsamples_ = input_.get_nsamples();
  integrator_ = input_.get_integrator();
  switch(integrator_)
  {
    case Integrator::SIMPSONS:
      app_log()<<"SIMPSONS\n";
      break;
    case Integrator::MONTECARLO:
      app_log()<<"MONTECARLO\n";
      break;
  }

  //Resize the data arrays.
  data_.resize(getFullDataSize());

}

MagnetizationDensity::MagnetizationDensity(const MagnetizationDensity& magdens,DataLocality dl): MagnetizationDensity(magdens)
{
 data_locality_ = dl;
}
void MagnetizationDensity::startBlock(int steps)
{
};
   
size_t MagnetizationDensity::getFullDataSize() { return npoints_ * DIM; }

void MagnetizationDensity::accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) 
{
  app_log()<<"!!!!!!!!!!!!!!!!! in accumulate 000000####\n"; 
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker     = walkers[iw];
    ParticleSet& pset     = psets[iw];
    TrialWaveFunction& wfn   = wfns[iw];

    QMCT::RealType weight = walker.Weight;
    std::vector<Value> sxgrid(nsamples_,0.0);
    std::vector<Value> sygrid(nsamples_,0.0);
    std::vector<Value> szgrid(nsamples_,0.0);

    //Temporary variables to hold the integrated sx, sy, sz estimates.  
    //These are overwritten as needed in the loop over particles.  
    Value sx(0.0);
    Value sy(0.0);
    Value sz(0.0);

    const int np = pset.getTotalNum();
    assert(weight >= 0);
    for (int p = 0; p < np; ++p)
    {
      size_t sxindex = computeBin(pset.R[p],0);
      generateSpinIntegrand(pset,wfn,p,sxgrid,sygrid,szgrid);
      sx=integrateMagnetizationDensity(sxgrid);
      sy=integrateMagnetizationDensity(sygrid); 
      sz=integrateMagnetizationDensity(szgrid); 
     
      data_[sxindex]+=std::real(sx*weight); 
      data_[sxindex+1]+=std::real(sy*weight); 
      data_[sxindex+2]+=std::real(sz*weight); 
      
      app_log()<<"p="<<p<<std::endl;
      app_log()<<"  sx="<<std::real(sx)<<" at "<<sxindex<<std::endl;
      app_log()<<"  sy="<<std::real(sy)<<" at "<<sxindex+1<<std::endl;
      app_log()<<"  sz="<<std::real(sz)<<" at "<<sxindex+2<<std::endl;
      for(int q=0; q< sxgrid.size(); q++)
      {
        app_log()<<"  "<<sxgrid[q]<<"  "<<sygrid[q]<<"  "<<szgrid[q]<<std::endl; 
      }
    }
  }


};

size_t MagnetizationDensity::computeBin(const QMCT::PosType& r, const unsigned int component)
{
  assert(component < QMCT::DIM);  
  QMCT::PosType u = lattice_.toUnit(r - rcorner_);
  size_t point = 0; //This establishes the real space grid point.  
  for (int d = 0; d < QMCT::DIM; ++d)
    point += gdims_[d] * ((int)(grid_[d] * (u[d] - std::floor(u[d])))); //periodic only
  
  //We now have the real space grid point.  For each grid point, we store sx[point],sy[point],sz[point].
  //Thus, the actual position in array is DIM*point+component.  
  return DIM*point+component; 
}

std::unique_ptr<OperatorEstBase> MagnetizationDensity::spawnCrowdClone() const 
{ 
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  //Everyone else has this attempt to set up a non-implemented memory saving optimization.
  //We won't rock the boat.  
  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for MagnetizationDensity");
  }

  UPtr<MagnetizationDensity> spawn(std::make_unique<MagnetizationDensity>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
};

void MagnetizationDensity::generateSpinIntegrand(ParticleSet& pset_target,
                                                  TrialWaveFunction& psi_target,
                                                  const int iat,
                                                  std::vector<Value>& sxgrid,
                                                  std::vector<Value>& sygrid,
                                                  std::vector<Value>& szgrid)
                                            
{
  VirtualParticleSet vp(pset_target, nsamples_);
  std::vector<Real> sgrid(nsamples_);
  std::vector<Real> ds(nsamples_,0.0);
  std::vector<Position> dV(nsamples_,0); 
  std::vector<Value> ratios(nsamples_,0);
  generateGrid(sgrid);
  for (int samp=0; samp<nsamples_; samp++)
    ds[samp]=sgrid[samp]-pset_target.spins[iat];
 
  //This is what I would like to do.  This gives some crazy ratio values though.  Debug.
  //
  //vp.makeMovesWithSpin(pset_target,iat,dV,ds);
  //psi_target.evaluateRatios(vp,ratios);
    
  const std::complex<Real> eye(0, 1.0);
 
  for (int samp=0; samp<nsamples_; samp++)
  {
    ///This is a hacky replacement for virtual particle moves.  
    pset_target.makeMoveWithSpin(iat, 0.0, ds[samp]);
    ratios[samp] = psi_target.calcRatio(pset_target, iat);
    pset_target.rejectMove(iat); //reject the move
    psi_target.resetPhaseDiff();
    /// end 
    //
    app_log()<<"  i="<<samp<<" dV="<<dV[samp]<<" ds="<<ds[samp]<<" ratio="<<ratios[samp]<<std::endl;
    sxgrid[samp] = std::real(2.0 * std::cos(sgrid[samp] + pset_target.spins[iat]) * ratios[samp]);
    sygrid[samp] = std::real(2.0 * std::sin(sgrid[samp] + pset_target.spins[iat]) * ratios[samp]);
    szgrid[samp] = std::real(-2.0 * eye * std::sin(ds[samp]) * ratios[samp]);
  } 
  
}

void MagnetizationDensity::generateUniformGrid(std::vector<Real>& sgrid, const Real start, const Real stop)
{
  size_t num_gridpoints = sgrid.size();
  Real delta = (stop - start) / (num_gridpoints - 1.0);
  for (int i = 0; i < num_gridpoints; i++)
    sgrid[i] = start + i * delta;
}

template<typename RAN_GEN>
inline void MagnetizationDensity::generateRandomGrid(std::vector<Real>& sgrid, RAN_GEN& rng, Real start, Real stop)
{
  size_t npoints=sgrid.size();
  for(int i=0; i<npoints;i++)
    sgrid[i] = (stop-start)*rng();
}

void MagnetizationDensity::generateGrid(std::vector<Real>& sgrid)
{
  sgrid.resize(nsamples_); 
  Real start=0.0;
  Real stop=2*M_PI;

  switch(integrator_)
  {
    case Integrator::SIMPSONS:
      generateUniformGrid(sgrid,start,stop);
      break;
    case Integrator::MONTECARLO:
      throw std::runtime_error("Monte Carlo sampling not implemented yet");
      break;
  }
}
template<class VAL>
VAL MagnetizationDensity::integrateBySimpsonsRule(const std::vector<VAL>& fgrid, Real gridDx)
{
  VAL sint(0.0);
  int gridsize = fgrid.size();
  for (int is = 1; is < gridsize - 1; is += 2)
  {
    sint += (4. / 3.) * gridDx * fgrid[is];
  }
  for (int is = 2; is < gridsize - 1; is += 2)
  {
    sint += (2. / 3.) * gridDx * fgrid[is];
  }
  sint += (1. / 3.) * gridDx * fgrid[0];
  sint += (1. / 3.) * gridDx * fgrid[gridsize - 1];
  return sint;
}

template<class VAL>
VAL MagnetizationDensity::integrateMagnetizationDensity(const std::vector<VAL>& fgrid)
{
  
  Real start=0.0;
  Real stop=2*M_PI;
  Real deltax=(stop-start)/(nsamples_ - 1.0);

  VAL val=0.0;
  switch(integrator_)
  {
    case Integrator::SIMPSONS:
      //There's a normalizing 2pi factor in this integral.  Take care of it here.  
      val=integrateBySimpsonsRule(fgrid,deltax)/(2.0*M_PI);
      break;
    case Integrator::MONTECARLO:
      //Integrand has form that can be monte carlo sampled.  This means the 2*PI
      //is taken care of if the integrand is uniformly sampled between [0,2PI),
      //which it is.  The following is just an average.  
      val=std::accumulate(fgrid.begin(),fgrid.end(),VAL(0));
      val/=Real(nsamples_);
      break;
  }

  return val;
}
void MagnetizationDensity::registerOperatorEstimator(hdf_archive& file)
{
  std::vector<size_t> my_indexes;

  std::vector<int> ng(DIM*npoints_);

  hdf_path hdf_name{my_name_};
//  for (int s = 0; s < species_.size(); ++s)
//  {
  h5desc_.emplace_back(hdf_name);
  auto& oh = h5desc_.back();
  oh.set_dimensions(ng, 0);
//  }
}

} //namespace qmcplusplus
