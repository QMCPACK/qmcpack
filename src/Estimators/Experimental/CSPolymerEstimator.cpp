//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
#include "Estimators/CSPolymerEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus
{

/** constructor
 * @param h QMCHamiltonian to define the components
 * @param hcopy number of copies of QMCHamiltonians
 */
CSPolymerEstimator::CSPolymerEstimator(QMCHamiltonian& h, int hcopy,
                                       MultiChain* polymer)
//   : Reptile(polymer)
{
  Reptile = polymer;
  NumCopies=hcopy;
  NumObservables = h.size();
  scalars.resize(NumCopies+NumCopies*(NumCopies-1)/2);
  scalars_saved=scalars;
  //d_data.resize(NumCopies*3+NumCopies*(NumCopies-1)/2);
}

CSPolymerEstimator::CSPolymerEstimator(const CSPolymerEstimator& mest):
  PolymerEstimator(mest)
//       Reptile(mest.Reptile)
{
  Reptile=mest.Reptile;
  NumCopies=mest.NumCopies;
  //d_data.resize(mest.d_data.size());
}

ScalarEstimatorBase* CSPolymerEstimator::clone()
{
  return new CSPolymerEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 *@param record storage of scalar records (name,value)
 */
void
CSPolymerEstimator::add2Record(RecordNamedProperty<RealType>& record)
{
  FirstIndex = record.add("LE0");
  int dummy=record.add("LESQ0");
  dummy=record.add("WPsi0");
  char aname[32];
  for(int i=1; i<NumCopies; i++)
  {
    sprintf(aname,"LE%i",i);
    dummy=record.add(aname);
    sprintf(aname,"LESQ%i",i);
    dummy=record.add(aname);
    sprintf(aname,"WPsi%i",i);
    dummy=record.add(aname);
  }
  for(int i=0; i<NumCopies-1; i++)
  {
    for(int j=i+1; j<NumCopies; j++)
    {
      sprintf(aname,"DiffS%iS%i",i,j);
      dummy=record.add(aname);
    }
  }
  //msg.add(d_data.begin(),d_data.end());
}

void
CSPolymerEstimator::accumulate(const MCWalkerConfiguration& W
                               , WalkerIterator first, WalkerIterator last, RealType wgt)
{
  //Directionless=2
  for(int i=0,ii=0; i<NumCopies; i++)
  {
    RealType uw(Reptile->UmbrellaWeight[i]);
    RealType ehead=Reptile->front()->Action(i,2);
    RealType etail=Reptile->back()->Action(i,2);
    scalars[ii]((ehead+etail)*OneOverTau,uw);
    //d_data[ii++] += uw*(ehead+etail)*OneOverTau;
    //d_data[ii++] += uw*(ehead*ehead+etail*etail)*OneOverTau;
    //d_data[ii++] += uw;
  }
  //d_wgt += 1.0;
  //TinyVector<RealType,4> e,uw;
  //for(int i=0; i<NumCopies; i++)
  //{
  //  //get the pointer to the i-th row
  //  const RealType* restrict prop=awalker.getPropertyBase(i);
  //  uw[i] = prop[UMBRELLAWEIGHT];
  //  e[i] = prop[LOCALENERGY];
  //  d_data[ii++]+=uw[i]*e[i];
  //  d_data[ii++]+=uw[i]*e[i]*e[i];
  //  d_data[ii++]+=uw[i];
  //}
  //for(int i=0; i<NumCopies-1; i++)
  //  for(int j=i+1; j<NumCopies; j++)
  //    d_data[ii++]+=uw[i]*e[i]-uw[j]*e[j];
}

void CSPolymerEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}
