//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Estimators/CSEnergyEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus {

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  CSEnergyEstimator::CSEnergyEstimator(QMCHamiltonian& h, int hcopy) 
  {
    int NumObservables = h.sizeOfObservables();

    NumCopies=hcopy;
    FirstHamiltonian = h.startIndex();
    LastHamiltonian = FirstHamiltonian+NumObservables;

    //add names
    h_components.push_back("LocEne");
    h_components.push_back("LocPot");
    for(int i=0; i<NumObservables; ++i) 
      h_components.push_back(h.getObservableName(i));

    scalars.resize(NumCopies  + 
		   h_components.size()*(NumCopies+NumCopies*(NumCopies-1)/2));
    scalars_saved.resize(scalars.size());
  }

  ScalarEstimatorBase* CSEnergyEstimator::clone()
  {
    return new CSEnergyEstimator(*this);
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void 
  CSEnergyEstimator::add2Record(RecordNamedProperty<RealType>& record) {
    char aname[80];
    FirstIndex = record.size();
    
    for(int i=0; i<NumCopies; ++i)
    {
      for(int k=0; k<h_components.size(); ++k)
      {
        sprintf(aname,"%s_%i",h_components[k].c_str(),i);   
        int dummy=record.add(aname);
      }
    }

    for(int i=0; i<NumCopies; ++i)
    {
      sprintf(aname,"wpsi_%i",i);   
      int dummy=record.add(aname);
    }

    for(int i=0; i<NumCopies; i++) {
      for(int j=i+1; j<NumCopies; j++) {
        for(int k=0; k<h_components.size(); ++k)
        {
          sprintf(aname,"d%s_%d_%d",h_components[k].c_str(),i,j); 
          int dummy=record.add(aname);
        }
      }
    }


    LastIndex=record.size();
    tmp_data.resize(NumCopies,h_components.size());
    uweights.resize(NumCopies);
    clear();

    //msg.add(d_data.begin(),d_data.end());
  }

  void CSEnergyEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
  {
    //NEED TO IMPLEMENT for hdf5
  }

  void 
  CSEnergyEstimator::accumulate(const Walker_t& awalker, RealType wgt) 
  {
	std::vector<double> weightaverage(NumCopies);
    //first copy data to tmp_dat to calculate differences
    for(int i=0; i<NumCopies; i++) 
    {
      const RealType* restrict prop=awalker.getPropertyBase(i);
      RealType* restrict prop_saved=tmp_data[i];
      uweights[i]=prop[UMBRELLAWEIGHT];
      *prop_saved++=prop[LOCALENERGY];
      *prop_saved++=prop[LOCALPOTENTIAL];
      std::copy(prop+FirstHamiltonian,prop+LastHamiltonian,prop_saved);
    }

    int ii=0;
    const RealType *hptr=tmp_data.data();
    for(int i=0; i<NumCopies; i++) 
    {
      RealType uw=uweights[i];
      for(int k=0; k<tmp_data.cols(); ++k) scalars[ii++](*hptr++,uw);
    }

    for(int i=0; i<NumCopies; i++) 
    {
      scalars[ii++](uweights[i],1.0);
    }

    int ii_i(0);  //index of observable ii for copy i
    int ii_j(0);  //index of observable ii for copy j
    for(int i=0; i<NumCopies; i++) 
    {
      for(int j=i+1; j<NumCopies; j++)
      {
        for(int k=0; k<tmp_data.cols(); ++k)
        {
          ii_i=i*tmp_data.cols()+k;
          ii_j=j*tmp_data.cols()+k; 
	  
          //reset is used here because we do no accumulating.  Just overwrite.
          scalars[ii++].reset(scalars[ii_i].mean()-scalars[ii_j].mean(),1.0);
	}
      }
    }
  }

}
