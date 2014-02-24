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
    int NumObservables = h.size();

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
    //    cerr << "scalars.size() = " << scalars.size() << endl;
    //d_data.resize(NumCopies*3+NumCopies*(NumCopies-1)/2);
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

  void CSEnergyEstimator::registerObservables(vector<observable_helper*>& h5dec, hid_t gid)
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
      //scalars.average()
      RealType ui_avg=scalars_saved[ii-1].mean();
      if (ui_avg > 0) weightaverage[i]=ui_avg;
      else weightaverage[i]=1;
     // app_log()<<"i="<<i<<" ii="<<ii<<" uweight[i]="<<uweights[i]<<" scalars[ii]="<<scalars_saved[ii-1].result()<<" wa="<<weightaverage[i]<<endl; 
    }

    for(int i=0; i<NumCopies; i++) 
    {
      RealType ui=uweights[i];
      for(int j=i+1; j<NumCopies; j++)
      {
        RealType uj=uweights[j];
	// cerr << "ui = " << ui << "  uj = " << uj << endl;
	// cerr << "diff        = " << tmp_data(j,0) - tmp_data(i,0) << endl;
	// cerr << "LOCALENERGY = " << uj*awalker.getPropertyBase(j)[LOCALENERGY] 
	//   - ui*awalker.getPropertyBase(i)[LOCALENERGY] << endl;

        for(int k=0; k<tmp_data.cols(); ++k){
			 
          scalars[ii++]((ui*tmp_data(i,k)/weightaverage[i]-uj*tmp_data(j,k)/weightaverage[j]),1.0);
         // app_log()<<"i="<<i<<" ii="<<ii<<" wai="<<weightaverage[i]<<" waj="<<weightaverage[j]<<endl;
	  }
	  // scalars[ii++](awalker.getPropertyBase(j)[LOCALENERGY] -
	  // 		awalker.getPropertyBase(i)[LOCALENERGY],1.0);
      }
    }
    //d_data[ii++]+=uw[i]*e[i]-uw[j]*e[j];
  }

  void 
  CSEnergyEstimator::evaluateDiff() 
  {
//    int ii=0;
//    for(int i=0; i<NumCopies; i++,ii+=3) 
//    {
//      RealType r= d_wgt/d_data[ii+2];
//      d_data[ii] *= r;
//      d_data[ii+1] *= r;
//    }
//
//    //d_wgt=1.0;
//    for(int i=0; i<NumCopies-1; i++) 
//      for(int j=i+1; j<NumCopies; j++)
//        d_data[ii++]+=d_data[j*3]-d_data[i*3];
  }

//  ///Set CurrentWalker to zero so that accumulation is done in a vectorized way
//  void CSEnergyEstimator::reset() 
//  {
//    //d_wgt=0.0;
//    //std::fill(d_data.begin(), d_data.end(),0.0);
//  }
//
//  void CSEnergyEstimator::report(RecordNamedProperty<RealType>& record, RealType wgtinv)
//  {
//  }
//
//  /** calculate the averages and reset to zero
//   * @param record a container class for storing scalar records (name,value)
//   * @param wgtinv the inverse weight
//   *
//   * Disable collection. CSEnergyEstimator does not need to communiate at all.
//   */
//  void CSEnergyEstimator::report(RecordNamedProperty<RealType>& record, 
//      RealType wgtinv, BufferType& msg) 
//  {
//  //  msg.get(d_data.begin(),d_data.end());
//  //  report(record,wgtinv);
//  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: CSEnergyEstimator.cpp 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
