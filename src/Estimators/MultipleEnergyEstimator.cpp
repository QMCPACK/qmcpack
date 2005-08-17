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
#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/DistanceTable.h"

namespace ohmmsqmc {

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  MultipleEnergyEstimator::MultipleEnergyEstimator(QMCHamiltonian& h, int hcopy) : CurrentWalker(0) {

    NumCopies=hcopy;
    NumOperators = h.size();
    FirstHamiltonian=h.startIndex();

    esum.resize(NumCopies,LE_INDEX);
    esum_name.resize(LE_INDEX,NumCopies);

    elocal.resize(NumCopies,NumOperators);
    elocal_name.resize(NumCopies,NumOperators);

    char aname[32];
    //make the name tables
    //(localenergy, variance, weight)* + (each hamiltonian term)*
    for(int i=0; i<NumCopies; i++) {
      sprintf(aname,"LE%i",i);   esum_name(ENERGY_INDEX,i)=aname;
      sprintf(aname,"Var%i",i);  esum_name(ENERGY_SQ_INDEX,i)=aname;
      sprintf(aname,"WPsi%i",i); esum_name(WEIGHT_INDEX,i)=aname;

      for(int j=0; j<NumOperators; j++) {
        sprintf(aname,"%s%i",h.getName(j).c_str(),i);
        elocal_name(i,j)=aname;
      }
    }
  }

  void 
  MultipleEnergyEstimator
  ::initialize(MCWalkerConfiguration& W, 
      vector<QMCHamiltonian*>& h, 
      vector<TrialWaveFunction*>& psi,
      RealType tau,
      bool require_register) {

    NumWalkers = W.getActiveWalkers();

    //allocate UmbrellaEnergy
    int numPtcls(W.getTotalNum());

    //UmbrellaEnergy.resize(NumWalkers,NumCopies);
    //UmbrellaWeight.resize(NumWalkers,NumCopies);
    RatioIJ.resize(NumWalkers,NumCopies*(NumCopies-1)/2);

    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 

    vector<RealType> sumratio(NumCopies), logpsi(NumCopies);

    int iw(0);
    while(it != it_end) {

      Walker_t& thisWalker(**it);

      if(require_register) {
        (*it)->DataSet.rewind();
        W.registerData(thisWalker,(*it)->DataSet);
        //W.DataSet[iw]->rewind();
        //W.registerData(thisWalker,*(W.DataSet[iw]));
      } else {
        W.R = thisWalker.R;
        DistanceTable::update(W);
      }

      //evalaute the wavefunction and hamiltonian
      for(int ipsi=0; ipsi< NumCopies;ipsi++) {			  
        psi[ipsi]->G.resize(numPtcls);
        psi[ipsi]->L.resize(numPtcls);
        //Need to modify the return value of OrbitalBase::registerData
        if(require_register) {
	  logpsi[ipsi]=psi[ipsi]->registerData(W,(*it)->DataSet);
	  //logpsi[ipsi]=psi[ipsi]->registerData(W,*(W.DataSet[iw]));
        } else {
	  logpsi[ipsi]=psi[ipsi]->evaluateLog(W); 		 
        }
        psi[ipsi]->G=W.G;
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=h[ipsi]->evaluate(W);
        h[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
        sumratio[ipsi]=1.0;
      } 							
     
      //Check SIMONE's note
      //Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
      int indexij(0);
      RealType *rPtr=RatioIJ[iw];
      for(int ipsi=0; ipsi< NumCopies-1; ipsi++) {			  
        for(int jpsi=ipsi+1; jpsi< NumCopies; jpsi++){     		 
          RealType r=exp(2.0*(logpsi[jpsi]-logpsi[ipsi])); 
          rPtr[indexij++]=r;
          sumratio[ipsi] += r;                            
          sumratio[jpsi] += 1.0/r;		
        }                                              
      }                                               

      //Re-use Multiplicity as the sumratio
      thisWalker.Multiplicity=sumratio[0];

      //DON't forget DRIFT!!!
      thisWalker.Drift=0.0;

      for(int ipsi=0; ipsi< NumCopies; ipsi++) {
        RealType wgt=1.0/sumratio[ipsi];
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=wgt;
        thisWalker.Drift += wgt*psi[ipsi]->G;
      }
      thisWalker.Drift *= tau;
      ++it;++iw;
    }
  }

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     *@param record storage of scalar records (name,value)
     */
  void 
  MultipleEnergyEstimator::add2Record(RecordNamedProperty<RealType>& record) {
    FirstColumnIndex = record.add(esum_name(0).c_str());
    for(int i=1; i<esum_name.size(); i++) record.add(esum_name(i).c_str());
    for(int i=0; i<elocal_name.size(); i++) record.add(elocal_name(i).c_str());
  }

  void 
  MultipleEnergyEstimator::accumulate(const Walker_t& awalker, RealType wgt) {

    //const RealType* restrict etot=UmbrellaEnergy[CurrentWalker];
    //const RealType* restrict wsum=UmbrellaWeight[CurrentWalker];
    int iloc=0;
    RealType *restrict esumPtr = esum.data();

    for(int i=0; i<NumCopies; i++) {
      //get the pointer to the i-th row
      const RealType* restrict prop=awalker.getPropertyBase(i);
      RealType invr = prop[UMBRELLAWEIGHT];
      RealType e = prop[LOCALENERGY];
      *esumPtr++ += invr*e;   //esum(i,ENERGY_INDEX)    += invr*e;
      *esumPtr++ += invr*e*e; //esum(i,ENERGY_SQ_INDEX) += invr*e*e; //how to variance
      *esumPtr++ += invr;     //esum(i,WEIGHT_INDEX)    += invr;

      //accumulate elocal(i,j) with the weight,
      //The content the Walker::DynProperties other than LOCALENERGY are NOT weighted.
      //const RealType *t= awalker.getEnergyBase(i);
      //for(int j=0; j<NumOperators; j++) {
      //  *elocPtr++ += invr*(*t++);
      //}
      
      //const RealType *t= awalker.getPropertyBase(i)+FirstHamiltonian;
      for(int j=0,h=FirstHamiltonian; j<NumOperators; j++,h++) {
        elocal(iloc++) += invr*prop[h];
      }
    }

    ++CurrentWalker;

    //reset to zero 
    if(CurrentWalker == NumWalkers) CurrentWalker=0;
  }

  ///Set CurrentWalker to zero so that accumulation is done in a vectorized way
  void MultipleEnergyEstimator::reset() { 
    CurrentWalker=0;
    elocal=0.0; esum=0.0;
  }

    /** calculate the averages and reset to zero
     *@param record a container class for storing scalar records (name,value)
     *@param wgtinv the inverse weight
     *
     * Disable collection. MultipleEnergyEstimator does not need to communiate at all.
     */
  void MultipleEnergyEstimator::report(RecordNamedProperty<RealType>& record, RealType wgtinv) {
//#ifdef HAVE_MPI
//    if(CollectSum) {
//      int esumSize(esum.size());
//      //pack the energy to be summed: v can be either static or data member
//      vector<RealType> v(esumSize+elocal.size());
//      std::copy(esum.begin(),esum.end(),v.begin());
//      std::copy(elocal.begin(),esum.end(),v.begin()+esumSize);
//
//      gsum(v,0);
//
//      std::copy(v.begin(),v.begin()+esumSize,esum.begin());
//      std::copy(v.begin()+esumSize,v.end(),elocal.begin());
//    }
//#endif
//
    //(localenergy, variance, weight)* 
    int ir(FirstColumnIndex);

    for(int i=0; i<NumCopies; i++) {
      RealType r = 1.0/esum(i,WEIGHT_INDEX);
      RealType e = esum(i,ENERGY_INDEX)*r;
      esum(i,ENERGY_INDEX)=e;
      esum(i,ENERGY_SQ_INDEX)=esum(i,ENERGY_SQ_INDEX)*r-e*e;
      esum(i,WEIGHT_INDEX)*=wgtinv; 
    }

    //swap the row/column indices for (localenergy*, variance*, weight*)
    for(int l=0; l<LE_INDEX;l++) {
      for(int i=0; i<NumCopies; i++) 
        record[ir++]=esum(i,l);
    }

    //(each hamiltonian term)*
    for(int i=0; i<elocal.size(); i++, ir++) {
      record[ir]=wgtinv*elocal(i);
    }
    //int n(elocal.size());
    //const RealType* restrict eit(elocal.data());
    //while(n) {record[ir++]=wgtinv*(*eit++);--n;}

    //set the ScalarEstimator<T>::b_average and b_variace
    b_average = esum(0,ENERGY_INDEX);
    b_variance = esum(0,ENERGY_SQ_INDEX);

    reset();
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
