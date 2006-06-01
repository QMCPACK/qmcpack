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
#ifndef QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#include <strstream>
#include "OhmmsData/libxmldefs.h"
#include "Estimators/ScalarEstimatorBase.h"

namespace qmcplusplus {

  /*** A class to evaluate the local energy 
   *
   *The LocalEnergyEstimator evaluates 
   <ul>
   <li> LocalEnergy
   <li> Variance of the LocalEnergy
   </ul>
   and the values of each QHCHamiltonianBase elements added by an application.
   Typical local energies are
   <li> Kinetic Energy
   <li> External Potential Energy
   <li> Electron-Electron Energy
   <li> Ion-Ion Energy (for molecules only)
   <li> Conserved Quantity (VMC only)
   </ul>
   The method of evaluating the estimators is as follows
   \f[
   \langle O \rangle = \frac{\sum_{i=1}^N w_i O_i}{\sum_{i=1}^N w_i}
   \f]
   where the sum runs over the total number of accumulated walkers 
   and \f$ w_i \f$ is the weight of the \f$ ith \f$ walker.
   *
   *The formula for the LocalEnergy 
   \f[ E_L({\bf R}) = \frac{\hat{H} \Psi({\bf R})}{ \Psi({\bf R})} \f]
  */
  template<class T>
  class LocalEnergyEstimator: public ScalarEstimatorBase<T> {

    enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};

    using ScalarEstimatorBase<T>::CollectSum;
    using ScalarEstimatorBase<T>::b_average;
    using ScalarEstimatorBase<T>::b_variance;

    ///locator of the first data this object handles
    int LocalEnergyIndex;
    int LocalPotentialIndex;
    int FirstHamiltonian;
    int SizeOfHamiltonians;
    QMCHamiltonian& Href;

    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;

    ///vector to contain all the constituents of the local energy
    std::vector<T>  elocal;

  public:

    typedef typename ScalarEstimatorBase<T>::Walker_t Walker_t;
    typedef typename ScalarEstimatorBase<T>::WalkerIterator WalkerIterator;
  
    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    LocalEnergyEstimator(QMCHamiltonian& h):Href(h) { 
      int hterms(h.size());
      SizeOfHamiltonians = hterms;
      FirstHamiltonian = h.startIndex();
      elocal.resize(SizeOfHamiltonians+LE_MAX);
      elocal_name.resize(SizeOfHamiltonians+LE_MAX);
      elocal_name[ENERGY_INDEX] = "LocalEnergy";
      elocal_name[ENERGY_SQ_INDEX] = "Variance";
      elocal_name[POTENTIAL_INDEX] = "LocalPotential";
      int ii(LE_MAX);
      for(int i=0; i < SizeOfHamiltonians; i++) elocal_name[ii++] = h.getName(i);
    }

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     *@param record storage of scalar records (name,value)
     */
    void add2Record(RecordNamedProperty<T>& record) {
      LocalEnergyIndex = record.add(elocal_name[ENERGY_INDEX].c_str());
      for(int i=1; i<elocal_name.size(); i++)
	record.add(elocal_name[i].c_str());
    }

    inline void accumulate(const Walker_t& awalker, T wgt) {
      const T* restrict ePtr = awalker.getPropertyBase();
      T e = ePtr[LOCALENERGY];
      elocal[ENERGY_INDEX] += wgt*e;
      elocal[ENERGY_SQ_INDEX] += wgt*e*e;
      elocal[POTENTIAL_INDEX] += wgt*ePtr[LOCALPOTENTIAL];
      for(int i=0, target=LE_MAX, source=FirstHamiltonian; i<SizeOfHamiltonians; 
          i++,target++,source++) {
	elocal[target] += wgt*ePtr[source];
      }
      /*
      T e = awalker.Properties(LOCALENERGY);
      const T* restrict e_ptr = awalker.getEnergyBase();
      //energy_sum += wgt*e; energy_sq_sum += wgt*e*e;
      elocal[ENERGY_INDEX] += wgt*e;
      elocal[ENERGY_SQ_INDEX] += wgt*e*e;
      elocal[POTENTIAL_INDEX] += wgt*awalker.Properties(LOCALPOTENTIAL);
      for(int ii=LE_MAX, i=0; i<SizeOfHamiltonians; ii++,i++) {
	elocal[ii] += wgt*e_ptr[i]; //wgt*(awalker.E[i]);
      }
      */
    }

    void accumulate(WalkerIterator first, WalkerIterator last) {
      T deltaE = Href.getEnsembleAverage();
      int wsum=0;
      while(first != last) {
        accumulate(**first,(*first)->Weight);
        ++first; wsum++;
      }
      elocal[ENERGY_INDEX] += static_cast<T>(wsum)*deltaE;
    }

    ///reset all the cumulative sums to zero
    inline void reset() { 
      for(int i=0; i<elocal.size(); i++) {
	elocal[i] = T();
      }
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    inline void report(RecordNamedProperty<T>& record, T wgtinv) {
      if(CollectSum) gsum(elocal,0);
      register int ir=LocalEnergyIndex;
      b_average =  elocal[ENERGY_INDEX]*wgtinv;
      b_variance = elocal[ENERGY_SQ_INDEX]*wgtinv-b_average*b_average;
      record[ir++] = b_average;
      record[ir++] = b_variance;
      record[ir++] = elocal[POTENTIAL_INDEX]*wgtinv;
      for(int i=0, ii=LE_MAX; i<SizeOfHamiltonians; i++,ii++) {
	record[ir++] = elocal[ii]*wgtinv;
      }
      reset();
    }

  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
