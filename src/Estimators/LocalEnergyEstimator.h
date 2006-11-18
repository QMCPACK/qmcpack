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
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

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
  class LocalEnergyEstimator: public ScalarEstimatorBase {

    //typedef PooledData<T>                            BufferType;
    enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};

    ///locator of the first data this object handles
    int LocalEnergyIndex;
    int LocalPotentialIndex;
    int FirstHamiltonian;
    int SizeOfHamiltonians;
    QMCHamiltonian& Href;

    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;

    ///vector to contain all the constituents of the local energy
    std::vector<RealType>  elocal;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    LocalEnergyEstimator(QMCHamiltonian& h);

    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     * @param msg buffer for message passing
     */
    void add2Record(RecordListType& record, BufferType& msg) {
      LocalEnergyIndex = record.add(elocal_name[ENERGY_INDEX].c_str());
      for(int i=1; i<elocal_name.size(); i++) record.add(elocal_name[i].c_str());
      //add elocal to the message buffer
      msg.add(elocal.begin(),elocal.end());
    }

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      const RealType* restrict ePtr = awalker.getPropertyBase();
      RealType e = ePtr[LOCALENERGY];
      elocal[ENERGY_INDEX] += wgt*e;
      elocal[ENERGY_SQ_INDEX] += wgt*e*e;
      elocal[POTENTIAL_INDEX] += wgt*ePtr[LOCALPOTENTIAL];
      for(int i=0, target=LE_MAX, source=FirstHamiltonian; i<SizeOfHamiltonians; 
          i++,target++,source++) {
        elocal[target] += wgt*ePtr[source];
      }
    }

    inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) {
      accumulate(awalker,awalker.Weight);
    }

    inline RealType accumulate(WalkerIterator first, WalkerIterator last) {
      //RealType deltaE = Href.getEnsembleAverage();
      int wsum=0;
      while(first != last) {
        accumulate(**first,(*first)->Weight);
        ++first; wsum++;
      }
      return wsum;
      //elocal[ENERGY_INDEX] += static_cast<RealType>(wsum)*deltaE;
    }

    ///reset all the cumulative sums to zero
    inline void reset() { 
      std::fill(elocal.begin(), elocal.end(),0.0);
    }

    ///copy the value to a message buffer
    inline void copy2Buffer(BufferType& msg) {
      msg.put(elocal.begin(),elocal.end());
    }

    /** calculate the averages and reset to zero
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     */
    void report(RecordListType& record, RealType wgtinv, BufferType& msg);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
