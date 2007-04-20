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
    //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};
    enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

    //int LocalPotentialIndex;
    int FirstHamiltonian;
    int SizeOfHamiltonians;

    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    LocalEnergyEstimator(QMCHamiltonian& h);
    LocalEnergyEstimator(const LocalEnergyEstimator& est);

    /** implement virtual function
     */
    ScalarEstimatorBase* clone();
    /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
     * @param record storage of scalar records (name,value)
     * @param msg buffer for message passing
     */
    void add2Record(RecordListType& record, BufferType& msg) {
      FirstIndex = record.add(elocal_name[0].c_str());
      for(int i=1; i<elocal_name.size(); i++) record.add(elocal_name[i].c_str());
      LastIndex=FirstIndex + elocal_name.size();
      //add elocal to the message buffer
      msg.add(d_data.begin(),d_data.end());
    }

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      const RealType* restrict ePtr = awalker.getPropertyBase();
      //RealType e = ePtr[LOCALENERGY];
      int target=0;
      d_sum=d_data[target++] += wgt*ePtr[LOCALENERGY];
      d_sumsq=d_data[target++] += wgt*ePtr[LOCALENERGY]*ePtr[LOCALENERGY];
      d_data[target++] += wgt*ePtr[LOCALPOTENTIAL];
      d_data[target++] += wgt*ePtr[LOCALPOTENTIAL]*ePtr[LOCALPOTENTIAL];
      for(int i=0, source=FirstHamiltonian; i<SizeOfHamiltonians; i++,source++) {
        d_data[target++] += wgt*ePtr[source];
        d_data[target++] += wgt*ePtr[source]*ePtr[source];
      }
      d_wgt+=wgt;
    }

    inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) {
      accumulate(awalker,awalker.Weight);
    }

    inline void accumulate(WalkerIterator first, WalkerIterator last) {
      while(first != last) {
        accumulate(**first,(*first)->Weight);
        ++first;
      }
      //elocal[ENERGY_INDEX] += static_cast<RealType>(wsum)*deltaE;
    }

    inline void report(RecordListType& record, RealType wgtinv)
    {
      //for(int i=0; i<d_data.size(); i++) d_data[i] *= wgtinv;
      //std::copy(d_data.begin(),d_data.end(),record.begin()+FirstIndex);
      d_average =  d_data[0]*wgtinv;
      d_variance = d_data[1]*wgtinv-d_average*d_average;
      cout << "Average = " << d_average << " variance = " << d_variance << endl;
      //record[ir++] = d_average;
      //record[ir++] = d_variance;
      //record[ir++] = d_data[POTENTIAL_INDEX]*wgtinv;
      //for(int i=0, ii=LE_MAX; i<SizeOfHamiltonians; i++,ii++) {
      //  record[ir++] = d_data[ii]*wgtinv;
      //}
      std::fill(d_data.begin(), d_data.end(),0.0);
    }

    void reset()
    {
      d_wgt=0.0;
      std::fill(d_data.begin(), d_data.end(),0.0);
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
