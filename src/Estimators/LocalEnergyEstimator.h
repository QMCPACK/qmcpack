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

    void add2Record(RecordListType& record);

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      const RealType* restrict ePtr = awalker.getPropertyBase();
      scalars[0](ePtr[LOCALENERGY],wgt);
      scalars[1](ePtr[LOCALPOTENTIAL],wgt);
      for(int target=2, source=FirstHamiltonian; target<scalars.size(); 
          ++target, ++source)
        scalars[target](ePtr[source],wgt);
    }

    inline void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt) {
      while(first != last) {
        accumulate(**first,wgt);
        ++first;
      }
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
