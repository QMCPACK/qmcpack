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
#ifndef QMCPLUSPLUS_WFMCONLYESTIMATOR_H
#define QMCPLUSPLUS_WFMCONLYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus {

  /*** A class to evaluate the local energy 
   *
   *The WFMCOnlyEstimator evaluates 
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
  class WFMCOnlyEstimator: public ScalarEstimatorBase {

    //typedef PooledData<T>                            BufferType;
    //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};
    enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

    //int LocalPotentialIndex;
    int FirstHamiltonian;
    int SizeOfHamiltonians;
    int PsiIndex;

    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    WFMCOnlyEstimator(QMCHamiltonian& h);

    /** implement virtual function
     */
    ScalarEstimatorBase* clone();

    void add2Record(RecordListType& record);

    inline void accumulate(const Walker_t& awalker, RealType wgt) {
      const RealType* restrict ePtr = awalker.getPropertyBase();
      ///weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
      RealType wwght= wgt* awalker.Weight;
//       RealType wwght= wgt;
      
      scalars[0](ePtr[LOCALENERGY],wwght);
      scalars[1](ePtr[LOCALPOTENTIAL],wwght);
      int target=2;
      int source;
      for(source=FirstHamiltonian; source<(FirstHamiltonian+SizeOfHamiltonians); 
          ++target, ++source)
        {
	  scalars[target](ePtr[source],wwght);
	}
        scalars[ target](wwght*ePtr[PsiIndex+FirstHamiltonian],1);
        scalars[ target+1](wwght,1);
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
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3449 $   $Date: 2009-01-11 08:47:20 -0600 (Sun, 11 Jan 2009) $
 * $Id: WFMCOnlyEstimator.h 3449 2009-01-11 14:47:20Z jnkim $ 
 ***************************************************************************/
