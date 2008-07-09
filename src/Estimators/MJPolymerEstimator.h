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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_FULL_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_FULL_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"

namespace qmcplusplus {

  struct MJPolymerEstimator: public PolymerEstimator {

    MJPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
    MJPolymerEstimator(const MJPolymerEstimator& mest);
    ScalarEstimatorBase* clone();

    void add2Record(RecordNamedProperty<RealType>& record);

    inline  void accumulate(const Walker_t& awalker, RealType wgt) {}

//     inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) { }

    void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt);

    void evaluateDiff();

    private:
          ///vector to contain the names of all the constituents of the local energy
      std::vector<string> elocal_name;
      int FirstHamiltonian;
      int SizeOfHamiltonians;
      double KEconst;
//       QMCHamiltonian* Hpointer;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: MJPolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $ 
 ***************************************************************************/
