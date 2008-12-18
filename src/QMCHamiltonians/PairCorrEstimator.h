//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_PAIRCOOR_HAMILTONIAN_H
#define QMCPLUSPLUS_PAIRCOOR_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus 
{

  class PairCorrEstimator: public QMCHamiltonianBase
  {
    public:

    PairCorrEstimator(ParticleSet& elns);

    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
    {
      return evaluate(P);
    }

    void addObservables(PropertySetType& plist);
    void setObservables(PropertySetType& plist);
    void setParticlePropertyList(PropertySetType& plist, int offset);
    bool put(xmlNodePtr cur);
    bool get(std::ostream& os) const;
    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

    private:
    /// maximum distance 
    RealType Dmax;
    /// bin size 
    RealType Delta;
    /// one of bin size 
    RealType DeltaInv;
    ///save pair indices
    vector<int> pair_ids;
    ///data
    Matrix<RealType> gof_r;
    ///prefix of each gof_r
    vector<string> gof_r_prefix;
    /** resize the internal data
     *
     * The argument list is not completed
     */
    void resize();
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
