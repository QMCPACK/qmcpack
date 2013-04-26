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
#ifndef QMCPLUSPLUS_LOCAL_MOMENT_HAMILTONIAN_H
#define QMCPLUSPLUS_LOCAL_MOMENT_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{

/** local moment estimator
 *
 * Compute local moment for up - down spins and all ion types in some particleset
 */
class LocalMomentEstimator: public QMCHamiltonianBase
{
public:

  /** constructor
   * @param elns target particle set
   * @param sources source particle set
   *
   * Use the elns.DistTables to evaluate the pair correlation functions.
   */
  LocalMomentEstimator(ParticleSet& elns, ParticleSet& srcs);

  void resetTargetParticleSet(ParticleSet& P);

  /* evaluate the pair correlation functions */
  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.size();
    int dm=std::floor(Dmax*100);
    for (int i=0; i<lm.size1()*lm.size2(); i++)
    {
      std::stringstream sstr;
      sstr << "lm" <<dm<<names[i];
      plist.add(sstr.str());
    }
  }
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(vector<observable_helper*>& h5list, hid_t gid) const;
  void setObservables(PropertySetType& plist)
  {
    int k(0);
    for(int i(0); i<lm.size1(); i++)
      for(int j(0); j<lm.size2(); j++,k++)
        plist[myIndex+k]=lm(i,j);
  }
  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    int k(0);
    for(int i(0); i<lm.size1(); i++)
      for(int j(0); j<lm.size2(); j++,k++)
        plist[myIndex+k+offset]=lm(i,j);
  }
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

private:
  /// maximum distance
  RealType Dmax;
  Matrix<RealType> lm;
  vector<int> ion_id, el_id;
  vector<RealType> el_nrm;
  int nag;
  DistanceTableData* d_table;
  vector<string> names;
  ParticleSet& ions;
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
