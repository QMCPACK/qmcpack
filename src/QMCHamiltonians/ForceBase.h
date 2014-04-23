//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

struct ForceBase
{
  /** cheat, need to use virtual inheriance to clean up*/
  typedef QMCTraits::RealType real_type;

  int FirstForceIndex;
  int myTableIndex;
  int Nnuc;
  int Nel;
  int tries;
  bool FirstTime;
  bool addionion;

  ParticleSet& Ions;
  ParticleSet::ParticlePos_t forces;
  ParticleSet::ParticlePos_t forces_IonIon;
  SymTensor<real_type,OHMMS_DIM> stress_IonIon;
  SymTensor<real_type,OHMMS_DIM> stress_ee;
  SymTensor<real_type,OHMMS_DIM> stress;
  
  string prefix;
  string pairName;

  // Data for variance reduction of Chiesa et al.
  // PRL 94, 036404 (2005)
  real_type Rcut;
  int m;
  vector<real_type> ck;
  inline real_type g(real_type r)
  {
    if (r > Rcut)
      return 1.0;
    real_type sum=0.0;
    real_type r2kplusm = r;
    for (int i=0; i<m; i++)
      r2kplusm *= r;
    for (int k=0; k<ck.size(); k++)
    {
      sum += ck[k] * r2kplusm;
      r2kplusm *= r;
    }
    return sum;
  }


  void InitVarReduction (real_type Rcut, int m, int numFuncs);


  ForceBase(ParticleSet& ions, ParticleSet& elns);
  virtual ~ForceBase() {}

  void registerObservablesF(vector<observable_helper*>& h5list, hid_t gid) const;

  void addObservablesF(QMCTraits::PropertySetType& plist);
  void addObservablesStress(QMCTraits::PropertySetType& plist);
  void setObservablesF(QMCTraits::PropertySetType& plist);
  void setObservablesStress(QMCTraits::PropertySetType& plist);
  void setParticleSetF(QMCTraits::PropertySetType& plist, int offset);
  void setParticleSetStress(QMCTraits::PropertySetType& plist, int offset);
};

struct BareForce: public QMCHamiltonianBase, public ForceBase
{
  BareForce(ParticleSet& ions, ParticleSet& elns);
  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(vector<observable_helper*>& h5list, hid_t gid) const
  {
    registerObservablesF(h5list,gid);
  }

  /** default implementation to add named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   */
  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    setParticleSetF(plist, offset);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "Force Base Hamiltonian: " << pairName;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
