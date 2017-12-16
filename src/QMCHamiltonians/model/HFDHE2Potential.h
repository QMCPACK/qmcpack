//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HFDHE2POTENTIAL_H
#define QMCPLUSPLUS_HFDHE2POTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

#include "QMCHamiltonians/model/HFDHE2Potential_tail.h"


namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HFDHE2Potential for the indentical source and target particle sets.
 */
struct HFDHE2Potential: public QMCHamiltonianBase
{

  Return_t tailcorr,TCValue,rc,A,alpha,c1,c2,c3,D;
  // remember that the default units are Hartree and Bohrs
  DistanceTableData* d_table;
  ParticleSet* PtclRef;
  HFDHE2Potential_tail* TCorr;

  // epsilon = 3.42016039e-5, rm = 5.607384357
  // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
  HFDHE2Potential(ParticleSet& P): PtclRef(&P)
  {
    Dependants=1;
    A = 18.63475757;
    alpha = -2.381392669;
    c1=1.460008056;
    c2=14.22016431;
    c3=187.2033646;
    D = 6.960524706;
    d_table = DistanceTable::add(P,DT_AOS);
//       Return_t rho = P.G.size()/P.Lattice.Volume;
//       Return_t N0 = P.G.size();
    rc = P.Lattice.WignerSeitzRadius;
//       tailcorr = 2.0*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
//       std::cout <<"  HFDHE2Potential tail correction is  "<<tailcorr<< std::endl;
    depName = "HFDHE2tail";
    Return_t r2 = (rc*rc);
    Return_t rm2 = 1.0/r2;
    Return_t rm6 = std::pow(rm2,3);
    Return_t rm8 = rm6*rm2;
    Return_t rm10 = rm8*rm2;
    tailcorr = (A*std::exp(alpha*rc) - (c1*rm6+c2*rm8+c3*rm10)*dampF(rc));
  }

  ~HFDHE2Potential() { }

  QMCHamiltonianBase* makeDependants(ParticleSet& qp )
  {
    TCorr = new HFDHE2Potential_tail(qp);
    return TCorr;
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(P,DT_AOS);
    PtclRef=&P;
//       Return_t rho = P.G.size()/P.Lattice.Volume;
//       Return_t N0 = P.G.size();
    rc = P.Lattice.WignerSeitzRadius;
//       tailcorr = 2.0*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value = 0.0;
    TCValue = 0.0;
    for(int i=0; i<d_table->getTotNadj(); i++)
    {
      Return_t r1 = d_table->r(i);
      if ( r1 < rc)
      {
        Return_t r2 = (r1*r1);
        Return_t rm2 = 1.0/r2;
        Return_t rm6 = std::pow(rm2,3);
        Return_t rm8 = rm6*rm2;
        Return_t rm10 = rm8*rm2;
        Value += (A*std::exp(alpha*r1) - (c1*rm6+c2*rm8+c3*rm10)*dampF(r1)) - tailcorr;
        TCValue +=tailcorr;
      }
    }
    TCorr->set_TC(TCValue);
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t dampF(Return_t r)
  {
    if (r < D)
    {
      Return_t t1=(D/r - 1.0);
      return std::exp(-t1*t1);
    }
    else
      return 1.0;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFDHE2Potential (T/S): " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HFDHE2Potential(qp);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.add("HFDHE2");
//       plist.add("HFDHE2tail");
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
//       plist[myIndex+1]=tailcorr;
  }
};
}
#endif
