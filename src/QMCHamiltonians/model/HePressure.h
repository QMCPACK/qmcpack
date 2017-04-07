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
    
    
#ifndef QMCPLUSPLUS_HEPRESSURE_H
#define QMCPLUSPLUS_HEPRESSURE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"



namespace qmcplusplus
{

/** @ingroup hamiltonian
 @brief Evaluate the He Pressure.

 where d is the dimension of space and /Omega is the volume.
**/

struct HePressure: public QMCHamiltonianBase
{
  double A,alpha,c0,c1,c2,c3,c1p,c2p,c3p,D,rdV,dV,dVnorm ;
  double pNorm,tailcorr,rcutoff,kNorm;
  bool ZV;
  bool ZB;

  ParticleSet* sourcePtcl;
  DistanceTableData* d_table;

  /** constructor
   *
   * HePressure operators need to be re-evaluated during optimization.
   */
  HePressure(ParticleSet& P ): sourcePtcl(&P)
  {
    UpdateMode.set(OPTIMIZABLE,1);
    d_table = DistanceTable::add(P,DT_AOS);
    Return_t rho = P.G.size()/P.Lattice.Volume;
    Return_t N0 = P.G.size();
    rcutoff = P.Lattice.WignerSeitzRadius;
    Return_t rs = std::pow(3.0/(4.0*M_PI*rho),1.0/3.0);
    pNorm = -1.0/(3.0*P.Lattice.Volume);
    kNorm = 2.0/(3.0*P.Lattice.Volume);
    dVnorm = 2.0*pNorm*rs;
    A = 18.63475757;
    alpha = -2.381392669;
    c1=1.460008056;
    c2=14.22016431;
    c3=187.2033646;
    c0 = A*alpha;
    c1p = -6*c1;
    c2p = -8*c2;
    c3p = -10*c3;
    D = 6.960524706;
    Return_t r1 = rcutoff;
    Return_t r2 = (r1*r1);
    Return_t rm3 = 1.0/(r2*r1);
    Return_t rm5 = rm3/r2;
    Return_t rm7 = rm5/r2;
    tailcorr =  (10*c3)/7*rm7 + (8*c2)/5*rm5 + (2*c1)*rm3 + ( A*std::exp(alpha*r1) * (6 - alpha*rcutoff*(6 - alpha*rcutoff* (3 - alpha* rcutoff))))/(alpha*alpha*alpha);
    tailcorr *= pNorm;
    tailcorr *= 2.0*M_PI*rho*N0;
    std::cout << std::setprecision(12) <<"  HePressure:: Tail Correction for "<<sourcePtcl->getName()<<" is: "<<tailcorr<< std::endl;
//       std::cout <<"  Cutoff is: "<<rcutoff<< std::endl;
  }
  ///destructor
  ~HePressure() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(P,DT_AOS);
    sourcePtcl= &P;
    Return_t rho = P.G.size()/P.Lattice.Volume;
    Return_t N0 = P.G.size();
    rcutoff = P.Lattice.WignerSeitzRadius;
    Return_t rs = std::pow(3.0/(4.0*M_PI*rho),1.0/3.0);
    pNorm = -1.0/(3.0*P.Lattice.Volume);
    kNorm = 2.0/(3.0*P.Lattice.Volume);
    dVnorm = 2.0*pNorm*rs;
    Return_t r1 = rcutoff;
    Return_t r2 = (r1*r1);
    Return_t rm3 = 1.0/(r2*r1);
    Return_t rm5 = rm3/r2;
    Return_t rm7 = rm5/r2;
    tailcorr =  (10*c3)/7*rm7 + (8*c2)/5*rm5 + (2*c1)*rm3 + ( A*std::exp(alpha*r1) * (6 - alpha*rcutoff*(6 - alpha*rcutoff* (3 - alpha* rcutoff))))/(alpha*alpha*alpha);
    tailcorr *= pNorm;
    tailcorr *= 2.0*M_PI*rho*N0;
//       std::cout <<"  Pressure Tail Correction for "<<sourcePtcl->getName()<<" is: "<<tailcorr<< std::endl;
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value = 0.0;
//       dV = 0.0;
    Return_t KE = P.PropertyList[LOCALENERGY]-P.PropertyList[LOCALPOTENTIAL];
    for(int i=0; i<d_table->getTotNadj(); i++)
    {
      Return_t RR = d_table->r(i);
      if ( RR < rcutoff)
      {
        Return_t rm1 = 1.0/RR;
        Return_t rm2 = rm1*rm1;
        Return_t rm6 = std::pow(rm2,3);
        Return_t rm8 = rm6*rm2;
        Return_t rm10 = rm8*rm2;
        Return_t TMP1;
        if (RR < D)
        {
          Return_t t1=(D*rm1 - 1.0);
          Return_t dampF = std::exp(-t1*t1);
          TMP1 = c0*RR*std::exp(alpha*RR) - (c1p*rm6 + c2p*rm8 + c3p*rm10 + 2.0* (c1*rm6+c2*rm8+c3*rm10)*D*rm1*t1 )*dampF;
        }
        else
          TMP1 = c0*RR*std::exp(alpha*RR) - (c1p*rm6 + c2p*rm8 + c3p*rm10);
        Value += TMP1;
        dV += TMP1/RR;
      }
    }
    dV *= dVnorm;
    Value *= pNorm;
    rdV = Value + tailcorr;
    Value += kNorm*KE;
    Value += tailcorr;
    return 0.0;
  }

//     inline Return_t dampF(Return_t r) {
//       if (r < D){
//         Return_t t1=(D/r - 1.0);
//         return std::exp(-t1*t1);
//       }
//       else
//         return 1.0;
//     }

  inline Return_t
  evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HePressure: "<< sourcePtcl->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HePressure(qp);
  }
  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.add("HePress");
    plist.add("rijdrV");
    plist.add("drV");
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
    plist[myIndex+1]=rdV;
    plist[myIndex+2]=dV;
  }
  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    plist[myIndex+offset]=Value;
    plist[myIndex+1+offset]=rdV;
    plist[myIndex+2+offset]=dV;
  }

};
}
#endif


