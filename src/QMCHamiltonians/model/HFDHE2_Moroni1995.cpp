//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/model/HFDHE2_Moroni1995.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
HFDHE2_Moroni1995_phy::HFDHE2_Moroni1995_phy(ParticleSet& P): PtclRef(&P)
{
  Dependants = 1;
  A = 18.63475757;
  alpha = -2.381392669;
  c1=1.460008056;
  c2=14.22016431;
  c3=187.2033646;
  D = 6.960524706;
  d_table = DistanceTable::add(P,DT_AOS);
  Return_t rho = P.G.size()/P.Lattice.Volume, N0 = P.G.size();
  rc = P.Lattice.WignerSeitzRadius;
  Return_t rcm2 = 1.0/(rc*rc), rcm6 = rcm2*rcm2*rcm2, rcm8 = rcm6*rcm2, rcm10 = rcm8*rcm2;
  mirrorshift = -2.0*(A*std::exp(alpha*rc) - (c1*rcm6+c2*rcm8+c3*rcm10)*dampF(rc));
}

void HFDHE2_Moroni1995_phy::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(P,DT_AOS);
  PtclRef=&P;
  Return_t rc = P.Lattice.WignerSeitzRadius;
  Return_t rcm2 = 1.0/(rc*rc), rcm6 = rcm2*rcm2*rcm2, rcm8 = rcm6*rcm2, rcm10 = rcm8*rcm2;
  mirrorshift = -2.0*(A*std::exp(alpha*rc) - (c1*rcm6+c2*rcm8+c3*rcm10)*dampF(rc));
}

HFDHE2_Moroni1995_phy::Return_t HFDHE2_Moroni1995_phy::evaluate(ParticleSet& P)
{
  Value = 0.0;
  smooth = 0.0;
  for(int i=0; i<d_table->getTotNadj(); i++)
  {
    Return_t r1 = d_table->r(i);
    if (r1 < rc)
    {
      Return_t rm2 = 1.0/(r1*r1), rm6 = rm2*rm2*rm2, rm8 = rm6*rm2, rm10 = rm8*rm2;
      Return_t rd1 = 2.0*rc - r1, rdm2 = 1.0/(rd1*rd1), rdm6 = rdm2*rdm2*rdm2, rdm8 = rdm6*rdm2, rdm10 = rdm8*rdm2;
      Value += (A*std::exp(alpha*r1) - (c1*rm6+c2*rm8+c3*rm10)*dampF(r1));
      smooth += A*std::exp(alpha*rd1) - (c1*rdm6+c2*rdm8+c3*rdm10)*dampF(rd1) + mirrorshift;
    }
  }
  Value += smooth;
//    dep->setValue(-smooth);
  return Value;
}

void HFDHE2_Moroni1995_phy::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  HFDHE2_Moroni1995_phy* myclone=new HFDHE2_Moroni1995_phy(qp);
  targetH.addOperator(myclone,myName,true);
  myclone->addCorrection(targetH);
}

void HFDHE2_Moroni1995_phy::addCorrection(QMCHamiltonian& targetH)
{
  HFDHE2_Moroni1995_aux* auxterm = new HFDHE2_Moroni1995_aux(this);
  std::string auxname=myName+"_aux";
  targetH.addOperator(auxterm,auxname,false);
}

/*
  void HFDHE2_Moroni1995_aux::registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const {
    int loc=h5list.size();
    std::vector<int> onedim(1,1);
    h5list.push_back(new observable_helper("HFDHE2_aux"));
    h5list[loc]->set_dimensions(onedim,myIndex);
    h5list[loc]->open(gid);

    h5list.push_back(new observable_helper("KperP"));
    ++loc;
    h5list[loc]->set_dimensions(onedim,myIndex);
    h5list[loc]->open(gid);
  }
*/

/*
void HFDHE2_Moroni1995_aux::addObservables(PropertySetType& plist)
{
  myIndex=plist.add("HFDHE2_aux");
  //    plist.add("KperP");
}

void HFDHE2_Moroni1995_aux::setObservables(PropertySetType& plist) {
  plist[myIndex]=Value;
  //      plist[myIndex+1]=KValue;
}

void HFDHE2_Moroni1995_aux::setParticlePropertyList(PropertySetType& plist, int offset)
{
  plist[myIndex+offset]=Value;
  //      plist[myIndex+1+offset]=KValue;
}
*/
}
