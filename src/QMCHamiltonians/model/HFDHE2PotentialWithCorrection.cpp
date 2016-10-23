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
    
    
#include "QMCHamiltonians/model/HFDHE2Potential.h"
#include "Particle/DistanceTableData.h"
namespace qmcplusplus
{

HFDHE2Potential::HFDHE2Potential(ParticleSet& P)
{
  Dependants=1;
  A = 18.63475757;
  alpha = -2.381392669;
  c1=1.460008056;
  c2=14.22016431;
  c3=187.2033646;
  D = 6.960524706;
  rc = P.Lattice.WignerSeitzRadius;
  kpre_factor = 3.157733e+5/static_cast<Return_t>(P.getTotalNum());
  tailcorr=mycorrection(rc);
  //Return_t r2 = (rc*rc);
  //Return_t rm2 = 1.0/r2;
  //Return_t rm6 = std::pow(rm2,3);
  //Return_t rm8 = rm6*rm2;
  //Return_t rm10 = rm8*rm2;
  //tailcorr = (A*std::exp(alpha*rc) - (c1*rm6+c2*rm8+c3*rm10)*dampF(rc));
}

void HFDHE2Potential::resetTargetParticleSet(ParticleSet& P)
{
  rc = P.Lattice.WignerSeitzRadius;
  tailcorr=mycorrection(rc);
}

HFDHE2Potential::Return_t HFDHE2Potential::evaluate(ParticleSet& P)
{
  Value = 0.0;
  TCValue = 0.0;
  const DistableTableData* d_table=P.DistTables[0];
  for(int i=0; i<d_table->getTotNadj(); ++i)
  {
    Return_t r1 = d_table->r(i);
    if ( r1 < rc)
    {
      Value += mycorrection(r1)-tailcorr;
      //Return_t r2 = (r1*r1);
      //Return_t rm2 = 1.0/r2;
      //Return_t rm6 = std::pow(rm2,3);
      //Return_t rm8 = rm6*rm2;
      //Return_t rm10 = rm8*rm2;
      //Value += (A*std::exp(alpha*r1) - (c1*rm6+c2*rm8+c3*rm10)*dampF(r1)) - tailcorr;
      TCValue +=tailcorr;
    }
  }
  //TCorr->set_TC(TCValue);
  return Value;
}

void  HFDHE2Potential::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi
                                       , QMCHamiltonian& targetH)
{
  HFDHE2Potential* myclone=new HFDHE2Potential(qp);
  targetH.addOperator(myclone,myName,true);
  myclone->addCorrection(targetH);
}

void  HFDHE2Potential::addCorrection(QMCHamiltonian& targetH)
{
  HFDHE2Potential_tail* tcrr = new HFDHE2Potential_tail(this);
  std::string corrected=myName+"_tail";
  targetH.addOperator(tcrr,corrected,false);
}

void HFDHE2Potential_tail::registerObservables(std::vector<observable_helper*>& h5list
    , hid_t gid) const
{
  int loc=h5list.size();
  std::vector<int> onedim(1,1);
  h5list.push_back(new observable_helper("HFDHE2tail"));
  h5list[loc]->set_dimensions(onedim,myIndex);
  h5list[loc]->open(gid);
  ++loc;
  h5list.push_back(new observable_helper("KperP"));
  h5list[loc]->set_dimensions(onedim,myIndex);
  h5list[loc]->open(gid);
}

void HFDHE2Potential_tail::addObservables(PropertySetType& plist)
{
  myIndex=plist.add("HFDHE2tail");
  plist.add("KperP");
}

}
