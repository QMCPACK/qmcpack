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
    
    
#include "QMCHamiltonians/model/HeSAPT_smoothed.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
HeSAPT_smoothed_phy::HeSAPT_smoothed_phy(ParticleSet& P): PtclRef(&P)
{
  Dependants = 1;
  d_table = DistanceTable::add(P,DT_AOS);
  Return_t rho = P.G.size()/P.Lattice.Volume, N0 = P.G.size();
  rc = P.Lattice.WignerSeitzRadius;
  mirrorshift = -2.0*((-4.27467067939-32.3451957988*rc-37.5775397337*rc*rc+4.0/rc)*std::exp(-5.72036885319*rc) + (18.2962387439-6.16632555293*rc+6.91482730781*rc*rc)*std::exp(-2.80857770752*rc) - damp(7,rc)*1.460977837725/std::pow(rc,6.) - damp(9,rc)*14.11785737/std::pow(rc,8.) - damp(11,rc)*183.691075/std::pow(rc,10.) + damp(12,rc)*76.70/std::pow(rc,11.) - damp(13,rc)*3372./std::pow(rc,12.) + damp(14,rc)*3806./std::pow(rc,13.) - damp(15,rc)*85340./std::pow(rc,14.) + damp(16,rc)*171000./std::pow(rc,15.) - damp(17,rc)*2860000./std::pow(rc,16.));
}

void HeSAPT_smoothed_phy::resetTargetParticleSet (ParticleSet& P)
{
  d_table = DistanceTable::add(P,DT_AOS);
  PtclRef=&P;
  rc = P.Lattice.WignerSeitzRadius;
  mirrorshift = -2.0*((-4.27467067939-32.3451957988*rc-37.5775397337*rc*rc+4.0/rc)*std::exp(-5.72036885319*rc) + (18.2962387439-6.16632555293*rc+6.91482730781*rc*rc)*std::exp(-2.80857770752*rc) - damp(7,rc)*1.460977837725/std::pow(rc,6.) - damp(9,rc)*14.11785737/std::pow(rc,8.) - damp(11,rc)*183.691075/std::pow(rc,10.) + damp(12,rc)*76.70/std::pow(rc,11.) - damp(13,rc)*3372./std::pow(rc,12.) + damp(14,rc)*3806./std::pow(rc,13.) - damp(15,rc)*85340./std::pow(rc,14.) + damp(16,rc)*171000./std::pow(rc,15.) - damp(17,rc)*2860000./std::pow(rc,16.));
}

HeSAPT_smoothed_phy::Return_t HeSAPT_smoothed_phy::evaluate(ParticleSet& P)
{
  Value = 0.0;
  smooth = 0.0;
  for(int i=0; i < d_table->getTotNadj(); i++)
  {
    Return_t r1 = d_table->r(i);
    if (r1 < rc)
    {
      Return_t r1i = 1.0/r1, r2i = r1i*r1i, r6i = r2i*r2i*r2i, r8i = r6i*r2i, r10i = r8i*r2i, r11i = r10i*r1i, r12i = r10i*r2i, r13i = r11i*r2i, r14i = r12i*r2i, r15i = r13i*r2i, r16i = r14i*r2i;
      Return_t rd1 = 2.0*rc - r1, rd1i = 1.0/rd1, rd2i = rd1i*rd1i, rd6i = rd2i*rd2i*rd2i, rd8i = rd6i*rd2i, rd10i = rd8i*rd2i, rd11i = rd10i*rd1i, rd12i = rd10i*rd2i, rd13i = rd11i*rd2i, rd14i = rd12i*rd2i, rd15i = rd13i*rd2i, rd16i = rd14i*rd2i;
      Value += ((-4.27467067939-32.3451957988*r1-37.5775397337*r1*r1+4.0*r1i)*std::exp(-5.72036885319*r1) + (18.2962387439-6.16632555293*r1+6.91482730781*r1*r1)*std::exp(-2.80857770752*r1) - damp(7,r1)*1.460977837725*r6i - damp(9,r1)*14.11785737*r8i - damp(11,r1)*183.691075*r10i + damp(12,r1)*76.70*r11i - damp(13,r1)*3372.*r12i + damp(14,r1)*3806.*r13i - damp(15,r1)*85340.*r14i + damp(16,r1)*171000.*r15i - damp(17,r1)*2860000.*r16i);
      smooth += ((-4.27467067939-32.3451957988*rd1-37.5775397337*rd1*rd1+4.0*rd1i)*std::exp(-5.72036885319*rd1) + (18.2962387439-6.16632555293*rd1+6.91482730781*rd1*rd1)*std::exp(-2.80857770752*rd1) - damp(7,rd1)*1.460977837725*rd6i - damp(9,rd1)*14.11785737*rd8i - damp(11,rd1)*183.691075*rd10i + damp(12,rd1)*76.70*rd11i - damp(13,rd1)*3372.*rd12i + damp(14,rd1)*3806.*rd13i - damp(15,rd1)*85340.*rd14i + damp(16,rd1)*171000.*rd15i - damp(17,rd1)*2860000.*rd16i + mirrorshift);
    }
  }
  Value += smooth;
  // temporary tail corr.  See Mathematica notebook.
  return Value;
}

void HeSAPT_smoothed_phy::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  HeSAPT_smoothed_phy* myClone = new HeSAPT_smoothed_phy(qp);
  targetH.addOperator(myClone, myName, true);
  myClone->addCorrection(targetH);
}

void HeSAPT_smoothed_phy::addCorrection(QMCHamiltonian& targetH)
{
  HeSAPT_smoothed_aux* auxTerm = new HeSAPT_smoothed_aux(this);
  std::string auxName = myName+"_aux";
  targetH.addOperator(auxTerm, auxName, false);
}
}
