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
    
    
#include "QMCHamiltonians/model/LennardJones_smoothed.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
LennardJones_smoothed_phy::LennardJones_smoothed_phy(ParticleSet& P, RealType e /* =3.23648511e-5 */, RealType s /* =4.830139998 */): d_table(NULL), PtclRef(&P), epsilon(e), sigma(s)
{
  Dependants = 1;
  s6 = sigma*sigma*sigma*sigma*sigma*sigma;
  d_table = DistanceTable::add(P,DT_AOS);
  rc = P.Lattice.WignerSeitzRadius;
  RealType rc6i = std::pow(rc,-6.0);
  mirrorshift = -2.0*(s6*rc6i - 1.0)*rc6i;
}


void LennardJones_smoothed_phy::resetTargetParticleSet(ParticleSet& P)
{
  PtclRef=&P;
  d_table = DistanceTable::add(P,DT_AOS);
  rc = P.Lattice.WignerSeitzRadius;
  RealType rc6i = std::pow(rc,-6.0);
  mirrorshift = -2.0*(s6*rc6i - 1.0)*rc6i;
}


LennardJones_smoothed_phy::Return_t LennardJones_smoothed_phy::evaluate(ParticleSet& P)
{
  smooth = 0.0;
  Value = 0.0;
  for(int i=0; i<d_table->getTotNadj(); i++)
  {
    Return_t r1 = d_table->r(i);
    if (r1 < rc)
    {
      Return_t r6i = 1.0/(r1*r1*r1*r1*r1*r1),
               rd1 = 2.0*rc - r1, rd6i = 1.0/(rd1*rd1*rd1*rd1*rd1*rd1);
      Value += (s6*r6i - 1.0)*r6i;
      smooth += ((s6*rd6i - 1.0)*rd6i + mirrorshift);
    }
  }
  Value += smooth;
  Value *= 4.0*epsilon*s6;
  smooth *= 4.0*epsilon*s6;
  return Value;
}


void LennardJones_smoothed_phy::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  LennardJones_smoothed_phy* myClone = new LennardJones_smoothed_phy(qp);
  targetH.addOperator(myClone, myName, true);
  myClone->addCorrection(targetH);
}


void LennardJones_smoothed_phy::addCorrection(QMCHamiltonian& targetH)
{
  LennardJones_smoothed_aux* auxTerm = new LennardJones_smoothed_aux(this);
  std::string auxName = myName+"_aux";
  targetH.addOperator(auxTerm, auxName, false);
}
}
