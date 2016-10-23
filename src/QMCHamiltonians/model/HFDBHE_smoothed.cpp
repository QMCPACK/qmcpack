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
    
    
#include "QMCHamiltonians/model/HFDBHE_smoothed.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HFDBHE_smoothed for the indentical source and target particle sets.
 */
HFDBHE_smoothed_phy::HFDBHE_smoothed_phy(ParticleSet& P): PtclRef(&P)
{
  Dependants = 1;
  d_table = DistanceTable::add(P);
  rc = P.Lattice.WignerSeitzRadius;
  mirrorshift = -2.0*(6.394277071*std::exp(-1.863335173*rc-0.072712207*rc*rc) - (1.461*std::pow(rc,-6.0)+14.11*std::pow(rc,-8.0)+183.5*std::pow(rc,-10.0))*dampF(rc));
}


void HFDBHE_smoothed_phy::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(P);
  PtclRef=&P;
  rc = P.Lattice.WignerSeitzRadius;
  mirrorshift = -2.0*(6.394277071*std::exp(-1.863335173*rc-0.072712207*rc*rc) - (1.461*std::pow(rc,-6.0)+14.11*std::pow(rc,-8.0)+183.5*std::pow(rc,-10.0))*dampF(rc));
}


HFDBHE_smoothed_phy::Return_t HFDBHE_smoothed_phy::evaluate(ParticleSet& P)
{
  Value = 0.0;
  smooth = 0.0;
  for(int i=0; i<d_table->getTotNadj(); i++)
  {
    Return_t r1 = d_table->r(i);
    if (r1 < rc)
    {
      Return_t r2i = 1.0/(r1*r1), r6i = r2i*r2i*r2i, r8i = r6i*r2i, r10i = r8i*r2i;
      Return_t rd1 = 2.0*rc - r1, rd2i = 1.0/(rd1*rd1), rd6i = rd2i*rd2i*rd2i, rd8i = rd6i*rd2i, rd10i = rd8i*rd2i;
      Value += (6.394277071*std::exp(-1.863335173*r1-0.072712207*r1*r1) - (1.461*r6i+14.11*r8i+183.5*r10i)*dampF(r1));
      smooth += (6.394277071*std::exp(-1.863335173*rd1-0.072712207*rd1*rd1) - (1.461*rd6i+14.11*rd8i+183.5*rd10i)*dampF(rd1) + mirrorshift);
    }
  }
  Value += smooth;
  return Value;
}


void HFDBHE_smoothed_phy::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  HFDBHE_smoothed_phy* myClone = new HFDBHE_smoothed_phy(qp);
  targetH.addOperator(myClone, myName, true);
  myClone->addCorrection(targetH);
}

void HFDBHE_smoothed_phy::addCorrection(QMCHamiltonian& targetH)
{
  HFDBHE_smoothed_aux* auxTerm = new HFDBHE_smoothed_aux(this);
  std::string auxName = myName+"_aux";
  targetH.addOperator(auxTerm, auxName, false);
}
}
