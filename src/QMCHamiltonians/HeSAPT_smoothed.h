//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HESAPT_SMOOTHED_H
#define QMCPLUSPLUS_HESAPT_SMOOTHED_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HeSAPT_smoothed_phy for the indentical source and target particle sets.
 */
struct HeSAPT_smoothed_phy: public QMCHamiltonianBase
{
  Return_t rc, mirrorshift, smooth;
  // remember that the default units are Hartrees and Bohrs
  DistanceTableData* d_table;
  ParticleSet* PtclRef;

  HeSAPT_smoothed_phy(ParticleSet& P);

  ~HeSAPT_smoothed_phy() { }

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t damp(int n, Return_t r)
  {
    Return_t wert, br, sum, term;
    int i;
    br = 2.41324077320*r;
    if (r >= 1.0)
    {
      sum = term = 1.0;
      for (i=1; i<=n; i++)
      {
        term *= (br/i);
        sum += term;
      }
      wert = 1.0 - std::exp(-br)*sum;
    }
    else
    {
      sum = 0.0;
      term = 1.0;
      for (i=1; i<=n; i++)
        term *= (br/i);
      for (i=n+1; i<=30; i++)
      {
        term *= (br/i);
        sum += term;
      }
      wert = std::exp(-br)*sum;
    }
    return wert;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HeSAPT_smoothed_phy: " << PtclRef->getName();
    return true;
  }

  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

  void addCorrection(QMCHamiltonian& targetH);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HeSAPT_smoothed_phy(qp);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex = plist.add(myName);
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex] = Value;
  }
};


/** @ingroup hamiltonian
 *@brief HeSAPT_smoothed_aux for the indentical source and target particle sets.
 */
struct HeSAPT_smoothed_aux: public QMCHamiltonianBase
{
  const HeSAPT_smoothed_phy* phyH;
  Return_t tailcorr;
  // remember that the default units are Hartrees and Bohrs

  HeSAPT_smoothed_aux(const HeSAPT_smoothed_phy* orig): phyH(orig)
  {
    Return_t rho = phyH->PtclRef->G.size()/phyH->PtclRef->Lattice.Volume, N0 = phyH->PtclRef->G.size(), rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    std::string warningmsg = std::string("  Warning: Tail correction is set to ZERO because ...\n")
                        + "  No pre-calculated tail correction value is found for the current\n"
                        + "  system size and density.  Please do a separate calculation\n"
                        + "  of the tail correction and add it onto the LocalEnergy.\n"
                        + "  And if necessary, edit QMCHamiltonians/HeSAPT_smoothed.h\n"
                        + "  to add the new tailcorr value to the database.";
    if (N0 == 32)
    {
      if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
        tailcorr = -0.0002705961159;	// RHO = 0.3648/sigma^3 : 01x
      else
        if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
          tailcorr = -0.0002178620081;	// RHO = 0.3280/sigma^3 : 21x
        else
          if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
            tailcorr = -0.0002443783556;	// RHO = 0.3470/sigma^3 : 31x
          else
            if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
              tailcorr = -0.0003280694504;	// RHO = 0.4009/sigma^3 : 41x
            else
              if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                tailcorr = -0.0003926875054;	// RHO = 0.4378/sigma^3 : 51x
              else
                if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                  tailcorr = -0.0004941981444;	// RHO = 0.4899/sigma^3 : 61x
                else
                {
                  tailcorr = 0.0;
                  app_log() << warningmsg << std::endl;
                }
    }
    else
      if (N0 == 64)
      {
        if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
          tailcorr = -0.0002651066072;	// RHO = 0.3648/sigma^3 : 02x
        else
          if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
            tailcorr = -0.0002137701388;	// RHO = 0.3280/sigma^3 : 22x
          else
            if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
              tailcorr = -0.0002395972888;	// RHO = 0.3470/sigma^3 : 32x
            else
              if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                tailcorr = -0.0003209410359;	// RHO = 0.4009/sigma^3 : 42x
              else
                if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                  tailcorr = -0.0003835898243;	// RHO = 0.4378/sigma^3 : 52x
                else
                  if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                    tailcorr = -0.0004817718362;	// RHO = 0.4899/sigma^3 : 62x
                  else
                  {
                    tailcorr = 0.0;
                    app_log() << warningmsg << std::endl;
                  }
      }
      else
        if (N0 == 128)
        {
          if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
            tailcorr = -0.0002618281770;	// RHO = 0.3648/sigma^3 : 03x
          else
            if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
              tailcorr = -0.0002113168481;	// RHO = 0.3280/sigma^3 : 23x
            else
              if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
                tailcorr = -0.0002367366424;	// RHO = 0.3470/sigma^3 : 33x
              else
                if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                  tailcorr = -0.0003166990570;	// RHO = 0.4009/sigma^3 : 43x
                else
                  if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                    tailcorr = -0.0003781940744;	// RHO = 0.4378/sigma^3 : 53x
                  else
                    if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                      tailcorr = -0.0004744316477;	// RHO = 0.4899/sigma^3 : 63x
                    else
                    {
                      tailcorr = 0.0;
                      app_log() << warningmsg << std::endl;
                    }
        }
        else
          if (N0 == 256)
          {
            if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
              tailcorr = -0.0002598285275;	// RHO = 0.3648/sigma^3 : 04x
            else
              if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
                tailcorr = -0.0002098169307;	// RHO = 0.3280/sigma^3 : 24x
              else
                if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
                  tailcorr = -0.0002349898273;	// RHO = 0.3470/sigma^3 : 34x
                else
                  if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                    tailcorr = -0.0003141175814;	// RHO = 0.4009/sigma^3 : 44x
                  else
                    if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                      tailcorr = -0.0003749179381;	// RHO = 0.4378/sigma^3 : 54x
                    else
                      if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                        tailcorr = -0.0004699888991;	// RHO = 0.4899/sigma^3 : 64x
                      else
                      {
                        tailcorr = 0.0;
                        app_log() << warningmsg << std::endl;
                      }
          }
          else
          {
            tailcorr = 0.0;
            app_log() << warningmsg << std::endl;
          }
//      tailcorr = -0.000274497151179;	// N = 32, RHO = 0.022 /Angstrom^3 : 11x
    // N = 32, rho = 0.022 Angstroms^-3 = 0.003260063604 a0^-3
    // Will be different for other densities.
    // Generalize later.
//      tailcorr = -0.0002705961159;	// N = 32, RHO = 0.3648 sigma^2 : 01x
//      tailcorr = -0.0002651066072;	// N = 64, RHO = 0.3648 sigma^2 : 02x
//      tailcorr = -0.0002618281770;	// N = 128, RHO = 0.3648 sigma^2 : 03x
//      tailcorr = -0.0002598285275;	// N = 256, RHO = 0.3648 sigma^2 : 04x
//      tailcorr = -0.0002178620081;	// N = 32, RHO = 0.3280 sigma^2 : 21x
//      tailcorr = -0.0002137701388;	// N = 64, RHO = 0.3280 sigma^2 : 22x
//      tailcorr = -0.0002113168481;	// N = 128, RHO = 0.3280 sigma^2 : 23x
//      tailcorr = -0.0002098169307;	// N = 256, RHO = 0.3280 sigma^2 : 24x
//      tailcorr = -0.0002443783556;	// N = 32, RHO = 0.3470 sigma^2 : 31x
//      tailcorr = -0.0002395972888;	// N = 64, RHO = 0.3470 sigma^2 : 32x
//      tailcorr = -0.0002367366424;	// N = 128, RHO = 0.3470 sigma^2 : 33x
//      tailcorr = -0.0002349898273;	// N = 256, RHO = 0.3470 sigma^2 : 34x
//      tailcorr = -0.0003280694504;	// N = 32, RHO = 0.4009 sigma^2 : 41x
//      tailcorr = -0.0003209410359;	// N = 64, RHO = 0.4009 sigma^2 : 42x
//      tailcorr = -0.0003166990570;	// N = 128, RHO = 0.4009 sigma^2 : 43x
//      tailcorr = -0.0003141175814;	// N = 256, RHO = 0.4009 sigma^2 : 44x
//      tailcorr = -0.0003926875054;	// N = 32, RHO = 0.4378 sigma^2 : 51x
//      tailcorr = -0.0003835898243;	// N = 64, RHO = 0.4378 sigma^2 : 52x
//      tailcorr = -0.0003781940744;	// N = 128, RHO = 0.4378 sigma^2 : 53x
//      tailcorr = -0.0003749179381;	// N = 256, RHO = 0.4378 sigma^2 : 54x
//      tailcorr = -0.0004941981444;	// N = 32, RHO = 0.4899 sigma^2 : 61x
//      tailcorr = -0.0004817718362;	// N = 64, RHO = 0.4899 sigma^2 : 62x
//      tailcorr = -0.0004744316477;	// N = 128, RHO = 0.4899 sigma^2 : 63x
//      tailcorr = -0.0004699888991;	// N = 256, RHO = 0.4899 sigma^2 : 64x
    app_log() << "  HeSAPT_smoothed_aux tail correction is " << tailcorr << std::endl;
  }

  ~HeSAPT_smoothed_aux() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    Return_t rho = phyH->PtclRef->G.size()/phyH->PtclRef->Lattice.Volume, N0 = phyH->PtclRef->G.size(), rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    std::string warningmsg = std::string("  Warning: Tail correction is set to ZERO because ...\n")
                        + "  No pre-calculated tail correction value is found for the current\n"
                        + "  system size and density.  Please do a separate calculation\n"
                        + "  of the tail correction and add it onto the LocalEnergy.\n"
                        + "  And if necessary, edit QMCHamiltonians/HeSAPT_smoothed.h\n"
                        + "  to add the new tailcorr value to the database.";
    if (N0 == 32)
    {
      if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
        tailcorr = -0.0002705961159;	// RHO = 0.3648/sigma^3 : 01x
      else
        if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
          tailcorr = -0.0002178620081;	// RHO = 0.3280/sigma^3 : 21x
        else
          if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
            tailcorr = -0.0002443783556;	// RHO = 0.3470/sigma^3 : 31x
          else
            if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
              tailcorr = -0.0003280694504;	// RHO = 0.4009/sigma^3 : 41x
            else
              if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                tailcorr = -0.0003926875054;	// RHO = 0.4378/sigma^3 : 51x
              else
                if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                  tailcorr = -0.0004941981444;	// RHO = 0.4899/sigma^3 : 61x
                else
                {
                  tailcorr = 0.0;
                  app_log() << warningmsg << std::endl;
                }
    }
    else
      if (N0 == 64)
      {
        if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
          tailcorr = -0.0002651066072;	// RHO = 0.3648/sigma^3 : 02x
        else
          if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
            tailcorr = -0.0002137701388;	// RHO = 0.3280/sigma^3 : 22x
          else
            if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
              tailcorr = -0.0002395972888;	// RHO = 0.3470/sigma^3 : 32x
            else
              if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                tailcorr = -0.0003209410359;	// RHO = 0.4009/sigma^3 : 42x
              else
                if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                  tailcorr = -0.0003835898243;	// RHO = 0.4378/sigma^3 : 52x
                else
                  if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                    tailcorr = -0.0004817718362;	// RHO = 0.4899/sigma^3 : 62x
                  else
                  {
                    tailcorr = 0.0;
                    app_log() << warningmsg << std::endl;
                  }
      }
      else
        if (N0 == 128)
        {
          if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
            tailcorr = -0.0002618281770;	// RHO = 0.3648/sigma^3 : 03x
          else
            if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
              tailcorr = -0.0002113168481;	// RHO = 0.3280/sigma^3 : 23x
            else
              if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
                tailcorr = -0.0002367366424;	// RHO = 0.3470/sigma^3 : 33x
              else
                if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                  tailcorr = -0.0003166990570;	// RHO = 0.4009/sigma^3 : 43x
                else
                  if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                    tailcorr = -0.0003781940744;	// RHO = 0.4378/sigma^3 : 53x
                  else
                    if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                      tailcorr = -0.0004744316477;	// RHO = 0.4899/sigma^3 : 63x
                    else
                    {
                      tailcorr = 0.0;
                      app_log() << warningmsg << std::endl;
                    }
        }
        else
          if (N0 == 256)
          {
            if (std::abs(rho - 3.237245782e-3) < 1.0e-7)
              tailcorr = -0.0002598285275;	// RHO = 0.3648/sigma^3 : 04x
            else
              if (std::abs(rho - 2.910681515e-3) < 1.0e-7)
                tailcorr = -0.0002098169307;	// RHO = 0.3280/sigma^3 : 24x
              else
                if (std::abs(rho - 3.079288066e-3) < 1.0e-7)
                  tailcorr = -0.0002349898273;	// RHO = 0.3470/sigma^3 : 34x
                else
                  if (std::abs(rho - 3.557598229e-3) < 1.0e-7)
                    tailcorr = -0.0003141175814;	// RHO = 0.4009/sigma^3 : 44x
                  else
                    if (std::abs(rho - 3.885049900e-3) < 1.0e-7)
                      tailcorr = -0.0003749179381;	// RHO = 0.4378/sigma^3 : 54x
                    else
                      if (std::abs(rho - 4.347386811e-3) < 1.0e-7)
                        tailcorr = -0.0004699888991;	// RHO = 0.4899/sigma^3 : 64x
                      else
                      {
                        tailcorr = 0.0;
                        app_log() << warningmsg << std::endl;
                      }
          }
          else
          {
            tailcorr = 0.0;
            app_log() << warningmsg << std::endl;
          }
    app_log() << "  HeSAPT_smoothed_aux tail correction is " << tailcorr << std::endl;
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value = -(phyH->smooth) + tailcorr;
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HeSAPT_smoothed_aux: " << phyH->PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return 0;
  }
};
}
#endif
