//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Configuration.h"
#include "CuspCorrection.h"
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaAtomicBasisSet.h"
#include "QMCWaveFunctions/lcao/MultiQuinticSpline1D.h"
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"


namespace qmcplusplus
{
bool readCuspInfo(const std::string& cuspInfoFile,
                  const std::string& objectName,
                  int OrbitalSetSize,
                  Matrix<CuspCorrectionParameters>& info)
{
  bool success = true;
  std::string cname;
  int ncenter = info.rows();
  int nOrbs   = info.cols();
  app_log() << "Reading cusp info from : " << cuspInfoFile << std::endl;
  Libxml2Document adoc;
  if (!adoc.parse(cuspInfoFile))
  {
    app_log() << "Could not find precomputed cusp data for spo set: " << objectName << std::endl;
    app_log() << "Recalculating data.\n";
    return false;
  }
  xmlNodePtr head = adoc.getRoot();
  head            = head->children;
  xmlNodePtr cur  = NULL, ctr;
  while (head != NULL)
  {
    getNodeName(cname, head);
    if (cname == "sposet")
    {
      std::string name;
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add(name, "name");
      spoAttrib.put(head);
      if (name == objectName)
      {
        cur = head;
        break;
      }
    }
    head = head->next;
  }
  if (cur == NULL)
  {
    app_log() << "Could not find precomputed cusp data for spo set: " << objectName << std::endl;
    app_log() << "Recalculating data.\n";
    return false;
  }
  else
  {
    app_log() << "Found precomputed cusp data for spo set: " << objectName << std::endl;
  }
  cur = cur->children;
  while (cur != NULL)
  {
    getNodeName(cname, cur);
    if (cname == "center")
    {
      int num = -1;
      OhmmsAttributeSet Attrib;
      Attrib.add(num, "num");
      Attrib.put(cur);
      if (num < 0 || num >= ncenter)
      {
        APP_ABORT("Error with cusp info xml block. incorrect center number. \n");
      }
      ctr = cur->children;
      while (ctr != NULL)
      {
        getNodeName(cname, ctr);
        if (cname == "orbital")
        {
          int orb = -1;
          OhmmsAttributeSet orbAttrib;
          QMCTraits::RealType a1(0.0), a2, a3, a4, a5, a6, a7, a8, a9;
          orbAttrib.add(orb, "num");
          orbAttrib.add(a1, "redo");
          orbAttrib.add(a2, "C");
          orbAttrib.add(a3, "sg");
          orbAttrib.add(a4, "rc");
          orbAttrib.add(a5, "a1");
          orbAttrib.add(a6, "a2");
          orbAttrib.add(a7, "a3");
          orbAttrib.add(a8, "a4");
          orbAttrib.add(a9, "a5");
          orbAttrib.put(ctr);
          if (orb < OrbitalSetSize)
          {
            info(num, orb).redo     = a1;
            info(num, orb).C        = a2;
            info(num, orb).sg       = a3;
            info(num, orb).Rc       = a4;
            info(num, orb).alpha[0] = a5;
            info(num, orb).alpha[1] = a6;
            info(num, orb).alpha[2] = a7;
            info(num, orb).alpha[3] = a8;
            info(num, orb).alpha[4] = a9;
          }
        }
        ctr = ctr->next;
      }
    }
    cur = cur->next;
  }
  return success;
}

void splitPhiEta(int center, const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi, LCAOrbitalSet& Eta)
{
  typedef QMCTraits::RealType RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);
  std::vector<bool> correct_this_center(corrCenter.size(), false);
  correct_this_center[center] = corrCenter[center];

  Phi.myBasisSet->queryOrbitalsForSType(correct_this_center, is_s_orbital);

  int nOrbs = Phi.OrbitalSetSize;
  int bss   = Phi.BasisSetSize;

  for (int i = 0; i < bss; i++)
  {
    if (is_s_orbital[i])
    {
      auto& cref(*(Eta.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0; //Eta->C(k,i) = 0.0;
    }
    else
    {
      auto& cref(*(Phi.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0; //Phi->C(k,i) = 0.0;
    }
  }
}

void removeSTypeOrbitals(const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi)
{
  typedef QMCTraits::RealType RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);

  Phi.myBasisSet->queryOrbitalsForSType(corrCenter, is_s_orbital);

  int nOrbs = Phi.OrbitalSetSize;
  int bss   = Phi.BasisSetSize;

  for (int i = 0; i < bss; i++)
  {
    if (is_s_orbital[i])
    {
      auto& cref(*(Phi.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0;
    }
  }
}


// Will be the corrected value for r < rc and the original wavefunction for r > rc
void computeRadialPhiBar(ParticleSet* targetP,
                         ParticleSet* sourceP,
                         int curOrb_,
                         int curCenter_,
                         SPOSet* Phi,
                         Vector<QMCTraits::RealType>& xgrid,
                         Vector<QMCTraits::RealType>& rad_orb,
                         const CuspCorrectionParameters& data)
{
  CuspCorrection cusp(targetP, sourceP);
  cusp.setPsi(Phi);
  cusp.cparam    = data;
  cusp.curOrb    = curOrb_;
  cusp.curCenter = curCenter_;

  for (int i = 0; i < xgrid.size(); i++)
  {
    rad_orb[i] = cusp.phiBar(xgrid[i]);
  }
}


}; // namespace qmcplusplus
