//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2029 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/LCAO/LCAOSpinorBuilder.h"

namespace qmcplusplus
{
SPOSet* LCAOSpinorBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  /*
  ReportEngine PRE(ClassName, "createSPO(xmlNodePtr)");
  std::string spo_name(""), id, cusp_file(""), optimize("no");
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(spo_name, "name");
  spoAttrib.add(id, "id");
  spoAttrib.add(cusp_file, "cuspInfo");
  spoAttrib.add(optimize, "optimize");
  spoAttrib.put(cur);

  if (myBasisSet == nullptr)
    PRE.error("Missing basisset.", true);

  if (optimize == "yes")
    app_log() << "  SPOSet " << spo_name << " is optimizable\n";

  LCAOrbitalSet* lcos = nullptr;
#if !defined(QMC_COMPLEX)
  LCAOrbitalSetWithCorrection* lcwc = nullptr;
  if (doCuspCorrection)
    lcos = lcwc = new LCAOrbitalSetWithCorrection(sourcePtcl, targetPtcl, myBasisSet, optimize == "yes");
  else
    lcos = new LCAOrbitalSet(myBasisSet, optimize == "yes");
#else
  lcos = new LCAOrbitalSet(myBasisSet, optimize == "yes");
#endif
  loadMO(*lcos, cur);

#if !defined(QMC_COMPLEX)
  if (doCuspCorrection)
  {
    int num_centers = sourcePtcl.getTotalNum();

    // Sometimes sposet attribute is 'name' and sometimes it is 'id'
    if (id == "")
      id = spo_name;

    int orbital_set_size = lcos->getOrbitalSetSize();
    Matrix<CuspCorrectionParameters> info(num_centers, orbital_set_size);

    bool valid = false;
    if (myComm->rank() == 0)
    {
      valid = readCuspInfo(cusp_file, id, orbital_set_size, info);
    }
#ifdef HAVE_MPI
    myComm->comm.broadcast_value(valid);
    if (valid)
    {
      for (int orb_idx = 0; orb_idx < orbital_set_size; orb_idx++)
      {
        for (int center_idx = 0; center_idx < num_centers; center_idx++)
        {
          broadcastCuspInfo(info(center_idx, orb_idx), *myComm, 0);
        }
      }
    }
#endif
    if (!valid)
    {
      generateCuspInfo(orbital_set_size, num_centers, info, targetPtcl, sourcePtcl, *lcwc, id, *myComm);
    }

    applyCuspCorrection(info, num_centers, orbital_set_size, targetPtcl, sourcePtcl, *lcwc, id);
  }
#endif

  return lcos;
  */
}

} // namespace qmcplusplus
