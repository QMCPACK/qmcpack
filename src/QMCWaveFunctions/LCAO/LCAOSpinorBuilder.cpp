//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "LCAOSpinorBuilder.h"
#include "QMCWaveFunctions/SpinorSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/ProgressReportEngine.h"
#include "hdf/hdf_archive.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{
LCAOSpinorBuilder::LCAOSpinorBuilder(ParticleSet& els, ParticleSet& ions, Communicate* comm, xmlNodePtr cur)
    : LCAOrbitalBuilder(els, ions, comm, cur)
{
  ClassName = "LCAOSpinorBuilder";

  if (h5_path == "")
  {
    myComm->barrier_and_abort("LCAOSpinorBuilder only works with href");
  }
}

SPOSet* LCAOSpinorBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createSPO(xmlNodePtr)");
  std::string spo_name(""), optimize("no");
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(spo_name, "name");
  spoAttrib.add(optimize, "optimize");
  spoAttrib.put(cur);

  if (myBasisSet == nullptr)
    PRE.error("Missing basisset.", true);

  if (optimize == "yes")
    app_log() << "  SPOSet " << spo_name << " is optimizable\n";

  std::unique_ptr<LCAOrbitalSet> upspo = std::make_unique<LCAOrbitalSet>(std::unique_ptr<BasisSet_t>(myBasisSet->makeClone()), optimize == "yes");
  std::unique_ptr<LCAOrbitalSet> dnspo = std::make_unique<LCAOrbitalSet>(std::unique_ptr<BasisSet_t>(myBasisSet->makeClone()), optimize == "yes");

  loadMO(*upspo, *dnspo, cur);

  //create spinor and register up/dn
  SpinorSet* spinor_set = new SpinorSet();
  spinor_set->set_spos(std::move(upspo), std::move(dnspo));

  return spinor_set;
}

bool LCAOSpinorBuilder::loadMO(LCAOrbitalSet& up, LCAOrbitalSet& dn, xmlNodePtr cur)
{
  bool PBC = false;
  int norb = up.getBasisSetSize();
  std::string debugc("no");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb, "size");
  aAttrib.add(debugc, "debug");
  aAttrib.put(cur);

  up.setOrbitalSetSize(norb);
  dn.setOrbitalSetSize(norb);

  xmlNodePtr occ_ptr = nullptr;
  cur                = cur->xmlChildrenNode;
  while (cur != nullptr)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "occupation")
    {
      occ_ptr = cur;
    }
    cur = cur->next;
  }

  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      myComm->barrier_and_abort("LCAOSpinorBuilder::loadMO missing or incorrect path to H5 file.");
    hin.push("PBC");
    PBC = false;
    hin.read(PBC, "PBC");
    hin.close();
  }
  myComm->bcast(PBC);
  if (PBC)
    myComm->barrier_and_abort("LCAOSpinorBuilder::loadMO lcao spinors not implemented in PBC");

  bool success = putFromH5(up, dn, occ_ptr);


  if (debugc == "yes")
  {
    app_log() << "UP:  Single-particle orbital coefficients dims=" << up.C->rows() << " x " << up.C->cols()
              << std::endl;
    app_log() << *up.C << std::endl;
    app_log() << "DN:  Single-particle orbital coefficients dims=" << dn.C->rows() << " x " << dn.C->cols()
              << std::endl;
    app_log() << *dn.C << std::endl;
  }
  return success;
}

bool LCAOSpinorBuilder::putFromH5(LCAOrbitalSet& up, LCAOrbitalSet& dn, xmlNodePtr occ_ptr)
{
#ifdef QMC_COMPLEX
#if defined(HAVE_LIBHDF5)
  if (up.getBasisSetSize() == 0 || dn.getBasisSetSize() == 0)
  {
    myComm->barrier_and_abort("LCASpinorBuilder::loadMO  detected ZERO BasisSetSize");
    return false;
  }

  int norbs = up.getOrbitalSetSize();

  bool success = true;
  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      myComm->barrier_and_abort("LCAOSpinorBuilder::putFromH5 missing or incorrect path to H5 file");

    std::string setname;
    Matrix<RealType> upReal;
    Matrix<RealType> upImag;
    setname = "/KPTS_0/eigenset_0";
    readRealMatrixFromH5(hin, setname, upReal);
    setname = "/KPTS_0/eigenset_0_imag";
    readRealMatrixFromH5(hin, setname, upImag);

    assert(upReal.rows() == upImag.rows());
    assert(upReal.cols() == upImag.cols());

    Matrix<ValueType> upTemp(upReal.rows(), upReal.cols());
    for (int i = 0; i < upTemp.rows(); i++)
    {
      for (int j = 0; j < upTemp.cols(); j++)
      {
        upTemp[i][j] = ValueType(upReal[i][j], upImag[i][j]);
      }
    }

    Matrix<RealType> dnReal;
    Matrix<RealType> dnImag;
    setname = "/KPTS_0/eigenset_1";
    readRealMatrixFromH5(hin, setname, dnReal);
    setname = "/KPTS_0/eigenset_1_imag";
    readRealMatrixFromH5(hin, setname, dnImag);

    assert(dnReal.rows() == dnImag.rows());
    assert(dnReal.cols() == dnImag.cols());

    Matrix<ValueType> dnTemp(dnReal.rows(), dnReal.cols());
    for (int i = 0; i < dnTemp.rows(); i++)
    {
      for (int j = 0; j < dnTemp.cols(); j++)
      {
        dnTemp[i][j] = ValueType(dnReal[i][j], dnImag[i][j]);
      }
    }

    assert(upReal.rows() == dnReal.rows());
    assert(upReal.cols() == dnReal.cols());

    Occ.resize(upReal.rows());
    success = putOccupation(norbs, occ_ptr);

    int n = 0, i = 0;
    while (i < norbs)
    {
      if (Occ[n] > 0.0)
      {
        std::copy(upTemp[n], upTemp[n + 1], (*up.C)[i]);
        std::copy(dnTemp[n], dnTemp[n + 1], (*dn.C)[i]);
        i++;
      }
      n++;
    }

    hin.close();
  }

#ifdef HAVE_MPI
  myComm->comm.broadcast_n(up.C->data(), up.C->size());
  myComm->comm.broadcast_n(dn.C->data(), dn.C->size());
#endif

#else
  myComm->barrier_and_abort("LCAOSpinorBuilder::putFromH5 HDF5 is disabled");
#endif

#else
  myComm->barrier_and_abort("LCAOSpinorBuilder::putFromH5 Must build with QMC_COMPLEX");
#endif

  return success;
}

bool LCAOSpinorBuilder::putOccupation(int norbs, xmlNodePtr occ_ptr)
{
  Occ = 0.0;
  for (int i = 0; i < norbs; i++)
    Occ[i] = 1.0;
  std::vector<int> occ_in;
  std::string occ_mode("table");
  if (occ_ptr == nullptr)
  {
    occ_mode = "ground";
  }
  else
  {
    const XMLAttrString o(occ_ptr, "mode");
    if (!o.empty())
      occ_mode = o;
  }
  if (occ_mode == "excited")
  {
    putContent(occ_in, occ_ptr);
    for (int k = 0; k < occ_in.size(); k++)
    {
      if (occ_in[k] < 0)
        Occ[-occ_in[k] - 1] = 0.0;
      else
        Occ[occ_in[k] - 1] = 1.0;
    }
  }
  else if (occ_mode == "table")
  {
    putContent(Occ, occ_ptr);
  }
  return true;
}

} // namespace qmcplusplus
