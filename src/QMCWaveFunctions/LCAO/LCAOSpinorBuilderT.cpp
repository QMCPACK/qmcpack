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

#include "LCAOSpinorBuilderT.h"

#include "Message/CommOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/SpinorSetT.h"
#include "Utilities/ProgressReportEngine.h"
#include "hdf/hdf_archive.h"

namespace qmcplusplus
{
template<class T>
LCAOSpinorBuilderT<T>::LCAOSpinorBuilderT(ParticleSetT<T>& els,
                                          ParticleSetT<T>& ions,
                                          Communicate* comm,
                                          xmlNodePtr cur)
    : LCAOrbitalBuilderT<T>(els, ions, comm, cur)
{
  this->ClassName = "LCAOSpinorBuilder";

  if (this->h5_path == "")
    this->myComm->barrier_and_abort("LCAOSpinorBuilder only works with href");
}

template<class T>
std::unique_ptr<SPOSetT<T>> LCAOSpinorBuilderT<T>::createSPOSetFromXML(xmlNodePtr cur)
{
  ReportEngine PRE(this->ClassName, "createSPO(xmlNodePtr)");
  std::string spo_name(""), optimize("no");
  std::string basisset_name("LCAOBSet");
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(spo_name, "name");
  spoAttrib.add(optimize, "optimize");
  spoAttrib.add(basisset_name, "basisset");
  spoAttrib.put(cur);

  BasisSet_t* myBasisSet = nullptr;
  if (this->basisset_map_.find(basisset_name) == this->basisset_map_.end())
    this->myComm->barrier_and_abort("basisset \"" + basisset_name + "\" cannot be found\n");
  else
    myBasisSet = this->basisset_map_[basisset_name].get();

  if (optimize == "yes")
    app_log() << "  SPOSet " << spo_name << " is optimizable\n";

  auto upspo =
      std::make_unique<LCAOrbitalSetT<T>>(spo_name + "_up", std::unique_ptr<BasisSet_t>(myBasisSet->makeClone()));
  auto dnspo =
      std::make_unique<LCAOrbitalSetT<T>>(spo_name + "_dn", std::unique_ptr<BasisSet_t>(myBasisSet->makeClone()));

  loadMO(*upspo, *dnspo, cur);

  // create spinor and register up/dn
  auto spinor_set = std::make_unique<SpinorSetT<T>>(spo_name);
  spinor_set->set_spos(std::move(upspo), std::move(dnspo));
  return spinor_set;
}

template<class T>
bool LCAOSpinorBuilderT<T>::loadMO(LCAOrbitalSetT<T>& up, LCAOrbitalSetT<T>& dn, xmlNodePtr cur)
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

  hdf_archive hin(this->myComm);
  if (this->myComm->rank() == 0)
  {
    if (!hin.open(this->h5_path, H5F_ACC_RDONLY))
      this->myComm->barrier_and_abort("LCAOSpinorBuilder::loadMO missing "
                                      "or incorrect path to H5 file.");
    hin.push("PBC");
    PBC = false;
    hin.read(PBC, "PBC");
    hin.close();
  }
  this->myComm->bcast(PBC);
  if (PBC)
    this->myComm->barrier_and_abort("LCAOSpinorBuilder::loadMO lcao spinors not implemented in PBC");

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

template<class T>
bool LCAOSpinorBuilderT<T>::putFromH5(LCAOrbitalSetT<T>& up, LCAOrbitalSetT<T>& dn, xmlNodePtr occ_ptr)
{
  if (up.getBasisSetSize() == 0 || dn.getBasisSetSize() == 0)
  {
    this->myComm->barrier_and_abort("LCASpinorBuilder::loadMO  detected ZERO BasisSetSize");
    return false;
  }

  bool success = true;
  hdf_archive hin(this->myComm);
  if (this->myComm->rank() == 0)
  {
    if (!hin.open(this->h5_path, H5F_ACC_RDONLY))
      this->myComm->barrier_and_abort("LCAOSpinorBuilder::putFromH5 missing or "
                                      "incorrect path to H5 file");

    Matrix<RealType> upReal;
    Matrix<RealType> upImag;
    std::string setname = "/Super_Twist/eigenset_0";
    this->readRealMatrixFromH5(hin, setname, upReal);
    setname += "_imag";
    this->readRealMatrixFromH5(hin, setname, upImag);


    assert(upReal.rows() == upImag.rows());
    assert(upReal.cols() == upImag.cols());

    Matrix<ValueType> upTemp(upReal.rows(), upReal.cols());
    for (int i = 0; i < upTemp.rows(); i++)
    {
      for (int j = 0; j < upTemp.cols(); j++)
      {
        upTemp[i][j] = ValueType{upReal[i][j], upImag[i][j]};
      }
    }

    Matrix<RealType> dnReal;
    Matrix<RealType> dnImag;
    setname = "/Super_Twist/eigenset_1";
    this->readRealMatrixFromH5(hin, setname, dnReal);
    setname += "_imag";
    this->readRealMatrixFromH5(hin, setname, dnImag);

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

    this->Occ.resize(upReal.rows());
    success = this->putOccupation(up, occ_ptr);

    int norbs = up.getOrbitalSetSize();

    int n = 0, i = 0;
    while (i < norbs)
    {
      if (this->Occ[n] > 0.0)
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
  this->myComm->comm.broadcast_n(up.C->data(), up.C->size());
  this->myComm->comm.broadcast_n(dn.C->data(), dn.C->size());
#endif

  return success;
}

#ifdef QMC_COMPLEX
#ifndef MIXED_PRECISION
template class LCAOSpinorBuilderT<std::complex<double>>;
#else
template class LCAOSpinorBuilderT<std::complex<float>>;
#endif
#endif
} // namespace qmcplusplus
