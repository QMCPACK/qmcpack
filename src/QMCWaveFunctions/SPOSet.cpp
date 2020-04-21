//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/SPOSet.h"
#include "Message/Communicate.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"
#include <simd/simd.hpp>
#include "Utilities/ProgressReportEngine.h"
#include <io/hdf_archive.h>
#include <limits>

namespace qmcplusplus
{
SPOSet::SPOSet(bool ion_deriv, bool optimizable)
    :
#if !defined(ENABLE_SOA)
      Identity(false),
      BasisSetSize(0),
      C(nullptr),
#endif
      ionDerivs(ion_deriv),
      Optimizable(optimizable),
      OrbitalSetSize(0)
{
  className = "invalid";
#if !defined(ENABLE_SOA)
  IsCloned = false;
  //default is false: LCOrbitalSet.h needs to set this true and recompute needs to check
  myComm = nullptr;
#endif
}

void SPOSet::evaluate(const ParticleSet& P, PosType& r, ValueVector_t& psi)
{
  APP_ABORT("Need specialization for SPOSet::evaluate(const ParticleSet& P, PosType &r)\n");
}

void SPOSet::mw_evaluateValue(const RefVector<SPOSet>& spo_list,
                              const RefVector<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector_t>& psi_v_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateValue(P_list[iw], iat, psi_v_list[iw]);
}

void SPOSet::evaluateDetRatios(const VirtualParticleSet& VP,
                               ValueVector_t& psi,
                               const ValueVector_t& psiinv,
                               std::vector<ValueType>& ratios)
{
  assert(psi.size() == psiinv.size());
  for (int iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    evaluateValue(VP, iat, psi);
    ratios[iat] = simd::dot(psi.data(), psiinv.data(), psi.size());
  }
}

void SPOSet::mw_evaluateDetRatios(const RefVector<SPOSet>& spo_list,
                                  const RefVector<const VirtualParticleSet>& vp_list,
                                  const RefVector<ValueVector_t>& psi_list,
                                  const RefVector<const ValueVector_t>& psiinv_list,
                                  std::vector<std::vector<ValueType>>& ratios_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateDetRatios(vp_list[iw], psi_list[iw], psiinv_list[iw], ratios_list[iw]);

}

void SPOSet::mw_evaluateVGL(const RefVector<SPOSet>& spo_list,
                            const RefVector<ParticleSet>& P_list,
                            int iat,
                            const RefVector<ValueVector_t>& psi_v_list,
                            const RefVector<GradVector_t>& dpsi_v_list,
                            const RefVector<ValueVector_t>& d2psi_v_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateVGL(P_list[iw], iat, psi_v_list[iw], dpsi_v_list[iw], d2psi_v_list[iw]);
}

void SPOSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluateThirdDeriv(). \n");
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix_t& logdet,
                                  GradMatrix_t& dlogdet,
                                  HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix_t& logdet,
                                  GradMatrix_t& dlogdet,
                                  HessMatrix_t& grad_grad_logdet,
                                  GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


SPOSet* SPOSet::makeClone() const
{
  APP_ABORT("Missing  SPOSet::makeClone for " + className);
  return 0;
}

#if !defined(ENABLE_SOA)
bool SPOSet::setIdentity(bool useIdentity)
{
  Identity = useIdentity;
  if (Identity)
    return true;

  if (C == nullptr && (OrbitalSetSize > 0) && (BasisSetSize > 0))
  {
    C = new ValueMatrix_t(OrbitalSetSize, BasisSetSize);
  }
  else
  {
    app_error() << "either OrbitalSetSize or BasisSetSize has an invalid value !!\n";
    app_error() << "OrbitalSetSize = " << OrbitalSetSize << std::endl;
    app_error() << "BasisSetSize = " << BasisSetSize << std::endl;
    APP_ABORT("SPOSet::setIdentiy ");
  }

  return true;
}

/** Parse the xml file for information on the Dirac determinants.
 *@param cur the current xmlNode
 */
bool SPOSet::put(xmlNodePtr cur)
{
#undef FunctionName
#define FunctionName                                      \
  printf("Calling FunctionName from %s\n", __FUNCTION__); \
  FunctionNameReal
  //Check if HDF5 present
  ReportEngine PRE("SPOSet", "put(xmlNodePtr)");

  //Special case for sposet hierarchy: go up only once.
  OhmmsAttributeSet locAttrib;
  std::string cur_name;
  locAttrib.add(cur_name, "name");
  locAttrib.put(cur);
  xmlNodePtr curtemp;
  if (cur_name == "spo-up" || cur_name == "spo-dn")
    curtemp = cur->parent;
  else
    curtemp = cur->parent->parent;

  std::string MOtype, MOhref;
  bool H5file = false;
  OhmmsAttributeSet H5checkAttrib;
  H5checkAttrib.add(MOtype, "type");
  H5checkAttrib.add(MOhref, "href");
  H5checkAttrib.put(curtemp);
  std::string MOhref2;
  if (MOtype == "MolecularOrbital" && MOhref != "")
  {
    MOhref2 = XMLAttrString(curtemp, "href");
    H5file  = true;
    PRE.echo(curtemp);
  }

  //initialize the number of orbital by the basis set size
  int norb = BasisSetSize;
  std::string debugc("no");
  bool PBC = false;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb, "orbitals");
  aAttrib.add(norb, "size");
  aAttrib.add(debugc, "debug");
  aAttrib.put(cur);
  setOrbitalSetSize(norb);
  xmlNodePtr occ_ptr   = NULL;
  xmlNodePtr coeff_ptr = NULL;
  cur                  = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "occupation")
    {
      occ_ptr = cur;
    }
    else if (cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var")
    {
      coeff_ptr = cur;
    }
    cur = cur->next;
  }
  if (coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << std::endl;
    return setIdentity(true);
  }
  bool success = putOccupation(occ_ptr);
  if (H5file == false)
    success = putFromXML(coeff_ptr);
  else
  {
    hdf_archive hin(myComm);

    if (myComm->rank() == 0)
    {
      if (!hin.open(MOhref2, H5F_ACC_RDONLY))
        APP_ABORT("SPOSet::putFromH5 missing or incorrect path to H5 file.");
      //TO REVIEWERS:: IDEAL BEHAVIOUR SHOULD BE:
      /*
       if(!hin.push("PBC")
           PBC=false;
       else
          if (!hin.read(PBC,"PBC"))
              APP_ABORT("Could not read PBC dataset in H5 file. Probably corrupt file!!!.");
      // However, it always succeeds to enter the if condition even if the group does not exists...
      */
      hin.push("PBC");
      PBC = false;
      hin.read(PBC, "PBC");
      hin.close();
      if (PBC)
        APP_ABORT("SPOSet::putFromH5 PBC is not supported by AoS builds");
    }
    myComm->bcast(PBC);
    success = putFromH5(MOhref2, coeff_ptr);
  }

  bool success2 = transformSPOSet();
  if (debugc == "yes")
  {
    app_log() << "   Single-particle orbital coefficients dims=" << C->rows() << " x " << C->cols() << std::endl;
    app_log() << C << std::endl;
  }

  return success && success2;
}

void SPOSet::checkObject()
{
  if (!(OrbitalSetSize == C->rows() && BasisSetSize == C->cols()))
  {
    app_error() << "   SPOSet::checkObject Linear coeffient for SPOSet is not consistent with the input." << std::endl;
    OHMMS::Controller->abort();
  }
}


bool SPOSet::putFromXML(xmlNodePtr coeff_ptr)
{
  Identity  = true;
  int norbs = 0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norbs, "size");
  aAttrib.add(norbs, "orbitals");
  aAttrib.put(coeff_ptr);
  if (norbs < OrbitalSetSize)
  {
    return false;
    APP_ABORT("SPOSet::putFromXML missing or incorrect size");
  }
  if (norbs)
  {
    Identity = false;
    std::vector<ValueType> Ctemp;
    Ctemp.resize(norbs * BasisSetSize);
    setIdentity(Identity);
    putContent(Ctemp, coeff_ptr);
    int n = 0, i = 0;
    std::vector<ValueType>::iterator cit(Ctemp.begin());
    while (i < OrbitalSetSize)
    {
      if (Occ[n] > std::numeric_limits<RealType>::epsilon())
      {
        std::copy(cit, cit + BasisSetSize, (*C)[i]);
        i++;
      }
      n++;
      cit += BasisSetSize;
    }
  }
  return true;
}

/** read data from a hdf5 file
 * @param norb number of orbitals to be initialized
 * @param fname hdf5 file name
 * @param coeff_ptr xmlnode for coefficients
 */
bool SPOSet::putFromH5(const std::string& fname, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  int norbs  = OrbitalSetSize;
  int neigs  = BasisSetSize;
  int setVal = -1;
  std::string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal, "spindataset");
  aAttrib.add(neigs, "size");
  aAttrib.add(neigs, "orbitals");
  aAttrib.put(coeff_ptr);
  setIdentity(false);
  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(fname, H5F_ACC_RDONLY))
      APP_ABORT("SPOSet::putFromH5 missing or incorrect path to H5 file.");

    Matrix<RealType> Ctemp(neigs, BasisSetSize);
    char name[72];
    sprintf(name, "%s%d", "/KPTS_0/eigenset_", setVal);
    setname = name;
    if (!hin.readEntry(Ctemp, setname))
    {
      setname = "SPOSet::putFromH5 Missing " + setname + " from HDF5 File.";
      APP_ABORT(setname.c_str());
    }
    hin.close();

    int n = 0, i = 0;
    while (i < norbs)
    {
      if (Occ[n] > 0.0)
      {
        std::copy(Ctemp[n], Ctemp[n + 1], (*C)[i]);
        i++;
      }
      n++;
    }
  }
  myComm->bcast(C->data(), C->size());
#else
  APP_ABORT("SPOSet::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


bool SPOSet::putOccupation(xmlNodePtr occ_ptr)
{
  //die??
  if (BasisSetSize == 0)
  {
    APP_ABORT("SPOSet::putOccupation detected ZERO BasisSetSize");
    return false;
  }
  Occ.resize(std::max(BasisSetSize, OrbitalSetSize));
  Occ = 0.0;
  for (int i = 0; i < OrbitalSetSize; i++)
    Occ[i] = 1.0;
  std::vector<int> occ_in;
  std::string occ_mode("table");
  if (occ_ptr == NULL)
  {
    occ_mode = "ground";
  }
  else
  {
    const XMLAttrString o(occ_ptr, "mode");
    if (!o.empty())
      occ_mode = o;
  }
  //Do nothing if mode == ground
  if (occ_mode == "excited")
  {
    putContent(occ_in, occ_ptr);
    for (int k = 0; k < occ_in.size(); k++)
    {
      if (occ_in[k] < 0) //remove this, -1 is to adjust the base
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
#endif

void SPOSet::basic_report(const std::string& pad)
{
  app_log() << pad << "size = " << size() << std::endl;
  app_log() << pad << "state info:" << std::endl;
  //states.report(pad+"  ");
  app_log().flush();
}

void SPOSet::evaluateVGH(const ParticleSet& P,
                         int iat,
                         ValueVector_t& psi,
                         GradVector_t& dpsi,
                         HessVector_t& grad_grad_psi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate(P,iat,psi,dpsi,dhpsi) (vector quantities)\n");
}

void SPOSet::evaluateVGHGH(const ParticleSet& P,
                           int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           HessVector_t& grad_grad_psi,
                           GGGVector_t& grad_grad_grad_psi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate(P,iat,psi,dpsi,dhpsi,dghpsi) (vector quantities)\n");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix_t& gradphi)
{
  APP_ABORT("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix_t& grad_phi,
                                HessMatrix_t& grad_grad_phi,
                                GradMatrix_t& grad_lapl_phi)
{
  APP_ABORT("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector_t& psi, ValueVector_t& dpsi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate_spin(P,iat,psi,dpsi) (vector quantities)\n");
}

#ifdef QMC_CUDA

void SPOSet::evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::ValueType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi,
                      gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                      int row_stride)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi,
                      gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                      int row_stride,
                      int k,
                      bool klinear)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

#endif
} // namespace qmcplusplus
