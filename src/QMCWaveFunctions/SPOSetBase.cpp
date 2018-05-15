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


#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"
#include <simd/simd.hpp>
#include "Utilities/ProgressReportEngine.h"
#include <io/hdf_archive.h>
#include <limits>

namespace qmcplusplus
{

SPOSetBase::SPOSetBase()
  :OrbitalSetSize(0),Optimizable(false),ionDerivs(false),builder_index(-1)
#if !defined(ENABLE_SOA)
  ,Identity(false),BasisSetSize(0),C(nullptr)
#endif
{
  CanUseGLCombo=false;
  className="invalid";
#if !defined(ENABLE_SOA)
  IsCloned=false;
  //default is false: LCOrbitalSet.h needs to set this true and recompute needs to check
  myComm=nullptr;
#endif
}

/** default implementation */
SPOSetBase::ValueType
SPOSetBase::RATIO(const ParticleSet& P, int iat, const ValueType* restrict arow)
{
  int ip=omp_get_thread_num();
  // YYYY to fix
  /*
  ValueVector_t psi(t_logpsi[ip],OrbitalSetSize);
  evaluate(P,iat,psi);
  return simd::dot(psi.data(),arow,OrbitalSetSize,ValueType());
  */
  return ValueType();
}

void SPOSetBase::evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl, bool newp)
{
  APP_ABORT("SPOSetBase::evaluateVGL not implemented.");
}

void SPOSetBase::evaluateValues(VirtualParticleSet& VP, ValueMatrix_t& psiM)
{
  for(int iat=0; iat<VP.getTotalNum(); ++iat)
  {
    ValueVector_t psi(psiM[iat],OrbitalSetSize);
    evaluate(VP,iat,psi);
  }
}

void SPOSetBase::evaluateThirdDeriv(const ParticleSet& P, int first, int last,
                                    GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateThirdDeriv(). \n");
}

void SPOSetBase::evaluate_notranspose(const ParticleSet& P, int first, int last
                                      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSetBase::evaluate_notranspose(const ParticleSet& P, int first, int last,
                                      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


SPOSetBase* SPOSetBase::makeClone() const
{
  APP_ABORT("Missing  SPOSetBase::makeClone for "+className);
  return 0;
}

#if !defined(ENABLE_SOA)
bool SPOSetBase::setIdentity(bool useIdentity)
{
  Identity = useIdentity;
  if(Identity) return true;

  if ( C== nullptr && (OrbitalSetSize > 0) && (BasisSetSize > 0) )
  {
    C=new ValueMatrix_t(OrbitalSetSize,BasisSetSize);
  }
  else
  {
    app_error() << "either OrbitalSetSize or BasisSetSize has an invalid value !!\n";
    app_error() << "OrbitalSetSize = " << OrbitalSetSize << std::endl;
    app_error() << "BasisSetSize = " << BasisSetSize << std::endl;
    APP_ABORT("SPOSetBase::setIdentiy ");
  }

  return true;
}

/** Parse the xml file for information on the Dirac determinants.
 *@param cur the current xmlNode
 */
bool SPOSetBase::put(xmlNodePtr cur)
{
  #undef FunctionName
#define FunctionName printf("Calling FunctionName from %s\n",__FUNCTION__);FunctionNameReal
  //Check if HDF5 present
  ReportEngine PRE("SPOSetBase","put(xmlNodePtr)");

  //Special case for sposet hierarchy: go up only once.
  OhmmsAttributeSet locAttrib;
  std::string cur_name;
  locAttrib.add (cur_name, "name");
  locAttrib.put(cur);
  xmlNodePtr curtemp;
  if (cur_name=="spo-up" || cur_name=="spo-dn")
     curtemp=cur->parent;
  else
     curtemp=cur->parent->parent;

  std::string MOtype,MOhref;
  bool H5file=false;
  OhmmsAttributeSet H5checkAttrib;
  H5checkAttrib.add(MOtype,"type");
  H5checkAttrib.add(MOhref,"href");
  H5checkAttrib.put(curtemp);
  xmlChar* MOhreftemp;
  if(MOtype=="MolecularOrbital" && MOhref!="")
  {
     MOhreftemp=xmlGetProp(curtemp, (xmlChar*)"href");
     H5file=true;
     PRE.echo(curtemp);
  }

  const char* MOhref2((const char*)MOhreftemp);

  //initialize the number of orbital by the basis set size
  int norb= BasisSetSize;
  std::string debugc("no");
  double orbital_mix_magnitude = 0.0;
  bool PBC=false;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(norb,"size");
  aAttrib.add(debugc,"debug");
  aAttrib.add(orbital_mix_magnitude, "orbital_mix_magnitude");
  aAttrib.put(cur);
  setOrbitalSetSize(norb);
  xmlNodePtr occ_ptr=NULL;
  xmlNodePtr coeff_ptr=NULL;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      occ_ptr=cur;
    }
    else if(cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var")
    {
      coeff_ptr=cur;
    }
    cur=cur->next;
  }
  if(coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << std::endl;
    return setIdentity(true);
  }
  bool success=putOccupation(occ_ptr);
  if(H5file==false)
    success = putFromXML(coeff_ptr);
  else
  {
    hdf_archive hin(myComm);

    if(myComm->rank()==0){
      if(!hin.open(MOhref2,H5F_ACC_RDONLY))
        APP_ABORT("SPOSetBase::putFromH5 missing or incorrect path to H5 file.");
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
      PBC=false;
      hin.read(PBC,"PBC");
      hin.close();
      if (PBC)
        APP_ABORT("SPOSetBase::putFromH5 PBC is not supported by AoS builds");
    }
    myComm->bcast(PBC);
    success = putFromH5(MOhref2, coeff_ptr);
  }

  bool success2 = transformSPOSet();
  if(debugc=="yes")
  {
    app_log() << "   Single-particle orbital coefficients dims="
      << C->rows() << " x " << C->cols() << std::endl;
    app_log() << C << std::endl;
  }

  init_LCOrbitalSetOpt(orbital_mix_magnitude);

  return success && success2;

 return true;
}

void SPOSetBase::checkObject()
{
  if(!(OrbitalSetSize == C->rows() && BasisSetSize == C->cols()))
  {
    app_error() << "   SPOSetBase::checkObject Linear coeffient for SPOSet is not consistent with the input." << std::endl;
    OHMMS::Controller->abort();
  }
}



bool SPOSetBase::putFromXML(xmlNodePtr coeff_ptr)
{
  Identity=true;
  int norbs=0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norbs,"size");
  aAttrib.add(norbs,"orbitals");
  aAttrib.put(coeff_ptr);
  if(norbs < OrbitalSetSize)
  {
    return false;
    APP_ABORT("SPOSetBase::putFromXML missing or incorrect size");
  }
  if(norbs)
  {
    Identity=false;
    std::vector<ValueType> Ctemp;
    Ctemp.resize(norbs*BasisSetSize);
    setIdentity(Identity);
    putContent(Ctemp,coeff_ptr);
    int n=0,i=0;
    std::vector<ValueType>::iterator cit(Ctemp.begin());
    while(i<OrbitalSetSize)
    {
      if(Occ[n]>std::numeric_limits<RealType>::epsilon())
      {
        std::copy(cit,cit+BasisSetSize,(*C)[i]);
        i++;
      }
      n++;
      cit+=BasisSetSize;
    }
  }
  return true;
}

/** read data from a hdf5 file
 * @param norb number of orbitals to be initialized
 * @param fname hdf5 file name
 * @param coeff_ptr xmlnode for coefficients
 */
bool SPOSetBase::putFromH5(const char* fname, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  int norbs=OrbitalSetSize;
  int neigs=BasisSetSize;
  int setVal=-1;
  std::string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal,"spindataset");
  aAttrib.add(neigs,"size");
  aAttrib.add(neigs,"orbitals");
  aAttrib.put(coeff_ptr);
  setIdentity(false);
  hdf_archive hin(myComm);
  if(myComm->rank()==0)
  {
    if(!hin.open(fname,H5F_ACC_RDONLY))
      APP_ABORT("SPOSetBase::putFromH5 missing or incorrect path to H5 file.");

    Matrix<RealType> Ctemp(neigs,BasisSetSize);
    char name[72];
    sprintf(name,"%s%d","/KPTS_0/eigenset_",setVal);
    setname=name;
    if(!hin.read(Ctemp,setname))
    {
       setname="SPOSetBase::putFromH5 Missing "+setname+" from HDF5 File.";
       APP_ABORT(setname.c_str());
    }
    hin.close();

    int n=0,i=0;
    while(i<norbs)
    {
      if(Occ[n]>0.0)
      {
        std::copy(Ctemp[n],Ctemp[n+1],(*C)[i]);
        i++;
      }
      n++;
    }
  }
  myComm->bcast(C->data(),C->size());
#else
  APP_ABORT("SPOSetBase::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


bool SPOSetBase::putOccupation(xmlNodePtr occ_ptr)
{
  //die??
  if(BasisSetSize ==0)
  {
    APP_ABORT("SPOSetBase::putOccupation detected ZERO BasisSetSize");
    return false;
  }
  Occ.resize(std::max(BasisSetSize,OrbitalSetSize));
  Occ=0.0;
  for(int i=0; i<OrbitalSetSize; i++)
    Occ[i]=1.0;
  std::vector<int> occ_in;
  std::string occ_mode("table");
  if(occ_ptr == NULL)
  {
    occ_mode="ground";
  }
  else
  {
    const xmlChar* o=xmlGetProp(occ_ptr,(const xmlChar*)"mode");
    if(o)
      occ_mode = (const char*)o;
  }
  //Do nothing if mode == ground
  if(occ_mode == "excited")
  {
    putContent(occ_in,occ_ptr);
    for(int k=0; k<occ_in.size(); k++)
    {
      if(occ_in[k]<0) //remove this, -1 is to adjust the base
        Occ[-occ_in[k]-1]=0.0;
      else
        Occ[occ_in[k]-1]=1.0;
    }
  }
  else
    if(occ_mode == "table")
    {
      putContent(Occ,occ_ptr);
    }
  return true;
}
#endif

void SPOSetBase::basic_report(const std::string& pad)
{
  app_log()<<pad<<"size = "<<size()<< std::endl;
  app_log()<<pad<<"state info:"<< std::endl;
  //states.report(pad+"  ");
  app_log().flush();
}


void SPOSetBase::evaluateBasis (const ParticleSet &P, int first, int last,
                                ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad,
                                ValueMatrix_t &basis_lapl)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateBasis.\n");
}

void SPOSetBase::evaluateForDeriv (const ParticleSet &P, int first, int last,
                                   ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad,
                                   ValueMatrix_t &basis_lapl)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateBasis.\n");
}

void SPOSetBase::copyParamsFromMatrix (const opt_variables_type& active,
                                       const ValueMatrix_t &mat, std::vector<RealType> &destVec)
{
  APP_ABORT("Need specialization of SPOSetBase::copyParamsFromMatrix.");
}

void SPOSetBase::evaluateGradSource (const ParticleSet &P
                                     , int first, int last, const ParticleSet &source
                                     , int iat_src, GradMatrix_t &gradphi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

void SPOSetBase::evaluateGradSource (const ParticleSet &P, int first, int last,
                                     const ParticleSet &source, int iat_src,
                                     GradMatrix_t &grad_phi,
                                     HessMatrix_t &grad_grad_phi,
                                     GradMatrix_t &grad_lapl_phi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

#ifdef QMC_CUDA

void SPOSetBase::evaluate(const ParticleSet& P, const PosType& r, std::vector<RealType> &psi)
{
  APP_ABORT("Not implemented.\n");
}


void SPOSetBase::evaluate (std::vector<Walker_t*> &walkers, int iat,
                           gpu::device_vector<CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSetBase.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSetBase::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &new_pos,
                           gpu::device_vector<CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSetBase.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSetBase::evaluate (std::vector<Walker_t*> &walkers,
                           std::vector<PosType> &new_pos,
                           gpu::device_vector<CudaValueType*> &phi,
                           gpu::device_vector<CudaValueType*> &grad_lapl_list,
                           int row_stride)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSetBase.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSetBase::evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSetBase.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSetBase::evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSetBase.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}
#endif
}


