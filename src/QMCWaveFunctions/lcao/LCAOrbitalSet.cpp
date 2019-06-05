//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include <Numerics/MatrixOperators.h>
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{
LCAOrbitalSet::LCAOrbitalSet(basis_type* bs)
    : myBasisSet(nullptr), C(nullptr), params_supplied(false), BasisSetSize(0), Identity(true), IsCloned(false)
{
  if (bs != nullptr)
    setBasisSet(bs);
}

LCAOrbitalSet::~LCAOrbitalSet()
{
  if (!IsCloned && C != nullptr)
    delete C;
}

void LCAOrbitalSet::setBasisSet(basis_type* bs)
{
  myBasisSet   = bs;
  BasisSetSize = myBasisSet->getBasisSetSize();
  Temp.resize(BasisSetSize);
  Temph.resize(BasisSetSize);
  Tempgh.resize(BasisSetSize);
}

bool LCAOrbitalSet::setIdentity(bool useIdentity)
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
    APP_ABORT("LCAOrbitalBuilder::setIdentiy ");
  }

  return true;
}

SPOSet* LCAOrbitalSet::makeClone() const
{
  LCAOrbitalSet* myclone = new LCAOrbitalSet(*this);
  myclone->myBasisSet    = myBasisSet->makeClone();
  myclone->IsCloned      = true;
  return myclone;
}

void LCAOrbitalSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  if (Identity)
  { //PAY ATTENTION TO COMPLEX
    myBasisSet->evaluateV(P, iat, psi.data());
  }
  else
  {
    Vector<ValueType> vTemp(Temp.data(0), BasisSetSize);
    myBasisSet->evaluateV(P, iat, vTemp.data());
    simd::gemv(*C, Temp.data(0), psi.data());
  }
}

/** Find a better place for other user classes, Matrix should be padded as well */
template<typename T, unsigned D>
inline void Product_ABt(const VectorSoaContainer<T, D>& A, const Matrix<T>& B, VectorSoaContainer<T, D>& C)
{
  constexpr char transa = 't';
  constexpr char transb = 'n';
  constexpr T zone(1);
  constexpr T zero(0);
  BLAS::gemm(transa, transb, B.rows(), D, B.cols(), zone, B.data(), B.cols(), A.data(), A.capacity(), zero, C.data(),
             C.capacity());
}

inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp,
                                             ValueVector_t& psi,
                                             GradVector_t& dpsi,
                                             ValueVector_t& d2psi) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, psi.data());
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dpsi[j][0] = gx[j];
    dpsi[j][1] = gy[j];
    dpsi[j][2] = gz[j];
  }
  std::copy_n(temp.data(4), OrbitalSetSize, d2psi.data());
}

inline void LCAOrbitalSet::evaluate_vgh_impl(const vgh_type& temp,
                                             ValueVector_t& psi,
                                             GradVector_t& dpsi,
                                             HessVector_t& d2psi) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, psi.data());
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  const ValueType* restrict hxx = temp.data(4);
  const ValueType* restrict hxy = temp.data(5);
  const ValueType* restrict hxz = temp.data(6);
  const ValueType* restrict hyy = temp.data(7);
  const ValueType* restrict hyz = temp.data(8);
  const ValueType* restrict hzz = temp.data(9);

  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dpsi[j][0] = gx[j];
    dpsi[j][1] = gy[j];
    dpsi[j][2] = gz[j];
    
    d2psi[j](0,0)                 = hxx[j]; 
    d2psi[j](0,1) = d2psi[j](1,0) = hxy[j]; 
    d2psi[j](0,2) = d2psi[j](2,0) = hxz[j]; 
    d2psi[j](1,1)                 = hyy[j]; 
    d2psi[j](2,1) = d2psi[j](1,2) = hyz[j]; 
    d2psi[j](2,2)                 = hzz[j]; 
  }
}

inline void LCAOrbitalSet::evaluate_vghgh_impl(const vghgh_type& temp,
                                             int i,
                                             ValueMatrix_t& psi,
                                             GradMatrix_t& dpsi,
                                             HessMatrix_t& d2psi,
                                             GGGMatrix_t& dghpsi) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, psi[i]);
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  const ValueType* restrict hxx = temp.data(4);
  const ValueType* restrict hxy = temp.data(5);
  const ValueType* restrict hxz = temp.data(6);
  const ValueType* restrict hyy = temp.data(7);
  const ValueType* restrict hyz = temp.data(8);
  const ValueType* restrict hzz = temp.data(9);
  const ValueType* restrict gh_xxx = temp.data(10);
  const ValueType* restrict gh_xxy = temp.data(11);
  const ValueType* restrict gh_xxz = temp.data(12);
  const ValueType* restrict gh_xyy = temp.data(13);
  const ValueType* restrict gh_xyz = temp.data(14);
  const ValueType* restrict gh_xzz = temp.data(15);
  const ValueType* restrict gh_yyy = temp.data(16);
  const ValueType* restrict gh_yyz = temp.data(17);
  const ValueType* restrict gh_yzz = temp.data(18);
  const ValueType* restrict gh_zzz = temp.data(19);

  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dpsi[i][j][0] = gx[j];
    dpsi[i][j][1] = gy[j];
    dpsi[i][j][2] = gz[j];
    
    d2psi[i][j](0,0)                 = hxx[j]; 
    d2psi[i][j](0,1) = d2psi[i][j](1,0) = hxy[j]; 
    d2psi[i][j](0,2) = d2psi[i][j](2,0) = hxz[j]; 
    d2psi[i][j](1,1)                 = hyy[j]; 
    d2psi[i][j](2,1) = d2psi[i][j](1,2) = hyz[j]; 
    d2psi[i][j](2,2)                 = hzz[j]; 
    
    dghpsi[i][j][0](0,0) = gh_xxx[j];  //x|xx  
    dghpsi[i][j][0](0,1) = gh_xxy[j];  //x|xy
    dghpsi[i][j][0](0,2) = gh_xxz[j];  //x|xz
    dghpsi[i][j][0](1,0) = gh_xxy[j];  //x|yx = xxy
    dghpsi[i][j][0](1,1) = gh_xyy[j];  //x|yy
    dghpsi[i][j][0](1,2) = gh_xyz[j];  //x|yz
    dghpsi[i][j][0](2,0) = gh_xxz[j];  //x|zx = xxz
    dghpsi[i][j][0](2,1) = gh_xyz[j];  //x|zy = xyz
    dghpsi[i][j][0](2,2) = gh_xzz[j];  //x|zz
    
    dghpsi[i][j][1](0,0) = gh_xxy[j];  //y|xx = xxy
    dghpsi[i][j][1](0,1) = gh_xyy[j];  //y|xy = xyy
    dghpsi[i][j][1](0,2) = gh_xyz[j];  //y|xz = xyz
    dghpsi[i][j][1](1,0) = gh_xyy[j];  //y|yx = xyy
    dghpsi[i][j][1](1,1) = gh_yyy[j];  //y|yy
    dghpsi[i][j][1](1,2) = gh_yyz[j];  //y|yz 
    dghpsi[i][j][1](2,0) = gh_xyz[j];  //y|zx = xyz
    dghpsi[i][j][1](2,1) = gh_yyz[j];  //y|zy = yyz
    dghpsi[i][j][1](2,2) = gh_yzz[j];  //y|zz 
   
    dghpsi[i][j][2](0,0) = gh_xxz[j];  //z|xx = xxz
    dghpsi[i][j][2](0,1) = gh_xyz[j];  //z|xy = xyz
    dghpsi[i][j][2](0,2) = gh_xzz[j];  //z|xz = xzz
    dghpsi[i][j][2](1,0) = gh_xyz[j];  //z|yx = xyz
    dghpsi[i][j][2](1,1) = gh_yyz[j];  //z|yy = yyz
    dghpsi[i][j][2](1,2) = gh_yzz[j];  //z|yz = yzz
    dghpsi[i][j][2](2,0) = gh_xzz[j];  //z|zx = xzz
    dghpsi[i][j][2](2,1) = gh_yzz[j];  //z|zy = yzz
    dghpsi[i][j][2](2,2) = gh_zzz[j];  //z|zz  
  }
}

inline void LCAOrbitalSet::evaluate_vghgh_impl(const vghgh_type& temp,
                                             ValueVector_t& psi,
                                             GradVector_t& dpsi,
                                             HessVector_t& d2psi,
                                             GGGVector_t& dghpsi) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, psi.data());
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  const ValueType* restrict hxx = temp.data(4);
  const ValueType* restrict hxy = temp.data(5);
  const ValueType* restrict hxz = temp.data(6);
  const ValueType* restrict hyy = temp.data(7);
  const ValueType* restrict hyz = temp.data(8);
  const ValueType* restrict hzz = temp.data(9);
  const ValueType* restrict gh_xxx = temp.data(10);
  const ValueType* restrict gh_xxy = temp.data(11);
  const ValueType* restrict gh_xxz = temp.data(12);
  const ValueType* restrict gh_xyy = temp.data(13);
  const ValueType* restrict gh_xyz = temp.data(14);
  const ValueType* restrict gh_xzz = temp.data(15);
  const ValueType* restrict gh_yyy = temp.data(16);
  const ValueType* restrict gh_yyz = temp.data(17);
  const ValueType* restrict gh_yzz = temp.data(18);
  const ValueType* restrict gh_zzz = temp.data(19);

  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dpsi[j][0] = gx[j];
    dpsi[j][1] = gy[j];
    dpsi[j][2] = gz[j];
    
    d2psi[j](0,0)                 = hxx[j]; 
    d2psi[j](0,1) = d2psi[j](1,0) = hxy[j]; 
    d2psi[j](0,2) = d2psi[j](2,0) = hxz[j]; 
    d2psi[j](1,1)                 = hyy[j]; 
    d2psi[j](2,1) = d2psi[j](1,2) = hyz[j]; 
    d2psi[j](2,2)                 = hzz[j]; 
    
    dghpsi[j][0](0,0) = gh_xxx[j];  //x|xx  
    dghpsi[j][0](0,1) = gh_xxy[j];  //x|xy
    dghpsi[j][0](0,2) = gh_xxz[j];  //x|xz
    dghpsi[j][0](1,0) = gh_xxy[j];  //x|yx = xxy
    dghpsi[j][0](1,1) = gh_xyy[j];  //x|yy
    dghpsi[j][0](1,2) = gh_xyz[j];  //x|yz
    dghpsi[j][0](2,0) = gh_xxz[j];  //x|zx = xxz
    dghpsi[j][0](2,1) = gh_xyz[j];  //x|zy = xyz
    dghpsi[j][0](2,2) = gh_xzz[j];  //x|zz
    
    dghpsi[j][1](0,0) = gh_xxy[j];  //y|xx = xxy
    dghpsi[j][1](0,1) = gh_xyy[j];  //y|xy = xyy
    dghpsi[j][1](0,2) = gh_xyz[j];  //y|xz = xyz
    dghpsi[j][1](1,0) = gh_xyy[j];  //y|yx = xyy
    dghpsi[j][1](1,1) = gh_yyy[j];  //y|yy
    dghpsi[j][1](1,2) = gh_yyz[j];  //y|yz 
    dghpsi[j][1](2,0) = gh_xyz[j];  //y|zx = xyz
    dghpsi[j][1](2,1) = gh_xyy[j];  //y|xy = xyy
    dghpsi[j][1](2,2) = gh_yzz[j];  //y|zz 
   
    dghpsi[j][2](0,0) = gh_xzz[j];  //z|xx = xzz
    dghpsi[j][2](0,1) = gh_xyz[j];  //z|xy = xyz
    dghpsi[j][2](0,2) = gh_xzz[j];  //z|xz = xzz
    dghpsi[j][2](1,0) = gh_xyz[j];  //z|yx = xyz
    dghpsi[j][2](1,1) = gh_yyz[j];  //z|yy = yyz
    dghpsi[j][2](1,2) = gh_yzz[j];  //z|yz = yzz
    dghpsi[j][2](2,0) = gh_xzz[j];  //z|zx = xzz
    dghpsi[j][2](2,1) = gh_yzz[j];  //z|zy = yzz
    dghpsi[j][2](2,2) = gh_zzz[j];  //z|zz
  }
}


void LCAOrbitalSet::evaluate(const ParticleSet& P,
                             int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             ValueVector_t& d2psi)
{
  //TAKE CARE OF IDENTITY
  myBasisSet->evaluateVGL(P, iat, Temp);
  if (Identity)
    evaluate_vgl_impl(Temp, psi, dpsi, d2psi);
  else
  {
    Product_ABt(Temp, *C, Tempv);
    evaluate_vgl_impl(Tempv, psi, dpsi, d2psi);
  }
}

void LCAOrbitalSet::evaluateDetRatios(const VirtualParticleSet& VP,
                                      ValueVector_t& psi,
                                      const ValueVector_t& psiinv,
                                      std::vector<ValueType>& ratios)
{
  Vector<ValueType> vTemp(Temp.data(0), BasisSetSize);
  Vector<ValueType> invTemp(Temp.data(1), BasisSetSize);

  MatrixOperators::product_Atx(*C, psiinv, invTemp.data());

  for (size_t j = 0; j < VP.getTotalNum(); j++)
  {
    myBasisSet->evaluateV(VP, j, vTemp.data());
    ratios[j] = simd::dot(vTemp.data(), invTemp.data(), BasisSetSize);
  }
}

void LCAOrbitalSet::evaluate(const ParticleSet& P,
                             int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             HessVector_t& dhpsi)
{
  //TAKE CARE OF IDENTITY
  myBasisSet->evaluateVGH(P, iat, Temph);
  if (Identity)
    evaluate_vgh_impl(Temph, psi, dpsi, dhpsi);
  else
  {
    Product_ABt(Temph, *C, Temphv);
    evaluate_vgh_impl(Temphv, psi, dpsi, dhpsi);
  }
}

void LCAOrbitalSet::evaluate(const ParticleSet& P,
                             int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             HessVector_t& dhpsi,
                             GGGVector_t& dghpsi)
{
// APP_ABORT("LCAORbitalSet::evaluate(psi,gpsi,hpsi,ghpsi) not implemented\n");
  
  //TAKE CARE OF IDENTITY
  myBasisSet->evaluateVGHGH(P, iat, Tempgh);
  if (Identity)
    evaluate_vghgh_impl(Tempgh, psi, dpsi, dhpsi, dghpsi);
  else
  {
    Product_ABt(Tempgh, *C, Tempghv);
    evaluate_vghgh_impl(Tempghv, psi, dpsi, dhpsi, dghpsi);
  }
}

/* implement using gemm algorithm */
inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp,
                                             int i,
                                             ValueMatrix_t& logdet,
                                             GradMatrix_t& dlogdet,
                                             ValueMatrix_t& d2logdet) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, logdet[i]);
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dlogdet[i][j][0] = gx[j];
    dlogdet[i][j][1] = gy[j];
    dlogdet[i][j][2] = gz[j];
  }
  std::copy_n(temp.data(4), OrbitalSetSize, d2logdet[i]);
}

inline void LCAOrbitalSet::evaluate_vgh_impl(const vgh_type& temp,
                                             int i,
                                             ValueMatrix_t& psi,
                                             GradMatrix_t& dpsi,
                                             HessMatrix_t& d2psi) const
{
  std::copy_n(temp.data(0), OrbitalSetSize, psi[i]);
  const ValueType* restrict gx = temp.data(1);
  const ValueType* restrict gy = temp.data(2);
  const ValueType* restrict gz = temp.data(3);
  const ValueType* restrict hxx = temp.data(4);
  const ValueType* restrict hxy = temp.data(5);
  const ValueType* restrict hxz = temp.data(6);
  const ValueType* restrict hyy = temp.data(7);
  const ValueType* restrict hyz = temp.data(8);
  const ValueType* restrict hzz = temp.data(9);

  for (size_t j = 0; j < OrbitalSetSize; j++)
  {
    dpsi[i][j][0] = gx[j];
    dpsi[i][j][1] = gy[j];
    dpsi[i][j][2] = gz[j];
    
    d2psi[i][j](0,0)                 = hxx[j]; 
    d2psi[i][j](0,1) = d2psi[i][j](1,0) = hxy[j]; 
    d2psi[i][j](0,2) = d2psi[i][j](2,0) = hxz[j]; 
    d2psi[i][j](1,1)                 = hyy[j]; 
    d2psi[i][j](2,1) = d2psi[i][j](1,2) = hyz[j]; 
    d2psi[i][j](2,2)                 = hzz[j]; 
  }
}

void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                         int first,
                                         int last,
                                         ValueMatrix_t& logdet,
                                         GradMatrix_t& dlogdet,
                                         ValueMatrix_t& d2logdet)
{
  if (Identity)
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGL(P, iat, Temp);
      evaluate_vgl_impl(Temp, i, logdet, dlogdet, d2logdet);
    }
  }
  else
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGL(P, iat, Temp);
      Product_ABt(Temp, *C, Tempv);
      evaluate_vgl_impl(Tempv, i, logdet, dlogdet, d2logdet);
    }
  }
}

void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                         int first,
                                         int last,
                                         ValueMatrix_t& logdet,
                                         GradMatrix_t& dlogdet,
                                         HessMatrix_t& grad_grad_logdet)
{
  if (Identity)
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGH(P, iat, Temph);
      evaluate_vgh_impl(Temph, i, logdet, dlogdet, grad_grad_logdet);
    }
  }
  else
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGH(P, iat, Temph);
      Product_ABt(Temph, *C, Temphv);
      evaluate_vgh_impl(Temphv, i, logdet, dlogdet, grad_grad_logdet);
    }
  }
}

void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                         int first,
                                         int last,
                                         ValueMatrix_t& logdet,
                                         GradMatrix_t& dlogdet,
                                         HessMatrix_t& grad_grad_logdet,
                                         GGGMatrix_t& grad_grad_grad_logdet)
{
  if (Identity)
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGHGH(P, iat, Tempgh);
      evaluate_vghgh_impl(Tempgh, i, logdet, dlogdet, grad_grad_logdet, grad_grad_grad_logdet);
    }
  }
  else
  {
    for (size_t i = 0, iat = first; iat < last; i++, iat++)
    {
      myBasisSet->evaluateVGHGH(P, iat, Tempgh);
      Product_ABt(Tempgh, *C, Tempghv);
      evaluate_vghgh_impl(Tempghv, i, logdet, dlogdet, grad_grad_logdet, grad_grad_grad_logdet);
    }
  }
}

void LCAOrbitalSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("LCAOrbitalSet::evaluateThirdDeriv(P,istart,istop,ggg_logdet) not implemented\n");
}

void LCAOrbitalSet::buildOptVariables(const std::vector<std::pair<int, int>>& rotations)
{
#if !defined(QMC_COMPLEX)
  const size_t nmo = OrbitalSetSize;
  const size_t nb  = BasisSetSize;
  m_init_B         = *C;

  // create active rotations
  m_act_rot_inds = rotations;

  // This will add the orbital rotation parameters to myVars
  // and will also read in initial parameter values supplied in input file
  int p, q;
  int nparams_active = m_act_rot_inds.size();
  if (params_supplied)
    if (nparams_active != params.size())
      APP_ABORT("The number of supplied orbital rotation parameters does not match number prdouced by the slater "
                "expansion. \n");

  for (int i = 0; i < nparams_active; i++)
  {
    p = m_act_rot_inds[i].first;
    q = m_act_rot_inds[i].second;
    std::stringstream sstr;
    sstr << objectName << "_orb_rot_" << (p < 10 ? "0" : "") << (p < 100 ? "0" : "") << (p < 1000 ? "0" : "") << p
         << "_" << (q < 10 ? "0" : "") << (q < 100 ? "0" : "") << (q < 1000 ? "0" : "") << q;

    // If the user input parameteres, use those. Otherwise, initialize the parameters to zero
    if (params_supplied)
    {
      myVars.insert(sstr.str(), params[i]);
    }
    else
    {
      myVars.insert(sstr.str(), 0.0);
    }
  }

  //Printing the parameters
  if (false)
  {
    app_log() << std::string(16, ' ') << "Parameter name" << std::string(15, ' ') << "Value\n";
    myVars.print(app_log());
  }

  std::vector<RealType> param(m_act_rot_inds.size());
  for (int i = 0; i < m_act_rot_inds.size(); i++)
  {
    param[i] = myVars[i];
  }
  apply_rotation(param);

  if (params_supplied)
  {
    app_log() << "PRINTING MO COEFFICIENTS CREATED BY SUPPLIED PARAMS\n";
    for (int j = 0; j < nmo; j++)
    {
      for (int i = 0; i < nb; i++)
      {
        app_log() << " " << std::right << std::fixed << std::setprecision(16) << std::setw(23) << std::scientific
                  << *(C->data() + j * nb + i);

        if ((j * nb + i + 1) % 4 == 0)
        {
          app_log() << std::endl;
        }
      }
    }
    app_log() << "\n done printing molecular orbital coefficients\n";
    app_log() << std::endl;
  }

  if (!Optimizable)
  {
    //THIS ALLOWS FOR ORBITAL PARAMETERS TO BE READ IN EVEN WHEN THOSE PARAMETERS ARE NOT BEING OPTIMIZED
    //this assumes there are only CI coefficients ahead of the M_orb_coefficients
    myVars.Index.erase(myVars.Index.begin(), myVars.Index.end());
    myVars.NameAndValue.erase(myVars.NameAndValue.begin(), myVars.NameAndValue.end());
    myVars.ParameterType.erase(myVars.ParameterType.begin(), myVars.ParameterType.end());
    myVars.Recompute.erase(myVars.Recompute.begin(), myVars.Recompute.end());
  }
#endif
}

void LCAOrbitalSet::evaluateDerivatives(ParticleSet& P,
                                        const opt_variables_type& optvars,
                                        std::vector<ValueType>& dlogpsi,
                                        std::vector<ValueType>& dhpsioverpsi,
                                        const ValueType& psiCurrent,
                                        const std::vector<ValueType>& Coeff,
                                        const std::vector<size_t>& C2node_up,
                                        const std::vector<size_t>& C2node_dn,
                                        const ValueVector_t& detValues_up,
                                        const ValueVector_t& detValues_dn,
                                        const GradMatrix_t& grads_up,
                                        const GradMatrix_t& grads_dn,
                                        const ValueMatrix_t& lapls_up,
                                        const ValueMatrix_t& lapls_dn,
                                        const ValueMatrix_t& M_up,
                                        const ValueMatrix_t& M_dn,
                                        const ValueMatrix_t& Minv_up,
                                        const ValueMatrix_t& Minv_dn,
                                        const GradMatrix_t& B_grad,
                                        const ValueMatrix_t& B_lapl,
                                        const std::vector<int>& detData_up,
                                        const size_t N1,
                                        const size_t N2,
                                        const size_t NP1,
                                        const size_t NP2,
                                        const std::vector<std::vector<int>>& lookup_tbl)
{
#if !defined(QMC_COMPLEX)
  bool recalculate(false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (optvars.recompute(kk))
      recalculate = true;
  }
  if (recalculate)
  {
    ParticleSet::ParticleGradient_t myG_temp, myG_J;
    ParticleSet::ParticleLaplacian_t myL_temp, myL_J;
    const int NP = P.getTotalNum();
    myG_temp.resize(NP);
    myG_temp = 0.0;
    myL_temp.resize(NP);
    myL_temp = 0.0;
    myG_J.resize(NP);
    myG_J = 0.0;
    myL_J.resize(NP);
    myL_J            = 0.0;
    const size_t nmo = OrbitalSetSize;
    const size_t nb  = BasisSetSize;
    const size_t nel = P.last(0) - P.first(0);

    const RealType* restrict C_p = Coeff.data();
    for (int i = 0; i < Coeff.size(); i++)
    {
      const size_t upC     = C2node_up[i];
      const size_t dnC     = C2node_dn[i];
      const ValueType tmp1 = C_p[i] * detValues_dn[dnC];
      const ValueType tmp2 = C_p[i] * detValues_up[upC];
      for (size_t k = 0, j = N1; k < NP1; k++, j++)
      {
        myG_temp[j] += tmp1 * grads_up(upC, k);
        myL_temp[j] += tmp1 * lapls_up(upC, k);
      }
      for (size_t k = 0, j = N2; k < NP2; k++, j++)
      {
        myG_temp[j] += tmp2 * grads_dn(dnC, k);
        myL_temp[j] += tmp2 * lapls_dn(dnC, k);
      }
    }

    myG_temp *= (1 / psiCurrent);
    myL_temp *= (1 / psiCurrent);

    // calculation of myG_J which will be used to represent \frac{\nabla\psi_{J}}{\psi_{J}}
    // calculation of myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}}
    // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold
    for (int iat = 0; iat < (myL_temp.size()); iat++)
    {
      myG_J[iat] = (P.G[iat] - myG_temp[iat]);
      myL_J[iat] = (P.L[iat] + dot(P.G[iat], P.G[iat]) - myL_temp[iat]);
    }


    table_method_eval(dlogpsi, dhpsioverpsi, myL_J, myG_J, nel, nmo, psiCurrent, Coeff, C2node_up, C2node_dn,
                      detValues_up, detValues_dn, grads_up, grads_dn, lapls_up, lapls_dn, M_up, M_dn, Minv_up, Minv_dn,
                      B_grad, B_lapl, detData_up, N1, N2, NP1, NP2, lookup_tbl);
  }
#endif
}

#if !defined(QMC_COMPLEX)
void LCAOrbitalSet::apply_rotation(const std::vector<RealType>& param)
{
  assert(param.size() == m_act_rot_inds.size());

  const size_t nmo = OrbitalSetSize;
  const size_t nb  = BasisSetSize;
  ValueMatrix_t rot_mat(nmo, nmo);
  rot_mat = ValueType(0);

  // read out the parameters that define the rotation into an antisymmetric matrix
  for (int i = 0; i < m_act_rot_inds.size(); i++)
  {
    const int p      = m_act_rot_inds[i].first;
    const int q      = m_act_rot_inds[i].second;
    const RealType x = param[i];

    rot_mat[q][p] = x;
    rot_mat[p][q] = -x;
  }

  exponentiate_antisym_matrix(rot_mat);

  BLAS::gemm('N', 'T', nb, nmo, nmo, RealType(1.0), m_init_B.data(), nb, rot_mat.data(), nmo, RealType(0.0), C->data(),
             nb);
}


// compute exponential of a real, antisymmetric matrix by diagonalizing and exponentiating eigenvalues
void LCAOrbitalSet::exponentiate_antisym_matrix(ValueMatrix_t& mat)
{
  const int n = mat.rows();
  std::vector<std::complex<RealType>> mat_h(n * n, 0);
  std::vector<RealType> eval(n, 0);
  std::vector<std::complex<RealType>> work(2 * n, 0);
  std::vector<RealType> rwork(3 * n, 0);
  std::vector<std::complex<RealType>> mat_d(n * n, 0);
  std::vector<std::complex<RealType>> mat_t(n * n, 0);
  // exponentiating e^X = e^iY (Y hermitian)
  // i(-iX) = X, so -iX is hermitian
  // diagonalize -iX = UDU^T, exponentiate e^iD, and return U e^iD U^T
  // construct hermitian analogue of mat by multiplying by -i
  for (int i = 0; i < n; ++i)
  {
    for (int j = i; j < n; ++j)
    {
      mat_h[i + n * j] = std::complex<RealType>(0, -1.0 * mat[j][i]);
      mat_h[j + n * i] = std::complex<RealType>(0, 1.0 * mat[j][i]);
    }
  }
  // diagonalize the matrix
  char JOBZ('V');
  char UPLO('U');
  int N(n);
  int LDA(n);
  int LWORK(2 * n);
  int info = 0;
  LAPACK::heev(JOBZ, UPLO, N, &mat_h.at(0), LDA, &eval.at(0), &work.at(0), LWORK, &rwork.at(0), info);
  if (info != 0)
  {
    std::ostringstream msg;
    msg << "heev failed with info = " << info << " in MultiSlaterDeterminantFast::exponentiate_antisym_matrix";
    app_log() << msg.str() << std::endl;
    APP_ABORT(msg.str());
  }
  // ierate through diagonal matrix, exponentiate terms
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      mat_d[i + j * n] = (i == j) ? std::exp(std::complex<RealType>(0.0, eval[i])) : std::complex<RealType>(0.0, 0.0);
    }
  }
  // perform matrix multiplication
  // assume row major
  BLAS::gemm('N', 'C', n, n, n, std::complex<RealType>(1.0, 0), &mat_d.at(0), n, &mat_h.at(0), n,
             std::complex<RealType>(0.0, 0.0), &mat_t.at(0), n);
  BLAS::gemm('N', 'N', n, n, n, std::complex<RealType>(1.0, 0), &mat_h.at(0), n, &mat_t.at(0), n,
             std::complex<RealType>(0.0, 0.0), &mat_d.at(0), n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    {
      if (mat_d[i + n * j].imag() > 1e-12)
      {
        app_log() << "warning: large imaginary value in orbital rotation matrix: (i,j) = (" << i << "," << j
                  << "), im = " << mat_d[i + n * j].imag() << std::endl;
      }
      mat[j][i] = mat_d[i + n * j].real();
    }
}

void LCAOrbitalSet::table_method_eval(std::vector<ValueType>& dlogpsi,
                                      std::vector<ValueType>& dhpsioverpsi,
                                      const ParticleSet::ParticleLaplacian_t& myL_J,
                                      const ParticleSet::ParticleGradient_t& myG_J,
                                      const size_t nel,
                                      const size_t nmo,
                                      const ValueType& psiCurrent,
                                      const std::vector<RealType>& Coeff,
                                      const std::vector<size_t>& C2node_up,
                                      const std::vector<size_t>& C2node_dn,
                                      const ValueVector_t& detValues_up,
                                      const ValueVector_t& detValues_dn,
                                      const GradMatrix_t& grads_up,
                                      const GradMatrix_t& grads_dn,
                                      const ValueMatrix_t& lapls_up,
                                      const ValueMatrix_t& lapls_dn,
                                      const ValueMatrix_t& M_up,
                                      const ValueMatrix_t& M_dn,
                                      const ValueMatrix_t& Minv_up,
                                      const ValueMatrix_t& Minv_dn,
                                      const GradMatrix_t& B_grad,
                                      const ValueMatrix_t& B_lapl,
                                      const std::vector<int>& detData_up,
                                      const size_t N1,
                                      const size_t N2,
                                      const size_t NP1,
                                      const size_t NP2,
                                      const std::vector<std::vector<int>>& lookup_tbl)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUIDE TO THE MATICES BEING BUILT
----------------------------------------------
The idea here is that there is a loop over all unique determinants. For each determiant the table method is employed to calculate the contributions to the parameter derivatives (dhpsioverpsi/dlogpsi)

  loop through unquie determinants 
    loop through parameters
      evaluate contributaion to dlogpsi and dhpsioverpsi
\noindent 

  BLAS GUIDE  for matrix multiplication of  [  alpha * A.B + beta * C = C ]
  Matrix A is of dimensions a1,a2 and Matrix B is b1,b2   in which a2=b1
  The BLAS command is as follows...

 BLAS::gemm('N','N', b2, a1, a2 ,alpha, B, b2, A, a2, beta, C, b2);

Below is a human readable format for the matrix multiplications performed below...

This notation is inspired by http://dx.doi.org/10.1063/1.4948778
\newline
\hfill\break
$
    A_{i,j}=\phi_j(r_{i}) \\
    T = A^{-1} \widetilde{A} \\
    B_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i}) \\
    \hat{O_{I}} = \hat{O}D_{I} \\
    D_{I}=det(A_{I}) \newline 
    \psi_{MS} = \sum_{I=0} C_{I} D_{I\uparrow}D_{I\downarrow} \\
    \Psi_{total} = \psi_{J}\psi_{MS} \\
    \alpha_{I} = P^{T}_{I}TQ_{I} \\
    M_{I} = P^{T}_{I} \widetilde{M} Q_{I} = P^{T}_{I} (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )Q_{I} \\
$
\newline
There are three constants I use in the expressions for dhpsioverpsi and dlogpsi
\newline
\hfill\break
$
  const0 = C_{0}*det(A_{0\downarrow})+\sum_{I=1} C_{I}*det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{I=1} C_{I}*\hat{O}det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const2 = \sum_{I=1} C_{I}*det(A_{I\downarrow})* Tr[\alpha_{I}^{-1}M_{I}]*det(\alpha_{I}) \\
$
\newline
Below is a translation of the shorthand I use to represent matrices independent of ``excitation matrix".
\newline
\hfill\break
$
    Y_{1} =  A^{-1}B   \\
    Y_{2} = A^{-1}BA^{-1}\widetilde{A} \\
    Y_{3} = A^{-1}\widetilde{B} \\
    Y_{4} = \widetilde{M} = (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )\\
$
\newline
Below is a translation of the shorthand I use to represent matrices dependent on ``excitation" with respect to the reference Matrix and sums of matrices. Above this line I have represented these excitation matrices with a subscript ``I" but from this point on The subscript will be omitted and it is clear that whenever a matrix depends on $P^{T}_I$ and $Q_{I}$ that this is an excitation matrix. The reference matrix is always $A_{0}$ and is always the Hartree Fock Matrix.
\newline
\hfill\break
$
    Y_{5} = TQ \\
    Y_{6} = (P^{T}TQ)^{-1} = \alpha_{I}^{-1}\\
    Y_{7} = \alpha_{I}^{-1} P^{T} \\
    Y_{11} = \widetilde{M}Q \\
    Y_{23} = P^{T}\widetilde{M}Q \\
    Y_{24} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q \\
    Y_{25} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1} \\
    Y_{26} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1}P^{T}\\
$
\newline
So far you will notice that I have not included up or down arrows to specify what spin the matrices are of. This is because we are calculating the derivative of all up or all down spin orbital rotation parameters at a time. If we are finding the up spin derivatives then any term that is down spin will be constant. The following assumes that we are taking up-spin MO rotation parameter derivatives. Of course the down spin expression can be retrieved by swapping the up and down arrows. I have dubbed any expression with lowercase p prefix as a "precursor" to an expression actually used...
\newline
\hfill\break
$
    \dot{C_{I}} = C_{I}*det(A_{I\downarrow})\\
    \ddot{C_{I}} = C_{I}*\hat{O}det(A_{I\downarrow}) \\
    pK1 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}) \\
    pK2 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK3 = \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK4 = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK5 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T}) \\
$
\newline
Now these p matrices will be used to make various expressions via BLAS commands.
\newline
\hfill\break
$
    K1T = const0^{-1}*pK1.T =const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}T) \\
    TK1T = T.K1T = const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    TK2AiB = T.K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    K2XA =  const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}X\widetilde{A})\\
    TK2XA = T.K2XA = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}X\widetilde{A})\\ \\
    K2T = \frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK2T = T.K2T =\frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\
    MK2T = \frac{const0}{const1} Y_{4}.K2T= const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (\widetilde{M}Q\alpha_{I}^{-1}P^{T}T)\\ \\
    K3T = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK3T = T.K3T  = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK4T = T.K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\ \\
    K5T =  const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
    TK5T = T.K5T  = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (T Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
$
\newline
Now with all these matrices and constants the expressions of dhpsioverpsi and dlogpsi can be created.




In addition I will be using a special generalization of the kinetic operator which I will denote as O. Our Slater matrix with the special O operator applied to each element will be called B_bar

$
``Bbar"_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i})
$
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
  ValueMatrix_t Table;
  ValueMatrix_t Bbar;
  ValueMatrix_t Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y11, Y23, Y24, Y25, Y26;
  ValueMatrix_t pK1, K1T, TK1T, pK2, K2AiB, TK2AiB, K2XA, TK2XA, K2T, TK2T, MK2T, pK3, K3T, TK3T, pK4, K4T, TK4T, pK5,
      K5T, TK5T;

  Table.resize(nel, nmo);

  Bbar.resize(nel, nmo);

  Y1.resize(nel, nel);
  Y2.resize(nel, nmo);
  Y3.resize(nel, nmo);
  Y4.resize(nel, nmo);

  pK1.resize(nmo, nel);
  K1T.resize(nmo, nmo);
  TK1T.resize(nel, nmo);

  pK2.resize(nmo, nel);
  K2AiB.resize(nmo, nmo);
  TK2AiB.resize(nel, nmo);
  K2XA.resize(nmo, nmo);
  TK2XA.resize(nel, nmo);
  K2T.resize(nmo, nmo);
  TK2T.resize(nel, nmo);
  MK2T.resize(nel, nmo);

  pK3.resize(nmo, nel);
  K3T.resize(nmo, nmo);
  TK3T.resize(nel, nmo);

  pK4.resize(nmo, nel);
  K4T.resize(nmo, nmo);
  TK4T.resize(nel, nmo);

  pK5.resize(nmo, nel);
  K5T.resize(nmo, nmo);
  TK5T.resize(nel, nmo);

  const int parameters_size(m_act_rot_inds.size());
  const int parameter_start_index(0);

  const size_t num_unique_up_dets(detValues_up.size());
  const size_t num_unique_dn_dets(detValues_dn.size());

  const RealType* restrict cptr = Coeff.data();
  const size_t nc               = Coeff.size();
  const size_t* restrict upC(C2node_up.data());
  const size_t* restrict dnC(C2node_dn.data());
  //B_grad holds the gardient operator
  //B_lapl holds the laplacian operator
  //B_bar will hold our special O operator

  const int offset1(N1);
  const int offset2(N2);
  const int NPother(NP2);

  RealType* T(Table.data());

  //possibly replace wit BLAS calls
  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nmo; j++)
      Bbar(i, j) = B_lapl(i, j) + 2 * dot(myG_J[i + offset1], B_grad(i, j)) + myL_J[i + offset1] * M_up(i, j);

  const RealType* restrict B(Bbar.data());
  const RealType* restrict A(M_up.data());
  const RealType* restrict Ainv(Minv_up.data());
  //IMPORTANT NOTE: THE Dets[0]->psiMinv OBJECT DOES NOT HOLD THE INVERSE IF THE MULTIDIRACDETERMINANTBASE ONLY CONTAINES ONE ELECTRON. NEED A FIX FOR THIS CASE
  // The T matrix should be calculated and stored for use
  // T = A^{-1} \widetilde A
  //REMINDER: that the ValueMatrix_t "matrix" stores data in a row major order and that BLAS commands assume column major
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), A, nmo, Ainv, nel, RealType(0.0), T, nmo);

  BLAS::gemm('N', 'N', nel, nel, nel, RealType(1.0), B, nmo, Ainv, nel, RealType(0.0), Y1.data(), nel);
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), T, nmo, Y1.data(), nel, RealType(0.0), Y2.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), B, nmo, Ainv, nel, RealType(0.0), Y3.data(), nmo);

  //possibly replace with BLAS call
  Y4 = Y3 - Y2;

  //Need to create the constants: (Oi, const0, const1, const2)to take advantage of minimal BLAS commands;
  //Oi is the special operator applied to the slater matrix "A subscript i" from the total CI expansion
  //\hat{O_{i}} = \hat{O}D_{i} with D_{i}=det(A_{i}) and Multi-Slater component defined as \sum_{i=0} C_{i} D_{i\uparrow}D_{i\downarrow}
  std::vector<RealType> Oi(num_unique_dn_dets);

  for (int index = 0; index < num_unique_dn_dets; index++)
    for (int iat = 0; iat < NPother; iat++)
      Oi[index] += lapls_dn(index, iat) + 2 * dot(grads_dn(index, iat), myG_J[offset2 + iat]) +
          myL_J[offset2 + iat] * detValues_dn[index];

  //const0 = C_{0}*det(A_{0\downarrow})+\sum_{i=1} C_{i}*det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{i=1} C_{i}*\hat{O}det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const2 = \sum_{i=1} C_{i}*det(A_{i\downarrow})* Tr[\alpha_{i}^{-1}M_{i}]*det(\alpha_{i})
  RealType const0(0.0), const1(0.0), const2(0.0);
  for (size_t i = 0; i < nc; ++i)
  {
    const RealType c  = cptr[i];
    const size_t up   = upC[i];
    const size_t down = dnC[i];

    const0 += c * detValues_dn[down] * (detValues_up[up] / detValues_up[0]);
    const1 += c * Oi[down] * (detValues_up[up] / detValues_up[0]);
  }

  std::fill(pK1.begin(), pK1.end(), 0.0);
  std::fill(pK2.begin(), pK2.end(), 0.0);
  std::fill(pK3.begin(), pK3.end(), 0.0);
  std::fill(pK4.begin(), pK4.end(), 0.0);
  std::fill(pK5.begin(), pK5.end(), 0.0);

  //Now we are going to loop through all unique determinants.
  //The few lines above are for the reference matrix contribution.
  //Although I start the loop below from index 0, the loop only performs actions when the index is >= 1
  //the detData object contains all the information about the P^T and Q matrices (projection matrices) needed in the table method
  const int* restrict data_it = detData_up.data();
  for (int index = 0, datum = 0; index < num_unique_up_dets; index++)
  {
    const int k = data_it[datum];

    if (k == 0)
    {
      datum += 3 * k + 1;
    }

    else
    {
      //Number of rows and cols of P^T
      const int prows = k;
      const int pcols = nel;
      //Number of rows and cols of Q
      const int qrows = nmo;
      const int qcols = k;

      Y5.resize(nel, k);
      Y6.resize(k, k);

      //Any matrix multiplication of P^T or Q is simply a projection
      //Explicit matrix multiplication can be avoided; instead column or row copying can be done
      //BlAS::copy(size of col/row being copied,
      //           Matrix pointer + place to begin copying,
      //           storage spacing (number of elements btw next row/col element),
      //           Pointer to resultant matrix + place to begin pasting,
      //           storage spacing of resultant matrix)
      //For example the next 4 lines is the matrix multiplication of T*Q = Y5
      std::fill(Y5.begin(), Y5.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(nel, T + data_it[datum + 1 + k + i], nmo, Y5.data() + i, k);
      }

      std::fill(Y6.begin(), Y6.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y5.data() + (data_it[datum + 1 + i]) * k, 1, (Y6.data() + i * k), 1);
      }


      Vector<ValueType> WS;
      Vector<IndexType> Piv;
      WS.resize(k);
      Piv.resize(k);
      RealType PhaseR = 0.0;
      InvertWithLog(Y6.data(), k, k, WS.data(), Piv.data(), PhaseR);

      Y11.resize(nel, k);
      Y23.resize(k, k);
      Y24.resize(k, k);
      Y25.resize(k, k);
      Y26.resize(k, nel);

      std::fill(Y11.begin(), Y11.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(nel, Y4.data() + (data_it[datum + 1 + k + i]), nmo, Y11.data() + i, k);
      }

      std::fill(Y23.begin(), Y23.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y11.data() + (data_it[datum + 1 + i]) * k, 1, (Y23.data() + i * k), 1);
      }

      BLAS::gemm('N', 'N', k, k, k, RealType(1.0), Y23.data(), k, Y6.data(), k, RealType(0.0), Y24.data(), k);
      BLAS::gemm('N', 'N', k, k, k, RealType(1.0), Y6.data(), k, Y24.data(), k, RealType(0.0), Y25.data(), k);


      Y26.resize(k, nel);

      std::fill(Y26.begin(), Y26.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y25.data() + i, k, Y26.data() + (data_it[datum + 1 + i]), nel);
      }


      Y7.resize(k, nel);

      std::fill(Y7.begin(), Y7.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y6.data() + i, k, Y7.data() + (data_it[datum + 1 + i]), nel);
      }

      // c_Tr_AlphaI_MI is a constant contributing to constant const2
      // c_Tr_AlphaI_MI = Tr[\alpha_{I}^{-1}(P^{T}\widetilde{M} Q)]
      RealType c_Tr_AlphaI_MI = 0.0;
      for (int i = 0; i < k; i++)
      {
        c_Tr_AlphaI_MI += Y24(i, i);
      }

      for (int p = 0; p < lookup_tbl[index].size(); p++)
      {
        //el_p is the element position that contains information about the CI coefficient, and det up/dn values associated with the current unique determinant
        const int el_p(lookup_tbl[index][p]);
        const RealType c  = cptr[el_p];
        const size_t up   = upC[el_p];
        const size_t down = dnC[el_p];

        const RealType alpha_1(c * detValues_dn[down] * detValues_up[up] / detValues_up[0] * c_Tr_AlphaI_MI);
        const RealType alpha_2(c * detValues_dn[down] * detValues_up[up] / detValues_up[0]);
        const RealType alpha_3(c * Oi[down] * detValues_up[up] / detValues_up[0]);
        const RealType alpha_4(c * detValues_dn[down] * detValues_up[up] * (1 / psiCurrent));

        const2 += alpha_1;

        for (int i = 0; i < k; i++)
        {
          BLAS::axpy(nel, alpha_1, Y7.data() + i * nel, 1, pK1.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_2, Y7.data() + i * nel, 1, pK2.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_3, Y7.data() + i * nel, 1, pK3.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_4, Y7.data() + i * nel, 1, pK4.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_2, Y26.data() + i * nel, 1, pK5.data() + (data_it[datum + 1 + k + i]) * nel, 1);
        }
      }
      datum += 3 * k + 1;
    }
  }


  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK1.data(), nel, RealType(0.0), K1T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K1T.data(), nmo, T, nmo, RealType(0.0), TK1T.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, Y3.data(), nmo, pK2.data(), nel, RealType(0.0), K2AiB.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2AiB.data(), nmo, T, nmo, RealType(0.0), TK2AiB.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, Y2.data(), nmo, pK2.data(), nel, RealType(0.0), K2XA.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2XA.data(), nmo, T, nmo, RealType(0.0), TK2XA.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, const1 / (const0 * const0), T, nmo, pK2.data(), nel, RealType(0.0), K2T.data(),
             nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2T.data(), nmo, T, nmo, RealType(0.0), TK2T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, const0 / const1, K2T.data(), nmo, Y4.data(), nmo, RealType(0.0), MK2T.data(),
             nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK3.data(), nel, RealType(0.0), K3T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K3T.data(), nmo, T, nmo, RealType(0.0), TK3T.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, RealType(1.0), T, nmo, pK4.data(), nel, RealType(0.0), K4T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K4T.data(), nmo, T, nmo, RealType(0.0), TK4T.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK5.data(), nel, RealType(0.0), K5T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K5T.data(), nmo, T, nmo, RealType(0.0), TK5T.data(), nmo);


  for (int mu = 0, k = parameter_start_index; k < (parameter_start_index + parameters_size); k++, mu++)
  {
    int kk = myVars.where(k);
    const int i(m_act_rot_inds[mu].first), j(m_act_rot_inds[mu].second);
    if (i <= nel - 1 && j > nel - 1)
    {
      dlogpsi[kk] += ValueType( detValues_up[0] * (Table(i, j)) * const0 * (1 / psiCurrent) + (K4T(i, j) - K4T(j, i) - TK4T(i, j)) );

      dhpsioverpsi[kk] += ValueType( -0.5 * Y4(i, j) -
          0.5 *
              (-K5T(i, j) + K5T(j, i) + TK5T(i, j) + K2AiB(i, j) - K2AiB(j, i) - TK2AiB(i, j) - K2XA(i, j) +
               K2XA(j, i) + TK2XA(i, j) - MK2T(i, j) + K1T(i, j) - K1T(j, i) - TK1T(i, j) -
               const2 / const1 * K2T(i, j) + const2 / const1 * K2T(j, i) + const2 / const1 * TK2T(i, j) + K3T(i, j) -
               K3T(j, i) - TK3T(i, j) - K2T(i, j) + K2T(j, i) + TK2T(i, j)) );
    }
    else if (i <= nel - 1 && j <= nel - 1)
    {
      dlogpsi[kk] += ValueType( detValues_up[0] * (Table(i, j) - Table(j, i)) * const0 * (1 / psiCurrent) +
          (K4T(i, j) - TK4T(i, j) - K4T(j, i) + TK4T(j, i)) );

      dhpsioverpsi[kk] += ValueType( -0.5 * (Y4(i, j) - Y4(j, i)) -
          0.5 *
              (-K5T(i, j) + K5T(j, i) + TK5T(i, j) - TK5T(j, i) + K2AiB(i, j) - K2AiB(j, i) - TK2AiB(i, j) +
               TK2AiB(j, i) - K2XA(i, j) + K2XA(j, i) + TK2XA(i, j) - TK2XA(j, i) - MK2T(i, j) + MK2T(j, i) +
               K1T(i, j) - K1T(j, i) - TK1T(i, j) + TK1T(j, i) - const2 / const1 * K2T(i, j) +
               const2 / const1 * K2T(j, i) + const2 / const1 * TK2T(i, j) - const2 / const1 * TK2T(j, i) + K3T(i, j) -
               K3T(j, i) - TK3T(i, j) + TK3T(j, i) - K2T(i, j) + K2T(j, i) + TK2T(i, j) - TK2T(j, i)) );
    }
    else
    {
      dlogpsi[kk] += ValueType( (K4T(i, j) - K4T(j, i)) );

      dhpsioverpsi[kk] += ValueType( -0.5 *
          (-K5T(i, j) + K5T(j, i) + K2AiB(i, j) - K2AiB(j, i) - K2XA(i, j) + K2XA(j, i)

           + K1T(i, j) - K1T(j, i) - const2 / const1 * K2T(i, j) + const2 / const1 * K2T(j, i) + K3T(i, j) - K3T(j, i) -
           K2T(i, j) + K2T(j, i)) );
    }
  }
}
#endif // !defined(QMC_COMPLEX)


} // namespace qmcplusplus
