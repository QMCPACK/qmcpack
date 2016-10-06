//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef ATOMIC_ORBITAL_H
#define ATOMIC_ORBITAL_H

#include <config/stdlib/math.h>
#include <einspline/multi_bspline.h>
#include <QMCWaveFunctions/SPOSetBase.h>
#include <Lattice/CrystalLattice.h>
#include <Configuration.h>
#include <Utilities/NewTimer.h>


namespace qmcplusplus
{

/******************************************************************
// This is just a template trick to avoid template specialization //
// in AtomicOrbital.                                              //
******************************************************************/

template<typename StorageType>  struct AtomicOrbitalTraits {};
template<> struct AtomicOrbitalTraits<double>
{
  typedef multi_UBspline_1d_d SplineType;
};
template<> struct AtomicOrbitalTraits<std::complex<double> >
{
  typedef multi_UBspline_1d_z SplineType;
};

inline void EinsplineMultiEval (multi_UBspline_1d_d *spline,
                                double x, double *val)
{
  eval_multi_UBspline_1d_d (spline, x, val);
}
inline void EinsplineMultiEval (multi_UBspline_1d_z *spline,
                                double x, std::complex<double> *val)
{
  eval_multi_UBspline_1d_z (spline, x, val);
}
inline void EinsplineMultiEval (multi_UBspline_1d_d *spline, double x,
                                double *val, double *grad, double *lapl)
{
  eval_multi_UBspline_1d_d_vgl (spline, x, val, grad, lapl);
}
inline void EinsplineMultiEval (multi_UBspline_1d_z *spline, double x,
                                std::complex<double> *val, std::complex<double> *grad,
                                std::complex<double> *lapl)
{
  eval_multi_UBspline_1d_z_vgl (spline, x, val, grad, lapl);
}


template<typename StorageType>
class AtomicOrbital
{
public:
  typedef QMCTraits::PosType                    PosType;
  typedef QMCTraits::RealType                   RealType;
  typedef CrystalLattice<RealType,OHMMS_DIM>    UnitCellType;
  typedef Vector<double>                        RealValueVector_t;
  typedef Vector<TinyVector<double,OHMMS_DIM> > RealGradVector_t;
  typedef Vector<std::complex<double> >              ComplexValueVector_t;
  typedef Vector<TinyVector<std::complex<double>,OHMMS_DIM> > ComplexGradVector_t;
  typedef Vector<Tensor<double,OHMMS_DIM> >     RealHessVector_t;
  typedef Vector<Tensor<std::complex<double>,OHMMS_DIM> >  ComplexHessVector_t;
  typedef typename AtomicOrbitalTraits<StorageType>::SplineType SplineType;

private:
  // Store in order
  // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
  std::vector<StorageType> YlmVec, dYlm_dthetaVec, dYlm_dphiVec, ulmVec, dulmVec, d2ulmVec;

  SplineType *RadialSpline;
  // The first index is n in r^n, the second is lm = l*(l+1)+m
  Array<StorageType,3> PolyCoefs;
  NewTimer YlmTimer, SplineTimer, SumTimer;
  RealType rmagLast;
  std::vector<PosType> TwistAngles;
public:
  PosType Pos;
  RealType CutoffRadius, SplineRadius, PolyRadius;
  int SplinePoints;
  int PolyOrder;
  int lMax, Numlm, NumBands;
  UnitCellType Lattice;

  inline void set_pos  (PosType pos)
  {
    Pos = pos;
  }
  inline void set_lmax (int lmax)
  {
    lMax = lmax;
  }
  inline void set_cutoff (RealType cutoff)
  {
    CutoffRadius = cutoff;
  }
  inline void set_spline (RealType radius, int points)
  {
    SplineRadius = radius;
    SplinePoints = points;
  }
  inline void set_polynomial (RealType radius, int order)
  {
    PolyRadius = radius;
    PolyOrder = order;
  }
  inline void set_num_bands (int num_bands)
  {
    NumBands = num_bands;
  }
  SplineType* get_radial_spline ()
  {
    return RadialSpline;
  }
  Array<StorageType,3>& get_poly_coefs()
  {
    return PolyCoefs;
  }

  inline void registerTimers()
  {
    YlmTimer.reset();
    SplineTimer.reset();
    TimerManager.addTimer (&YlmTimer);
    TimerManager.addTimer (&SplineTimer);
    TimerManager.addTimer (&SumTimer);
  }

  void allocate();

  void set_band (int band, Array<std::complex<double>,2> &spline_data,
                 Array<std::complex<double>,2> &poly_coefs,
                 PosType twist);
  inline void CalcYlm(PosType rhat,
                      std::vector<std::complex<double> > &Ylm,
                      std::vector<std::complex<double> > &dYlm_dtheta,
                      std::vector<std::complex<double> > &dYlm_dphi);

  inline void CalcYlm(PosType rhat,
                      std::vector<double> &Ylm,
                      std::vector<double> &dYlm_dtheta,
                      std::vector<double> &dYlm_dphi);

  inline bool evaluate (PosType r, ComplexValueVector_t &vals);
  inline bool evaluate (PosType r, ComplexValueVector_t &val,
                        ComplexGradVector_t &grad,
                        ComplexValueVector_t &lapl);
  inline bool evaluate (PosType r, ComplexValueVector_t &val,
                        ComplexGradVector_t &grad,
                        ComplexHessVector_t &lapl);
  inline bool evaluate (PosType r, RealValueVector_t &vals);
  inline bool evaluate (PosType r, RealValueVector_t &val,
                        RealGradVector_t &grad,
                        RealValueVector_t &lapl);
  inline bool evaluate (PosType r, RealValueVector_t &val,
                        RealGradVector_t &grad,
                        RealHessVector_t &lapl);


  AtomicOrbital() : RadialSpline(NULL),
    YlmTimer("AtomicOrbital::CalcYlm"),
    SplineTimer("AtomicOrbital::1D spline"),
    SumTimer("AtomicOrbital::Summation"),
    rmagLast(std::numeric_limits<RealType>::max())
  {
    // Nothing else for now
  }
};


template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r, ComplexValueVector_t &vals)
{
  PosType dr = r - Pos;
  PosType u = Lattice.toUnit(dr);
  PosType img;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    img[i] = round(u[i]);
    u[i] -= img[i];
  }
  dr = Lattice.toCart(u);
  double r2 = dot(dr,dr);
  if (r2 > CutoffRadius*CutoffRadius)
    return false;
  double rmag = std::sqrt(r2);
  PosType rhat = (1.0/rmag)*dr;
  // Evaluate Ylm's
  CalcYlm (rhat, YlmVec, dYlm_dthetaVec, dYlm_dphiVec);
  if (std::abs(rmag - rmagLast) > 1.0e-6)
  {
    // Evaluate radial functions
    if (rmag > PolyRadius)
      EinsplineMultiEval (RadialSpline, rmag, &(ulmVec[0]));
    else
    {
      for (int index=0; index<ulmVec.size(); index++)
        ulmVec[index]  = StorageType();
      double r2n = 1.0;
      for (int n=0; n <= PolyOrder; n++)
      {
        int index = 0;
        for (int i=0; i<vals.size(); i++)
          for (int lm=0; lm<Numlm; lm++)
            ulmVec[index++] += r2n*PolyCoefs(n,i,lm);
        r2n *= rmag;
      }
    }
    rmagLast = rmag;
  }
  SumTimer.start();
  int index = 0;
  for (int i=0; i<vals.size(); i++)
  {
    vals[i] = std::complex<double>();
    for (int lm=0; lm < Numlm; lm++)
      vals[i] += ulmVec[index++] * YlmVec[lm];
    double phase = -2.0*M_PI*dot(TwistAngles[i],img);
    // fprintf (stderr, "phase[%d] = %1.2f pi\n", i, phase/M_PI);
    // fprintf (stderr, "img = [%f,%f,%f]\n", img[0], img[1], img[2]);
    double s,c;
    sincos(phase,&s,&c);
    vals[i] *= std::complex<double>(c,s);
  }
  SumTimer.stop();
  return true;
}


template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r, RealValueVector_t &vals)
{
  PosType dr = r - Pos;
  PosType u = Lattice.toUnit(dr);
  PosType img;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    img[i] = round(u[i]);
    u[i] -= img[i];
  }
  dr = Lattice.toCart(u);
  double r2 = dot(dr,dr);
  if (r2 > CutoffRadius*CutoffRadius)
    return false;
  double rmag = std::sqrt(r2);
  PosType rhat = (1.0/rmag)*dr;
  // Evaluate Ylm's
  CalcYlm (rhat, YlmVec, dYlm_dthetaVec, dYlm_dphiVec);
  if (std::abs(rmag - rmagLast) > 1.0e-6)
  {
    // Evaluate radial functions
    if (rmag > PolyRadius)
    {
      SplineTimer.start();
      EinsplineMultiEval (RadialSpline, rmag, &(ulmVec[0]));
      SplineTimer.stop();
    }
    else
    {
      for (int index=0; index<ulmVec.size(); index++)
        ulmVec[index] = StorageType();
      double r2n = 1.0;
      for (int n=0; n <= PolyOrder; n++)
      {
        int index = 0;
        for (int i=0; i<vals.size(); i++)
          for (int lm=0; lm<Numlm; lm++)
            ulmVec[index++] += r2n*PolyCoefs(n,i,lm);
        r2n *= rmag;
      }
    }
    rmagLast = rmag;
  }
  SumTimer.start();
  int index = 0;
  for (int i=0; i<vals.size(); i++)
  {
    vals[i] = 0.0;
    StorageType tmp = 0.0;
    for (int lm=0; lm < Numlm; lm++, index++)
      tmp += ulmVec[index] * YlmVec[lm];
    //vals[i] += real(ulmVec[index++] * YlmVec[lm]);
    // vals[i] += (ulmVec[index].real() * YlmVec[lm].real() -
    // 	    ulmVec[index].imag() * YlmVec[lm].imag());
    double phase = -2.0*M_PI*dot(TwistAngles[i],img);
    double s,c;
    sincos(phase,&s,&c);
    vals[i] = real(std::complex<double>(c,s) * tmp);
  }
  SumTimer.stop();
  return true;
}

template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r,
                                      RealValueVector_t &vals,
                                      RealGradVector_t  &grads,
                                      RealHessVector_t &hess)
{
  APP_ABORT(" AtomicOrbital<StorageType>::evaluate not implemented for Hess. \n");
  return true;
}


template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r,
                                      RealValueVector_t &vals,
                                      RealGradVector_t  &grads,
                                      RealValueVector_t &lapl)
{
  PosType dr = r - Pos;
  PosType u = Lattice.toUnit(dr);
  PosType img;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    img[i] = round(u[i]);
    u[i] -= img[i];
  }
  dr = Lattice.toCart(u);
  double r2 = dot(dr,dr);
  if (r2 > CutoffRadius*CutoffRadius)
    return false;
  double rmag = std::sqrt(r2);
  double rInv = 1.0/rmag;
  PosType rhat = rInv*dr;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  PosType thetahat = PosType(costheta*cosphi,
                             costheta*sinphi,
                             -sintheta);
  PosType phihat   = PosType(-sinphi, cosphi, 0.0 );
  // Evaluate Ylm's
  CalcYlm (rhat, YlmVec, dYlm_dthetaVec, dYlm_dphiVec);
  // Evaluate radial functions
  if (rmag > PolyRadius)
  {
    SplineTimer.start();
    EinsplineMultiEval
    (RadialSpline, rmag, &(ulmVec[0]), &(dulmVec[0]), &(d2ulmVec[0]));
    SplineTimer.stop();
  }
  else
  {
    for (int index=0; index<ulmVec.size(); index++)
    {
      ulmVec[index]   = StorageType();
      dulmVec[index]  = StorageType();
      d2ulmVec[index] = StorageType();
    }
    double r2n = 1.0, r2nm1=0.0, r2nm2=0.0;
    double dn = 0.0;
    double dnm1 = -1.0;
    for (int n=0; n <= PolyOrder; n++)
    {
      int index = 0;
      for (int i=0; i<vals.size(); i++)
        for (int lm=0; lm<Numlm; lm++,index++)
        {
          StorageType c = PolyCoefs(n,i,lm);
          ulmVec[index]   += r2n*c;
          dulmVec[index]  += dn * r2nm1 * c;
          d2ulmVec[index] += dn*dnm1 * r2nm2 * c;
        }
      dn += 1.0;
      dnm1 += 1.0;
      r2nm2 = r2nm1;
      r2nm1 = r2n;
      r2n *= rmag;
    }
  }
  SumTimer.start();
  int index = 0;
  for (int i=0; i<vals.size(); i++)
  {
    vals[i] = 0.0;
    for (int j=0; j<OHMMS_DIM; j++)
      grads[i][j] = 0.0;
    lapl[i] = 0.0;
    // Compute e^{-i k.L} phase factor
    double phase = -2.0*M_PI*dot(TwistAngles[i],img);
    double s,c;
    sincos(phase,&s,&c);
    std::complex<double> e2mikr(c,s);
    StorageType tmp_val, tmp_lapl,
                grad_rhat, grad_thetahat, grad_phihat;
    tmp_val = tmp_lapl = grad_rhat = grad_thetahat = grad_phihat =
                                       StorageType();
    int lm=0;
    for (int l=0; l<= lMax; l++)
      for (int m=-l; m<=l; m++,lm++,index++)
      {
        std::complex<double> im(0.0,(double)m);
        tmp_val       += ulmVec[index]  * YlmVec[lm];
        grad_rhat     += dulmVec[index] * YlmVec[lm];
        grad_thetahat += ulmVec[index]  * rInv * dYlm_dthetaVec[lm];
        grad_phihat   += (ulmVec[index] * dYlm_dphiVec[lm])/(rmag*sintheta);
        //grad_phihat += (ulmVec[index] * im *YlmVec[lm])/(rmag*sintheta);
        tmp_lapl += YlmVec[lm] *
                    (-(double)(l*(l+1))*rInv*rInv * ulmVec[index]
                     + d2ulmVec[index] + 2.0*rInv *dulmVec[index]);
      }
    vals[i]  = real(e2mikr*tmp_val);
    lapl[i]  = real(e2mikr*tmp_lapl);
    grads[i] = (real(e2mikr*grad_rhat    ) * rhat     +
                real(e2mikr*grad_thetahat) * thetahat +
                real(e2mikr*grad_phihat  ) * phihat);
  }
  SumTimer.stop();
  rmagLast = rmag;
  return true;
}

template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r,
                                      ComplexValueVector_t &vals,
                                      ComplexGradVector_t  &grads,
                                      ComplexHessVector_t &hess)
{
  APP_ABORT(" AtomicOrbital<StorageType>::evaluate not implemented for Hess. \n");
  return true;
}

template<typename StorageType> inline bool
AtomicOrbital<StorageType>::evaluate (PosType r, ComplexValueVector_t &vals,
                                      ComplexGradVector_t &grads,
                                      ComplexValueVector_t &lapl)
{
  PosType dr = r - Pos;
  PosType u = Lattice.toUnit(dr);
  PosType img;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    img[i] = round(u[i]);
    u[i] -= img[i];
  }
  dr = Lattice.toCart(u);
  double r2 = dot(dr,dr);
  if (r2 > CutoffRadius*CutoffRadius)
    return false;
  double rmag = std::sqrt(r2);
  double rInv = 1.0/rmag;
  PosType rhat = rInv*dr;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  PosType thetahat = PosType(costheta*cosphi,
                             costheta*sinphi,
                             -sintheta);
  PosType phihat   = PosType(-sinphi, cosphi, 0.0 );
  // Evaluate Ylm's
  CalcYlm (rhat, YlmVec, dYlm_dthetaVec, dYlm_dphiVec);
  // Evaluate radial functions
  if (rmag > PolyRadius)
  {
    SplineTimer.start();
    EinsplineMultiEval
    (RadialSpline, rmag, &(ulmVec[0]), &(dulmVec[0]), &(d2ulmVec[0]));
    SplineTimer.stop();
  }
  else
  {
    for (int index=0; index<ulmVec.size(); index++)
    {
      ulmVec[index]   = StorageType();
      dulmVec[index]  = StorageType();
      d2ulmVec[index] = StorageType();
    }
    double r2n = 1.0, r2nm1=0.0, r2nm2=0.0;
    double dn = 0.0;
    double dnm1 = -1.0;
    for (int n=0; n <= PolyOrder; n++)
    {
      int index = 0;
      for (int i=0; i<vals.size(); i++)
        for (int lm=0; lm<Numlm; lm++,index++)
        {
          StorageType   c = PolyCoefs(n,i,lm);
          ulmVec[index]   += r2n*c;
          dulmVec[index]  += dn * r2nm1 * c;
          d2ulmVec[index] += dn*dnm1 * r2nm2 * c;
        }
      dn += 1.0;
      dnm1 += 1.0;
      r2nm2 = r2nm1;
      r2nm1 = r2n;
      r2n *= rmag;
    }
  }
  SumTimer.start();
  int index = 0;
  for (int i=0; i<vals.size(); i++)
  {
    vals[i] = 0.0;
    for (int j=0; j<OHMMS_DIM; j++)
      grads[i][j] = 0.0;
    lapl[i] = 0.0;
    int lm=0;
    StorageType grad_rhat, grad_thetahat, grad_phihat;
    // Compute e^{-i k.L} phase factor
    double phase = -2.0*M_PI*dot(TwistAngles[i],img);
    double s,c;
    sincos(phase,&s,&c);
    std::complex<double> e2mikr(c,s);
    for (int l=0; l<= lMax; l++)
      for (int m=-l; m<=l; m++,lm++,index++)
      {
        std::complex<double> im(0.0,(double)m);
        vals[i]  += ulmVec[index] * YlmVec[lm];
        grad_rhat     += dulmVec[index] * YlmVec[lm];
        grad_thetahat += ulmVec[index] * rInv * dYlm_dthetaVec[lm];
        grad_phihat   += (ulmVec[index] * im*YlmVec[lm])/(rmag*sintheta);
        lapl[i] += YlmVec[lm] *
                   (-(double)(l*(l+1))*rInv*rInv * ulmVec[index]
                    + d2ulmVec[index] + 2.0*rInv *dulmVec[index]);
      }
    vals[i] *= e2mikr;
    lapl[i] *= e2mikr;
    for (int j=0; j<OHMMS_DIM; j++)
    {
      grads[i][j] = e2mikr*(grad_rhat*rhat[j] + grad_thetahat*thetahat[j]
                            + grad_phihat * phihat[j]);
    }
  }
  SumTimer.stop();
  rmagLast = rmag;
  return true;
}





// Fast implementation
// See Geophys. J. Int. (1998) 135,pp.307-309
template<typename StorageType> inline void
AtomicOrbital<StorageType>::CalcYlm (PosType rhat,
                                     std::vector<std::complex<double> > &Ylm,
                                     std::vector<std::complex<double> > &dYlm_dtheta,
                                     std::vector<std::complex<double> > &dYlm_dphi)
{
  YlmTimer.start();
  const double fourPiInv = 0.0795774715459477;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cottheta = costheta/sintheta;
  double cosphi, sinphi;
  cosphi=rhat[0]/sintheta;
  sinphi=rhat[1]/sintheta;
  std::complex<double> e2iphi(cosphi, sinphi);
  double lsign = 1.0;
  double dl = 0.0;
  double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
  for (int l=0; l<=lMax; l++)
  {
    XlmVec[2*l]  = lsign;
    dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
    XlmVec[0]    = lsign*XlmVec[2*l];
    dXlmVec[0]   = lsign*dXlmVec[2*l];
    double dm = dl;
    double msign = lsign;
    for (int m=l; m>0; m--)
    {
      double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
      XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
      dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
      // Copy to negative m
      XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
      dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
      msign *= -1.0;
      dm -= 1.0;
    }
    double sum = 0.0;
    for (int m=-l; m<=l; m++)
      sum += XlmVec[l+m]*XlmVec[l+m];
    // Now, renormalize the Ylms for this l
    double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
    for (int m=-l; m<=l; m++)
    {
      XlmVec[l+m]  *= norm;
      dXlmVec[l+m] *= norm;
    }
    // Multiply by azimuthal phase and store in Ylm
    std::complex<double> e2imphi (1.0, 0.0);
    std::complex<double> eye(0.0, 1.0);
    for (int m=0; m<=l; m++)
    {
      Ylm[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
      Ylm[l*(l+1)-m]  =  XlmVec[l-m]*qmcplusplus::conj(e2imphi);
      dYlm_dphi[l*(l+1)+m ]  =  (double)m * eye *XlmVec[l+m]*e2imphi;
      dYlm_dphi[l*(l+1)-m ]  = -(double)m * eye *XlmVec[l-m]*qmcplusplus::conj(e2imphi);
      dYlm_dtheta[l*(l+1)+m] = dXlmVec[l+m]*e2imphi;
      dYlm_dtheta[l*(l+1)-m] = dXlmVec[l-m]*qmcplusplus::conj(e2imphi);
      e2imphi *= e2iphi;
    }
    dl += 1.0;
    lsign *= -1.0;
  }
  YlmTimer.stop();
}

// Fast implementation
// See Geophys. J. Int. (1998) 135,pp.307-309
template<typename StorageType> inline void
AtomicOrbital<StorageType>::CalcYlm (PosType rhat,
                                     std::vector<double> &Ylm,
                                     std::vector<double> &dYlm_dtheta,
                                     std::vector<double> &dYlm_dphi)
{
  YlmTimer.start();
  const double fourPiInv = 0.0795774715459477;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cottheta = costheta/sintheta;
  double cosphi, sinphi;
  cosphi=rhat[0]/sintheta;
  sinphi=rhat[1]/sintheta;
  std::complex<double> e2iphi(cosphi, sinphi);
  double lsign = 1.0;
  double dl = 0.0;
  double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
  for (int l=0; l<=lMax; l++)
  {
    XlmVec[2*l]  = lsign;
    dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
    XlmVec[0]    = lsign*XlmVec[2*l];
    dXlmVec[0]   = lsign*dXlmVec[2*l];
    double dm = dl;
    double msign = lsign;
    for (int m=l; m>0; m--)
    {
      double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
      XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
      dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
      // Copy to negative m
      XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
      dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
      msign *= -1.0;
      dm -= 1.0;
    }
    double sum = 0.0;
    for (int m=-l; m<=l; m++)
      sum += XlmVec[l+m]*XlmVec[l+m];
    // Now, renormalize the Ylms for this l
    double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
    for (int m=-l; m<=l; m++)
    {
      XlmVec[l+m]  *= norm;
      dXlmVec[l+m] *= norm;
    }
    // Multiply by azimuthal phase and store in Ylm
    Ylm[l*(l+1)]         =  XlmVec[l];
    dYlm_dphi[l*(l+1) ]  = 0.0;
    dYlm_dtheta[l*(l+1)] = dXlmVec[l];
    std::complex<double> e2imphi = e2iphi;
    for (int m=1; m<=l; m++)
    {
      Ylm[l*(l+1)+m]         =  XlmVec[l+m]*e2imphi.real();
      Ylm[l*(l+1)-m]         =  XlmVec[l+m]*e2imphi.imag();
      dYlm_dphi[l*(l+1)+m ]  = -(double)m * XlmVec[l+m] *e2imphi.imag();
      dYlm_dphi[l*(l+1)-m ]  =  (double)m * XlmVec[l+m] *e2imphi.real();
      dYlm_dtheta[l*(l+1)+m] = dXlmVec[l+m]*e2imphi.real();
      dYlm_dtheta[l*(l+1)-m] = dXlmVec[l+m]*e2imphi.imag();
      e2imphi *= e2iphi;
    }
    dl += 1.0;
    lsign *= -1.0;
  }
  YlmTimer.stop();
}




}
#endif
