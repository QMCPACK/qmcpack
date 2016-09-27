//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef GUARD_QDWF_H
#define GUARD_QDWF_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/uGrid1D.h"
#include "Numerics/Spline3D/CubicSpline.h"
#include "Optimize/VarList.h"
#include <gsl/gsl_sf_airy.h>
#include <math.h>
#include <iostream>


class QDwf
{


  std::vector<double> X;
  std::vector<double> sigx;
  std::vector<double> a;

  double sigy;//var
  double Y0;//var

  double z0;//var
  double Ez;//var

  double VB;
  double zB;
  double kappa;
  double EL;
  double psiB;

  double onethird;

  double fac1,fac2;


  uGrid1D* m_grid;
  CubicSpline* m_spline;

  gsl_mode_t mode;

  int funcz(double z)
  {
    double denom = gsl_sf_airy_Ai(z,mode);
    double numer = fac2*gsl_sf_airy_Ai_deriv(z,mode);
    double xsign = numer-denom*kappa;
    if(xsign < 0 )
    {
      return -1;
    }
    else
    {
      return 1;
    }
  }

  inline double zprime(double z)
  {
    return fac2*(z-z0-EL/Ez);
  }

public:

  typedef double value_type;

  QDwf()
  {
    VB = 24.535;
    zB = 9.1854;
    //    z0 = 10.6858;
    //    Ez = 2.2;
    Ez = 2.04866;
    z0 = 10.8049;
    //    sigxr = 0.08;
    //    xr = 80.72;
    sigy = 0.107;
    Y0 = 25.3;
    /// set the x parameters
    int nxp = 1;
    a.resize(nxp);
    sigx.resize(nxp);
    X.resize(nxp);
    a[0] = 1.00000;
    sigx[0] = 0.08;
    X[0] = 80.628;
    //    a[1] = 1.00e-4; sigx[1] = 0.08;  X[1] = 62;
    /*
    a[0] = 1.00000; sigx[0] = 0.08; X[0] = 80.628;
    a[1] = 0.08260; sigx[1] = 0.08; X[1] = 76.000;
    a[2] = 2.17e-4; sigx[2] = 0.04; X[2] = 67.000;
    */
    /// optimised parameters
    /*
    a[0] = 1.00000; sigx[0] = 0.092403; X[0] = 80.569;
    a[1] = 0.071703; sigx[1] = 0.18445; X[1] = 76.064;
    a[2] = 2.86e-3; sigx[2] = 0.052224; X[2] = 67.000;
    sigy = 0.11112;
    Y0 = 25.167;
    */
    onethird = 1.0/3.0;
    reset();
    std::cout << "Airy Solved: " << EL << '\t' << kappa << '\t'
         << psiB << std::endl;
    /*
    /// create and initialise the uniform one-dimensional uGrid1D
    if( !m_grid ) m_grid = new uGrid1D;
    double zstart = zB; double zend = 16.00; int npts = 1000;
    m_grid->init(zstart,zend,npts);


    /// create and initialise the CubicSpline
    if(!m_spline) m_spline = new CubicSpline(m_grid);
    for(int i = 0; i < npts; i++){
      double z = m_grid->m_coord[i]; double zz = zprime(z);
      //(*m_spline)(i) = gsl_sf_airy_Ai(zz,mode);//evaluate(z);
      m_spline->F(i)[0] = gsl_sf_airy_Ai(zz,mode);
      m_spline->F(i)[1] = fac2*gsl_sf_airy_Ai_deriv(zz,mode);
    }
    */
  }

  template<class T1>
  void put(xmlNodePtr cur, VarRegistry<T1>& vlist)
  {
    //    vlist.add("C_sigxr",&sigxr,1);
    //    vlist.add("C_xr",&xr,1);
    vlist.add("C_X",&X[0],1);
    vlist.add("C_sigx",&sigx[0],1);
    vlist.add("C_a",&a[0],1);
    vlist.add("C_sigy",&sigy,1);
    vlist.add("C_Y0",&Y0,1);
    vlist.add("C_Ez",&Ez,1);
    vlist.add("C_z0",&z0,1);
  }

  void reset();

  inline void set_point(const posvec_t& r) {}

  inline double evaluate(const posvec_t& r,
                         posvec_t& gradf,
                         double& lapf)
  {
    int nsize = sigx.size();
    /// \f$\Psi_{x}\f$ sum of Gaussians
    double psix = 0.0,gfx=0.0,lfx=0.0;
    for(int i = 0; i < nsize; i++)
    {
      double x = r[0] - X[i];
      double sx = sigx[i] * x;
      double sxx = sx * x;
      double phix = a[i] * exp(-sxx);
      psix += phix;
      gfx += -2 * sx * phix;
      lfx += -2 * sigx[i] * ( 1.0 - 2.0 * sxx ) * phix;
    }
    /// \f$\Psi_{y}\f$ single Gaussian
    double y = r[1] - Y0;
    double sy = sigy * y;
    double syy = sy * y;
    double psiy = exp(-syy);
    double gfy = - 2 * sy * psiy;
    double lfy = - 2 * sigy * ( 1.0 - 2.0 * syy ) * psiy;
    /// \f$\Psi_{z}\f$ AiryAi with exponential tail
    double z = r[2];
    double psiz,gfz,lfz;
    if( z <= zB )
    {
      psiz = psiB * exp ( kappa * ( z - zB ) );
      gfz = kappa * psiz;
      lfz = kappa * gfz;
    }
    else
    {
      psiz = m_spline->evaluate(z,gfz,lfz);
    }
    double psi = psix * psiy * psiz;
    gradf[0] = gfx * psiy * psiz;
    gradf[1] = psix * gfy * psiz;
    gradf[2] = psix * gfy * psiz;
    lapf = lfx * psiy * psiz + psix * lfy * psiz + psix * psiy * lfz;
    return psi;
  }


};
#endif
