#include "QMCWaveFunctions/QDwf.h"

void QDwf::reset()
{
  //cout << "Calling QDwf::reset:: " << Ez << " " << z0 << endl;
  fac1 = zB - z0;
  fac2 = pow(2.0*Ez,onethird);
  kappa = sqrt(2.0*VB);
  double azp0;
  for(int iter = 0; iter < 4; iter++)
  {
    double a = -2.33810;
    double b = -1.0;
    while( b - a > 1.e-5 )
    {
      double m = a + 0.5 * ( b - a );
      int fa = funcz(a);
      int fm = funcz(m);
      if( fa == fm )
      {
        a = m;
      }
      else
      {
        b = m;
      }
    }
    azp0 = 0.5 * ( a + b );
    EL = Ez * ( fac1 - azp0 / fac2 );
    kappa = sqrt( 2.0 * ( VB - EL ) );
    //      cout << Ez << '\t' << fac1 << '\t' << azp0 << '\t'
    //	   << fac2 << '\t' << EL << '\t' << VB-EL << '\t' << kappa << endl;
  }
  psiB = gsl_sf_airy_Ai( azp0, mode );
  /// create and initialise the uniform one-dimensional uGrid1D
  int npts = 1000;
  if( !m_grid )
  {
    m_grid = new uGrid1D;
    double zstart = zB;
    double zend = 16.00;
    m_grid->init(zstart,zend,npts);
  }
  /// create and initialise the CubicSpline
  if(!m_spline)
    m_spline = new CubicSpline(m_grid);
  for(int i = 0; i < npts; i++)
  {
    double z = m_grid->m_coord[i];
    double zz = zprime(z);
    //(*m_spline)(i) = gsl_sf_airy_Ai(zz,mode);//evaluate(z);
    m_spline->F(i)[0] = gsl_sf_airy_Ai(zz,mode);
    m_spline->F(i)[1] = fac2*gsl_sf_airy_Ai_deriv(zz,mode);
  }
  //    cout << z0 << '\t' << Ez << '\t' << kappa << '\t' << EL << endl;
  //    exit(-1);
  return;
}

