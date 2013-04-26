#include "SQD/HFConfiguration.h"
#include "PseudoGen/PseudoGen.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimGridBase.h"
#include <fstream>
#include <iostream>


namespace ohmmshf
{

void PseudoGen::plot_ascii()
{
  char* fname = "Ge.pp.ASCII";
  ofstream gplot(fname);
  gplot.precision(10);
  gplot.setf(ios::scientific,ios::floatfield);
  cout << "Writing Pseudopotential to file " << fname << endl;
  for(int i=0; i<Psi.m_grid->size(); i++)
  {
    value_type r =  Psi.m_grid->r(i);
    value_type SJ_num = 1.0-exp(-Params(0)*r);
    value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
    gplot << setw(20) << r << setw(20) <<  (-1.0*4.0/r)*(SJ_num/SJ_den) << endl;
  }
  gplot.close();
}

void PseudoGen::plot_siesta_grid()
{
  double A = 0.125e-01;
  double B = 0.774610055208e-04;
  int npts = 1141;
  RadialGrid_t* siesta_grid = new LogGridZero<double>;
  siesta_grid->set(A,B,npts);
  RadialOrbital_t VCharge_siesta_grid(siesta_grid);
  RadialOrbital_t VCharge(Psi(0));
  for(int ob=0; ob<Psi.size(); ob++)
  {
    for(int j=0; j<VCharge.size(); j++)
    {
      VCharge(j)+= Psi(ob,j)*Psi(ob,j);
      // VCharge(j)+= AEorbitals[ob](j)*AEorbitals[ob](j);
    }
  }
  int imin = 0;
  value_type deriv = (VCharge(imin+1)-VCharge(imin))/VCharge.dr(imin);
  VCharge.spline(imin,deriv,VCharge.size()-1,0.0);
  Transform2GridFunctor<OneDimGridFunctor<value_type>,
                        OneDimGridFunctor<value_type> > transform(VCharge,VCharge_siesta_grid);
  transform.generate();
  cout << "Total valence charge = " << integrate_RK2(VCharge) << endl;
  cout << "Total valence charge = " << integrate_RK2(VCharge_siesta_grid) << endl;
  RadialOrbital_t PP(siesta_grid);
  value_type prefactor = -2.0*4.0;
  for(int i=0; i<npts; i++)
  {
    value_type r =  siesta_grid->r(i);
    value_type SJ_num = 1.0-exp(-Params(0)*r);
    value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
    PP(i) = prefactor*(SJ_num/SJ_den);
  }
  char* fname = "Ge.siesta_grid.psf";
  ofstream siesta(fname);
  cout << "Writing Pseudopential to file " << fname << endl;
  siesta << " Radial grid follows" << endl;
  int nlines = (npts-1)/4;
  int remainder = npts-1-nlines*4;
  siesta.precision(12);
  siesta.setf(ios::scientific,ios::floatfield);
  siesta << setiosflags(ios::uppercase);
  int j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << siesta_grid->r(j)
           << setw(20) << siesta_grid->r(j+1)
           << setw(20) << siesta_grid->r(j+2)
           << setw(20) << siesta_grid->r(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << siesta_grid->r(i);
    siesta << endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << endl;
  siesta << "  0" << endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << PP(j)
           << setw(20) << PP(j+1)
           << setw(20) << PP(j+2)
           << setw(20) << PP(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << PP(i);
    siesta << endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << endl;
  siesta << "  1" << endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << PP(j)
           << setw(20) << PP(j+1)
           << setw(20) << PP(j+2)
           << setw(20) << PP(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << PP(i);
    siesta << endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << endl;
  siesta << "  2" << endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << PP(j)
           << setw(20) << PP(j+1)
           << setw(20) << PP(j+2)
           << setw(20) << PP(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << PP(i);
    siesta << endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << endl;
  siesta << "  3" << endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << PP(j)
           << setw(20) << PP(j+1)
           << setw(20) << PP(j+2)
           << setw(20) << PP(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=0; i<remainder; i++)
      siesta << setw(20) << PP(i);
    siesta << endl;
  }
  double CC = 0.0;
  siesta << " Core charge follows" << endl;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << CC
           << setw(20) << CC
           << setw(20) << CC
           << setw(20) << CC
           << endl;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << CC;
    siesta << endl;
  }
  siesta << " Valence charge follows" << endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << setw(20) << VCharge_siesta_grid(j)
           << setw(20) << VCharge_siesta_grid(j+1)
           << setw(20) << VCharge_siesta_grid(j+2)
           << setw(20) << VCharge_siesta_grid(j+3)
           << endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << setw(20) << VCharge_siesta_grid(i);
    siesta << endl;
  }
  siesta.close();
}


}
