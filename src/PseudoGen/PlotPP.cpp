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
  std::ofstream gplot(fname);
  gplot.precision(10);
  gplot.setf(std::ios::scientific,std::ios::floatfield);
  std::cout << "Writing Pseudopotential to file " << fname << std::endl;
  for(int i=0; i<Psi.m_grid->size(); i++)
  {
    value_type r =  Psi.m_grid->r(i);
    value_type SJ_num = 1.0-exp(-Params(0)*r);
    value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
    gplot << std::setw(20) << r << std::setw(20) <<  (-1.0*4.0/r)*(SJ_num/SJ_den) << std::endl;
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
  std::cout << "Total valence charge = " << integrate_RK2(VCharge) << std::endl;
  std::cout << "Total valence charge = " << integrate_RK2(VCharge_siesta_grid) << std::endl;
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
  std::ofstream siesta(fname);
  std::cout << "Writing Pseudopential to file " << fname << std::endl;
  siesta << " Radial grid follows" << std::endl;
  int nlines = (npts-1)/4;
  int remainder = npts-1-nlines*4;
  siesta.precision(12);
  siesta.setf(std::ios::scientific,std::ios::floatfield);
  siesta << setiosflags(std::ios::uppercase);
  int j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << siesta_grid->r(j)
           << std::setw(20) << siesta_grid->r(j+1)
           << std::setw(20) << siesta_grid->r(j+2)
           << std::setw(20) << siesta_grid->r(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << siesta_grid->r(i);
    siesta << std::endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << std::endl;
  siesta << "  0" << std::endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << PP(j)
           << std::setw(20) << PP(j+1)
           << std::setw(20) << PP(j+2)
           << std::setw(20) << PP(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << PP(i);
    siesta << std::endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << std::endl;
  siesta << "  1" << std::endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << PP(j)
           << std::setw(20) << PP(j+1)
           << std::setw(20) << PP(j+2)
           << std::setw(20) << PP(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << PP(i);
    siesta << std::endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << std::endl;
  siesta << "  2" << std::endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << PP(j)
           << std::setw(20) << PP(j+1)
           << std::setw(20) << PP(j+2)
           << std::setw(20) << PP(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << PP(i);
    siesta << std::endl;
  }
  siesta << " Down Pseudopotential follows (l on next line)" << std::endl;
  siesta << "  3" << std::endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << PP(j)
           << std::setw(20) << PP(j+1)
           << std::setw(20) << PP(j+2)
           << std::setw(20) << PP(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=0; i<remainder; i++)
      siesta << std::setw(20) << PP(i);
    siesta << std::endl;
  }
  double CC = 0.0;
  siesta << " Core charge follows" << std::endl;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << CC
           << std::setw(20) << CC
           << std::setw(20) << CC
           << std::setw(20) << CC
           << std::endl;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << CC;
    siesta << std::endl;
  }
  siesta << " Valence charge follows" << std::endl;
  j=1;
  for(int i=0; i<nlines; i++)
  {
    siesta << std::setw(20) << VCharge_siesta_grid(j)
           << std::setw(20) << VCharge_siesta_grid(j+1)
           << std::setw(20) << VCharge_siesta_grid(j+2)
           << std::setw(20) << VCharge_siesta_grid(j+3)
           << std::endl;
    j+=4;
  }
  if(remainder)
  {
    for(int i=j; i<npts; i++)
      siesta << std::setw(20) << VCharge_siesta_grid(i);
    siesta << std::endl;
  }
  siesta.close();
}


}
