//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "OhmmsPETE/TinyVector.h"
#include "Numerics/SphericalTensor.h"
int main(int argc, char** argv)
{
  int l=5;
  TinyVector<double,3> pos(0.1,0.3,-0.45), deltax(0.0001,0.0,0.0), deltay(0.0,0.0001,0.0), deltaz(0.0,0.0,0.0001), tpos, g;
  if(argc>1)
    l = atoi(argv[1]);
  if(argc>4)
  {
    pos[0] = atof(argv[2]);
    pos[1] = atof(argv[3]);
    pos[2] = atof(argv[4]);
  }
  SphericalTensor<double,TinyVector<double,3> > Y1(l,true), Y2(l,true), Yp(l,true), Ym(l,true);
  Y1.evaluate(pos);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  for(int lm=0; lm<Y1.size(); lm++)
  {
    std::cout << lm << std::endl;
    std::cout << std::setw(20) << std::setprecision(12) << Y1.Ylm[lm]
              << std::setprecision(12) << Y1.gradYlm[lm] << std::endl;
    //std::cout << std::setw(20) << std::setprecision(12) << Y2.Ylm[lm]
    //  << std::setprecision(12) << Y2.gradYlm[lm] << std::endl;
  }
//  for(int lm=0; lm<Y1.size(); lm++) {
//    tpos = pos+deltax; Yp.evaluate(tpos);
//    tpos = pos-deltax; Ym.evaluate(tpos);
//    g[0] = (Yp.Ylm[lm]-Ym.Ylm[lm])/0.0002;
//    tpos = pos+deltay; Yp.evaluate(tpos);
//    tpos = pos-deltay; Ym.evaluate(tpos);
//    g[1] = (Yp.Ylm[lm]-Ym.Ylm[lm])/0.0002;
//    tpos = pos+deltaz; Yp.evaluate(tpos);
//    tpos = pos-deltaz; Ym.evaluate(tpos);
//    g[2] = (Yp.Ylm[lm]-Ym.Ylm[lm])/0.0002;
//    std::cout << lm << std::endl;
//    std::cout << std::setw(20) << std::setprecision(12) << Y1.Ylm[lm]
//      << std::setprecision(12) << Y1.gradYlm[lm] - g << std::endl;
//  }
}
