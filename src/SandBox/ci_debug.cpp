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
    
    



/**@file excitation.cpp
 * @brief Test code for multideterminant tree
 */
#include <Utilities/OhmmsInfo.h>
#include <Utilities/RandomGenerator.h>
#include <Utilities/Timer.h>
#include <Message/Communicate.h>
#include <Message/OpenMP.h>
#include <Numerics/MatrixOperators.h>
#include <Numerics/DeterminantOperators.h>
#include <QMCWaveFunctions/Fermion/ci_node.h>
#include <QMCWaveFunctions/Fermion/ci_builder.h>
#include <Utilities/IteratorUtility.h>
#include <fstream>
#include <limits>
using namespace qmcplusplus;

template<typename CT, typename T>
inline void check_ratios(const CT& a, const CT& b, T eps)
{
  T del;
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  for(int i=0; i<a.size(); ++i)
    if(std::abs(del=(a[i]/b[i]-1.0))>eps)
      std::cout << i << " recursive=" << a[i] << " error=" << del << std::endl;
}


int main(int argc, char** argv)
{
  const int max_states=8;
  const int vmax=4;
  const int cmax=8;
  typedef ci_node_proxy node_type;
  std::vector<node_type> excitations;
  int multiplet=3;
  int nparent=0;
  int nchildren=0;
  std::vector<int> det_offset(multiplet+2,1);
  //add zero
  excitations.push_back(node_type());
  {
    ci_builder<max_states> ci_maker(vmax,cmax);
    //single is treated differently
    det_offset[1]=ci_maker.singles(excitations);
    for(int level=0,k=1; level<multiplet; ++level,++k)
      det_offset[k+1]=ci_maker.promote(excitations,level);
    for(int i=0; i<excitations.size(); ++i)
    {
      nchildren+=excitations[i].children.size();
      if(excitations[i].closed())
        continue;
      nparent++;
    }
  }
  std::ofstream fout("tree.xml");
  int npeers=0;
  int count=0;
  excitations[0].write_node<max_states>(fout,0,count,excitations);
  count=0;
  const int M=8;
  ci_node<double> CI;
  CI.build_tree(M,excitations);
  CI.set_peers(M,vmax);
  std::ofstream fout1("tree_1.xml");
  CI.write_node(fout1);
  typedef Matrix<double> mat_t;
  mat_t psiv(M,M), psic(cmax,M), Identity(M,M);
  //build up big  psi
  for(int i=0; i<psiv.size(); ++i)
    psiv(i)=Random();
  for(int i=0; i<psic.size(); ++i)
    psic(i)=Random();
  const double eps=1e-12;
  std::vector<double> dets(excitations.size()),ratios(excitations.size());
  {
    double det_0=CI.getRatios(psiv,psic,ratios,true);
    CI.debugRatios(psiv,psic,dets,true);
    std::cout << "Checking row replacement " << std::endl;
    for(int i=0; i<ratios.size(); ++i)
      ratios[i]*=det_0;
    check_ratios(ratios,dets,eps);
  }
  {
    double det_0=CI.getRatios(psiv,psic,ratios,false);
    CI.debugRatios(psiv,psic,dets,false);
    std::cout << "Checking column replacement " << std::endl;
    for(int i=0; i<ratios.size(); ++i)
      ratios[i]*=det_0;
    check_ratios(ratios,dets,eps);
  }
  return 0;
}
