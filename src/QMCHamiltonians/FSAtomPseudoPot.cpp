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
    
    
/** @file FSAtomPseudoPot.h
 * @brief Xml Parser Definition for FSAtom Pseudopotential Standard
 */
#include "QMCHamiltonians/FSAtomPseudoPot.h"
#include <fstream>
namespace qmcplusplus
{

void FSAtomPseudoPot::convert2HartreeBohr(RealType sc, bool is_r_times_v)
{
  if(is_r_times_v)
    for(int i=0; i<myFunc.size(); i++)
      myFunc(i) *= sc;
  else
    for(int i=0; i<myFunc.size(); i++)
      myFunc(i) *= sc*myFunc.r(i);
}

FSAtomPseudoPot::RealType FSAtomPseudoPot::getCutOff(RealType v0)
{
  bool ignore=true;
  int ng=myFunc.size()-2;
  RealType r=myFunc.r(ng);
  while(ignore&& ng)
  {
    std::cout << r << " " << myFunc(ng) << std::endl;
    ignore=(std::abs(myFunc(ng)-v0)<1e-12);
    r=myFunc.r(--ng);
  }
  return r;
}

/** create a LinearSpline<RealType>
 * @param sc scaling factor
 * @return a LinearSpline<RealType> on a LinearGrid
 */
FSAtomPseudoPot::return_type*
FSAtomPseudoPot::getLocalPot(RealType zeff)
{
  myFunc.spline();
  const RealType del=1.0e-3;
  //int ng=static_cast<int>(2*Rcut/del)+1;
  int ng=static_cast<int>(Rcut/del)+1;
  app_log() << "  FSAtomPseudoPot::getLocalPot grid = [0," << Rcut << "] npts=" << ng << std::endl;
  LinearGrid<RealType> *agrid=new LinearGrid<RealType>;
  agrid->set(0.0,Rcut,ng);
  RealType sc=-1.0/zeff;
  return_type* newFunc=new return_type(agrid);
  (*newFunc)(0) = 0.0;
  for(int i=1; i<ng-1; i++)
  {
    (*newFunc)(i)=sc*myFunc.splintNG((*agrid)[i]);
  }
  (*newFunc)(ng-1)=1.0;
  RealType y_prime=((*newFunc)(1)-(*newFunc)(0))/del;
  newFunc->spline(0,y_prime,ng-1,0.0);
  std::ofstream fout("pp.loc.dat");
  fout << "#  Local pseudopotential -rV(r)/Zeff 1 beyond 2*Rcut " << std::endl;
  fout.setf(std::ios::scientific, std::ios::floatfield);
  fout.precision(12);
  for(int i=1; i<ng; i++)
    fout << (*agrid)[i] << " " << (*newFunc)(i) << std::endl;
  fout << std::endl;
  for(int i=0; i<myFunc.size(); i++)
    fout << myFunc.r(i) << " " <<sc*myFunc.m_Y[i] << std::endl;
  return newFunc;
}

FSAtomPseudoPot::return_type*
FSAtomPseudoPot::getNonLocalPot(FSAtomPseudoPot& vloc)
{
  //remove local part
  myFunc.m_Y -= vloc.myFunc.m_Y;
  //RealType y_prime=(m_Y[1]-m_Y[0])/m_grid->dr(0);
  RealType y_prime=(myFunc(1)-myFunc(0))/myFunc.dr(0);
  myFunc.spline(0,y_prime,myFunc.size()-1,0.0);
  const RealType del=1.0e-3;
  int ng=static_cast<int>(Rcut/del)+1;
  LinearGrid<RealType> *agrid=new LinearGrid<RealType>;
  agrid->set(0.0,Rcut,ng);
  return_type* newFunc=new return_type(agrid);
  for(int i=1; i<ng-1; i++)
  {
    (*newFunc)(i)=myFunc.splintNG((*agrid)[i])/(*agrid)[i];
  }
  //force the boudary conditions
  (*newFunc)(0) =(*newFunc)(1);
  (*newFunc)(ng-1)=0.0;
  newFunc->spline();
//    char fname[32];
//    sprintf(fname,"pp.L%d.dat",AngL);
//    std::ofstream fout(fname);
//    fout.setf(std::ios::scientific, std::ios::floatfield);
//    fout.precision(12);
//    for(int i=0; i<ng; i++)
//      fout << (*agrid)[i] << " " << (*newFunc)(i) << std::endl;
//    fout <<  std::endl;
//    for(int i=0; i<myFunc.size(); i++)
//      fout << myFunc.r(i) << " " <<vFac*myFunc.m_Y[i]/myFunc.r(i) << std::endl;
  return newFunc;
}

bool FSAtomPseudoPot::put(xmlNodePtr cur)
{
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname ((const char*)cur->name);
    if(cname == "radfunc")
    {
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1=="data")
        {
          myFunc.m_Y.resize(myFunc.grid().size());
          putContent(myFunc.m_Y,cur1);
        }
        //else if(cname1 =="grid")
        //{
        //  app_warning() << "    FSAtomPseudoPot::vps/grid is ignored " << std::endl;
        //}
        cur1=cur1->next;
      }
    }
    cur=cur->next;
  }
  return true;
}
} // namespace qmcPlusPlus
