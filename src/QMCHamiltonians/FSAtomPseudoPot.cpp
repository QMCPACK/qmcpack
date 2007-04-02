//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and  Kenneth Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file FSAtomPseudoPot.h
 * @brief Xml Parser Definition for FSAtom Pseudopotential Standard
 */
#include "QMCHamiltonians/FSAtomPseudoPot.h"
namespace qmcplusplus {

  FSAtomPseudoPot::RealType FSAtomPseudoPot::getCutOff(RealType v0)
  {
    bool ignore=true;
    int ng=myFunc.size()-2;
    RealType r=myFunc.r(ng);
    while(ignore&& ng)
    {
      cout << r << " " << myFunc(ng) << endl;
      ignore=(std::abs(myFunc(ng)-v0)<1e-12);
      r=myFunc.r(--ng);
    }
    return r;
  }

  /** create a LinearSpline<RealType>
   * @param sc scaling factor
   * @return a LinearSpline<RealType> on a LinearGrid
   */
  FSAtomPseudoPot::return_type* FSAtomPseudoPot::getLocalPot(RealType sc)
  {
    myFunc.spline();

    const RealType del=1.0e-3;
    int ng=static_cast<int>(2*Rcut/del)+1;
    app_log() << "  FSAtomPseudoPot::getLocalPot grid = [0," << 2*Rcut << "] npts=" << ng << endl;

    LinearGrid<RealType> *agrid=new LinearGrid<RealType>;
    agrid->set(0.0,2*Rcut,ng);

    return_type* newFunc=new return_type(agrid);
    (*newFunc)(0) = 0.0;
    for(int i=1; i<ng-1; i++)
    {
      (*newFunc)(i)=sc*myFunc.splintNG((*agrid)[i]);
    }
    (*newFunc)(ng-1)=1.0;

    RealType y_prime=((*newFunc)(1)-(*newFunc)(0))/del;
    newFunc->spline(0,y_prime,ng-1,0.0);
    return newFunc;
  }

  FSAtomPseudoPot::return_type* FSAtomPseudoPot::getNonLocalPot(FSAtomPseudoPot& vloc, RealType vFac)
  {
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
      (*newFunc)(i)=vFac*myFunc.splintNG((*agrid)[i])/(*agrid)[i];
    }
    (*newFunc)(0) =(*newFunc)(1);
    //(*newFunc)(0) = 2.0*(*newFunc)(1) - (*newFunc)(2);
    (*newFunc)(ng-1)=0.0;

    //for(int i=0; i<ng; i++)
    //  cout << (*agrid)[i] << " " << (*newFunc)(i) << endl;
    //for(int i=1; i<m_Y.size(); i++)
    //  cout << (*m_grid)[i] << " " << m_Y[i]/(*m_grid)[i] << endl;
    newFunc->spline();
    return newFunc;
  }

  bool FSAtomPseudoPot::put(xmlNodePtr cur)
  {
    cur=cur->children;
    while(cur != NULL)
    {
      string cname ((const char*)cur->name);
      if(cname == "radfunc")
      {
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL)
        {
          string cname1((const char*)cur1->name);
          if(cname1=="data") 
          {
            myFunc.m_Y.resize(myFunc.grid().size());
            putContent(myFunc.m_Y,cur1);
          }
          //else if(cname1 =="grid")
          //{
          //  app_warning() << "    FSAtomPseudoPot::vps/grid is ignored " << endl;
          //}
          cur1=cur1->next;
        }
      }
      cur=cur->next;
    }
    return true;
  }
} // namespace qmcPlusPlus
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1737 $   $Date: 2007-02-12 10:42:06 -0600 (Mon, 12 Feb 2007) $
 * $Id: FSAtomPseudoPot.h 1737 2007-02-12 16:42:06Z jnkim $ 
 ***************************************************************************/
