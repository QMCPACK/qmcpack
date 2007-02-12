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
#ifndef QMCPLUSPLUS_FSATOMPSEDUOPOTENTIAL_H 
#define QMCPLUSPLUS_FSATOMPSEDUOPOTENTIAL_H 
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimLinearSpline.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus {

  template<typename T>
    struct FSAtomPseudoPot//: public OneDimCubicSpline<T>
    {
      typedef OneDimLinearSpline<T> return_type;
      typedef typename OneDimLinearSpline<T>::grid_type grid_type;
      return_type myFunc;
      int AngL;
      T Rcut;
      
      FSAtomPseudoPot(int l, T rc, grid_type* agrid): 
        myFunc(agrid), AngL(l), Rcut(rc)
      {
      }

      ~FSAtomPseudoPot()
      {
      }

      void convert2RV()
      {
        for(int i=0; i<myFunc.size(); i++) myFunc(i) *= myFunc.r(i);
      }

      T getCutOff(T v0)
      {
        bool ignore=true;
        int ng=myFunc.size()-2;
        T r=myFunc.r(ng);
        while(ignore&& ng)
        {
          cout << r << " " << myFunc(ng) << endl;
          ignore=(std::abs(myFunc(ng)-v0)<1e-12);
          r=myFunc.r(--ng);
        }
        return r;
      }

      /** create a LinearSpline<T>
       * @param sc scaling factor
       * @return a LinearSpline<T> on a LinearGrid
       */
      return_type* getLocalPot(T sc)
      {
        myFunc.spline();

        const T del=1.0e-3;
        int ng=static_cast<int>(2*Rcut/del)+1;
        app_log() << "  FSAtomPseudoPot::getLocalPot grid = [0," << 2*Rcut << "] npts=" << ng << endl;

        LinearGrid<T> *agrid=new LinearGrid<T>;
        agrid->set(0.0,2*Rcut,ng);

        return_type* newFunc=new return_type(agrid);
        (*newFunc)(0) = 0.0;
        for(int i=1; i<ng-1; i++)
        {
          (*newFunc)(i)=sc*myFunc.splintNG((*agrid)[i]);
        }
        (*newFunc)(ng-1)=1.0;

        T y_prime=((*newFunc)(1)-(*newFunc)(0))/del;
        newFunc->spline(0,y_prime,ng-1,0.0);
        return newFunc;
      }

      return_type* getNonLocalPot(FSAtomPseudoPot<T>& vloc, T vFac)
      {
        myFunc.m_Y -= vloc.myFunc.m_Y;
        //T y_prime=(m_Y[1]-m_Y[0])/m_grid->dr(0);
        T y_prime=(myFunc(1)-myFunc(0))/myFunc.dr(0);
        myFunc.spline(0,y_prime,myFunc.size()-1,0.0);

        const T del=1.0e-3;
        int ng=static_cast<int>(Rcut/del)+1;
        LinearGrid<T> *agrid=new LinearGrid<T>;
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

      bool put(xmlNodePtr cur)
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

    };
} // namespace qmcPlusPlus
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
