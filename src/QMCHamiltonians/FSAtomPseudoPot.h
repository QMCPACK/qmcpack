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
#include "Numerics/OneDimCubicSpline.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus {

  template<typename T>
    struct FSAtomPseudoPot:public OneDimCubicSpline<T>
    {
      typedef OneDimCubicSpline<T> mybase_type;
      typedef OneDimLinearSpline<T> return_type;
      typedef typename mybase_type::value_type value_type;
      typedef typename mybase_type::point_type point_type;
      typedef typename mybase_type::data_type  data_type;
      typedef typename mybase_type::grid_type  grid_type;

      using mybase_type::m_grid;
      using mybase_type::m_Y;

      bool OwnGrid;
      int AngL;
      
      FSAtomPseudoPot(): 
        OneDimCubicSpline<T>(0), OwnGrid(false), AngL(-1)
      {
      }

      ~FSAtomPseudoPot()
      {
        if(OwnGrid) delete m_grid;
      }

      T getCutOff(T v0)
      {
        bool ignore=true;
        int ng=m_grid->size()-2;
        T r=(*m_grid)(ng);
        while(ignore&& ng)
        {
          ignore=abs(m_Y[ng]-v0)<1e-12;
          r=(*m_grid)(--ng);
        }
        return r;
      }

      /** create a LinearSpline<T>
       * @param sc scaling factor
       * @return a LinearSpline<T> on a LinearGrid
       */
      return_type* getLocalPot(T sc)
      {
        m_Y *= sc;
        this->spline();

        T rcut=getCutOff(1.0);
        int ng=static_cast<int>(rcut*100)+1;

        cout << "  FSAtomPseudoPot::getLocalPot grid = [0," 
          << rcut << "] npts=" << ng << endl;

        LinearGrid<T> *agrid=new LinearGrid<T>;
        agrid->set(0.0,rcut,ng);

        return_type* newFunc=new return_type(agrid);
        for(int i=1; i<ng-1; i++)
        {
          (*newFunc)(i)=this->splint((*agrid)[i]);
        }
        (*newFunc)(0) =  2.0*(*newFunc)(1) - (*newFunc)(2);
        (*newFunc)(ng-1)=1.0;
        newFunc->spline();
        return newFunc;
      }

      return_type* getNonLocalPot(FSAtomPseudoPot<T>& vloc, T vFac)
      {
        m_Y -= vloc.m_Y;
        m_Y *= vFac;
        this->spline();

        T rcut=getCutOff(0.0);
        int ng=static_cast<int>(m_grid->rmax()*100)+1;

        cout << "  FSAtomPseudoPot::getNonLocalPot grid = [0," 
          << rcut << "] npts=" << ng << endl;
        LinearGrid<T> *agrid=new LinearGrid<T>;
        agrid->set(0.0,m_grid->rmax(),ng);

        return_type* newFunc=new return_type(agrid);
        for(int i=1; i<ng-1; i++)
        {
          (*newFunc)(i)=this->splint((*agrid)[i])/(*agrid)[i];
        }

        (*newFunc)(0) =  2.0*(*newFunc)(1) - (*newFunc)(2);
        (*newFunc)(ng-1)=0.0;
        newFunc->spline();
        return newFunc;
      }

      bool put(xmlNodePtr cur)
      {
        T rc=-1.0;
        OhmmsAttributeSet aAttrib;
        aAttrib.add(AngL,"l");
        aAttrib.add(rc,"cutoff");
        aAttrib.put(cur);

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
              if(cname1 =="grid")
                putGrid(cur1);
              else if(cname1=="data")
                putContent(m_Y,cur1);
              cur1=cur1->next;
            }
          }
          cur=cur->next;
        }
        return true;
      }

      void putGrid(xmlNodePtr cur)
      {
        if(m_grid == 0)
        {
          T ri = 1e-5;
          T rf = 100.0;
          T ascale = -1.0e0;
          T astep = 1.25e-2;
          int npts = 1001;

          string gridType("log");
          OhmmsAttributeSet radAttrib;
          radAttrib.add(gridType,"type"); 
          radAttrib.add(npts,"npts"); 
          radAttrib.add(ri,"ri"); radAttrib.add(rf,"rf");
          radAttrib.add(ascale,"ascale"); radAttrib.add(astep,"astep");
          radAttrib.add(ascale,"scale"); radAttrib.add(astep,"step");
          radAttrib.put(cur);

          if(gridType == "log") {
            if(ascale>0.0) {
              m_grid = new LogGridZero<T>;
              m_grid->set(astep,ascale,npts);
            } else {
              if(ri<numeric_limits<T>::epsilon())
              {
                ri=numeric_limits<T>::epsilon();
              }
              m_grid = new LogGrid<T>;
              m_grid->set(ri,rf,npts);
            }
          } else if(gridType == "linear") {
            m_grid = new LinearGrid<T>;
            m_grid->set(ri,rf,npts);
          }
        }
        this->resize(m_grid->size());
      }
    };
} // namespace qmcPlusPlus
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
