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
    
    



#ifndef GUARD_GRID3D_H
#define GUARD_GRID3D_H

#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid1D.h"
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

class Grid3D
{

public:

  /// number of points along each direction
  int n_x, n_y, n_z;

  /// total number of points in the Grid3D : 0 to m_size-1
  int m_size;

  /// array of the axes Grid1D
  std::vector<Grid1D> m_axis;

  /// Constructor :: create the three axes
  Grid3D()
  {
    m_axis.resize(3);
  }

  /// initialise an axis with entire data from the input
  void init(int idir,
            int nsections,
            int npts,
            double xi,
            const std::vector<int>& nrho,
            const std::vector<double>& dh)
  {
    m_axis[idir].init(nsections,npts,xi,nrho,dh);
    return;
  }


  ///  initialise the Grid3D from another Grid3D
  void init(const gridvec_t& ri,
            const gridvec_t& rf,
            const Grid3D* agrid)
  {
    for(int idir = 0; idir < 3; idir++)
      m_axis[idir].init(ri[idir],rf[idir],agrid->m_axis[idir]);
    set_size();
    return;
  }

  void set_size()
  {
    n_x = m_axis[0].Size();
    n_y = m_axis[1].Size();
    n_z = m_axis[2].Size();
    m_size = n_x * n_y * n_z;
    return;
  }



  /// index of a Grid3D point
  inline int index(const gridvec_t& ir)
  {
    return ir[0] + n_x * ( ir[1] + n_y * ir[2] );
  }

  /// index of a Grid3D point
  inline int index(int ix, int iy, int iz)
  {
    return ix + n_x * ( iy + n_y * iz );
  }

  /// Grid3D point corresponding to index
  inline gridvec_t ipt(int index)
  {
    gridvec_t ir;
    int nxy = n_x * n_y;
    ir[2] = index/nxy;
    index -= nxy*ir[2];
    ir[1] = index/n_x;
    index -= n_x*ir[1];
    ir[0] = index;
    return ir;
  }

  /// lowest Grid3D point
  inline gridvec_t ptl(const posvec_t& r)
  {
    gridvec_t ir;
    for(int i = 0; i < 3; i++)
      ir[i] = m_axis[i].xl(r[i]);
    return ir;
  }

  /// nearest Grid3D point
  inline gridvec_t ptn(const posvec_t& r)
  {
    gridvec_t ir;
    for(int i = 0; i < 3; i++)
      ir[i] = m_axis[i].xn(r[i]);
    return ir;
  }

  /// r coordinate of a Grid3D point
  inline posvec_t ptr(const gridvec_t& ir)
  {
    posvec_t r;
    for(int i = 0; i < 3; i++)
      r[i] = m_axis[i].m_coord[ir[i]];
    return r;
  }

  inline posvec_t ptr(int ix, int iy, int iz)
  {
    posvec_t r;
    r[0] = m_axis[0].m_coord[ix];
    r[1] = m_axis[1].m_coord[iy];
    r[2] = m_axis[2].m_coord[iz];
    return r;
  }

  double put(xmlNodePtr cur, double epsbym)
  {
    double units;                        /// conversion unit
    const double aB = 0.5291772108e-10;  /// Bohr radius in meters
    const double aBeff = epsbym * aB;
    const double u0 = 1.0/aBeff;
    /// read in the units
    std::string cname ((const char*)(cur->name));
    if( cname == "Grid3D")
    {
      std::string cunit = (char*)xmlGetProp(cur,(xmlChar*)"unit");
      if(cunit == "nm")
      {
        units = u0 * 1.0e-9;
        std::cout << "The input data are expressed in nano-meters. " << std::endl;
        std::cout << "Using a.u. :: 1nm = " << units << " a_B*, and 1 a_B* =  "
             << aBeff << " m" << std::endl;
      }
      else
        if(cunit == "A")
        {
          units = u0 * 1.0e-10;
          std::cout << "The input data are expressed in Angstroms. " << std::endl;
          std::cout << "Using a.u. :: 1A = " << units << " a_B* " << std::endl;
        }
        else
          if(cunit == "au")
          {
            units = u0 * aB;
            std::cout << "The input data are expressed in atomic units. " << std::endl;
            std::cout << "Using a.u. :: 1 a_B = " << units << " a_B* "  << std::endl;
          }
          else
          {
            units = 1.0;
            std::cout << "The input data are expressed in effective au. " << std::endl;
            std::cout << "Using a.u. :: 1nm = 1 a_B = 5.291772083(19)e-11m" << std::endl;
          }
      std::vector<int> n_rho;
      std::vector<double> d_h;
      xmlNodePtr node1 = cur->xmlChildrenNode;     /// node1->name = Grid1D
      while( node1 != NULL )
      {
        std::string name1 ((const char*)(node1->name));
        if( name1 == "Grid1D" )
        {
          /// axis number
          int dir = atoi((char*)xmlGetProp(node1,(xmlChar*)"dir"));
          /// number of sections
          int nsecs = atoi((char*)xmlGetProp(node1,(xmlChar*)"nsecs"));
          /// total number of points
          int npts = atoi((char*)xmlGetProp(node1,(xmlChar*)"npts"));
          /// coordinate origin
          double x0 = units * atof((char*)xmlGetProp(node1,(xmlChar*)"x0"));
          /// assign the Grid::initialisation-vectors
          n_rho.resize(nsecs);
          d_h.resize(nsecs);
          /// initilaise section to be read in.
          int isec = 0;
          xmlNodePtr node2 = node1->xmlChildrenNode;
          while( node2 != NULL )
          {
            std::string name2 ((const char*)(node2->name));
            if( name2 == "section" )
            {
              n_rho[isec]=atoi((char*)xmlGetProp(node2,(xmlChar*)"nx"));
              d_h[isec]=units*atof((char*)xmlGetProp(node2,(xmlChar*)"dx"));
              isec++;
            }
            node2 = node2->next;   /// next section if any
          }
          /// initialise axis
          m_axis[dir].init(nsecs,npts,x0,n_rho,d_h);
        }
        node1 = node1->next;      /// next Grid1D if any
      }
    }
    else
    {
      std::cout << "Grid3D not properly supplied. Exiting ... " << std::endl;
      exit(-1);
    }
    set_size();     /// set all dimension sizes
    return units;
  }


};
#endif
