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
    
    




#ifndef CUBICSPLINE3D_REGULARGRID_H
#define CUBICSPLINE3D_REGULARGRID_H

#include <bitset>


/*! \class template<class T> CubicSplineDiff
 * \brief spline functions on regular grid [Rmin, Rmax)
 */

template<class grid_type, class pos_type>
class CubicSpline3D
{
public:

  typedef typename grid_type::value_type value_type;
  pos_type Rmin, Rmax, dr, drinv, drinv2, dr2over6, drover6;
  std::bitset<4> xyz[8];
  grid_type Y;

  //!< default constructor
  CubicSpline3D():Y(NULL)
  {
    setxyz();
  }

  //!< default constructor
  CubicSpline3D(const CubicSpline3D<grid_type,pos_type>& a)
  {
    setxyz();
    Rmin = a.Rmin;
    Rmax = a.Rmax;
    dr = a.dr;
    drinv = a.drinv;
    drover6 = a.drover6;
    dr2over6 = a.dr2over6;
    //Y = a.Y;
  }

  void setxyz()
  {
    //4-bit to calculate the coefficient indices and sign
    //xyz[i](x,y,z,sign): initially, all the bits are set to 0
    //xyz[0] = std::bitset<4>(0,0,0,0); // -1(0) 3*0
    xyz[1].flip(2);
    xyz[1].flip(3); //xyz[1] = std::bitset<4>(0,0,1,1); // +1(1) 2*0
    xyz[2].flip(1);
    xyz[2].flip(2); //xyz[2] = std::bitset<4>(0,1,1,0); // -1(0) 1*0
    xyz[3].flip(1);
    xyz[3].flip(3); //xyz[3] = std::bitset<4>(0,1,0,1); // +1(1) 2*0
    xyz[4].flip(0);
    xyz[4].flip(3); //xyz[4] = std::bitset<4>(1,0,0,1); // +1(1) 2*0
    xyz[5].flip(0);
    xyz[5].flip(2); //xyz[5] = std::bitset<4>(1,0,1,0); // -1(0) 1*0
    xyz[6].flip();                  //xyz[6] = std::bitset<4>(1,1,1,1); // +1(1) 0*0
    xyz[7].flip(0);
    xyz[7].flip(1); //xyz[7] = std::bitset<4>(1,1,0,0); // -1(0) 1*0
  }

  value_type& operator()(int i, int j, int k)
  {
    // all shifted by 1 for the boundary
    return Y(i,j,k);
  }

  value_type operator()(int i, int j, int k)  const
  {
    // all shifted by 1 for the boundary
    return Y(i,j,k);
  }

  //!< resize the arrays
  void resize(const int l, const int m, const int n)
  {
    Y.resize(l+1,m+1,n+1);
    Rmin[0] = 0.0;
    Rmin[1] = 0.0;
    Rmin[2] = 0.0;
    Rmax[0] = 1.0;
    Rmax[1] = 1.0;
    Rmax[2] = 1.0;
    drinv[0] = static_cast<value_type>(l);
    drinv[1] = static_cast<value_type>(m);
    drinv[2] = static_cast<value_type>(n);
    drinv2[0] = drinv[0]*drinv[0];
    drinv2[1] = drinv[1]*drinv[1];
    drinv2[2] = drinv[2]*drinv[2];
    dr[0] = 1.0/drinv[0];
    dr[1] = 1.0/drinv[1];
    dr[2] = 1.0/drinv[2];
    drover6[0] = dr[0]/6.0;
    drover6[1] = dr[1]/6.0;
    drover6[2] = dr[2]/6.0;
    dr2over6[0] = dr[0]*drover6[0];
    dr2over6[1] = dr[1]*drover6[1];
    dr2over6[2] = dr[2]*drover6[2];
  }


  value_type evaluate(const pos_type& pos)
  {
    // find the location of pos
    value_type di = (pos[0]-Rmin[0])*drinv[0];
    value_type dj = (pos[1]-Rmin[1])*drinv[1];
    value_type dk = (pos[2]-Rmin[2])*drinv[2];
    // find the index of pos
    int im = static_cast<int>(di);
    int jm = static_cast<int>(dj);
    int km = static_cast<int>(dk);
    // find the index+1 of pos
    int ip = im+1;
    int jp = jm+1;
    int kp = km+1;
    // determin L_{nx}, L_{ny}, L_{nz}
    //\latex{$$L_{nx} = {R_{min}+nx*dx - x\over dx}  = -{x-R_{min}\over dx}+nx$$}
    pos_type Lm(static_cast<value_type>(im)-di,
                static_cast<value_type>(jm)-dj,
                static_cast<value_type>(km)-dk);
    // determin L_{nx+1}, L_{ny+1}, L_{nz+1}
    //\latex{$$L_{nx+1} = {R_{min}+(nx+1)*dx - x\over dx}
    // = -{x-R_{min}\over dx}+nx+1 = L_{nx}+1$$}
    pos_type Lp = Lm + 1.0;
    // determin C_{nx}, C_{ny}, C_{nz}
    pos_type Cm(dr2over6[0]*Lm[0]*(Lm[0]*Lm[0]-1),
                dr2over6[1]*Lm[1]*(Lm[1]*Lm[1]-1),
                dr2over6[2]*Lm[2]*(Lm[2]*Lm[2]-1));
    // determin C_{nx+1}, C_{ny+1}, C_{nz+1}
    pos_type Cp(dr2over6[0]*Lp[0]*(Lp[0]*Lp[0]-1),
                dr2over6[1]*Lp[1]*(Lp[1]*Lp[1]-1),
                dr2over6[2]*Lp[2]*(Lp[2]*Lp[2]-1));
    int i, j, k;
    value_type Lx, Ly, Lz;
    value_type Cx, Cy, Cz;
    value_type y, y2x, y2y, y2z;
    value_type res = 0.0;
    // precalculate d^2 phi/dx^2, d^2 phi/dy^2 and d^2 phi/dz^2
    for(int ixyz=0; ixyz<8; ixyz++)
    {
      // get Lx, Cx, Cvx, and index i to get the value
      if(xyz[ixyz][0])
        // (1,y,z)
      {
        Lx = Lp[0];
        Cx = Cp[0];
        i = im;
      }
      else
        // (0,y,z)
      {
        Lx = Lm[0];
        Cx = Cm[0];
        i = ip;
      }
      // get Ly, Cy, Cvy, and index j to get the value
      if(xyz[ixyz][1])
        // (x,1,z)
      {
        Ly = Lp[1];
        Cy = Cp[1];
        j = jm;
      }
      else
        // (x,0,z)
      {
        Ly = Lm[1];
        Cy = Cm[1];
        j = jp;
      }
      // get Lz, Cz, Cvz, and index k to get the value
      if(xyz[ixyz][2])
        // (x,y,1)
      {
        Lz = Lp[2];
        Cz = Cp[2];
        k = km;
      }
      else
        // (x,y,0)
      {
        Lz = Lm[2];
        Cz = Cm[2];
        k = kp;
      }
      y = Y(i,j,k);
      //\warning: not protected for negative indices
      y2x = (Y(i+1,j,k)+Y(i-1,j,k)-2.0*y)*drinv[0]*drinv[0];
      y2y = (Y(i,j+1,k)+Y(i,j-1,k)-2.0*y)*drinv[1]*drinv[1];
      y2z = (Y(i,j,k+1)+Y(i,j,k-1)-2.0*y)*drinv[2]*drinv[2];
      if(xyz[ixyz][3])
        //positive
      {
        res += (Lx*Ly*Lz*y +Cx*Ly*Lz*y2x +Lx*Cy*Lz*y2y+Lx*Ly*Cz*y2z);
      }
      else
      {
        res -= (Lx*Ly*Lz*y +Cx*Ly*Lz*y2x +Lx*Cy*Lz*y2y+Lx*Ly*Cz*y2z);
      }
    }
    return res;
  }

  value_type laplacian(const pos_type& pos, value_type& value, pos_type& gr)
  {
    // find the location of pos
    value_type di = (pos[0]-Rmin[0])*drinv[0];
    value_type dj = (pos[1]-Rmin[1])*drinv[1];
    value_type dk = (pos[2]-Rmin[2])*drinv[2];
    // find the index of pos
    int im = static_cast<int>(di);
    int jm = static_cast<int>(dj);
    int km = static_cast<int>(dk);
    // find the index+1 of pos
    int ip = im+1;
    int jp = jm+1;
    int kp = km+1;
    // determin L_{nx}, L_{ny}, L_{nz}
    //\latex{$$L_{nx} = {R_{min}+nx*dx - x\over dx}  = -{x-R_{min}\over dx}+nx$$}
    pos_type Lm(static_cast<value_type>(im)-di,
                static_cast<value_type>(jm)-dj,
                static_cast<value_type>(km)-dk);
    // determin L_{nx+1}, L_{ny+1}, L_{nz+1}
    //\latex{$$L_{nx+1} = {R_{min}+(nx+1)*dx - x\over dx}
    //= -{x-R_{min}\over dx}+nx+1 = L_{nx}+1$$}
    pos_type Lp = Lm + 1.0;
    // determin C_{nx}, C_{ny}, C_{nz}
    pos_type Cm(dr2over6[0]*Lm[0]*(Lm[0]*Lm[0]-1),
                dr2over6[1]*Lm[1]*(Lm[1]*Lm[1]-1),
                dr2over6[2]*Lm[2]*(Lm[2]*Lm[2]-1));
    // determin C_{nx+1}, C_{ny+1}, C_{nz+1}
    pos_type Cp(dr2over6[0]*Lp[0]*(Lp[0]*Lp[0]-1),
                dr2over6[1]*Lp[1]*(Lp[1]*Lp[1]-1),
                dr2over6[2]*Lp[2]*(Lp[2]*Lp[2]-1));
    // determine -dC_{nx}/dx, -dC_{ny}/dy, -dC_{nz}/dz
    pos_type Cvm(drover6[0]*(3*Lm[0]*Lm[0]-1),
                 drover6[1]*(3*Lm[1]*Lm[1]-1),
                 drover6[2]*(3*Lm[2]*Lm[2]-1));
    // determin -dC_{nx+1}/dx, -dC_{ny+1}/dy, -dC_{nz+1}/dz
    pos_type Cvp(drover6[0]*(3*Lp[0]*Lp[0]-1),
                 drover6[1]*(3*Lp[1]*Lp[1]-1),
                 drover6[2]*(3*Lp[2]*Lp[2]-1));
    int i,j,k;
    value_type Lx, Ly, Lz;
    value_type Cx, Cy, Cz;
    value_type Cvx, Cvy, Cvz;
    value_type lap = 0.0;
    value = 0;
    // precalculate d^2 phi/dx^2, d^2 phi/dy^2 and d^2 phi/dz^2
    for(int ixyz=0; ixyz<8; ixyz++)
    {
      // get Lx, Cx, Cvx, and index i to get the value
      if(xyz[ixyz][0])
        // (1,y,z)
      {
        Lx = Lp[0];
        Cx = Cp[0];
        Cvx = Cvp[0];
        i = im;
      }
      else
        // (0,y,z)
      {
        Lx = Lm[0];
        Cx = Cm[0];
        Cvx = Cvm[0];
        i = ip;
      }
      // get Ly, Cy, Cvy, and index j to get the value
      if(xyz[ixyz][1])
        // (x,1,z)
      {
        Ly = Lp[1];
        Cy = Cp[1];
        Cvy = Cvp[1];
        j = jm;
      }
      else
        // (x,0,z)
      {
        Ly = Lm[1];
        Cy = Cm[1];
        Cvy = Cvm[1];
        j = jp;
      }
      // get Lz, Cz, Cvz, and index k to get the value
      if(xyz[ixyz][2])
        // (x,y,1)
      {
        Lz = Lp[2];
        Cz = Cp[2];
        Cvz = Cvp[2];
        k = km;
      }
      else
        // (x,y,0)
      {
        Lz = Lm[2];
        Cz = Cm[2];
        Cvz = Cvm[2];
        k = kp;
      }
      value_type y = Y(i,j,k);
      //\warning: not protected for out-of-bound indices
      value_type y2x = (Y(i+1,j,  k)  +Y(i-1,j,  k)  -2.0*y)*drinv2[0];
      value_type y2y = (Y(i,  j+1,k)  +Y(i,  j-1,k)  -2.0*y)*drinv2[1];
      value_type y2z = (Y(i,  j,  k+1)+Y(i,  j,  k-1)-2.0*y)*drinv2[2];
      if(xyz[ixyz][3])
        //positive
      {
        value += (Lx*Ly*Lz*y + Cx*Ly*Lz*y2x + Lx*Cy*Lz*y2y + Lx*Ly*Cz*y2z);
        lap   += Lx*Ly*Lz*(y2x+y2y+y2z);
        // negative signs are taken care of here
        // Cvx = -dC/dx, Cvy = -dC/dy, and Cvz = -dC/dz
        // dL/dx = -1/dx
        gr[0] -= (drinv[0]*(Ly*Lz*y+Cy*Lz*y2y+Ly*Cz*y2z) + Cvx*Ly *Lz *y2x);
        gr[1] -= (drinv[1]*(Lz*Lx*y+Cz*Lx*y2z+Lz*Cx*y2x) + Lx *Cvy*Lz *y2y);
        gr[2] -= (drinv[2]*(Lx*Ly*y+Cx*Ly*y2x+Lx*Cy*y2y) + Lx *Ly *Cvz*y2z);
      }
      else
      {
        value -= (Lx*Ly*Lz*y +Cx*Ly*Lz*y2x + Lx*Cy*Lz*y2y + Lx*Ly*Cz*y2z);
        lap   -= Lx*Ly*Lz*(y2x+y2y+y2z);
        // negative signs are taken care of here
        gr[0] += (drinv[0]*(Ly*Lz*y +Cy*Lz*y2y + Ly*Cz*y2z) + Cvx*Ly *Lz *y2x);
        gr[1] += (drinv[1]*(Lz*Lx*y +Cz*Lx*y2z + Lz*Cx*y2x) + Lx *Cvy*Lz *y2y);
        gr[2] += (drinv[2]*(Lx*Ly*y +Cx*Ly*y2x + Lx*Cy*y2y) + Lx *Ly *Cvz*y2z);
      }
    }
    return lap;
  }

};

#endif
