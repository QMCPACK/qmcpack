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
    
    




#ifndef CUBICSPLINE3DCOLLECTION_REGULARGRID_H
#define CUBICSPLINE3DCOLLECTION_REGULARGRID_H

#include <bitset>

/*! \class template<class T> CubicSpline3D
 * \brief spline functions on 3D regular grid [Rmin, Rmax)^3
 * The number of grids per axis can vary.
 * A set of orbitals are associated with a CubicSpline3D.
 * This class cannot destroy the orbitals. They have to be provided externally.
 */

template<class grid_type, class pos_type>
struct CubicSpline3D
{

  typedef typename grid_type::value_type value_type;

  int Npt[3];
  int OrbitalID;
  pos_type Rmin, Rmax, dr, drinv, drinv2, dr2over6, drover6;
  std::bitset<4> xyz[8];     // \note possible to make it "static const"
  pos_type Lm, Lp, Cm, Cp, Cvm, Cvp;
  int I[3], J[3];

  std::vector<grid_type*> Orb; // a set of orbitals
  //!< default constructor
  CubicSpline3D(): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    // default bounds [0,1)x[0,1)x[0,1) with 1x1x1 grid
    Npt[0] = 1;
    Npt[1] = 1;
    Npt[2] = 1;
    Rmin[0] = 0.0;
    Rmin[1] = 0.0;
    Rmin[2] = 0.0;
    Rmax[0] = 1.0;
    Rmax[1] = 1.0;
    Rmax[2] = 1.0;
    resize(1,1,1); // set the default variables
  }

  //!< constructor using grid
  CubicSpline3D(int l, int m, int n): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    Rmin[0] = 0.0;
    Rmin[1] = 0.0;
    Rmin[2] = 0.0;
    Rmax[0] = 1.0;
    Rmax[1] = 1.0;
    Rmax[2] = 1.0;
    resize(l,m,n);
  }

  //!< copy constructor
  CubicSpline3D(const CubicSpline3D<grid_type,pos_type>& a): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    Rmin = a.Rmin;
    Rmax = a.Rmax;
    resize(a.Npt[0], a.Npt[1], a.Npt[2]);
    Orb.reserve(a.Orb.size());
    for(int iorb=0; iorb< a.Orb.size(); iorb++)
      Orb.push_back(a.Orb[iorb]);
  }

  ~CubicSpline3D() { }

  void setxyz();

  // assign value at (i,j,k) of iorb-th orbital
  value_type& operator()(int iorb, int i, int j, int k)
  {
    return Orb[iorb]->operator()(i,j,k);
  }

  // return value at (i,j,k) of iorb-th orbital
  value_type operator()(int iorb, int i, int j, int k)  const
  {
    return Orb[iorb]->operator()(i,j,k);
  }


  //!\fn   int add(grid_type* orbital)
  //!\param orbital a grid_type to be added to Orb
  //!\return the index of the added orbital
  //!\warning Assuming that the new orbital is the same kind, i.e., size and grid
  int add(grid_type* orbital)
  {
    int id = Orb.size();
    Orb.push_back(orbital);
    return id;
  }

  //!< resize the arrays
  void resize(const int l, const int m, const int n);

  //!< function to calculate value, gradient and laplacian
  //\return value[iorb] = value of iorb-th orbial
  //        gr[iorb] = gradient of iorb-th orbial
  //        lap[iorb] = laplacian of iorb-th orbial
  template<class ValArray, class PosArray>
  void laplacian(const pos_type& pos, ValArray& value, PosArray& gr, ValArray& lap);

  inline void get_grid(const pos_type& pos)
  {
    // find the location of pos
    value_type di = (pos[0]-Rmin[0])*drinv[0];
    value_type dj = (pos[1]-Rmin[1])*drinv[1];
    value_type dk = (pos[2]-Rmin[2])*drinv[2];
    // find the index of pos
    I[0] = static_cast<int>(di);
    I[1] = static_cast<int>(dj);
    I[2] = static_cast<int>(dk);
    // find the index+1 of pos
    J[0] = I[0]+1;
    J[1] = I[1]+1;
    J[2] = I[2]+1;
    // determin L_{nx}, L_{ny}, L_{nz}
    //\latex{$$L_{nx} = {R_{min}+nx*dx - x\over dx}  = -{x-R_{min}\over dx}+nx$$}
    Lm[0] = static_cast<value_type>(I[0])-di;
    Lm[1] = static_cast<value_type>(I[1])-dj;
    Lm[2] = static_cast<value_type>(I[2])-dk;
    // determin L_{nx+1}, L_{ny+1}, L_{nz+1}
    //\latex{$$L_{nx+1} = {R_{min}+(nx+1)*dx - x\over dx}
    //= -{x-R_{min}\over dx}+nx+1 = L_{nx}+1$$}
    Lp = Lm + 1.0;
    // determin C_{nx}, C_{ny}, C_{nz}
    Cm[0] = dr2over6[0]*Lm[0]*(Lm[0]*Lm[0]-1);
    Cm[1] = dr2over6[1]*Lm[1]*(Lm[1]*Lm[1]-1);
    Cm[2] = dr2over6[2]*Lm[2]*(Lm[2]*Lm[2]-1);
    // determin C_{nx+1}, C_{ny+1}, C_{nz+1}
    Cp[0] = dr2over6[0]*Lp[0]*(Lp[0]*Lp[0]-1);
    Cp[1] = dr2over6[1]*Lp[1]*(Lp[1]*Lp[1]-1);
    Cp[2] = dr2over6[2]*Lp[2]*(Lp[2]*Lp[2]-1);
    // determine -dC_{nx}/dx, -dC_{ny}/dy, -dC_{nz}/dz
    Cvm[0] = drover6[0]*(3*Lm[0]*Lm[0]-1);
    Cvm[1] = drover6[1]*(3*Lm[1]*Lm[1]-1);
    Cvm[2] = drover6[2]*(3*Lm[2]*Lm[2]-1);
    // determin -dC_{nx+1}/dx, -dC_{ny+1}/dy, -dC_{nz+1}/dz
    Cvp[0] = drover6[0]*(3*Lp[0]*Lp[0]-1);
    Cvp[1] = drover6[1]*(3*Lp[1]*Lp[1]-1);
    Cvp[2] = drover6[2]*(3*Lp[2]*Lp[2]-1);
  }
};


template<class grid_type, class pos_type>
void
CubicSpline3D<grid_type,pos_type>::resize(const int l, const int m, const int n)
{
  Npt[0] = l;
  Npt[1] = m;
  Npt[2] = n;
  dr[0] = (Rmax[0]-Rmin[0])/static_cast<value_type>(l);
  dr[1] = (Rmax[1]-Rmin[1])/static_cast<value_type>(m);
  dr[2] = (Rmax[2]-Rmin[2])/static_cast<value_type>(n);
  drinv[0] = 1.0/dr[0];
  drinv[1] = 1.0/dr[1];
  drinv[2] = 1.0/dr[2];
  drinv2[0] = drinv[0]*drinv[0];
  drinv2[1] = drinv[1]*drinv[1];
  drinv2[2] = drinv[2]*drinv[2];
  drover6[0] = dr[0]/6.0;
  drover6[1] = dr[1]/6.0;
  drover6[2] = dr[2]/6.0;
  dr2over6[0] = dr[0]*drover6[0];
  dr2over6[1] = dr[1]*drover6[1];
  dr2over6[2] = dr[2]*drover6[2];
}


template<class grid_type, class pos_type>
template<class ValArray, class PosArray>
void
CubicSpline3D<grid_type,pos_type>::laplacian(const pos_type& pos,
    ValArray& value, PosArray& gr,
    ValArray& lap)
{
  get_grid(pos);
//   // find the location of pos
//   value_type di = (pos[0]-Rmin[0])*drinv[0];
//   value_type dj = (pos[1]-Rmin[1])*drinv[1];
//   value_type dk = (pos[2]-Rmin[2])*drinv[2];
//   // find the index of pos
//   int im = static_cast<int>(di);
//   int jm = static_cast<int>(dj);
//   int km = static_cast<int>(dk);
//   // find the index+1 of pos
//   int ip = im+1;
//   int jp = jm+1;
//   int kp = km+1;
//   // determin L_{nx}, L_{ny}, L_{nz}
//   //\latex{$$L_{nx} = {R_{min}+nx*dx - x\over dx}  = -{x-R_{min}\over dx}+nx$$}
//   pos_type Lm(static_cast<value_type>(im)-di,
// 	      static_cast<value_type>(jm)-dj,
// 	      static_cast<value_type>(km)-dk);
//   // determin L_{nx+1}, L_{ny+1}, L_{nz+1}
//   //\latex{$$L_{nx+1} = {R_{min}+(nx+1)*dx - x\over dx}
//   //= -{x-R_{min}\over dx}+nx+1 = L_{nx}+1$$}
//   pos_type Lp = Lm + 1.0;
//   // determin C_{nx}, C_{ny}, C_{nz}
//   pos_type Cm(dr2over6[0]*Lm[0]*(Lm[0]*Lm[0]-1),
// 	      dr2over6[1]*Lm[1]*(Lm[1]*Lm[1]-1),
// 	      dr2over6[2]*Lm[2]*(Lm[2]*Lm[2]-1));
//   // determin C_{nx+1}, C_{ny+1}, C_{nz+1}
//   pos_type Cp(dr2over6[0]*Lp[0]*(Lp[0]*Lp[0]-1),
// 	      dr2over6[1]*Lp[1]*(Lp[1]*Lp[1]-1),
// 	      dr2over6[2]*Lp[2]*(Lp[2]*Lp[2]-1));
//   // determine -dC_{nx}/dx, -dC_{ny}/dy, -dC_{nz}/dz
//   pos_type Cvm(drover6[0]*(3*Lm[0]*Lm[0]-1),
// 	       drover6[1]*(3*Lm[1]*Lm[1]-1),
// 	       drover6[2]*(3*Lm[2]*Lm[2]-1));
//   // determin -dC_{nx+1}/dx, -dC_{ny+1}/dy, -dC_{nz+1}/dz
//   pos_type Cvp(drover6[0]*(3*Lp[0]*Lp[0]-1),
// 	       drover6[1]*(3*Lp[1]*Lp[1]-1),
// 	       drover6[2]*(3*Lp[2]*Lp[2]-1));
  int i,j,k;
  value_type Lx, Ly, Lz;
  value_type Cx, Cy, Cz;
  value_type Cvx, Cvy, Cvz;
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
      i = I[0];
    }
    else
      // (0,y,z)
    {
      Lx = Lm[0];
      Cx = Cm[0];
      Cvx = Cvm[0];
      i = J[0];
    }
    // get Ly, Cy, Cvy, and index j to get the value
    if(xyz[ixyz][1])
      // (x,1,z)
    {
      Ly = Lp[1];
      Cy = Cp[1];
      Cvy = Cvp[1];
      j = I[1];
    }
    else
      // (x,0,z)
    {
      Ly = Lm[1];
      Cy = Cm[1];
      Cvy = Cvm[1];
      j = J[1];
    }
    // get Lz, Cz, Cvz, and index k to get the value
    if(xyz[ixyz][2])
      // (x,y,1)
    {
      Lz = Lp[2];
      Cz = Cp[2];
      Cvz = Cvp[2];
      k = I[2];
    }
    else
      // (x,y,0)
    {
      Lz = Lm[2];
      Cz = Cm[2];
      Cvz = Cvm[2];
      k = J[2];
    }
    for(int iorb=0; iorb<Orb.size(); iorb++)
    {
      value_type y = (*Orb[iorb])(i,j,k);
      //\warning: not protected for negative indices
      value_type y2x = ((*Orb[iorb])(i+1,j,k)+(*Orb[iorb])(i-1,j,k)-2.0*y)*drinv2[0];
      value_type y2y = ((*Orb[iorb])(i,j+1,k)+(*Orb[iorb])(i,j-1,k)-2.0*y)*drinv2[1];
      value_type y2z = ((*Orb[iorb])(i,j,k+1)+(*Orb[iorb])(i,j,k-1)-2.0*y)*drinv2[2];
      if(xyz[ixyz][3])
        //positive
      {
        value[iorb] += (Lx*Ly*Lz*y + Cx*Ly*Lz*y2x + Lx*Cy*Lz*y2y + Lx*Ly*Cz*y2z);
        lap[iorb]   += Lx*Ly*Lz*(y2x+y2y+y2z);
        // negative signs are taken care of here
        // Cvx = -dC/dx, Cvy = -dC/dy, and Cvz = -dC/dz
        // dL/dx = -1/dx
        gr[iorb][0] -= (drinv[0]*(Ly*Lz*y+Cy*Lz*y2y+Ly*Cz*y2z) + Cvx*Ly *Lz *y2x);
        gr[iorb][1] -= (drinv[1]*(Lz*Lx*y+Cz*Lx*y2z+Lz*Cx*y2x) + Lx *Cvy*Lz *y2y);
        gr[iorb][2] -= (drinv[2]*(Lx*Ly*y+Cx*Ly*y2x+Lx*Cy*y2y) + Lx *Ly *Cvz*y2z);
      }
      else
      {
        value[iorb] -= (Lx*Ly*Lz*y +Cx*Ly*Lz*y2x + Lx*Cy*Lz*y2y + Lx*Ly*Cz*y2z);
        lap[iorb]   -= Lx*Ly*Lz*(y2x+y2y+y2z);
        // negative signs are taken care of here
        gr[iorb][0] += (drinv[0]*(Ly*Lz*y +Cy*Lz*y2y + Ly*Cz*y2z) + Cvx*Ly *Lz *y2x);
        gr[iorb][1] += (drinv[1]*(Lz*Lx*y +Cz*Lx*y2z + Lz*Cx*y2x) + Lx *Cvy*Lz *y2y);
        gr[iorb][2] += (drinv[2]*(Lx*Ly*y +Cx*Ly*y2x + Lx*Cy*y2y) + Lx *Ly *Cvz*y2z);
      }
    }
  }
}

template<class grid_type, class pos_type>
void
CubicSpline3D<grid_type,pos_type>::setxyz()
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

#endif
