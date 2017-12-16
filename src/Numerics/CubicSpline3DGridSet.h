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

template<class PosA>
struct CubicSpline3DGrid
{

  typedef PosA pos_type;
  typedef typename PosA::Type_t value_type;

  int Npt[3];
  value_type rcut2;
  pos_type CurrentPos;
  pos_type Rmin, Rmax, dr, drinv, drinv2, dr2over6, drover6;
  pos_type Lm, Lp, Cm, Cp, Cvm, Cvp;
  int I[3], Ip[3];

  CubicSpline3DGrid()
  {
    Rmin[0] = 0.0;
    Rmin[1] = 0.0;
    Rmin[2] = 0.0;
    Rmax[0] = 1.0;
    Rmax[1] = 1.0;
    Rmax[2] = 1.0;
    resize(1,1,1);
  }

  CubicSpline3DGrid(const CubicSpline3DGrid<PosA>& g)
  {
    Rmin = g.Rmin;
    Rmax = g.Rmax;
    resize(g.Npt[0], g.Npt[1], g.Npt[2]);
  }

  CubicSpline3DGrid<PosA>& operator=(const CubicSpline3DGrid<PosA>& g)
  {
    Rmin = g.Rmin;
    Rmax = g.Rmax;
    resize(g.Npt[0], g.Npt[1], g.Npt[2]);
    return *this;
  }

  // return 1/dr^2 for the 2nd-derivatives
  inline value_type c2(int i) const
  {
    return drinv2[i];
  }

  // return 1/dr for the 2nd-derivatives
  inline value_type c1(int i) const
  {
    return drinv[i];
  }

  void resize(int l, int m, int n)
  {
    Npt[0] = l;
    Npt[1] =m;
    Npt[2] = n;
    dr[0] = (Rmax[0]-Rmin[0])/static_cast<value_type>(l);
    dr[1] = (Rmax[1]-Rmin[1])/static_cast<value_type>(m);
    dr[2] = (Rmax[2]-Rmin[2])/static_cast<value_type>(n);
    // minimim distance 1/100 th of the distance, to be over written
    rcut2 = dot(dr,dr)/1000000;
    CurrentPos = PosA(1000, 1000, 1000);
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


  inline bool get(const pos_type& pos)
  {
    //how much the gain will be???
    //if(dot(CurrentPos,pos) > rcut2) {
    if(std::abs(CurrentPos[0]-pos[0])>1e-6
        && std::abs(CurrentPos[1]-pos[1])>1e-6 && std::abs(CurrentPos[2]-pos[2])>1e-6)
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
      Ip[0] = I[0]+1;
      Ip[1] = I[1]+1;
      Ip[2] = I[2]+1;
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
      CurrentPos = pos; // save input position
      return true; // new values are evaluated
    }
    else
      return false; // position is identical
  }
};

#include <bitset>

/*! \class template<class T> CubicSpline3D
 * \brief spline functions on 3D regular grid [Rmin, Rmax)^3
 * The number of grids per axis can vary.
 * A set of orbitals are associated with a CubicSpline3D.
 * This class cannot destroy the orbitals. They have to be provided externally.
 */

template<class grid_type,class orb_type>
struct CubicSpline3D
{

  typedef typename orb_type::value_type value_type;
  typedef typename grid_type::pos_type pos_type;

  int OrbitalID;
  std::bitset<4> xyz[8];      // \note possible to make it "static const"
  grid_type myGrid;      // grid
  std::vector<orb_type*> Orb; // a set of orbitals
  Vector<value_type> Val,Lap;
  Vector<pos_type>   Grad;

  //!< default constructor
  CubicSpline3D(): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    myGrid.resize(1,1,1); // set the default variables
  }

  //!< constructor using grid
  CubicSpline3D(int l, int m, int n): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    myGrid.resize(l,m,n);
  }

  //!< copy constructor
  CubicSpline3D(const CubicSpline3D<grid_type,orb_type>& a): OrbitalID(-1)
  {
    setxyz(); // set bit operators
    myGrid = a.myGrid;
    Val = a.Val;
    Lap = a.Lap;
    Grad = a.Grad;
    Orb.reserve(a.Orb.size());
    for(int iorb=0; iorb< a.Orb.size(); iorb++)
      Orb.push_back(a.Orb[iorb]);
  }

  ~CubicSpline3D() { }

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
  int add(orb_type* orbital)
  {
    int id = Orb.size();
    Orb.push_back(orbital);
    // simply resize it
    Val.resize(id+1);
    Lap.resize(id+1);
    Grad.resize(id+1);
    return id;
  }

  //!< function to calculate value, gradient and laplacian
  //\return value[iorb] = value of iorb-th orbial
  //        gr[iorb] = gradient of iorb-th orbial
  //        lap[iorb] = laplacian of iorb-th orbial
  inline void laplacian(const pos_type& pos)
  {
    //grid returns true, if pos != myGrid.CurrentPos;
    if(myGrid.get(pos))
    {
      get_laplacian(pos, Val, Grad, Lap);// evalulate everything
      for(int iorb=0; iorb<Orb.size(); iorb++)
        //assign the values, not sure if needed
      {
        Orb[iorb]->Val = Val[iorb];
        Orb[iorb]->Lap = Lap[iorb];
        Orb[iorb]->Grad =Grad[iorb];
      }
    }
  }

  template<class ValArray, class PosArray>
  inline void laplacian(const pos_type& pos, ValArray& value, PosArray& gr, ValArray& lap)
  {
    if(myGrid.get(pos))
    {
      get_laplacian(pos, value,gr,lap);
      Val = value;
      Grad = gr;
      Lap = lap;
    }
  }


  //!< function to calculate value, gradient and laplacian
  //\return value[iorb] = value of iorb-th orbial
  //        gr[iorb] = gradient of iorb-th orbial
  //        lap[iorb] = laplacian of iorb-th orbial
  template<class ValArray, class PosArray>
  inline void get_laplacian(const pos_type& pos, ValArray& value, PosArray& gr, ValArray& lap)
  {
    // using myGrid data members pre-calculated by myGrid.get(pos)
    // \li Linear coefficients: Lm, Lp
    // \li Cubic coefficients: Cm, Cp
    // \li Gradient of Cubic coefficients: Cvm, Cvp
    // \li I : index of the (nx,ny,nz)
    // \li Ip : index of the (nx+1,ny+1,nz+1)
    value =0.0;
    lap = 0.0;
    gr = 0.0;
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
        Lx = myGrid.Lp[0];
        Cx = myGrid.Cp[0];
        Cvx = myGrid.Cvp[0];
        i = myGrid.I[0];
      }
      else
        // (0,y,z)
      {
        Lx = myGrid.Lm[0];
        Cx = myGrid.Cm[0];
        Cvx = myGrid.Cvm[0];
        i = myGrid.Ip[0];
      }
      // get Ly, Cy, Cvy, and index j to get the value
      if(xyz[ixyz][1])
        // (x,1,z)
      {
        Ly = myGrid.Lp[1];
        Cy = myGrid.Cp[1];
        Cvy = myGrid.Cvp[1];
        j = myGrid.I[1];
      }
      else
        // (x,0,z)
      {
        Ly = myGrid.Lm[1];
        Cy = myGrid.Cm[1];
        Cvy = myGrid.Cvm[1];
        j = myGrid.Ip[1];
      }
      // get Lz, Cz, Cvz, and index k to get the value
      if(xyz[ixyz][2])
        // (x,y,1)
      {
        Lz = myGrid.Lp[2];
        Cz = myGrid.Cp[2];
        Cvz = myGrid.Cvp[2];
        k = myGrid.I[2];
      }
      else
        // (x,y,0)
      {
        Lz = myGrid.Lm[2];
        Cz = myGrid.Cm[2];
        Cvz = myGrid.Cvm[2];
        k = myGrid.Ip[2];
      }
      value_type lxy = Lx*Ly;
      value_type lyz = Ly*Lz;
      value_type lzx = Lz*Lx;
      value_type lxyz = lxy*Lz;
      for(int iorb=0; iorb<Orb.size(); iorb++)
      {
        value_type y = (*Orb[iorb])(i,j,k);
        //\warning: not protected for negative indices
        value_type y2x=((*Orb[iorb])(i+1,j,k)+(*Orb[iorb])(i-1,j,k)-2.0*y)*myGrid.c2(0);
        value_type y2y=((*Orb[iorb])(i,j+1,k)+(*Orb[iorb])(i,j-1,k)-2.0*y)*myGrid.c2(1);
        value_type y2z=((*Orb[iorb])(i,j,k+1)+(*Orb[iorb])(i,j,k-1)-2.0*y)*myGrid.c2(2);
        if(xyz[ixyz][3])
          //positive
        {
          value[iorb] += (lxyz*y + Cx*lyz*y2x + Cy*lzx*y2y + Cz*lxy*y2z);
          lap[iorb]   += lxyz*(y2x+y2y+y2z);
          // negative signs are taken care of here
          // Cvx = -dC/dx, Cvy = -dC/dy, and Cvz = -dC/dz
          // dL/dx = -1/dx
          gr[iorb][0] -= (myGrid.c1(0)*(lyz*y+Cy*Lz*y2y+Ly*Cz*y2z) + Cvx*lyz*y2x);
          gr[iorb][1] -= (myGrid.c1(1)*(lzx*y+Cz*Lx*y2z+Lz*Cx*y2x) + Cvy*lzx*y2y);
          gr[iorb][2] -= (myGrid.c1(2)*(lxy*y+Cx*Ly*y2x+Lx*Cy*y2y) + Cvz*lxy*y2z);
        }
        else
        {
          value[iorb] -= (lxyz*y + Cx*lyz*y2x + Cy*lzx*y2y + Cz*lxy*y2z);
          lap[iorb]   -= lxyz*(y2x+y2y+y2z);
          // negative signs are taken care of here
          // Cvx = -dC/dx, Cvy = -dC/dy, and Cvz = -dC/dz
          // dL/dx = -1/dx
          gr[iorb][0] += (myGrid.c1(0)*(lyz*y+Cy*Lz*y2y+Ly*Cz*y2z) + Cvx*lyz*y2x);
          gr[iorb][1] += (myGrid.c1(1)*(lzx*y+Cz*Lx*y2z+Lz*Cx*y2x) + Cvy*lzx*y2y);
          gr[iorb][2] += (myGrid.c1(2)*(lxy*y+Cx*Ly*y2x+Lx*Cy*y2y) + Cvz*lxy*y2z);
        }
      }
    }
  }

};

#endif
