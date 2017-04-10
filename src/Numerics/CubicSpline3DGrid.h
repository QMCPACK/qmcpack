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


template<class PosA>
struct CubicSpline3DGrid
{

  typedef PosA pos_type;
  typedef typename PosA::Type_t value_type;

  int Npt[3];
  value_type reps;
  pos_type CurrentPos;
  pos_type Rmin, Rmax, dr, drinv, drinv2, dr2over6, drover6;
  std::bitset<4> xyz[8];
  pos_type L[8], C[8], dC[8];
  int I[8][3];

  CubicSpline3DGrid()
  {
    setxyz();
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
    setxyz();
    Rmin = g.Rmin;
    Rmax = g.Rmax;
    resize(g.Npt[0], g.Npt[1], g.Npt[2]);
  }

  CubicSpline3DGrid<PosA>& operator=(const CubicSpline3DGrid<PosA>& g)
  {
    setxyz();
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
    reset();
  }

  void reset()
  {
    dr[0] = (Rmax[0]-Rmin[0])/static_cast<value_type>(Npt[0]);
    dr[1] = (Rmax[1]-Rmin[1])/static_cast<value_type>(Npt[1]);
    dr[2] = (Rmax[2]-Rmin[2])/static_cast<value_type>(Npt[2]);
    // minimim distance 1/100 th of the distance, to be over written
    reps = 1e-6; //rcut2 = dot(dr,dr)/1000000;
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

  inline bool get(const pos_type& pos)
  {
    //how much the gain will be???
    //if(dot(CurrentPos,pos) > rcut2) {
    if(std::abs(CurrentPos[0]-pos[0])>reps
        && std::abs(CurrentPos[1]-pos[1])>reps
        && std::abs(CurrentPos[2]-pos[2])>reps)
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
      CurrentPos = pos; // save input position
      for(int ixyz=0; ixyz<8; ixyz++)
      {
        //get Lx, Cx, Cvx, and index i to get the value
        if(xyz[ixyz][0])
          //(1,y,z)
        {
          L[ixyz][0]  = Lp[0];
          C[ixyz][0]  = Cp[0];
          dC[ixyz][0] = Cvp[0];
          I[ixyz][0]  = im;
        }
        else
          //(0,y,z)
        {
          L[ixyz][0]  = Lm[0];
          C[ixyz][0]  = Cm[0];
          dC[ixyz][0] = Cvm[0];
          I[ixyz][0]  = ip;
        }
        //get Ly, Cy, Cvy, and index j to get the value
        if(xyz[ixyz][1])
          //(x,1,z)
        {
          L[ixyz][1]  = Lp[1];
          C[ixyz][1]  = Cp[1];
          dC[ixyz][1] = Cvp[1];
          I[ixyz][1]  = jm;
        }
        else
          //(x,0,z)
        {
          L[ixyz][1]  = Lm[1];
          C[ixyz][1]  = Cm[1];
          dC[ixyz][1] = Cvm[1];
          I[ixyz][1]  = jp;
        }
        //get Lz, Cz, Cvz, and index k to get the value
        if(xyz[ixyz][2])
          //(x,y,1)
        {
          L[ixyz][2]  = Lp[2];
          C[ixyz][2]  = Cp[2];
          dC[ixyz][2] = Cvp[2];
          I[ixyz][2]  = km;
        }
        else
          //(x,y,0)
        {
          L[ixyz][2]  = Lm[2];
          C[ixyz][2]  = Cm[2];
          dC[ixyz][2] = Cvm[2];
          I[ixyz][2]  = kp;
        }
      }
      return true; // new values are evaluated
    }
    else
      return false; // position is identical
  }
};

/*! \class template<class T> CubicSpline3D
 * \brief spline functions on 3D regular grid [Rmin, Rmax)^3
 * The number of grids per axis can vary.
 * A set of orbitals are associated with a CubicSpline3D.
 * This class cannot destroy the orbitals. They have to be provided externally.
 */

template<class grid_type,class orb_type>
struct CubicSpline3D
{

  typedef typename grid_type::pos_type pos_type;
  typedef typename orb_type::value_type value_type;

  bool OwnGrid;
  value_type Val, Lap;
  pos_type Grad;

  grid_type* myGrid; // grid
  orb_type  Orb;

  //!< default constructor
  CubicSpline3D(): OwnGrid(false) {}

  //!< constructor using grid
  CubicSpline3D(int l, int m, int n): OwnGrid(true)
  {
    myGrid = new grid_type;
    myGrid->resize(l,m,n);
    resize(l, m, n);
  }

  //!< copy constructor
  CubicSpline3D(const CubicSpline3D<grid_type,orb_type>& a)
  {
    OwnGrid = false;
    myGrid = a.myGrid;
    resize(myGrid->Npt[0], myGrid->Npt[1], myGrid->Npt[2]);
  }

  void resize(int l, int m, int n)
  {
    Orb.resize(l+1,m+1,n+1);
  }

  ~CubicSpline3D()
  {
    if(OwnGrid)
    {
      std::cerr << "This CubicSpline3D owns grid. Delete it" << std::endl;
      delete myGrid;
    }
  }
  // assign value at (i,j,k) of iorb-th orbital
  value_type& operator()(int i, int j, int k)
  {
    return Orb(i,j,k);
  }

  // return value at (i,j,k) of iorb-th orbital
  value_type operator()(int i, int j, int k)  const
  {
    return Orb(i,j,k);
  }

  // always return the current value
  value_type laplacian(const pos_type& pos)
  {
    myGrid->get(pos);
    do_evaluate(pos);
    return Val;
  }

  // always return the current value
  value_type evaluate(const pos_type& pos)
  {
    myGrid->get(pos);
    do_evaluate(pos);
    return Val;
  }

  //!< function to calculate value, gradient and laplacian
  //update Val, Grad and Lap
  void do_evaluate(const pos_type& pos)
  {
    // using myGrid data members pre-calculated by myGrid.get(pos)
    // \li Linear coefficients: Lm, Lp
    // \li Cubic coefficients: Cm, Cp
    // \li Gradient of Cubic coefficients: Cvm, Cvp
    // \li I : index of the (nx,ny,nz)
    // \li Ip : index of the (nx+1,ny+1,nz+1)
    Val = 0.0;
    Lap = 0.0;
    Grad = 0.0;
    value_type val, lap;
    pos_type gr;
    // precalculate d^2 phi/dx^2, d^2 phi/dy^2 and d^2 phi/dz^2
    for(int ixyz=0; ixyz<8; ixyz++)
    {
      int i = myGrid->I[ixyz][0];
      int j = myGrid->I[ixyz][1];
      int k = myGrid->I[ixyz][2];
      value_type y = Orb(i,j,k);
      //\warning: not protected for negative indices
      value_type y2x=(Orb(i+1,j,k)+Orb(i-1,j,k)-2.0*y)*myGrid->c2(0);
      value_type y2y=(Orb(i,j+1,k)+Orb(i,j-1,k)-2.0*y)*myGrid->c2(1);
      value_type y2z=(Orb(i,j,k+1)+Orb(i,j,k-1)-2.0*y)*myGrid->c2(2);
      value_type lxy = myGrid->L[ixyz][0]*myGrid->L[ixyz][1];
      value_type lyz = myGrid->L[ixyz][1]*myGrid->L[ixyz][2];
      value_type lzx = myGrid->L[ixyz][2]*myGrid->L[ixyz][0];
      value_type lxyz = lxy*myGrid->L[ixyz][2];
      value_type cxx = myGrid->C[ixyz][0]*y2x;
      value_type cyy = myGrid->C[ixyz][1]*y2y;
      value_type czz = myGrid->C[ixyz][2]*y2z;
      val = lxyz*y + cxx*lyz + cyy*lzx + czz*lxy;
      lap = lxyz*(y2x+y2y+y2z);
      gr[0] = lyz*(myGrid->c1(0)*y+myGrid->dC[ixyz][0]*y2x)
              + myGrid->c1(0)*(cyy*myGrid->L[ixyz][2]+czz*myGrid->L[ixyz][1]);
      gr[1] = lzx*(myGrid->c1(1)*y+myGrid->dC[ixyz][1]*y2y)
              + myGrid->c1(1)*(czz*myGrid->L[ixyz][0]+cxx*myGrid->L[ixyz][2]);
      gr[2] = lxy*(myGrid->c1(2)*y+myGrid->dC[ixyz][2]*y2z)
              + myGrid->c1(2)*(cxx*myGrid->L[ixyz][1]+cyy*myGrid->L[ixyz][0]);
      if(myGrid->xyz[ixyz][3])
        //positive
      {
        Val += val;
        Lap += lap;
        Grad -= gr;
      }
      else
      {
        Val -= val;
        Lap -= lap;
        Grad += gr;
      }
    }
  }

};

#endif
