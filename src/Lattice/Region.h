//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_REGINON_H
#define OHMMS_REGION_H

/* \class Region defines a spatial region bound by [Ri,Rf)
   \brief Defined in unit vectors 0 <= Ri, Rf < 1
*/
template<class T, unsigned D>
struct Region
{

  typedef T Scalar_t;
  enum {DIM = D};
  T Ri[D], Rf[D];

  Region() { }

  Region(const T* r0, const T* dr)
  {
    set(r0,dr);
  }

  Region(const Region<T,D>& rg)
  {
    for(int i=0; i<D; i++)
      Ri[i] = rg.Ri[i];
    for(int i=0; i<D; i++)
      Rf[i] = rg.Rf[i];
  }

  Region<T,D>& operator=(const Region<T,D>& rg)
  {
    for(int i=0; i<D; i++)
      Ri[i] = rg.Ri[i];
    for(int i=0; i<D; i++)
      Rf[i] = rg.Rf[i];
    return *this;
  }

  ~Region() { }

  inline void set(const T* r0, const T*  dr)
  {
    for(int i=0; i<D; i++)
      Ri[i] = r0[i];
    for(int i=0; i<D; i++)
      Rf[i] = r0[i]+dr[i];
  }

  template<class Pos_t>
  inline bool inside(const Pos_t& r) const
  {
    for(int i=0; i<DIM; i++)
    {
      if(r[i] < Ri[i] || r[i] >= Rf[i])
        return false;
    }
    return true;
  }
};
#endif



