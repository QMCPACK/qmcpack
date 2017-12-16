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
    
    


#ifndef OHMMS_BITSET_IBMGCC_H
#define OHMMS_BITSET_IBMGCC_H

template<unsigned N>
struct bitset
{

  bool data[N];

  bitset()
  {
    reset();
  }
  inline void reset()
  {
    for(int i=0; i<N; i++)
      data[i] = false;
  }
  inline void set(int i)
  {
    data[i] = true;
  }
  inline void flip(int i)
  {
    data[i] =  (data[i])? false: true;
  }

  inline bool operator[](int i) const
  {
    return data[i];
  }
  inline bool& operator[](int i)
  {
    return data[i];
  }

  inline bool any() const
  {
    int i=0;
    while(i<N)
    {
      if(data[i++])
        return true;
    }
    return false;
  }
};
#endif // OHMMS_BITSET_IBMGCC_H

