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
    
    


#ifndef OHMMS_ITERATOR_UTILITIES_H
#define OHMMS_ITERATOR_UTILITIES_H

namespace qmcplusplus
{
/** delete the pointers in [first,last)
*/
template<class IT>
inline void delete_iter(IT first, IT last)
{
  while(first != last)
  {
    if(*first)
      delete *first;
    ++first;
  }
}


template<typename IT1, typename IT2>
inline void accumulate_elements(IT1 first, IT1 last, IT2 res)
{
  while(first != last)
    *res++ += *first++;
}

//  template<typename IT1, typename IT2, typename INT>
//    inline void accumulate_elements(IT1 first, IT2 res, INT n)
//    {
//      for(;n>0; n--) *res++ += *first++;
//    }

}
#endif
