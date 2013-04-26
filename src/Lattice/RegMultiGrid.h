//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_REGULARMULTIGRID_H
#define OHMMS_REGULARMULTIGRID_H
#include <vector>
using namespace std;

/*! \class RegMultiGrid<class T, unsigned D>
 *  \brief A class for grid layout using octree
 */
template<class T, unsigned D>
struct RegMultiGrid
{

  typedef RegMultiGrid<T,D> this_t;
  Region<T,D> R;

  RegMultiGrid():Up(NULL) { }
  ~RegMultiGrid();

  inline bool root() const
  {
    return Up == NULL;
  }
  inline bool empty() const
  {
    return Down.empty();
  }

  inline const RegMultiGrid* node() const
  {
    return Up;
  }
  inline const RegMultiGrid* child(int i) const
  {
    if(Down.empty())
      return NULL;
    else
      return Down[i];
  }

  void refine();

  void remove()
  {
    for(int i=0; i<Down.size(); i++)
      delete Down[i];
  }

  template<class pos_t>
  inline int node(const pos_t& pos)
  {
    if(R.inside(pos))
      return R.inside(pos);
  }

  this_t* Up;
  vector<this_t* > Down;
};

template<class T, unsigned D>
RegMultiGrid<T,D>::~RegMultiGrid()
{
  for(int i=0; i<Down.size(); i++)
    delete Down[i];
}

template<class T, unsigned D>
void RegMultiGrid<T,D>::refine() {}

template<class T>
void RegMultiGrid<T,3>::refine()
{
  if(Down.empty())
  {
    for(int i=0; i<8; i++)
      Down.push_back(new RegMultGrid<T,D>);
    T dx[D];
    dx[0] = R.Rf[0]-R.Ri[0];
    dx[1] = R.Rf[1]-R.Ri[1];
    dx[2] = R.Rf[2]-R.Ri[2];
    for(int idim=0; idim<D; idim++)
      dx[idim] = (R.Rf[idim]-R.Ri[dim])*0.5;
  }
  if(Down.empty())
  {
    for(int i=0; i<nsub; i++)
      Down.push_back(new RegMultGrid<T,D>);
  }
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
