//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Kenneth P. Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PARTICLE_BCONDS_H
#define QMCPLUSPLUS_PARTICLE_BCONDS_H

#include <config.h>
#include <Lattice/CrystalLattice.h>

namespace APPNAMESPACE 
{
  template<int N,unsigned D>
    struct PowerOfN
    {
      enum {value=N*PowerOfN<N,D-1>::value};
    };

  template<int N>
    struct PowerOfN<N,0>
    {
      enum {value=1};
    };

  const int TwoPowerD=PowerOfN<2,OHMMS_DIM>::value;

  /** generic Boundary condition handler
   *
   * @tparam T real data type
   * @tparam D physical dimension
   * @tparm SC supercell type
   *
   * Default method for any dimension with OPEN boundary condition.
   * \htmlonly
   <table>
   <th>
   <td>SC</td><td>3D</td><td>2D</td><td>1D</td><td>comment</td>
   </th>
   <tr><td>SUPERCELL_OPEN           </td><td>n n n</td><td>n n</td><td>n</td><td>open boudary conditions </td> </tr>
   <tr><td>SUPERCELL_BULK           </td><td>p p p</td><td>NA   </td><td>NA</td><td>periodic boundry conditions in 3 dimensions, general cell</td> </tr>
   <tr><td>SUPERCELL_BULK+TwoPowerD </td><td>p p p</td><td>NA   </td><td>NA</td><td>periodic boundry conditions in 3 dimensions, orthorombic cell</td></tr>
   <tr><td>SUPERCELL_SLAB           </td><td>p p n</td><td>p p</td><td>NA</td><td>periodic boundry conditions in 2 dimensions, general cell</td></tr>
   <tr><td>SUPERCELL_SLAB+TwoPowerD </td><td>p p n</td><td>p p</td><td>NA</td><td>periodic boundry conditions in 2 dimensions, orthorombic cell</td></tr>
   <tr><td>SUPERCELL_WIRE           </td><td>p n n</td><td>p n</td><td>p</td><td>periodic boundry conditions in 1 dimension</td></tr>
   </table>
   * \endhtmlonly
   * Specialization of DTD_BConds should implement
   * - apply_bc(TinyVector<T,D>& displ): apply BC on displ, Cartesian displacement vector, and returns |displ|^2
   * - apply_bc(dr,r,rinv): apply BC on displacements
   * - apply_bc(dr,r): apply BC without inversion calculations
   */
  template<class T, unsigned D, int SC>
    struct DTD_BConds 
    {

      /** constructor: doing nothing */
      inline DTD_BConds(const CrystalLattice<T,D>& lat){}

      /** apply BC on displ and return |displ|^2 
       * @param displ a displacement vector in the Cartesian coordinate
       * @return \f$|displ|^2\f$
       */
      inline T apply_bc(TinyVector<T,D>& displ) const
      {
        return dot(displ,displ);
      }

      /** apply BC on dr and evaluate r and rinv
       * @param dr vector of displacements, in and out
       * @param r vector of distances
       * @param rinv vector of 1/r
       *
       * The input displacement vectors are not modified with the open boundary conditions.
       */
      inline void apply_bc(std::vector<TinyVector<T,D> >& dr
          , std::vector<T>& r
          , std::vector<T>& rinv) const
      {
        const int n=dr.size();
        for(int i=0; i<n; ++i) rinv[i]=dot(dr[i],dr[i]);
        vec_sqrt(n,&rinv[0],&r[0]);
        vec_inv(n,&r[0],&rinv[0]);
      }

      inline void apply_bc(std::vector<TinyVector<T,D> >& dr
          , std::vector<T>& r) const
      {
        for(int i=0;i<dr.size();++i) r[i]=std::sqrt(dot(dr[i],dr[i]));
      }
    };

}

#if OHMMS_DIM == 3
#include <Lattice/ParticleBConds3D.h>
#elif OHMMS_DIM == 2
#include <Lattice/ParticleBConds2D.h>
#endif

#endif // OHMMS_PARTICLE_BCONDS_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
