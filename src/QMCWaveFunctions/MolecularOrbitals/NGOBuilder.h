//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H
#define QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H

#include "Configuration.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"
#include "io/hdf_archive.h"

namespace qmcplusplus
{

struct NGOrbital: public OptimizableFunctorBase
{
  typedef real_type                    value_type;
  typedef real_type                    point_type;
  typedef OneDimGridBase<real_type>    grid_type;
#if QMC_BUILD_LEVEL>2
  typedef OneDimQuinticSpline<real_type> functor_type;
#else
  typedef OneDimCubicSpline<real_type> functor_type;
#endif
  functor_type myFunc;
  real_type Y, dY, d2Y, d3Y;

  NGOrbital(grid_type* agrid):myFunc(agrid) { }

  template<typename VV>
  NGOrbital(grid_type* agrid, const VV& nv):myFunc(agrid,nv) { }

  void checkInVariables(opt_variables_type& active) {}
  void checkOutVariables(const opt_variables_type& active) {}
  void resetParameters(const opt_variables_type& active) {}
  void reset() {}
  inline real_type f(real_type r)
  {
    return myFunc.f(r);
  }
  inline real_type df(real_type r)
  {
    return myFunc.df(r);
  }
  bool put(xmlNodePtr cur)
  {
    return true;
  }
  OptimizableFunctorBase* makeClone() const;

  inline real_type evaluate(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r);
  }
  inline value_type evaluateAll(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r,dY,d2Y);
  }

  inline value_type evaluateWithThirdDeriv(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r,dY,d2Y,d3Y);
  }

  inline value_type operator()(int i) const
  {
    return myFunc(i);
  }
  inline value_type& operator()(int i)
  {
    return myFunc(i);
  }
  inline grid_type& grid()
  {
    return myFunc.grid();
  }
  inline void setGridManager(bool willmanage)
  {
    myFunc.setGridManager(willmanage);
  }

  inline void spline(int imin, value_type yp1, int imax, value_type ypn)
  {
    myFunc.spline(imin,yp1,imax,ypn);
  }

  inline void resize(int n)
  {
    myFunc.resize(n);
  }
};

/**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 */
class NGOBuilder: public QMCTraits
{

public:
  //typedef OneDimGridBase<RealType>                        GridType;
  //typedef OneDimGridFunctor<RealType>                     RadialOrbitalType;
  typedef NGOrbital                                     RadialOrbitalType;
  typedef NGOrbital::grid_type                          GridType;
  typedef SphericalBasisSet<RadialOrbitalType,GridType> CenteredOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///input grid in case transform is needed
  GridType *input_grid;
  ///maximum cutoff
  RealType m_rcut;
  ///the quantum number of this node
  QuantumNumberType m_nlms;
  ///the species
  std::string m_species;
  ///type of input function
  std::string m_infunctype;

  ///constructor
  NGOBuilder(xmlNodePtr cur=NULL);
  ///destructor
  ~NGOBuilder();

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter);

  ///add a grid
  bool addGrid(xmlNodePtr cur);
  bool addGridH5(hdf_archive &hin);

  /** add a radial functor
   * @param cur xml element
   * @param nlms quantum number
   */
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
  bool addRadialOrbitalH5(hdf_archive &hin, const QuantumNumberType& nlms);

  /** put common element
   * @param cur xml element
   */
  bool putCommon(xmlNodePtr cur);

private:
  void addGaussian(xmlNodePtr cur);
  void addGaussianH5(hdf_archive &hin);
  void addSlater(xmlNodePtr cur);
  void addNumerical(xmlNodePtr cur, const std::string& dsname);
  void addPade(xmlNodePtr cur);
  hid_t m_fileid;
};

}
#endif
