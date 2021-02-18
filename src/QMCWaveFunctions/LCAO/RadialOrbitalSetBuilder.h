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


#ifndef QMCPLUSPLUS_RADIAL_NUMERICALGRIDORBITALBUILDER_H
#define QMCPLUSPLUS_RADIAL_NUMERICALGRIDORBITALBUILDER_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "hdf/HDFVersion.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/OneDimGridFactory.h"
#include "Message/MPIObjectBase.h"


namespace qmcplusplus
{
/** abstract class to defer generation of spline functors
   * @tparam T precision of the final result
  */
template<typename T>
struct TransformerBase
{
  ///temporary grid in double precision
  typedef OneDimGridBase<double> grid_type;
  ///the multiple set
  typedef MultiQuinticSpline1D<T> FnOut;
  /** convert input 1D functor to the multi set
       * @param agrid  original grid
       * @param multiset the object that should be populated
       * @param ispline index of the this analytic function
       * @param int order quintic (or cubic) only quintic is used
       */
  virtual void convert(grid_type* agrid, FnOut* multiset, int ispline, int order) = 0;
  virtual ~TransformerBase() {}
};

template<typename T, typename FnIn>
struct A2NTransformer : TransformerBase<T>
{
  typedef typename TransformerBase<T>::grid_type grid_type;
  typedef typename TransformerBase<T>::FnOut FnOut;

  FnIn* m_ref; //candidate for unique_ptr
  A2NTransformer(FnIn* in) : m_ref(in) {}
  ~A2NTransformer() { delete m_ref; }

  void convert(grid_type* agrid, FnOut* multiset, int ispline, int order)
  {
    typedef OneDimQuinticSpline<OHMMS_PRECISION_FULL> spline_type;
    spline_type radorb(agrid);
    Transform2GridFunctor<FnIn, spline_type> transform(*m_ref, radorb);
    transform.generate(agrid->rmin(), agrid->rmax(), agrid->size());
    multiset->add_spline(ispline, radorb);
  }
};


/** Build a set of radial orbitals at the origin
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 *   - using MultiQuinticSpline1D
 *   - Gaussian and Slater mixing allowed
 *   - Slater's are assumed to fit in a 100 rmax grid
 */
template<typename COT>
class RadialOrbitalSetBuilder : public MPIObjectBase
{
public:
  //Now we count on COT::ValueType being real
  //since COT::ValueType == ROT::ValueType && ROT==MultiQuiniticSpline1D
  //And MultiQuiniticSpline1D does not distringuish between real_type of r
  //And the value_type of 'Psi(r)' this is safe for now
  //I think COT::gridtype needs to go to a template<typename RT, typename VT>
  //For this all to become 'correct'.
  using RealType = typename COT::RealType;
  typedef typename COT::RadialOrbital_t RadialOrbitalType;
  using GridType = typename COT::GridType;


  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  COT* m_orbitals;
  ///input grid in case transform is needed
  GridType* input_grid;
  ///the quantum number of this node
  QuantumNumberType m_nlms;

  ///constructor
  RadialOrbitalSetBuilder(Communicate* comm,
                          int radial_grid_size = 1001); //radial_grid_size is 1001 just magic?
  ///destructor: cleanup gtoTemp, stoTemp
  ~RadialOrbitalSetBuilder();

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(COT* oset, const std::string& acenter);

  ///add a grid
  bool addGrid(xmlNodePtr cur, const std::string& rad_type);
  bool addGridH5(hdf_archive& hin);

  /** add a radial functor
   * @param cur xml element
   * @param nlms quantum number
   */
  bool addRadialOrbital(xmlNodePtr cur, const std::string& m_infunctype, const QuantumNumberType& nlms);
  bool addRadialOrbitalH5(hdf_archive& hin, const std::string& radtype, const QuantumNumberType& nlms);

  /** This is when the radial orbitals are actually created */
  void finalize();

private:
  //only check the cutoff
  void addGaussian(xmlNodePtr cur);
  void addGaussianH5(hdf_archive& hin);
  void addSlater(xmlNodePtr cur);

  template<typename Fin, typename T>
  T find_cutoff(Fin& in, T rmax);

  /// hdf file only for numerical basis h5 file generated by SQD
  hdf_archive hin;

  /// Number of grid points
  int radial_grid_size_;

  ///maximum cutoff for an orbital, could come from XML always set to -1 by constructor
  RealType m_rcut;

  ///safe common cutoff radius
  RealType m_rcut_safe;

  /** radial functors to be finalized
   */
  std::vector<TransformerBase<RealType>*> radTemp;
  ///store the temporary analytic data
  std::vector<GaussianCombo<RealType>*> gtoTemp;
  std::vector<SlaterCombo<RealType>*> stoTemp;

  std::tuple<int, double, double> grid_param_in;
};

template<typename COT>
RadialOrbitalSetBuilder<COT>::RadialOrbitalSetBuilder(Communicate* comm, int radial_grid_size)
    : MPIObjectBase(comm),
      Normalized(true),
      m_orbitals(nullptr),
      input_grid(nullptr),
      radial_grid_size_(radial_grid_size),
      m_rcut(-1.0)
{}

template<typename COT>
RadialOrbitalSetBuilder<COT>::~RadialOrbitalSetBuilder()
{
  if (input_grid != nullptr)
    delete input_grid;
}

template<typename COT>
void RadialOrbitalSetBuilder<COT>::setOrbitalSet(COT* oset, const std::string& acenter)
{
  m_orbitals = oset;
}

template<typename COT>
bool RadialOrbitalSetBuilder<COT>::addGrid(xmlNodePtr cur, const std::string& rad_type)
{
  if (!m_orbitals)
  {
    APP_ABORT("NGOBuilder::addGrid SphericalOrbitals<ROT,GT>*, is not initialized");
  }

  if (rad_type == "Numerical")
  {
    hin.push("grid");
    addGridH5(hin);
    hin.pop();
  }
  else
  {
    GridType* agrid = OneDimGridFactory::createGrid(cur);
    m_orbitals->Grids.push_back(agrid);
  }

  //set zero to use std::max
  m_rcut_safe = 0;

  return true;
}

template<typename COT>
bool RadialOrbitalSetBuilder<COT>::addGridH5(hdf_archive& hin)
{
  if (!m_orbitals)
  {
    APP_ABORT("NGOBuilder::addGrid SphericalOrbitals<ROT,GT>*, is not initialized");
  }

  app_log() << "   Grid is created by the input paremters in h5" << std::endl;

  std::string gridtype;
  if (myComm->rank() == 0)
  {
    if (!hin.readEntry(gridtype, "grid_type"))
    {
      std::cerr << "Could not read grid_type in H5; Probably Corrupt H5 file" << std::endl;
      exit(0);
    }
  }
  myComm->bcast(gridtype);

  int npts    = 0;
  RealType ri = 0.0, rf = 10.0, rmax_safe = 10;
  if (myComm->rank() == 0)
  {
    double tt = 0;
    hin.read(tt, "grid_ri");
    ri = tt;
    hin.read(tt, "grid_rf");
    rf = tt;
    // Ye TODO: grid handling will all moved to XML.
    //hin.read(tt, "rmax_safe");
    //rmax_safe = tt;
    hin.read(npts, "grid_npts");
  }
  myComm->bcast(ri);
  myComm->bcast(rf);
  myComm->bcast(rmax_safe);
  myComm->bcast(npts);

  if (gridtype.empty())
  {
    APP_ABORT("Grid type is not specified.");
  }
  if (gridtype == "log")
  {
    app_log() << "    Using log grid ri = " << ri << " rf = " << rf << " npts = " << npts << std::endl;
    input_grid = new LogGrid<RealType>;
    input_grid->set(ri, rf, npts);
    m_orbitals->Grids.push_back(input_grid);
    input_grid = 0;
  }
  else if (gridtype == "linear")
  {
    app_log() << "    Using linear grid ri = " << ri << " rf = " << rf << " npts = " << npts << std::endl;
    input_grid = new LinearGrid<RealType>;
    input_grid->set(ri, rf, npts);
    m_orbitals->Grids.push_back(input_grid);
    input_grid = 0;
  }
  //set zero to use std::max
  m_rcut_safe = 0;
  return true;
}

/** Add a new Slater Type Orbital with quantum numbers \f$(n,l,m,s)\f$
   * \param cur  the current xmlNode to be processed
   * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
   * \return true is succeeds
   *
   */
template<typename COT>
bool RadialOrbitalSetBuilder<COT>::addRadialOrbital(xmlNodePtr cur,
                                                    const std::string& m_infunctype,
                                                    const QuantumNumberType& nlms)
{
  if (!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  std::string radtype(m_infunctype);
  std::string dsname("0");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(radtype, "type");
  aAttrib.add(m_rcut, "rmax");
  aAttrib.add(dsname, "ds");
  aAttrib.put(cur);
  m_nlms = nlms;
  if (radtype == "Gaussian" || radtype == "GTO")
  {
    addGaussian(cur);
  }
  else if (radtype == "Slater" || radtype == "STO")
  {
    addSlater(cur);
  }
  else
  {
    myComm->barrier_and_abort("Purely numerical atomic orbitals are not supported any longer.");
  }
  return true;
}

template<typename COT>
bool RadialOrbitalSetBuilder<COT>::addRadialOrbitalH5(hdf_archive& hin,
                                                      const std::string& radtype_atomicBasisSet,
                                                      const QuantumNumberType& nlms)
{
  if (!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  std::string dsname("0");
  std::string radtype(radtype_atomicBasisSet);
  if (myComm->rank() == 0)
    hin.read(radtype, "type");
  myComm->bcast(radtype);

  m_nlms = nlms;
  if (radtype == "Gaussian" || radtype == "GTO")
  {
    addGaussianH5(hin);
  }
  else if (radtype == "Slater" || radtype == "STO")
  {
    // addSlaterH5(hin);
    APP_ABORT(
        " RadType: Slater. Any type other than Gaussian not implemented in H5 format. Please contact developers.");
  }
  else
  {
    APP_ABORT(
        " RadType: Numerical. Any type other than Gaussian not implemented in H5 format. Please contact developers.");
  }
  return true;
}


template<typename COT>
void RadialOrbitalSetBuilder<COT>::addGaussian(xmlNodePtr cur)
{
  int L          = m_nlms[1];
  using gto_type = GaussianCombo<OHMMS_PRECISION_FULL>;
  gto_type* gset = new gto_type(L, Normalized);
  gset->putBasisGroup(cur);
  //Warning::Magic Number for max rmax of gaussians
  RealType r0 = find_cutoff(*gset, 100.);
  m_rcut_safe = std::max(m_rcut_safe, r0);
  radTemp.push_back(new A2NTransformer<RealType, gto_type>(gset));
  m_orbitals->RnlID.push_back(m_nlms);
}


template<typename COT>
void RadialOrbitalSetBuilder<COT>::addGaussianH5(hdf_archive& hin)
{
  int L          = m_nlms[1];
  using gto_type = GaussianCombo<OHMMS_PRECISION_FULL>;
  gto_type* gset = new gto_type(L, Normalized);
  gset->putBasisGroupH5(hin);
  //at least gamess derived xml seems to provide the max its grid goes to
  //So in priniciple this 100 should be coming in from input
  //m_rcut seems like it once served this purpose but is somehow
  //a class global variable even though it should apply here and
  //similar locations on a function by function basis.
  RealType r0 = find_cutoff(*gset, 100.);
  m_rcut_safe = 6 * std::max(m_rcut_safe, r0);
  radTemp.push_back(new A2NTransformer<RealType, gto_type>(gset));
  m_orbitals->RnlID.push_back(m_nlms);
}


/* Finalize this set using the common grid
 *
 * This function puts the All RadialOrbitals 
 * on a logarithmic grid with r_max of matching the largest r_max found
 * filling radTemp. The derivatives at the endpoint
 * are assumed to be all zero.  Note: for the radial orbital we use
 * \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
 * radial orbital and \f$ l \f$ is the angular momentum.
 *
 */
template<typename COT>
void RadialOrbitalSetBuilder<COT>::finalize()
{
  MultiQuinticSpline1D<RealType>* multiset = new MultiQuinticSpline1D<RealType>;
  int norbs                                = radTemp.size();

  // This is a temporary grid used in conversion, at the full precision of the calculation
  // to reduce error. It doesn't need to be on the heap but this is the result of a
  // series of design decisions requiring a base class pointer here.
  OneDimGridBase<OHMMS_PRECISION_FULL>* grid_prec;
  grid_prec = new LogGrid<OHMMS_PRECISION_FULL>;
  grid_prec->set(1.e-6, m_rcut_safe, 1001);
  multiset->initialize(*grid_prec, norbs);

  for (int ib = 0; ib < norbs; ++ib)
    radTemp[ib]->convert(grid_prec, multiset, ib, 5);

  m_orbitals->MultiRnl = multiset;

  app_log() << "  Setting cutoff radius " << m_rcut_safe << std::endl << std::endl;
  m_orbitals->setRmax(static_cast<RealType>(m_rcut_safe));

  delete grid_prec;
}

template<typename COT>
void RadialOrbitalSetBuilder<COT>::addSlater(xmlNodePtr cur)
{
  using sto_type = SlaterCombo<OHMMS_PRECISION_FULL>;
  sto_type* gset = new sto_type(m_nlms[1], Normalized);

  gset->putBasisGroup(cur);

  //need a find_cutoff for STO's, but this was previously in finalize and wiping out GTO's m_rcut_safe
  m_rcut_safe = std::max(m_rcut_safe, static_cast<RealType>(100));
  radTemp.push_back(new A2NTransformer<RealType, sto_type>(gset));
  m_orbitals->RnlID.push_back(m_nlms);
}


/** compute the safe cutoff radius of a radial functor
   */
/** temporary function to compute the cutoff without constructing NGFunctor */
template<typename COT>
template<typename Fin, typename T>
T RadialOrbitalSetBuilder<COT>::find_cutoff(Fin& in, T rmax)
{
  LogGridLight<OHMMS_PRECISION_FULL> agrid;
  //WARNING Magic number, should come from input or be set somewhere more cnetral.
  const OHMMS_PRECISION_FULL eps = 1e-6;
  bool too_small                 = true;
  agrid.set(eps, rmax, RadialOrbitalSetBuilder<COT>::radial_grid_size_);
  int i = radial_grid_size_ - 1;
  T r   = rmax;
  while (too_small && i > 0)
  {
    r         = agrid(i--);
    T x       = in.f(r);
    too_small = (std::abs(x) < eps);
  }
  return static_cast<OHMMS_PRECISION_FULL>(r);
}

} // namespace qmcplusplus
#endif
