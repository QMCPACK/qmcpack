//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
#ifndef OHMMS_QMC_RADIALGRIDFUNCTOR_BUILDERBASE_H
#define OHMMS_QMC_RADIALGRIDFUNCTOR_BUILDERBASE_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"

namespace ohmmsqmc {

  /**An abstract base class to manage the radial grids and orbitals that are associated with a center.
   *
   */
  struct RGFBuilderBase: public QMCTraits {

    typedef OneDimGridBase<RealType>                        GridType;
    typedef OneDimGridFunctor<ValueType>                    RadialOrbitalType;
    typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;

    ///the species
    string m_species;
    ///the radial orbitals
    CenteredOrbitalType* m_orbitals;

    RGFBuilderBase(): m_orbitals(NULL) { }

    virtual ~RGFBuilderBase(){}

    ///assign a CenteredOrbitalType to work on
    void setOrbitalSet(CenteredOrbitalType* oset, const string& acenter) { 
      m_orbitals = oset;
      m_species = acenter;
    }

    ///add a radial grid
    virtual bool addGrid(xmlNodePtr cur) = 0;
    ///add a radial orbital with quantum numbers nlms
    virtual bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms) = 0;
  };
  

  /**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
   *
   * For a center,
   *   - only one grid is used
   *   - any number of radial orbitals 
   */
  struct STO2GridBuilder: public RGFBuilderBase {

    ///constructor
    STO2GridBuilder(){}
    bool addGrid(xmlNodePtr cur);
    bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
  };

  /**Class to create a set of radial orbitals on a grid (e.g., AtomHF/Siesta)
   *
   * The grid and orbitals are stored in HDF5 format.
   */
  class RGFBuilder: public RGFBuilderBase {

  public:
    ///constructor
    RGFBuilder():Counter(0){}
    ~RGFBuilder();
    bool addGrid(xmlNodePtr cur);
    bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

  private:
    ///keeps track of number of times HDF5 file has been accessed
    int Counter;
    ///id for HDF5 file
    hid_t m_file_id;
    ///id for main group of HDF5 file
    hid_t m_group_id;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
