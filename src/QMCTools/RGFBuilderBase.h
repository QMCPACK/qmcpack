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
#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_BUILDERBASE_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_BUILDERBASE_H

#include "Configuration.h"
#include "Numerics/OneDimGridFunctor.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus {

  /**An abstract base class to manage the radial grids and orbitals that are associated with a center.
   *
   */
  struct RGFBuilderBase: public QMCTraits {

    typedef OneDimGridBase<RealType>                        GridType;
    typedef OneDimGridFunctor<RealType>                     RadialOrbitalType;
    typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;

    ///the species
    std::string m_species;
    ///the radial orbitals
    CenteredOrbitalType* m_orbitals;

    RGFBuilderBase(): m_orbitals(NULL), Counter(0), m_file_id(-1), m_group_id(-1){ }

    virtual ~RGFBuilderBase();

    ///assign a CenteredOrbitalType to work on
    void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter) { 
      m_orbitals = oset;
      m_species = acenter;
    }

    ///add a radial grid
    virtual bool addGrid(xmlNodePtr cur);

    ///add a radial orbital with quantum numbers nlms
    virtual bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms) = 0;

    ///add common parameters
    virtual bool putCommon(xmlNodePtr cur) = 0;

    void print(const std::string& aroot, int omode);
  
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
