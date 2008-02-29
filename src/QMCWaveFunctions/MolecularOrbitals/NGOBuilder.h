//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H
#define QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/OneDimGridFunctor.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"

namespace qmcplusplus {

  /**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
   *
   * For a center,
   *   - only one grid is used
   *   - any number of radial orbitals 
   */
  class NGOBuilder: public QMCTraits {

    public:
    typedef OneDimGridBase<RealType>                        GridType;
    typedef OneDimGridFunctor<RealType>                     RadialOrbitalType;
    typedef SphericalBasisSet<RadialOrbitalType,GridType> CenteredOrbitalType;

    ///true, if the RadialOrbitalType is normalized
    bool Normalized;
    ///the radial orbitals
    CenteredOrbitalType* m_orbitals;
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

    /** add a radial functor
     * @param cur xml element
     * @param nlms quantum number
     */
    bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

    /** put common element
     * @param cur xml element
     */
    bool putCommon(xmlNodePtr cur);

    private:
    void addGaussian(xmlNodePtr cur);
    void addSlater(xmlNodePtr cur);
    void addNumerical(xmlNodePtr cur, const string& dsname);
    void addPade(xmlNodePtr cur);
    hid_t m_fileid;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
