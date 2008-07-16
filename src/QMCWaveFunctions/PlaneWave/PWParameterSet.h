/////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
// -*- C++ -*-
/** @file PWParameterSet.cpp
 * @brief Utility class to handle hdf5
 */
#ifndef QMCPLUSPLUS_PWPARAMETERSET_H
#define QMCPLUSPLUS_PWPARAMETERSET_H
#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/MPIObjectBase.h"

namespace qmcplusplus {
  /** class to handle various name conventions for hdf5 file
   */
  struct PWParameterSet: public MPIObjectBase
  {
    ///true if spin channel exists
    bool hasSpin;
    ///version
    TinyVector<int,2> version;
    ///index of the twist angle
    int twistIndex;
    ///number of input bands
    int numBands;
    ///energy cutoff for QMC wavefunction
    double Ecut;
    ///cutoff radius for truncated orbitals
    double Rcut;
    ///radius of buffer layer for truncated orbitals
    double BufferRadius;
    ///cell multiplications
    TinyVector<int,OHMMS_DIM> BoxDup;
    ///tag for the parameters
    std::string paramTag;
    ///tag for the basis
    std::string basisTag;
    ///tag for the planewaves
    std::string pwTag;
    ///tag for the multipliers of the planewaves
    std::string pwMultTag;
    ///tag for eigentstates
    std::string eigTag;
    ///tag for twist angles
    std::string twistTag;
    ///tag for the band
    std::string bandTag;
    ///tag for the spin
    std::string spinTag;
    ///tag for eigvector
    std::string eigvecTag;
    ///xml processor
    ParameterSet m_param;

    PWParameterSet();

    bool put(xmlNodePtr cur)
    {
      return m_param.put(cur);
    }

    void checkVersion(hid_t h);

    double getEcut(double ecut);

    /** get the dimensions of the eigenvectors
     * @param h fileid
     * @return true, if the data is complex
     */
    bool getEigVectorType(hid_t h);

    bool hasComplexData(hid_t h);

    string getTwistAngleName();

    string getTwistName();

    string getTwistName(int i);

    string getBandName(int ib, int ispin);

    string getBandName(int ib);

    string getSpinName(int ispin);

    string getEigVectorName(const string& hg, int ib, int ispin);
    string getEigVectorName(int ib, int ispin);
    string getCenterName(const string& hg,int ib);
    string getOriginName(const string& hg,int ib);

    void writeParameters(hid_t gid);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
