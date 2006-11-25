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

namespace qmcplusplus {
  /** class to handle various name conventions for hdf5 file
   */
  struct PWParameterSet {
    ///version
    TinyVector<int,2> version;
    ///index of the twist angle
    int twistIndex;
    ///tag for the parameters
    std::string paramTag;
    ///tag for the basis
    std::string basisTag;
    ///tag for the planewaves
    std::string pwTag;
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

    string getTwistAngleName();

    string getTwistName();

    string getTwistName(int i);

    string getBandName(int i);

    string getSpinName(int i);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
