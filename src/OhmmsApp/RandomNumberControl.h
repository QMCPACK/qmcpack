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
#ifndef OHMMS_RANDOMNUMBERCONTROL_H__
#define OHMMS_RANDOMNUMBERCONTROL_H__

#include "OhmmsData/OhmmsElementBase.h"

namespace OHMMS {

  /**class RandomNumberControl
   *\brief Encapsulate data to initialize and save the status of the random number generator
   *
   * Default:  myName = "random"
   *
   */
  class RandomNumberControl: public OhmmsElementBase {

  public:

    /// constructors and destructors
    RandomNumberControl(const char* aname="random");
    
    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    xmlNodePtr initialize(xmlXPathContextPtr);

  private:

    bool NeverBeenInitialized;
    xmlNodePtr myCur;
  };
}

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
