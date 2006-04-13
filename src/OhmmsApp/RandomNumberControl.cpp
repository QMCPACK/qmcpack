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
#include "RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"

namespace OHMMS{

  /// constructors and destructors
  RandomNumberControl::RandomNumberControl(const char* aname)
    :OhmmsElementBase(aname), NeverBeenInitialized(true), myCur(NULL)
  { }
    
  /// generic output
  bool RandomNumberControl::get(std::ostream& os) const {
    return true;
  }


  /// generic input
  bool RandomNumberControl::put(std::istream& is) {
    return true;
  }

  /// reset the generator
  void RandomNumberControl::reset() {
    if(NeverBeenInitialized) {
      LOGMSG("\tInitialize random number generator: seed will be generated")
      qmcplusplus::Random.init(Controller->mycontext(), Controller->ncontexts(),-1);
      NeverBeenInitialized = false;
    }
  }

  xmlNodePtr 
  RandomNumberControl::initialize(xmlXPathContextPtr acontext) {

    ///initialize the random number generators
    xmlXPathObjectPtr rg_request 
      = xmlXPathEvalExpression((const xmlChar*)"//random",acontext);
    
    ///local xmlNodePtr that will be processed by put
    xmlNodePtr random_ptr = NULL;
    
    std::string parallel_random("true");
    //When no random node is given, use default
    //if(xmlXPathNodeSetIsEmpty(rg_request->nodesetval)) {
    //  add <random/>
    //  myCur = xmlNewNode(NULL,(const xmlChar*)"random");
    //  xmlNewProp(random_ptr,(const xmlChar*)"parallel",
    //    	 (const xmlChar*)parallel_random.c_str());
    //  xmlNewProp(random_ptr,(const xmlChar*)"seed",(const xmlChar*)("-1"));
    //  random_ptr = myCur;
    //} else {
    //  random_ptr = rg_request->nodesetval->nodeTab[0];
    //}
    if(!xmlXPathNodeSetIsEmpty(rg_request->nodesetval)) {
      random_ptr = rg_request->nodesetval->nodeTab[0];
    }
    
    put(random_ptr);
    xmlXPathFreeObject(rg_request);
    return myCur;
  }

  bool RandomNumberControl::put(xmlNodePtr cur){

    if(NeverBeenInitialized) {

      bool init_mpi = true;
      int iseed = -1; // default is to generate by Wall-clock
      if(cur != NULL) {
        xmlAttrPtr att = cur->properties;
        while(att != NULL) {
  	  std::string aname((const char*)(att->name));
  	  std::string vname((const char*)(att->children->content));
	  if(aname == "parallel") {
	    if(vname == "0" || vname == "false") init_mpi = false;
	  } else if(aname == "seed") {
	    iseed = atoi(vname.c_str());
	  }
	  att = att->next;
        }
      }

      int nproc = 1;
      int id = 0;
      if(init_mpi) {
	id = Controller->mycontext();
	nproc = Controller->ncontexts();
      }
      qmcplusplus::Random.init(id,nproc,iseed);
      NeverBeenInitialized = false; 

	//\todo write back to the DOM for restart
//       stringstream s;
//       s << iseed;
//       ///reset random seed: to use random generator to handle reassignment
//       xmlSetProp(cur,(const xmlChar*)("seed"),(const xmlChar*)(s.str().c_str()));
    }
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
