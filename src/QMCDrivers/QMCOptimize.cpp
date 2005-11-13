//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/QMCOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Optimize/CGOptimization.h"
#include "Optimize/DampedDynamics.h"
#if defined(HAVE_LIBGSL)
#include "Optimize/GSLMinimize.h"
#endif

namespace qmcplusplus {

  QMCOptimize::QMCOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), PartID(0), NumParts(1),optTarget(0),
    wfNode(NULL), optNode(NULL)
  { 
    //set the optimization flag
    QMCDriverMode.set(QMC_OPTIMIZE,1);
    //read to use vmc output (just in case)
    RootName = "vmc";
    QMCType ="opt";
    //default method is cg
    optmethod = "cg";
  }

  /** Clean up the vector */
  QMCOptimize::~QMCOptimize() { }

  /** Add configuration files for the optimization
   * @param a root of a hdf5 configuration file
   */
  void QMCOptimize::addConfiguration(const string& a) {
    if(a.size()) ConfigFile.push_back(a);
  }

  /** Reimplement QMCDriver::run
   */
  bool
  QMCOptimize::run() {

    //reset the rootname
    optTarget->setRootName(RootName);
    optTarget->setWaveFunctionNode(wfNode);

    //set the data members to start a new run
    optTarget->checkConfigurations(ConfigFile,PartID,NumParts);

    //estimator has to collect the data over mpi nodes
    Estimators->setCollectionMode(OHMMS::Controller->ncontexts()>1);

    //overwrite the Etarget by E_T if E_T is zero
    if(abs(branchEngine->E_T)>numeric_limits<RealType>::epsilon()) {
      optTarget->setTargetEnergy(branchEngine->E_T);
    }

    MinimizerBase<RealType>* solver;

    if(optmethod == "cg") {
      app_log() << " Conjugate-gradient optimization using CGOptimization"<<endl;
      solver = new CGOptimization<RealType>;
    } 
    else if(optmethod == "anneal") {
      app_log() << " Annealing optimization using DampedDynamics"<<endl;
      solver = new DampedDynamics<RealType>;
    } 
#if defined(HAVE_LIBGSL)
    else if(optmethod == "gslcg") {
      app_log() << " Conjugate-gradient optimization using GSLConjugateGradient"<<endl;
      solver = new GSLConjugateGradient;
    }
#endif

    //set the stream
    solver->setOstream(&app_log());

    if(optNode == NULL) 
      solver->put(qmcNode);
    else
      solver->put(optNode);

    bool success = solver->optimize(optTarget);

    delete solver;
    return success;
  }

  /** Parses the xml input file for parameter definitions for the wavefunction optimization.
   * @param q current xmlNode 
   * @return true if successful
   */
  bool
  QMCOptimize::put(xmlNodePtr q) {

    xmlNodePtr qsave=q;
    xmlNodePtr cur=qsave->children;
    int pid=OHMMS::Controller->mycontext();
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "mcwalkerset") {
        mcwalkerNodePtr.push_back(cur);
      } else if(cname == "optimizer") {
        xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
        if(att) { optmethod = (const char*)att; }
        optNode=cur;
      } else if(cname == "optimize") {
        xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
        if(att) { optmethod = (const char*)att; }
      }
      cur=cur->next;
    }  
    //will check the method

    int nfile=mcwalkerNodePtr.size();
    if(nfile) {
      int ng=OHMMS::Controller->ncontexts()/nfile;
      if(ng>=1) {
        NumParts=ng;
        PartID=pid%ng;
        int mygroup=pid/ng;
        string fname("invalid");
        OhmmsAttributeSet pAttrib;
        pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
        pAttrib.put(mcwalkerNodePtr[mygroup]);
        if(fname != "invalid") ConfigFile.push_back(fname);
      } else {//more files than the number of processor
        for(int ifile=0; ifile<nfile; ifile++) {
          int pid_target=pid;
          string fname("invalid");
          OhmmsAttributeSet pAttrib;
          pAttrib.add(pid_target,"node"); 
          pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
          pAttrib.put(mcwalkerNodePtr[ifile]);
          if(pid_target == pid && fname != "invalid") {
            ConfigFile.push_back(fname);
          }
        }
      }
    }

    app_log() 
      << " list of the configuration files used for optimization" << endl;
    for(int i=0; i<ConfigFile.size(); i++) 
      app_log() << "# " << ConfigFile[i] << " part " << PartID << "/" << NumParts << endl;

    if(optTarget == 0) {
      optTarget = new QMCCostFunction(W,Psi,H);
      //THIS IS TEST ONLY
      optTarget->setStream(&app_log());
    }
    return optTarget->put(q);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
