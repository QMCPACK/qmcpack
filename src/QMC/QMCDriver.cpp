//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#include "QMC/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  int QMCDriver::Counter = -1;
  
  QMCDriver::QMCDriver(MCWalkerConfiguration& w, 
		       TrialWaveFunction& psi, 
		       QMCHamiltonian& h, 
		       xmlNodePtr q):
    nAccept(0),  nReject(0), nTargetWalkers(0),
    Estimators(h), W(w), Psi(psi), H(h),
    Tau(0.01), FirstStep(0.2), nBlocks(100), nSteps(1000), 
    pStride(false), qmc_node(q), LogOut(NULL),
    QMCType("invalid") { 
    
    Counter++; 
    m_param.add(nSteps,"steps","int");
    m_param.add(nBlocks,"blocks","int");
    m_param.add(nTargetWalkers,"walkers","int");
    m_param.add(Tau,"Tau","AU");
    m_param.add(FirstStep,"FirstStep","AU");
  }

  QMCDriver::~QMCDriver() { 
    
    if(Estimators.size())
      W.setLocalEnergy(Estimators.average(0));
    
    if(LogOut) delete LogOut;
  }

  /**Sets the root file name for all output files.!
   * \param aname the root file name
   *
   * All output files will be of
   * the form "aname.s00X.suffix", where "X" is number
   * of previous QMC runs for the simulation and "suffix"
   * is the suffix for the output file. 
   */
  void QMCDriver::setFileRoot(const string& aname) {
    RootName = aname;
    
    char logfile[128];
    sprintf(logfile,"%s.%s",RootName.c_str(),QMCType.c_str());
    
    if(LogOut) delete LogOut;
    LogOut = new OhmmsInform(" ",logfile);
    
    LogOut->getStream() << "Starting a " << QMCType << " run " << endl;
  }

  /** initialize estimators and other internal data */  
  void QMCDriver::getReady() {
    
    Estimators.resetReportSettings(RootName);
    AcceptIndex = Estimators.addColumn("AcceptRatio");
    Estimators.reportHeader();
  }
  

  /** Add walkers to the end of the ensemble of walkers.  
   *@param nwalkers number of walkers to add
   *@return true, if the walker configuration is not empty.
   *
   * Assign positions to any new 
   * walkers \f[ {\bf R}[i] = {\bf R}[i-1] + g{\bf \chi}, \f]
   * where \f$ g \f$ is a constant and \f$ {\bf \chi} \f$
   * is a 3N-dimensional gaussian.
   * As a last step, for each walker calculate 
   * the properties given the new configuration
   <ul>
   <li> Local Energy \f$ E_L({\bf R} \f$
   <li> wavefunction \f$ \Psi({\bf R}) \f$
   <li> wavefunction squared \f$ \Psi^2({\bf R}) \f$
   <li> weight\f$ w({\bf R}) = 1.0\f$  
   <li> drift velocity \f$ {\bf v_{drift}}({\bf R})) \f$
   </ul>
  */
  void 
  QMCDriver::addWalkers(int nwalkers) {

    if(nwalkers>0) {
      //add nwalkers walkers to the end of the ensemble
      int nold = W.getActiveWalkers();

      LOGMSG("Adding " << nwalkers << " walkers to " << nold << " existing sets")

      W.createWalkers(nwalkers);
      LogOut->getStream() <<"Added " << nwalkers << " walkers" << endl;
      
      ParticleSet::ParticlePos_t rv(W.getTotalNum());
      RealType g = FirstStep;
      
      MCWalkerConfiguration::iterator it = W.begin(), itprev;
      int iw = 0;
      while(it != W.end()) {
	if(iw>=nold) {
	  makeGaussRandom(rv);
	  if(iw)
	    (*it)->R = (*itprev)->R+g*rv;
	  else
	    (*it)->R = W.R+g*rv;
	}
	itprev = it;
	it++;iw++;
      }
    } else {
      LOGMSG("Using the existing walkers. Numer of walkers =" << W.getActiveWalkers())
    }

    //calculate local energies and wave functions:
    //can be redundant but the overhead is small
    LOGMSG("Evaluate all the walkers before starting")

    MCWalkerConfiguration::iterator it = W.begin();
    while(it != W.end()) {
      W.R = (*it)->R;
      //evaluate the distance table
      DistanceTable::update(W);
      //evaluate the wavefunction
      ValueType psi = Psi.evaluate(W);
      RealType psi2 = psi*psi;
      //update the properties
      (*it)->Properties(PsiSq) = psi2;
      (*it)->Properties(Sign) = psi;
      (*it)->Properties(Weight) = 1.0;
      (*it)->Properties(Multiplicity) = 1.0;
      //evaluate the Hamiltonian
      (*it)->Properties(LocalEnergy) = H.evaluate(W);
      (*it)->E.resize(H.size()+1,0.0);
      //copy the data from the Hamiltonian to the walker
      H.get((*it)->E);
      ValueType vsq = Dot(W.G,W.G);
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      (*it)->Drift = scale*W.G;
      //(*it)->Drift = Tau*W.G;
      it++;
    }
  }
  
  /** Parses the xml input file for parameter definitions for a single qmc simulation.
   * \param q the current xmlNode
   *
   Available parameters added to the ParameterSeet
   <ul>
   <li> blocks: number of blocks, default 100
   <li> steps: number of steps per block, default 1000
   <li> walkers: target population of walkers, default 100
   <li> Tau: the timestep, default 0.001
   <li> stride: flag for printing the ensemble of walkers,  default false
   <ul>
   <li> true: print every block
   <li> false: print at the end of the run
   </ul>
   </ul>
   In addition, sets the stride for the scalar estimators
   such that the scalar estimators flush and print to
   file every block and calls the function to initialize
   the walkers.
   *Derived classes can add their parameters.
   */
  bool 
  QMCDriver::putQMCInfo(xmlNodePtr cur) {
    
    int defaultw = 100;
    int targetw = 0;
    
    nTargetWalkers=0;
    if(cur) {

      xmlAttrPtr att = cur->properties;
      while(att != NULL) {
	string aname((const char*)(att->name));
	const char* vname=(const char*)(att->children->content);
	if(aname == "blocks") nBlocks = atoi(vname);
	else if(aname == "steps") nSteps = atoi(vname);
	att=att->next;
      }
      
      xmlNodePtr tcur=cur->children;
      m_param.put(cur);
      while(tcur != NULL) {
	string cname((const char*)(tcur->name));
	if(cname == "record") {
	  int stemp;
	  att = tcur->properties;
	  while(att != NULL) {
	    string aname((const char*)(att->name));
	    if(aname == "stride") stemp=atoi((const char*)(att->children->content));
	    att=att->next;
	  }
	  LogOut->getStream() << 
	    "Simulation: print walker ensemble every block." << endl;
	  if(stemp >= 0) pStride = true;
	}
	tcur=tcur->next;
      }
    }
    
    //set the stride for the scalar estimators 
    Estimators.setStride(nSteps);
    
    int nw  = W.getActiveWalkers();
    int ndiff = 0;
    if(nw) {
      ndiff = nTargetWalkers-nw;
    } else {
      ndiff=(nTargetWalkers)? nTargetWalkers:defaultw;
    }

    addWalkers(ndiff);

    //always true
    return (W.getActiveWalkers()>0);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
