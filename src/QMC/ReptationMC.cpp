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
#include "QMC/ReptationMC.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMC/MolecuFixedNodeBranch.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"
namespace ohmmsqmc {

  struct PolymerChain: public std::deque<MCWalkerConfiguration::Walker_t*> {

    typedef MCWalkerConfiguration::Walker_t Walker_t;

    ///Boolean for the direction of a move
    bool MoveHead;

    ///The number of beads to be cut
    int  NumCuts;

    double NumMoves;
    double AcceptedMoves;

    ///The index of the middle
    int Middle;    
    int Last;

    ///The walkers for the proposed move(s).
    std::vector<Walker_t*> heads; 

    ///The walkers that belong to the tails.
    std::vector<Walker_t*> tails;

    ///The walker repository to reuse Walkers.
    std::vector<Walker_t*> repository;

    /** constructor    
     *@param awalker the walker which is cloned to form a chain
     *@param len the number of chains in a polymer
     *@param movables the number of heads and tails that are moved
     */
    PolymerChain(Walker_t* awalker,int len, int movables): 
      MoveHead(true), NumCuts(movables),NumMoves(0.0),AcceptedMoves(0.0) {

      //always add number of beads
      if(len%2 == 0) len++;
      Middle = len/2;

      Last = len-1;
      for(int i=0; i<len; i++) {
	Walker_t* acopy=new Walker_t(*awalker);
	acopy->ID = i;
	push_back(acopy);
      }
     
      tails.resize(movables+1,0);
      heads.resize(movables,0);
      for(int i=0; i<2*movables; i++) {
	Walker_t* acopy=new Walker_t(*awalker);
        acopy->ID = i+len;
	repository.push_back(acopy);
      }
    }
    
    /** destructor
     *
     * Need to clean up the walkers in the repository and the polymer chain
     */
    ~PolymerChain() {
      for(int i=0; i<repository.size(); i++) delete repository[i];
      for(int i=0; i<size(); i++) delete (*this)[i];
    }

    inline size_t getID() const { return (*this)[Middle]->ID;}
    inline void subCuts() {
      if(NumCuts>1) NumCuts--;
    }

    inline void addCuts() {
      NumCuts++;
      if(heads.size()<NumCuts) {
	heads.resize(2*NumCuts);
	tails.resize(2*NumCuts);
	for(int i=heads.size(); i<2*NumCuts; i++) 
	  repository.push_back(new Walker_t(*(repository[0])));
      }
    }

    /** make tails and heads to make moves
     *@return the pointer to the anchor Walker
     *
     *\if MoveHead == true, 
     *tails are built from the end of the chain. The anchor is the first walker.
     *\else
     *tails are built from the start of the chain. The anchor is the last walker.
     *\endif
     *The heads are copied from the repository and the heads walkers will contain 
     *the data with the new configuration by the drift-and-diffusion.
     */
    inline Walker_t* makeEnds() {
      NumMoves+=1.0;
      for(int i=0; i<NumCuts; i++) heads[i]=repository[i];
      Walker_t* anchor = 0;
      if(MoveHead) {
	anchor=(*this)[0];
	for(int i=0, j=Last; i<NumCuts+1; i++,j--) {
	  tails[i]=(*this)[j];
	}
      } else {
	anchor=(*this)[Last];
	for(int i=0; i<NumCuts+1; i++) {
	  tails[i]=(*this)[i];
	}
      }
      return anchor;
    }
    
    inline void updateEnds() {
      AcceptedMoves+=1.0;
      if(MoveHead){
	for(int i=0; i<NumCuts; i++) {
	  push_front(heads[i]);
	  pop_back();
	}
      }else {
	for(int i=0; i<NumCuts; i++) {
	  push_back(heads[i]);
	  pop_front();
	}
      }
      //copy NumCuts of the tails to the repository for next step
      for(int i=0; i<NumCuts; i++) {repository[i] = tails[i];}
    }
    
    /** rather stupid average **/
    inline void average(Walker_t& center) {
      center.ID = (*this)[Middle]->ID;
      center.R = (*this)[Middle]->R;
      center.Properties = ((*this)[0]->Properties + (*this)[Last]->Properties)*0.5;
      center.Properties(LOCALPOTENTIAL) = (*this)[Middle]->Properties(LOCALPOTENTIAL);
    }

    inline void flip() {
      MoveHead = !MoveHead; //flip the direction
    }
  };

  class LocalPotentialEstimator {

    int Middle;
    int Counter;
    ofstream* fout;
    std::vector<PolymerChain*>* Polymers;
    std::vector<double> PEavg;
    std::vector<double> PE2;

  public:
  
    LocalPotentialEstimator(std::vector<PolymerChain*>* polymers):
      Middle(0), Counter(0), fout(0),Polymers(polymers) {  
      Middle=(*polymers)[0]->Middle;
      PEavg.resize(Middle+1,0.0);
      PE2.resize(Middle+1,0.0);
    }

    ~LocalPotentialEstimator() {
      if(fout) {delete fout;}
    }

    inline void resetReportSettings(const string& aname) {
      if(fout) {
	delete fout;
      }
      std::string filename(aname);
      filename.append(".pe.dat");
      fout = new ofstream(filename.c_str());
    }

    inline void reset() { 
      Counter = 0;
      for(int i=0; i<PEavg.size(); i++) PEavg[i]=0.0;
      for(int i=0; i<PE2.size(); i++) PE2[i]=0.0;
    }

    inline void report(int iter) {
      (*fout) << iter;
      double wgtinv = 1.0/static_cast<double>(Counter);
      for(int i=0; i<PEavg.size(); i++) {
	double eavg = PEavg[i]*wgtinv;
	(*fout)  <<  setw(15) << eavg;
      }
      (*fout) << endl;
      reset();
    }

    inline void accumulate() {
      Counter+=Polymers->size();
      register double e;
      for(int i=0; i<Polymers->size(); i++) {
	PolymerChain& curLink = *((*Polymers)[i]);
	for(int k=0,j=curLink.size()-1; k<Middle; k++,j--) {
	  e =0.5*(curLink[k]->Properties(LOCALPOTENTIAL)+curLink[j]->Properties(LOCALPOTENTIAL));
	  PEavg[k] += e;
	  PE2[k] += e*e;
	}
        e= curLink[Middle]->Properties(LOCALPOTENTIAL);
	PEavg[Middle] += e;
	PE2[Middle] += e*e;
      }
    }
  };


  ReptationMC::ReptationMC(MCWalkerConfiguration& w, 
			   TrialWaveFunction& psi, 
			   QMCHamiltonian& h, 
			   xmlNodePtr q): 
    QMCDriver(w,psi,h,q), 
    UseBounce(false),
    ClonePolymer(true),
    PolymerLength(21),
    NumCuts(1),
    NumTurns(0)
    { 
      RootName = "rmc";
      QMCType ="rmc";
      m_param.add(PolymerLength,"chains","int");
      m_param.add(NumCuts,"cuts","int");
      m_param.add(UseBounce,"bounce","int");
      m_param.add(ClonePolymer,"clone","int");
    }

  ReptationMC::~ReptationMC() {
    for(int i=0; i<Polymers.size(); i++) delete Polymers[i];
  }
  
  ///initialize polymers
  void ReptationMC::initPolymers() {
    
    //overwrite the number of cuts for Bounce algorithm
    if(UseBounce) NumCuts = 1;

    RealType g = sqrt(Tau);
    LOGMSG("Moving " << NumCuts << " for each reptation step")
    MCWalkerConfiguration::iterator it=W.begin();
    MCWalkerConfiguration::iterator it_end=W.end();
    while(it != it_end) {
      (*it)->Properties(WEIGHT)=1.0;     
      PolymerChain* achain = new PolymerChain((*it),PolymerLength,NumCuts);
      Polymers.push_back(achain);it++;

      Walker_t* cur=(*achain)[0];
      W.R = cur->R;
      DistanceTable::update(W);
      ValueType logpsi(Psi.evaluateLog(W));
      cur->Properties(LOCALENERGY) = H.evaluate(W);
      H.copy(cur->getEnergyBase());
      cur->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
      cur->Drift = W.G;

      if(!ClonePolymer) {
	
	for(int i=0; i<NumCuts-1; i++ ) {
	  //create a 3N-Dimensional Gaussian with variance=1
	  makeGaussRandom(deltaR);
	  W.R = cur->R + g*deltaR + Tau*cur->Drift;
	  
	  //update the distance table associated with W
	  DistanceTable::update(W);
	  
	  //evaluate wave function
	  ValueType logpsic(Psi.evaluateLog(W));
	  cur = (*achain)[i+1];	  
	  cur->Properties(LOCALENERGY) = H.evaluate(W);
	  H.copy(cur->getEnergyBase());
	  cur->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
	  cur->R = W.R;
	  cur->Drift = W.G;
	}
      }
      ++it;
    }
  }

  bool ReptationMC::run() { 

    //create a distance table for one walker
    DistanceTable::create(1);
    
    if(put(qmc_node)){
      
      //set the data members to start a new run
      //    getReady();
      int PopIndex, E_TIndex;
      Estimators.resetReportSettings(RootName);
      AcceptIndex = Estimators.addColumn("AcceptRatio");
      Estimators.reportHeader();

      //resize the temporary container
      deltaR.resize(W.getTotalNum());

      initPolymers();

      IndexType block = 0;
      Pooma::Clock timer;
      IndexType accstep=0;
      IndexType nAcceptTot = 0;
      IndexType nRejectTot = 0;
      
      LocalPotentialEstimator pe(&Polymers);
      pe.resetReportSettings(RootName);

      //accumulate configuration: probably need to reorder
      HDFWalkerOutput WO(RootName);
      do {

	IndexType step = 0;
	timer.start();
	NumTurns = 0;

	do {
          movePolymers();
	  step++; accstep++;

	  MCWalkerConfiguration::iterator it=W.begin();
	  MCWalkerConfiguration::iterator it_end=W.end();
	  int ilink=0;
	  while(it != it_end) {
	    Polymers[ilink]->average(**it);
	    ++ilink; ++it;
	  }

	  Estimators.accumulate(W);
	  pe.accumulate();

// 	  cout << step << " ";
// 	  for(int i=Polymers[0]->Middle-1; i>=0 ; i--)
// 	    cout << Polymers[0]->PEavg[i] << " ";
// 	  cout << endl;
	} while(step<nSteps);
	timer.stop();
	
	nAcceptTot += nAccept;
	nRejectTot += nReject;

        RealType acceptedR = static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject); 
	Estimators.flush();
	Estimators.setColumn(AcceptIndex,acceptedR);
	Estimators.report(accstep);
	pe.report(accstep);

        //change NumCuts to make accstep ~ 50%
	LogOut->getStream() 
	  << "Block " << block << " " 
	  << timer.cpu_time() << " " << NumTurns << " " << Polymers[0]->getID() << endl;

	nAccept = 0; nReject = 0;
	block++;

	if(pStride) WO.get(W);

      } while(block<nBlocks);
      
      LogOut->getStream() 
	<< "ratio = " 
	<< static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
	<< endl;

      Estimators.finalize();
      return true;
    } else 
      return false;
  }

  bool 
  ReptationMC::put(xmlNodePtr q){
    xmlNodePtr qsave=q;
    bool success = putQMCInfo(q);
    success = Estimators.put(qsave);
    return success;
  }
  
  void 
  ReptationMC::movePolymers(){
    
    //Pooma::Clock timer;
    //RealType oneovertau = 1.0/Tau;
    //RealType oneover2tau = 0.5*oneovertau;
    RealType tauover2 = 0.5*Tau;
    RealType g = sqrt(Tau);
    
    typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;
    typedef MCWalkerConfiguration::Walker_t Walker_t;

    for(int ilink=0; ilink<Polymers.size(); ilink++) {

      PolymerChain& polymer = *(Polymers[ilink]);

      if(!UseBounce && Random()<0.5) {
	polymer.flip(); 	  
	NumTurns++;
      }

      Walker_t* anchor = polymer.makeEnds();
      
      //save the local energies of the anchor and tails
      //eloc_xp = the energy of the front
      //eloc_yp = the energy of the proposed move
      //eloc_x = the energy of the tail
      //eloc_y = the energy of the tail-1
      RealType eloc_xp=anchor->Properties(LOCALENERGY);
      RealType eloc_x = polymer.tails[0]->Properties(LOCALENERGY);
      RealType eloc_y = polymer.tails[1]->Properties(LOCALENERGY);

      NumCuts = polymer.NumCuts;
      RealType Wpolymer=0.0;

      for(int i=0; i<NumCuts; ) {

	Walker_t* head=polymer.heads[i];

	//create a 3N-Dimensional Gaussian with variance=1
	makeGaussRandom(deltaR);
	W.R = anchor->R + g*deltaR + Tau* anchor->Drift;

	//update the distance table associated with W
	DistanceTable::update(W);

	//evaluate wave function
	ValueType logpsi(Psi.evaluateLog(W));
	
	//update the properties of the front chain
	RealType eloc_yp = head->Properties(LOCALENERGY) = H.evaluate(W);
	H.copy(head->getEnergyBase());
	head->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
	head->R = W.R;

	//ValueType vsq = Dot(W.G,W.G);
	//ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
	//head->Drift = scale*W.G;
	head->Drift = W.G;

	//\f${x-y-\tau\nabla \ln \Psi_{T}(y))\f$
	//deltaR = anchor->R - W.R - heads[i]->Drift;
	//Gdrift *= exp(-oneover2tau*Dot(deltaR,deltaR));
	/* 
	   \f$ X= \{R_0, R_1, ... , R_M\}\f$
	   \f$ X' = \{R_1, .., R_M, R_{M+1}\}\f$
	   \f[ G_B(R_{M+1}\leftarrow R_{M}, \tau)/G_B(R_{0}\leftarrow R_{1}, \tau)
	   = exp\(-\tau/2[E_L(R_{M+1})+E_L(R_M)-E_L(R_1)-E_L(R_0)]\)\f]
	   *
	   -  eloc_yp = \f$E_L(R_{M+1})\f$
	   -  eloc_xp = \f$E_L(R_{M})\f$
	   -  eloc_y = \f$E_L(R_{1})\f$
	   -  eloc_x = \f$E_L(R_{0})\f$
	*/
	//Wpolymer *= exp(-oneover2tau*(eloc_yp+eloc_xp-eloc_x-eloc_y));
	Wpolymer +=(eloc_yp+eloc_xp-eloc_x-eloc_y);

	//move the anchor and swap the local energies for Wpolymer
	anchor=head;

	//increment the index
	i++;
	if(i<NumCuts) {
	  eloc_xp  = eloc_yp;
	  eloc_x = eloc_y;
	  eloc_y = polymer.tails[i+1]->Properties(LOCALENERGY);
	}
      }
      
      Wpolymer = exp(-tauover2*Wpolymer);
      double accept = std::min(1.0,Wpolymer);
      if(Random() < accept){//move accepted
	polymer.updateEnds();
	++nAccept;
      } else {
	++nReject; 
	if(UseBounce && Random()>accept) {
	  NumTurns++;
	  polymer.flip();
	}
      }

//       RealType Bounce =  UseBounce ? 1.0-accept: 0.5;
//       if(Random()<Bounce) {
// 	polymer.flip();
// 	LogOut->getStream() << "Bounce = " << Bounce << " " << NumTurns << " " << polymer.MoveHead << endl;
// 	NumTurns++;//increase the number of turns
//       }
    }
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
