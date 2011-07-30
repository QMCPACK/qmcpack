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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file InitMolecularSystem.cpp
 * @brief Implements InitMolecuarSystem operators.
 */
#include "QMCApp/InitMolecularSystem.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus {

  InitMolecularSystem::InitMolecularSystem(ParticleSetPool* pset, 
      const char* aname):
    OhmmsElementBase(aname), ptclPool(pset){ }

  bool InitMolecularSystem::put(xmlNodePtr cur) {

    string target("e"), source("i");

    OhmmsAttributeSet hAttrib;
    hAttrib.add(target,"target");
    hAttrib.add(source,"source");
    hAttrib.put(cur);

    ParticleSet* els=ptclPool->getParticleSet(target);
    if(els == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    ParticleSet* ions=ptclPool->getParticleSet(source);
    if(ions == 0) {
      ERRORMSG("No target particle " << target << " exists.")
      return false;
    }

    if(ions->getTotalNum()>1) {
      initMolecule(ions,els);
    } else {
      initAtom(ions,els);
    }
     
    return true;
  }

  void InitMolecularSystem::initAtom(ParticleSet* ions, ParticleSet* els) 
  {

    //3N-dimensional Gaussian
    ParticleSet::ParticlePos_t chi(els->getTotalNum());
    makeGaussRandom(chi);

    double q=std::sqrt(static_cast<double>(els->getTotalNum()))*0.5;
    int nel(els->getTotalNum()), items(0);
    while(nel) {
      els->R[items]=ions->R[0]+q*chi[items]; 
      --nel; ++items;
    }
  }

  struct LoneElectron 
  {
    int ID;
    double BondLength;
    inline LoneElectron(int id, double bl):ID(id),BondLength(bl){}
  };
  
  void InitMolecularSystem::initMolecule(ParticleSet* ions, ParticleSet* els) 
  {
    
    if(ions->getTotalNum()==1) return initAtom(ions,els);

    DistanceTableData* d_ii = DistanceTable::add(*ions);
    d_ii->create(1);
    d_ii->evaluate(*ions);

    const ParticleSet::ParticleIndex_t& grID(ions->GroupID);
    SpeciesSet& Species(ions->getSpeciesSet());

    int Centers = ions->getTotalNum();
    vector<int> Qtot(Centers), Qcore(Centers), Qval(Centers,0);
 
    //use charge as the core electrons first
    int icharge = Species.addAttribute("charge");

    //Assign default core charge
    for(int iat=0; iat<Centers; iat++) 
      Qtot[iat] = static_cast<int>(Species(icharge,grID[iat]));

    //cutoff radius (Bohr) this a random choice    
    double cutoff = 4.0;
    ParticleSet::ParticlePos_t chi(els->getTotalNum());
    //makeGaussRandom(chi);
    makeSphereRandom(chi);

    int numUp =els->last(0);
    int numDown = els->last(1)-els->first(0);
    int item = 0;

    int nup_tot=0, ndown_tot=numUp;

    vector<LoneElectron> loneQ;

    double rmin=cutoff; 
    for(int iat=0; iat<Centers; iat++) { 
      for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++){
        rmin = std::min(rmin,d_ii->r(nn));
      }
      //use 40% of the minimum bond
      double sep=rmin*0.4;
      int v2=Qtot[iat]/2;
      if(Qtot[iat]>v2*2) {
        loneQ.push_back(LoneElectron(iat,sep));
      }
      for(int k=0; k<v2; k++) {
        els->R[nup_tot++]=ions->R[iat]+sep*chi[item++];
        els->R[ndown_tot++]=ions->R[iat]+sep*chi[item++];
      }
    }

    vector<LoneElectron>::iterator it(loneQ.begin());
    vector<LoneElectron>::iterator it_end(loneQ.end());
    while(nup_tot<numUp && it != it_end) {
      els->R[nup_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
      ++it;
    }
    while(it != it_end) {
      els->R[ndown_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
      ++it;
    }

    //put all the electrons in a unit box
    if(els->Lattice.SuperCellEnum != SUPERCELL_OPEN) els->applyBC(els->R);

    /*
    //Overwrite the valence charge
    int iz = Species.addAttribute("atomicnumber");
    int nattrib=Species.numAttributes();
    int ival=Species.addAttribute("valence");
    //store the max distance from atom
    if(ival<nattrib) {
      for(int iat=0; iat<Centers; iat++) {
        Qval[iat]=static_cast<int>(Species(ival,grID[iat]));
	//Probably charge is missing: valence becomes charge
	if(Qval[iat]>Qtot[iat]) Qtot[iat]=Qval[iat];
	Qcore[iat]=Qtot[iat]-Qval[iat];
      }
    } else {
      for(int iat=0; iat<Centers; iat++) {
	Qcore[iat] = Qtot[iat];
	Qval[iat] = 0;
      }
    }

    //3N-dimensional Gaussian
    //assign the core
    int ncoreUp(0), items(0);
    for(int iat=0; iat<Centers; iat++) {
      double sep=sqrt(static_cast<double>(Qcore[iat]))*0.5;
      for(int iel=0; iel<Qcore[iat]/2; iel++) {
        els->R[ncoreUp]=ions->R[iat]+sep*chi[items++];
        els->R[ncoreUp+numUp]=ions->R[iat]+sep*chi[items++];
        ++ncoreUp;
      }
    }

    int vup=numUp-ncoreUp;
    int vdown=numDown-ncoreUp;
    int vtot=vup+vdown;
    int valIndex = ncoreUp;
    int ic=0;
    while(vtot && ic<Centers) {
      for(int nn=d_ii->M[ic]; nn<d_ii->M[ic+1]; nn++){
        double bl = d_ii->r(nn);
        int jc = d_ii->J[nn];
        //only assign if the half bond-length < cutoff
	//assign pairs of electrons (up and down)
        if(vtot && bl < 2.0*cutoff){
          ParticleSet::SingleParticlePos_t displ= ions->R[ic]+0.5*d_ii->dr(nn);
          bl*=0.1;
          if(vup) {
            els->R[ncoreUp] = displ+bl*chi[items];
	    --vup;--vtot;++items;
          }
          if(vdown) {
            els->R[ncoreUp+numUp]=displ+bl*chi[items];
            --vdown;--vtot;++items;
          }
	  ++ncoreUp;
        }
      }
      ++ic;
    }
    */
  }

  bool InitMolecularSystem::put(std::istream& is) {
    return true;
  }

  bool InitMolecularSystem::get(std::ostream& os) const {
    return true;
  }

  void InitMolecularSystem::reset() {
 
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
