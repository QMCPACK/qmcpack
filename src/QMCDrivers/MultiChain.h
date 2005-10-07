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
#ifndef OHMMS_QMC_MULTICHAIN_H
#define OHMMS_QMC_MULTICHAIN_H
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include <deque>
namespace ohmmsqmc {

  template<class IT>
  inline void delete_iter(IT first, IT last) {
    while(first != last) { delete *first; ++first;}
  }

  struct Bead: public MCWalkerConfiguration::Walker_t{

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::RealType RealType;
    typedef MCWalkerConfiguration::ParticlePos_t ParticlePos_t;
 
    Vector<int> BeadSignWgt;    
    //Vector<ParticlePos_t*> Gradients;
    vector<ParticlePos_t*> Gradients;
    Matrix<RealType> Action;
    RealType TransProb[2];

    inline Bead(const Bead& a) {
      makeCopyBead(a); 
    }

    inline Bead(const Walker_t& a){
      makeCopy(a);
      int rows=Properties.rows();
      Resize_Grad_and_Action(rows,R.size());
      BeadSignWgt.resize(rows);
    }

    inline Bead& operator=(const Bead& a) {
      makeCopyBead(a);
      return *this;
    }

    inline void makeCopyBead(const Bead& a){
      makeCopy(a);
      int rows=a.Gradients.size();
      Resize_Grad_and_Action(rows,a.size());
      Action=a.Action;
      for(int i=0; i<rows; i++) *Gradients[i] = *(a.Gradients[i]);
      BeadSignWgt.resize(rows);
      BeadSignWgt=a.BeadSignWgt;
      TransProb[0]=a.TransProb[0];
      TransProb[1]=a.TransProb[1];
    }

    inline void Resize_Grad_and_Action(int n, int m){
      int curg=Gradients.size();
      while(curg<n) {
        Gradients.push_back(new ParticlePos_t(m));
        ++curg;
      }
      Action.resize(n,3);
    }
  };

  struct MultiChain: public std::deque<Bead*> {

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::RealType RealType;

    /// Direction of growth
    int GrowthDirection;

    ///The index of the middle
    int Middle;    

    ///The index of the Last
    int Last;

    ///The number of H/Psi pairs
    int nPsi;

    RealType GlobalWgt;

    Vector<RealType> GlobalAction,UmbrellaWeight;
    Vector<int>      GlobalSignWgt,RefSign;

    /**  constructor    
     * @param abead Bead used to generate len Beads to form a chain
     * @param len size of the chain
     * @param direction initial growth direction
     * @param npsi number of Psi/H pairs
     */
    //MultiChain(Walker_t* awalker,int len, int direction, int npsi): 
    MultiChain(Bead* abead,int len, int direction, int npsi): 
      GrowthDirection(direction), nPsi(npsi){
      //always add number of beads
      if(len%2 == 0) len++;
      Middle = len/2;
      Last = len-1;
      for(int i=0; i<len; i++) {
        push_back(new Bead(*abead));
      }
      GlobalAction.resize(npsi);   GlobalAction=0.0;
      UmbrellaWeight.resize(npsi); UmbrellaWeight=1.0;
      GlobalSignWgt.resize(npsi);  GlobalSignWgt=0;
      RefSign.resize(npsi); RefSign=0;
    }
    


    /** destructor
     *
     * Need to clean up the walkers in the repository and the polymer chain
     */
    ~MultiChain() {
      delete_iter(this->begin(),this->end());
    }



    inline void flip(){ 
     GrowthDirection = abs(GrowthDirection-1); //flip the direction
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
