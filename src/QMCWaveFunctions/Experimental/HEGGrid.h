/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "Lattice/CrystalLattice.h"
#include <map>

namespace qmcplusplus {

  template<class T, unsigned D>
  struct HEGGrid { };

  //three-d specialization
  template<class T>
  struct HEGGrid<T,3> {
    typedef CrystalLattice<T,3> PL_t;
    typedef typename PL_t::SingleParticlePos_t PosType;

    ///number of kpoints of a half sphere excluding gamma
    int NumKptsHalf;
    ///maxmim ksq
    T MaxKsq;
    PL_t& Lattice;
    map<int,vector<PosType>*> rs;
    vector<PosType> kpt;
    vector<T>  mk2;

    HEGGrid(PL_t& lat):Lattice(lat) { }

    ~HEGGrid() { 
      typename map<int,vector<PosType>*>::iterator it(rs.begin()),
        it_end(rs.end());
      while(it != it_end) { delete (*it).second; ++it;}
    }

    /** return the estimated number of grid in each direction */
    inline int getNC(int nup) {
      return static_cast<int>(std::pow(static_cast<T>(nup),1.0/3.0))/2+1;
    }

    void sortGrid(int nc) {
      int first_ix2, first_ix3; 
      for(int ix1=0; ix1<=nc; ix1++) {
        if(ix1 == 0) first_ix2=0;
        else         first_ix2=-nc;
        for(int ix2=first_ix2; ix2<=nc; ix2++) {
          if(ix1 == 0 && ix2 == 0) first_ix3=1;
          else                     first_ix3=-nc;
          for(int ix3=first_ix3; ix3<=nc; ix3++) {
            int ih=ix1*ix1+ix2*ix2+ix3*ix3;
            typename std::map<int,vector<PosType>*>::iterator it = rs.find(ih);
            if(it == rs.end()) {
              vector<PosType>* ns = new vector<PosType>;
              ns->push_back(PosType(ix1,ix2,ix3));
              rs[ih] = ns;
            } else {
              (*it).second->push_back(PosType(ix1,ix2,ix3));
            }
          }
        }
      }
    }

    void createGrid(int nc, int nkpts) {

      if(rs.empty()) sortGrid(nc);

      NumKptsHalf=nkpts;

      kpt.resize(nkpts);
      mk2.resize(nkpts);

      int ikpt=0;
      //int checkNum=0;
      //int ke=0;
      MaxKsq=0.0;

      int rsbin=0;
      typename map<int, vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
      while(ikpt<nkpts && rs_it != rs_end) {
        typename vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
        T minus_ksq=-Lattice.ksq(*ns_it);
        while(ikpt<nkpts && ns_it!=ns_end) {
          kpt[ikpt]=Lattice.k_cart(*ns_it);
          mk2[ikpt]=minus_ksq;
          ++ikpt;
          ++ns_it;
        }
        ++rsbin;
        ++rs_it;
      }

      MaxKsq = Lattice.ksq(*((*rs_it).second->begin()));

      app_log() << "List of kpoints (half-sphere) " << endl;
      for(int ik=0; ik<kpt.size(); ik++) {
        app_log() << ik << " " << kpt[ik] << " " << mk2[ik] << endl;
      }
    }


    void createGrid(const PosType& twistAngle) {

      //unfold and add gamma
      int nkpts = 2*NumKptsHalf+1;
      kpt.resize(nkpts);
      mk2.resize(nkpts);

      //add gamma
      int ikpt=0;
      kpt[ikpt]=Lattice.k_cart(twistAngle);
      mk2[ikpt]=-Lattice.ksq(twistAngle);
      
      ++ikpt;
      typename map<int, vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
      while(ikpt<nkpts && rs_it != rs_end) {
        typename vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
        while(ikpt<nkpts && ns_it!=ns_end) {
          //add twist+k
          PosType k0(twistAngle+(*ns_it));
          T ksq=Lattice.ksq(k0);
          kpt[ikpt]=Lattice.k_cart(k0);
          mk2[ikpt]=-ksq;
          ++ikpt;
          //add twist-k
          k0=twistAngle-(*ns_it);
          ksq=Lattice.ksq(k0);
          kpt[ikpt]=Lattice.k_cart(k0);
          mk2[ikpt]=-ksq;
          ++ikpt;
          ++ns_it;
        }
        ++rs_it;
      }

      app_log() << "List of kpoints with twist = " << twistAngle << endl;
      for(int ik=0; ik<kpt.size(); ik++) {
        app_log() << ik << " " << kpt[ik] << " " <<-mk2[ik] << endl;
      }
    }
  };
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
