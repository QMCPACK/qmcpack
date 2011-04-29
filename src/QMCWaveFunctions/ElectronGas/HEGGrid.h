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
    vector<int> n_within_shell;

    HEGGrid(PL_t& lat):Lattice(lat) 
    { 
      n_within_shell.resize(31);
      n_within_shell[0]=1; n_within_shell[1]=7; n_within_shell[2]=19;
      n_within_shell[3]=27; n_within_shell[4]=33; n_within_shell[5]=57;
      n_within_shell[6]=81; n_within_shell[7]=93; n_within_shell[8]=123;
      n_within_shell[9]=147; n_within_shell[10]=171; n_within_shell[11]=179;
      n_within_shell[12]=203; n_within_shell[13]=251; n_within_shell[14]=257;
      n_within_shell[15]=305; n_within_shell[16]=341; n_within_shell[17]=365;
      n_within_shell[18]=389; n_within_shell[19]=437; n_within_shell[20]=461;
      n_within_shell[21]=485; n_within_shell[22]=515; n_within_shell[23]=587;
      n_within_shell[24]=619; n_within_shell[25]=691; n_within_shell[26]=739;
      n_within_shell[27]=751; n_within_shell[28]=799; n_within_shell[29]=847;
      n_within_shell[30]=895;
    }

    ~HEGGrid() { 
      typename map<int,vector<PosType>*>::iterator it(rs.begin()),
        it_end(rs.end());
      while(it != it_end) { delete (*it).second; ++it;}
    }

    /** return the estimated number of grid in each direction */
    inline int getNC(int nup) const {
      return static_cast<int>(std::pow(static_cast<T>(nup),1.0/3.0))/2+1;
    }

    //return the number of k-points upto nsh-shell
    inline int getNumberOfKpoints(int nsh) const
    {
      if(nsh<n_within_shell.size())
        return n_within_shell[nsh];
      else
        return -1;
    }

    //return the shell index for nkpt k-points
    inline int getShellIndex(int nkpt) const
    {
      vector<int>::const_iterator loc=std::upper_bound(n_within_shell.begin(),n_within_shell.end(),nkpt);
      if(loc<n_within_shell.end()) 
        return loc-n_within_shell.begin()-1;
      else
        return getNC(nkpt);
    }

    /** return the cell size  for the number of particles and rs
     * @param nptcl number of particles
     * @param rs_in rs
     */
    inline T getCellLength(int nptcl, T rs_in) const
    {
      return std::pow(4.0*M_PI*nptcl/3.0,1.0/3.0)*rs_in;
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


    void createGrid(int nc, int nkpts, const PosType& twistAngle) 
    {
      std::map<int,vector<PosType>*> rs_big;
      for(int ix1=-nc-1; ix1<=nc+1; ++ix1) 
      {
        for(int ix2=-nc-1; ix2<=nc+1; ++ix2) 
        {
          for(int ix3=-nc-1; ix3<=nc+1; ++ix3) 
          {
            PosType k0(ix1+twistAngle[0],ix2+twistAngle[1],ix3+twistAngle[2]);
            int ih=static_cast<int>(1e4*Lattice.ksq(k0));
            typename std::map<int,vector<PosType>*>::iterator it = rs_big.find(ih);
            if(it == rs_big.end()) 
            {
              vector<PosType>* ns = new vector<PosType>;
              ns->push_back(k0);
              rs_big[ih] = ns;
            } else {
              (*it).second->push_back(k0);
            }
          }
        }
      }

      kpt.resize(nkpts);
      mk2.resize(nkpts);

      typename map<int, vector<PosType>*>::iterator rs_it(rs_big.begin()), rs_end(rs_big.end());
      int ikpt=0;
      while(ikpt<nkpts && rs_it != rs_end) 
      {
        typename vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
        while(ikpt<nkpts && ns_it!=ns_end) 
        {
          //add twist+k
          PosType k0(*ns_it);
          T ksq=Lattice.ksq(k0);
          kpt[ikpt]=Lattice.k_cart(k0);
          mk2[ikpt]=-ksq;
          ++ikpt;
          ++ns_it;
        }
        ++rs_it;
      }

      app_log() << "List of kpoints with twist = " << twistAngle << endl;
      for(int ik=0; ik<kpt.size(); ik++) 
      {
        app_log() << ik << " " << kpt[ik] << " " <<-mk2[ik] << endl;
      }
    }
  };

  //three-d specialization
  template<class T>
  struct HEGGrid<T,2> {
    typedef CrystalLattice<T,2> PL_t;
    typedef typename PL_t::SingleParticlePos_t PosType;

    ///number of kpoints of a half sphere excluding gamma
    int NumKptsHalf;
    ///maxmim ksq
    T MaxKsq;
    PL_t& Lattice;
    map<int,vector<PosType>*> rs;
    vector<PosType> kpt;
    vector<T>  mk2;
    vector<int> n_within_shell;

    HEGGrid(PL_t& lat):Lattice(lat) 
    { 
      n_within_shell.resize(31);
      //fill this in
      n_within_shell[0]=1; 
    }

    ~HEGGrid() { 
      typename map<int,vector<PosType>*>::iterator it(rs.begin()),
        it_end(rs.end());
      while(it != it_end) { delete (*it).second; ++it;}
    }

    /** return the estimated number of grid in each direction */
    inline int getNC(int nup) {
      return static_cast<int>(std::pow(static_cast<T>(nup),1.0/2))/2+1;
    }

    //return the number of k-points upto nsh-shell
    inline int getNumberOfKpoints(int nsh) const
    {
      return -1;
    }

    //return the shell index for nkpt k-points
    inline int getShellIndex(int nkpt) const
    {
      return -1;
    }

    inline T getCellLength(int nptcl, T rs_in)
    {
      return std::sqrt(M_PI*nptcl)*rs_in;
    }

    void sortGrid(int nc) {
      int first_ix2, first_ix3; 
      for(int ix1=0; ix1<=nc; ix1++) {
        if(ix1 == 0) first_ix2=0;
        else         first_ix2=-nc;
        for(int ix2=first_ix2; ix2<=nc; ix2++) {
          int ih=ix1*ix1+ix2*ix2;
          typename std::map<int,vector<PosType>*>::iterator it = rs.find(ih);
          if(it == rs.end()) {
            vector<PosType>* ns = new vector<PosType>;
            ns->push_back(PosType(ix1,ix2));
            rs[ih] = ns;
          } else {
            (*it).second->push_back(PosType(ix1,ix2));
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

      app_log() << "Check this " << NumKptsHalf << " " << nkpts << endl;
      abort();

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
