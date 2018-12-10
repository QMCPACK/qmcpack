//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_HEGGRID_H
#define QMCPLUSPLUS_HEGGRID_H

#include "Lattice/CrystalLattice.h"
#include <map>

namespace qmcplusplus
{

template<typename T, unsigned D>
struct kpdata{
  TinyVector<T,D> k;
  T k2;
  int g;
};


template<typename T, unsigned D>
bool kpdata_comp(const kpdata<T,D>& left, const kpdata<T,D>& right)
{
  return left.k2 < right.k2;
}


template<class T, unsigned D>
struct HEGGrid { };



//three-d specialization
template<class T>
struct HEGGrid<T,3>
{
  typedef CrystalLattice<T,3> PL_t;
  typedef typename PL_t::SingleParticlePos_t PosType;
  typedef typename PL_t::Scalar_t            RealType;

  ///number of kpoints of a half sphere excluding gamma
  int NumKptsHalf;
  ///maxmim ksq
  T MaxKsq;
  PL_t& Lattice;
  std::map<int,std::vector<PosType>*> rs;


  std::vector<PosType> kpt;
  std::vector<T>       mk2;
  std::vector<int>     deg;
  std::vector<int> n_within_shell;
  PosType twist;


  typedef kpdata<T,3> kpdata_t;
  typedef std::vector<kpdata_t > kpoints_t;

  kpoints_t* kpoints_grid;
  int nctmp;


  HEGGrid(PL_t& lat):Lattice(lat),twist(0.0), nctmp(-1), kpoints_grid(0)
  {
    n_within_shell.resize(31);
    n_within_shell[0]=1;
    n_within_shell[1]=7;
    n_within_shell[2]=19;
    n_within_shell[3]=27;
    n_within_shell[4]=33;
    n_within_shell[5]=57;
    n_within_shell[6]=81;
    n_within_shell[7]=93;
    n_within_shell[8]=123;
    n_within_shell[9]=147;
    n_within_shell[10]=171;
    n_within_shell[11]=179;
    n_within_shell[12]=203;
    n_within_shell[13]=251;
    n_within_shell[14]=257;
    n_within_shell[15]=305;
    n_within_shell[16]=341;
    n_within_shell[17]=365;
    n_within_shell[18]=389;
    n_within_shell[19]=437;
    n_within_shell[20]=461;
    n_within_shell[21]=485;
    n_within_shell[22]=515;
    n_within_shell[23]=587;
    n_within_shell[24]=619;
    n_within_shell[25]=691;
    n_within_shell[26]=739;
    n_within_shell[27]=751;
    n_within_shell[28]=799;
    n_within_shell[29]=847;
    n_within_shell[30]=895;
  }

  ~HEGGrid()
  {
    typename std::map<int,std::vector<PosType>*>::iterator it(rs.begin()),
             it_end(rs.end());
    while(it != it_end)
    {
      delete (*it).second;
      ++it;
    }
    clear_kpoints();
  }

  /** return the estimated number of grid in each direction */
  inline int getNC(int nup) const
  {
    return static_cast<int>(std::pow(static_cast<T>(nup),1.0/3.0))/2+1;
  }

  /** return the estimated number of grid in each direction (upper bound) */
  inline int get_nc(int nstates) const
  {
    return static_cast<int>(std::pow(static_cast<T>(nstates),1.0/3.0)*.7)+1;
  }

  //return the number of k-points upto nsh-shell
  inline int getNumberOfKpoints(int nsh) const
  {
    if(nsh<n_within_shell.size())
      return n_within_shell[nsh];
    else
      return -1;
  }


  inline int getShellFromStates(int nst)
  {
    for(int i=0; i<n_within_shell.size(); i++)
      if (n_within_shell[i]==nst)
        return i;
    return -1;
  }

  //return the shell index for nkpt k-points
  inline int getShellIndex(int nkpt) const
  {
    std::vector<int>::const_iterator loc=std::upper_bound(n_within_shell.begin(),n_within_shell.end(),nkpt);
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

  void sortGrid(int nc)
  {
    int first_ix2, first_ix3;
    for(int ix1=0; ix1<=nc; ix1++)
    {
      if(ix1 == 0)
        first_ix2=0;
      else
        first_ix2=-nc;
      for(int ix2=first_ix2; ix2<=nc; ix2++)
      {
        if(ix1 == 0 && ix2 == 0)
          first_ix3=1;
        else
          first_ix3=-nc;
        for(int ix3=first_ix3; ix3<=nc; ix3++)
        {
          int ih=ix1*ix1+ix2*ix2+ix3*ix3;
          typename std::map<int,std::vector<PosType>*>::iterator it = rs.find(ih);
          if(it == rs.end())
          {
            std::vector<PosType>* ns = new std::vector<PosType>;
            ns->push_back(PosType(ix1,ix2,ix3));
            rs[ih] = ns;
          }
          else
          {
            (*it).second->push_back(PosType(ix1,ix2,ix3));
          }
        }
      }
    }
  }

  void createGrid(int nc, int nkpts)
  {
    if(rs.empty())
      sortGrid(nc);
    NumKptsHalf=nkpts;
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    int ikpt=0;
    //int checkNum=0;
    //int ke=0;
    MaxKsq=0.0;
    int rsbin=0;
    typename std::map<int, std::vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
    while(ikpt<nkpts && rs_it != rs_end)
    {
      typename std::vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
      T minus_ksq=-Lattice.ksq(*ns_it);
      while(ikpt<nkpts && ns_it!=ns_end)
      {
        kpt[ikpt]=Lattice.k_cart(*ns_it);
        mk2[ikpt]=minus_ksq;
        ++ikpt;
        ++ns_it;
      }
      ++rsbin;
      ++rs_it;
    }
    MaxKsq = Lattice.ksq(*((*rs_it).second->begin()));
    app_log() << "List of kpoints (half-sphere) " << std::endl;
    for(int ik=0; ik<kpt.size(); ik++)
    {
      app_log() << ik << " " << kpt[ik] << " " << mk2[ik] << std::endl;
    }
  }


  void clear_kpoints()
  {
    if(kpoints_grid!=0)
      delete kpoints_grid;
  }


  void create_kpoints(int nc, const PosType& tw, T tol=1e-6)
  {
    if(kpoints_grid==0)
      kpoints_grid = new kpoints_t;
    else if(nc<=nctmp)
      return;
    nctmp = nc;
    kpoints_t& kpoints = *kpoints_grid;
    app_log()<<"  resizing kpoint grid"<< std::endl;
    app_log()<<"  current size = "<<kpoints.size()<< std::endl;
    // make space for the kpoint grid
    int nkpoints = pow( 2*(nc+1)+1 , 3 );
    kpoints.resize(nkpoints);
    app_log()<<"  cubic size = "<<kpoints.size()<< std::endl;
    typename kpoints_t::iterator kptmp,kp=kpoints.begin(),kp_end=kpoints.end();
    // make the kpoint grid
    T k2max = std::numeric_limits<RealType>::max();
    for(int i0=-nc-1; i0<=nc+1; ++i0)
      for(int i1=-nc-1; i1<=nc+1; ++i1)
        for(int i2=-nc-1; i2<=nc+1; ++i2)
        {
          PosType k(i0+tw[0],i1+tw[1],i2+tw[2]);
          kp->k  = Lattice.k_cart(k);
          kp->k2 = Lattice.ksq(k);
          if(std::abs(i0)==(nc+1) || std::abs(i1)==(nc+1) || std::abs(i2)==(nc+1))
            k2max = std::min(k2max,kp->k2);
          ++kp;
        }
    // sort kpoints by magnitude
    sort(kpoints.begin(),kpoints.end(),kpdata_comp<T,3>);
    // eliminate kpoints outside of inscribing sphere
    int nkp = 0;
    kp = kpoints.begin();
    while(kp!=kp_end && kp->k2<k2max+1e-3)
    {
      nkp++;
      ++kp;
    }
    kpoints.resize(nkp);
    app_log()<<"  new spherical size = "<<kpoints.size()<< std::endl;
    kp_end = kpoints.end();
    // count degeneracies
    kp = kpoints.begin();
    while(kp!=kp_end)
    {
      T k2 = kp->k2;
      kptmp=kp;
      int g=1;
      ++kptmp;
      // look ahead to count
      while(kptmp!=kp_end && std::abs(kptmp->k2-k2)<tol)
      {
        g++;
        ++kptmp;
      }
      kp->g = g;
      // run over degenerate states to assign
      for(int n=0;n<g-1;++n)
        (++kp)->g = g;
      ++kp;
    }
    //app_log()<<"create_kpoints"<< std::endl;
    //app_log()<<"  nkpoints = "<<nkpoints<< std::endl;
    //app_log()<<"  kpoints"<< std::endl;
    //for(kp=kpoints.begin();kp!=kp_end;++kp)
    //  app_log()<<"    "<<kp->k2<<" "<<kp->g<<" "<<kp->k<< std::endl;
    //APP_ABORT("end create_kpoints");
  }


  void createGrid(int nc, int nkpts, const PosType& twistAngle)
  {
    twist = twistAngle;
    create_kpoints(nc,twistAngle);
    kpoints_t& kpoints = *kpoints_grid;
    if(nkpts>kpoints.size())
      APP_ABORT("HEGGrid::createGrid  requested more kpoints than created");
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    deg.resize(nkpts);
    for(int i=0;i<nkpts;++i)
    {
      const kpdata_t& kp = kpoints[i];
      kpt[i] =  kp.k;
      mk2[i] = -kp.k2;
      deg[i] =  kp.g;
    }
    app_log() << "List of kpoints with twist = " << twistAngle << std::endl;
    for(int ik=0; ik<kpt.size(); ik++)
      app_log() << ik << " " << kpt[ik] << " " <<-mk2[ik] << std::endl;
  }



  void createGrid(const std::vector<int>& states, T tol=1e-6)
  {
    createGrid(states,twist,tol);
  }


  void createGrid(const std::vector<int>& states,const PosType& twistAngle, T tol=1e-6)
  {
    int smax=0;
    for(int i=0;i<states.size();++i)
      smax = std::max(smax,states[i]);
    smax++;
    create_kpoints(get_nc(smax),twistAngle,tol);
    kpoints_t& kpoints = *kpoints_grid;
    if(smax>kpoints.size())
      APP_ABORT("HEGGrid::createGrid(states)  requested more kpoints than created");
    int nkpts = states.size(); 
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    deg.resize(nkpts);
    for(int i=0;i<states.size();++i)
    {
      const kpdata_t& kp = kpoints[states[i]];
      kpt[i] =  kp.k;
      mk2[i] = -kp.k2;
      deg[i] =  kp.g;
    }
  }
};



//two-d specialization
template<class T>
struct HEGGrid<T,2>
{
  typedef CrystalLattice<T,2> PL_t;
  typedef typename PL_t::SingleParticlePos_t PosType;

  ///number of kpoints of a half sphere excluding gamma
  int NumKptsHalf;
  ///maxmim ksq
  T MaxKsq;
  PL_t& Lattice;
  std::map<int,std::vector<PosType>*> rs;
  std::vector<PosType> kpt;
  std::vector<T>       mk2;
  std::vector<int>     deg;
  std::vector<int> n_within_shell;
  PosType twist;

  typedef kpdata<T,2> kpdata_t;
  typedef std::vector<kpdata_t > kpoints_t;

  kpoints_t* kpoints_grid;

  HEGGrid(PL_t& lat):Lattice(lat),kpoints_grid(0)
  {
    n_within_shell.resize(220);
    //fill this in
    n_within_shell[ 0 ]= 1 ;
    n_within_shell[ 1 ]= 5 ;
    n_within_shell[ 2 ]= 9 ;
    n_within_shell[ 3 ]= 13 ;
    n_within_shell[ 4 ]= 21 ;
    n_within_shell[ 5 ]= 25 ;
    n_within_shell[ 6 ]= 29 ;
    n_within_shell[ 7 ]= 37 ;
    n_within_shell[ 8 ]= 45 ;
    n_within_shell[ 9 ]= 49 ;
    n_within_shell[ 10 ]= 57 ;
    n_within_shell[ 11 ]= 61 ;
    n_within_shell[ 12 ]= 69 ;
    n_within_shell[ 13 ]= 81 ;
    n_within_shell[ 14 ]= 89 ;
    n_within_shell[ 15 ]= 97 ;
    n_within_shell[ 16 ]= 101 ;
    n_within_shell[ 17 ]= 109 ;
    n_within_shell[ 18 ]= 113 ;
    n_within_shell[ 19 ]= 121 ;
    n_within_shell[ 20 ]= 129 ;
    n_within_shell[ 21 ]= 137 ;
    n_within_shell[ 22 ]= 145 ;
    n_within_shell[ 23 ]= 149 ;
    n_within_shell[ 24 ]= 161 ;
    n_within_shell[ 25 ]= 169 ;
    n_within_shell[ 26 ]= 177 ;
    n_within_shell[ 27 ]= 185 ;
    n_within_shell[ 28 ]= 193 ;
    n_within_shell[ 29 ]= 197 ;
    n_within_shell[ 30 ]= 213 ;
    n_within_shell[ 31 ]= 221 ;
    n_within_shell[ 32 ]= 225 ;
    n_within_shell[ 33 ]= 233 ;
    n_within_shell[ 34 ]= 241 ;
    n_within_shell[ 35 ]= 249 ;
    n_within_shell[ 36 ]= 253 ;
    n_within_shell[ 37 ]= 261 ;
    n_within_shell[ 38 ]= 277 ;
    n_within_shell[ 39 ]= 285 ;
    n_within_shell[ 40 ]= 293 ;
    n_within_shell[ 41 ]= 301 ;
    n_within_shell[ 42 ]= 305 ;
    n_within_shell[ 43 ]= 317 ;
    n_within_shell[ 44 ]= 325 ;
    n_within_shell[ 45 ]= 333 ;
    n_within_shell[ 46 ]= 341 ;
    n_within_shell[ 47 ]= 349 ;
    n_within_shell[ 48 ]= 357 ;
    n_within_shell[ 49 ]= 365 ;
    n_within_shell[ 50 ]= 373 ;
    n_within_shell[ 51 ]= 377 ;
    n_within_shell[ 52 ]= 385 ;
    n_within_shell[ 53 ]= 401 ;
    n_within_shell[ 54 ]= 405 ;
    n_within_shell[ 55 ]= 421 ;
    n_within_shell[ 56 ]= 429 ;
    n_within_shell[ 57 ]= 437 ;
    n_within_shell[ 58 ]= 441 ;
    n_within_shell[ 59 ]= 457 ;
    n_within_shell[ 60 ]= 465 ;
    n_within_shell[ 61 ]= 473 ;
    n_within_shell[ 62 ]= 481 ;
    n_within_shell[ 63 ]= 489 ;
    n_within_shell[ 64 ]= 497 ;
    n_within_shell[ 65 ]= 505 ;
    n_within_shell[ 66 ]= 509 ;
    n_within_shell[ 67 ]= 517 ;
    n_within_shell[ 68 ]= 529 ;
    n_within_shell[ 69 ]= 545 ;
    n_within_shell[ 70 ]= 553 ;
    n_within_shell[ 71 ]= 561 ;
    n_within_shell[ 72 ]= 569 ;
    n_within_shell[ 73 ]= 577 ;
    n_within_shell[ 74 ]= 593 ;
    n_within_shell[ 75 ]= 601 ;
    n_within_shell[ 76 ]= 609 ;
    n_within_shell[ 77 ]= 613 ;
    n_within_shell[ 78 ]= 621 ;
    n_within_shell[ 79 ]= 633 ;
    n_within_shell[ 80 ]= 641 ;
    n_within_shell[ 81 ]= 657 ;
    n_within_shell[ 82 ]= 665 ;
    n_within_shell[ 83 ]= 673 ;
    n_within_shell[ 84 ]= 681 ;
    n_within_shell[ 85 ]= 697 ;
    n_within_shell[ 86 ]= 709 ;
    n_within_shell[ 87 ]= 717 ;
    n_within_shell[ 88 ]= 725 ;
    n_within_shell[ 89 ]= 733 ;
    n_within_shell[ 90 ]= 741 ;
    n_within_shell[ 91 ]= 749 ;
    n_within_shell[ 92 ]= 757 ;
    n_within_shell[ 93 ]= 761 ;
    n_within_shell[ 94 ]= 769 ;
    n_within_shell[ 95 ]= 777 ;
    n_within_shell[ 96 ]= 793 ;
    n_within_shell[ 97 ]= 797 ;
    n_within_shell[ 98 ]= 805 ;
    n_within_shell[ 99 ]= 821 ;
    n_within_shell[ 100 ]= 829 ;
    n_within_shell[ 101 ]= 845 ;
    n_within_shell[ 102 ]= 853 ;
    n_within_shell[ 103 ]= 861 ;
    n_within_shell[ 104 ]= 869 ;
    n_within_shell[ 105 ]= 877 ;
    n_within_shell[ 106 ]= 885 ;
    n_within_shell[ 107 ]= 889 ;
    n_within_shell[ 108 ]= 901 ;
    n_within_shell[ 109 ]= 917 ;
    n_within_shell[ 110 ]= 925 ;
    n_within_shell[ 111 ]= 933 ;
    n_within_shell[ 112 ]= 941 ;
    n_within_shell[ 113 ]= 949 ;
    n_within_shell[ 114 ]= 965 ;
    n_within_shell[ 115 ]= 973 ;
    n_within_shell[ 116 ]= 981 ;
    n_within_shell[ 117 ]= 989 ;
    n_within_shell[ 118 ]= 997 ;
    n_within_shell[ 119 ]= 1005 ;
    n_within_shell[ 120 ]= 1009 ;
    n_within_shell[ 121 ]= 1033 ;
    n_within_shell[ 122 ]= 1041 ;
    n_within_shell[ 123 ]= 1049 ;
    n_within_shell[ 124 ]= 1057 ;
    n_within_shell[ 125 ]= 1069 ;
    n_within_shell[ 126 ]= 1085 ;
    n_within_shell[ 127 ]= 1093 ;
    n_within_shell[ 128 ]= 1101 ;
    n_within_shell[ 129 ]= 1109 ;
    n_within_shell[ 130 ]= 1117 ;
    n_within_shell[ 131 ]= 1125 ;
    n_within_shell[ 132 ]= 1129 ;
    n_within_shell[ 133 ]= 1137 ;
    n_within_shell[ 134 ]= 1153 ;
    n_within_shell[ 135 ]= 1161 ;
    n_within_shell[ 136 ]= 1177 ;
    n_within_shell[ 137 ]= 1185 ;
    n_within_shell[ 138 ]= 1201 ;
    n_within_shell[ 139 ]= 1209 ;
    n_within_shell[ 140 ]= 1217 ;
    n_within_shell[ 141 ]= 1225 ;
    n_within_shell[ 142 ]= 1229 ;
    n_within_shell[ 143 ]= 1237 ;
    n_within_shell[ 144 ]= 1245 ;
    n_within_shell[ 145 ]= 1257 ;
    n_within_shell[ 146 ]= 1265 ;
    n_within_shell[ 147 ]= 1273 ;
    n_within_shell[ 148 ]= 1281 ;
    n_within_shell[ 149 ]= 1289 ;
    n_within_shell[ 150 ]= 1305 ;
    n_within_shell[ 151 ]= 1313 ;
    n_within_shell[ 152 ]= 1321 ;
    n_within_shell[ 153 ]= 1329 ;
    n_within_shell[ 154 ]= 1353 ;
    n_within_shell[ 155 ]= 1361 ;
    n_within_shell[ 156 ]= 1369 ;
    n_within_shell[ 157 ]= 1373 ;
    n_within_shell[ 158 ]= 1389 ;
    n_within_shell[ 159 ]= 1405 ;
    n_within_shell[ 160 ]= 1413 ;
    n_within_shell[ 161 ]= 1425 ;
    n_within_shell[ 162 ]= 1433 ;
    n_within_shell[ 163 ]= 1441 ;
    n_within_shell[ 164 ]= 1449 ;
    n_within_shell[ 165 ]= 1457 ;
    n_within_shell[ 166 ]= 1465 ;
    n_within_shell[ 167 ]= 1473 ;
    n_within_shell[ 168 ]= 1481 ;
    n_within_shell[ 169 ]= 1489 ;
    n_within_shell[ 170 ]= 1505 ;
    n_within_shell[ 171 ]= 1513 ;
    n_within_shell[ 172 ]= 1517 ;
    n_within_shell[ 173 ]= 1533 ;
    n_within_shell[ 174 ]= 1541 ;
    n_within_shell[ 175 ]= 1549 ;
    n_within_shell[ 176 ]= 1565 ;
    n_within_shell[ 177 ]= 1581 ;
    n_within_shell[ 178 ]= 1597 ;
    n_within_shell[ 179 ]= 1605 ;
    n_within_shell[ 180 ]= 1609 ;
    n_within_shell[ 181 ]= 1617 ;
    n_within_shell[ 182 ]= 1633 ;
    n_within_shell[ 183 ]= 1641 ;
    n_within_shell[ 184 ]= 1649 ;
    n_within_shell[ 185 ]= 1653 ;
    n_within_shell[ 186 ]= 1669 ;
    n_within_shell[ 187 ]= 1685 ;
    n_within_shell[ 188 ]= 1693 ;
    n_within_shell[ 189 ]= 1701 ;
    n_within_shell[ 190 ]= 1709 ;
    n_within_shell[ 191 ]= 1725 ;
    n_within_shell[ 192 ]= 1733 ;
    n_within_shell[ 193 ]= 1741 ;
    n_within_shell[ 194 ]= 1749 ;
    n_within_shell[ 195 ]= 1757 ;
    n_within_shell[ 196 ]= 1765 ;
    n_within_shell[ 197 ]= 1781 ;
    n_within_shell[ 198 ]= 1789 ;
    n_within_shell[ 199 ]= 1793 ;
    n_within_shell[ 200 ]= 1801 ;
    n_within_shell[ 201 ]= 1813 ;
    n_within_shell[ 202 ]= 1829 ;
    n_within_shell[ 203 ]= 1837 ;
    n_within_shell[ 204 ]= 1853 ;
    n_within_shell[ 205 ]= 1861 ;
    n_within_shell[ 206 ]= 1869 ;
    n_within_shell[ 207 ]= 1877 ;
    n_within_shell[ 208 ]= 1885 ;
    n_within_shell[ 209 ]= 1893 ;
    n_within_shell[ 210 ]= 1901 ;
    n_within_shell[ 211 ]= 1917 ;
    n_within_shell[ 212 ]= 1925 ;
    n_within_shell[ 213 ]= 1933 ;
    n_within_shell[ 214 ]= 1941 ;
    n_within_shell[ 215 ]= 1961 ;
    n_within_shell[ 216 ]= 1969 ;
    n_within_shell[ 217 ]= 1977 ;
    n_within_shell[ 218 ]= 1993 ;
    n_within_shell[ 219 ]= 2001 ;
  }

  ~HEGGrid()
  {
    typename std::map<int,std::vector<PosType>*>::iterator it(rs.begin()),
             it_end(rs.end());
    while(it != it_end)
    {
      delete (*it).second;
      ++it;
    }
  }

  /** return the estimated number of grid in each direction */
  inline int getNC(int nup) const
  {
    return static_cast<int>(std::pow(static_cast<T>(nup),1.0/2))/2+1;
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
    std::vector<int>::const_iterator loc=std::upper_bound(n_within_shell.begin(),n_within_shell.end(),nkpt);
    if(loc<n_within_shell.end())
      return loc-n_within_shell.begin()-1;
    else
      return getNC(nkpt);
  }

  inline int getShellFromStates(int nst)
  {
    for(int i=0; i<n_within_shell.size(); i++)
      if (n_within_shell[i]==nst)
        return i;
    return -1;
  }

  inline T getCellLength(int nptcl, T rs_in) const
  {
    return std::sqrt(M_PI*nptcl)*rs_in;
  }

  void sortGrid(int nc)
  {
    int first_ix2, first_ix3;
    for(int ix1=0; ix1<=nc; ix1++)
    {
      if(ix1 == 0)
        first_ix2=1;
      else
        first_ix2=-nc;
      for(int ix2=first_ix2; ix2<=nc; ix2++)
      {
        int ih=ix1*ix1+ix2*ix2;
        typename std::map<int,std::vector<PosType>*>::iterator it = rs.find(ih);
        if(it == rs.end())
        {
          std::vector<PosType>* ns = new std::vector<PosType>;
          ns->push_back(PosType(ix1,ix2));
          rs[ih] = ns;
        }
        else
        {
          (*it).second->push_back(PosType(ix1,ix2));
        }
      }
    }
  }

  void createGrid(int nc, int nkpts)
  {
    if(rs.empty())
      sortGrid(nc);
    NumKptsHalf=nkpts;
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    int ikpt=0;
    //int checkNum=0;
    //int ke=0;
    MaxKsq=0.0;
    int rsbin=0;
    typename std::map<int, std::vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
    while(ikpt<nkpts && rs_it != rs_end)
    {
      typename std::vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
      T minus_ksq=-Lattice.ksq(*ns_it);
      while(ikpt<nkpts && ns_it!=ns_end)
      {
        kpt[ikpt]=Lattice.k_cart(*ns_it);
        mk2[ikpt]=minus_ksq;
        ++ikpt;
        ++ns_it;
      }
      ++rsbin;
      ++rs_it;
    }
    MaxKsq = Lattice.ksq(*((*rs_it).second->begin()));
    app_log() << "List of kpoints (half-sphere) " << std::endl;
    for(int ik=0; ik<kpt.size(); ik++)
    {
      app_log() << ik << " " << kpt[ik] << " " << mk2[ik] << std::endl;
    }
  }


  void createGrid(const PosType& twistAngle)
  {
    twist = twistAngle;
    //unfold and add gamma
    int nkpts = 2*NumKptsHalf+1;
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    app_log() << "Check this " << NumKptsHalf << " " << nkpts << std::endl;
    abort();
    //add gamma
    int ikpt=0;
    kpt[ikpt]=Lattice.k_cart(twistAngle);
    mk2[ikpt]=-Lattice.ksq(twistAngle);
    ++ikpt;
    typename std::map<int, std::vector<PosType>*>::iterator rs_it(rs.begin()), rs_end(rs.end());
    while(ikpt<nkpts && rs_it != rs_end)
    {
      typename std::vector<PosType>::iterator ns_it((*rs_it).second->begin()), ns_end((*rs_it).second->end());
      while(ikpt<nkpts && ns_it!=ns_end)
      {
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
    app_log() << "List of kpoints with twist = " << twistAngle << std::endl;
    for(int ik=0; ik<kpt.size(); ik++)
    {
      app_log() << ik << " " << kpt[ik] << " " <<-mk2[ik] << std::endl;
    }
  }


  void createGrid(const std::vector<int>& states, T tol=1e-6)
  {
    APP_ABORT("HEGGrid::createGrid(states) has not been implemented");
  }

};
}

#endif

