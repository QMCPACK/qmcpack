//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef OHMMS_LINKCEDCELL_H
#define OHMMS_LINKCEDCELL_H
///////////////////////////////////////////////////////////////////////////
// Templated class for LCNode
///////////////////////////////////////////////////////////////////////////
/*! \class LCNode
 *  \author Jeongnim Kim
 *  \brief A node class.
 */
template<class PL>
struct LCNode
{

  typedef typename PL::SingleParticlePos_t SingleParticlePos_t;

  int ID; //!< Unique ID for this node
  int LocalNum; //!< Number of particles belonging to this node
  std::vector<int> Next;//!< ID's of adjacent nodes
  std::vector<SingleParticlePos_t> Del;//!< Displacements of adjacent nodes
  std::vector<int> PtclID; //!< Particle ID's for the particles of this node

  LCNode() { }
  ~LCNode() { }

  inline void init(int id, int nctot)
  {
    ID = id;
    Next.reserve(nctot);
    Del.reserve(nctot);
    PtclID.reserve(100);
  }//!< Sets the ID and reserve memory

  void reset(int nat)
  {
    LocalNum = 0;
    if(PtclID.size() < nat)
      PtclID = std::vector<int>(nat);
  }//!< Resets the number of particles of this node

  inline void add(int ijk, SingleParticlePos_t& disp)
  {
    Next.push_back(ijk);
    Del.push_back(disp);
  }//!< Add an adjacent node with a displacement vector

  inline void add(int i)
  {
    if(LocalNum >= PtclID.size())
      PtclID.push_back(i);
    else
      PtclID[LocalNum] = i;
    LocalNum++;
  }//!< Add a particle with the index i
};

/////////////////////////////////////////////////////////////////////////
/*! \class LCNodeSet
 *  \author Jeongnim Kim
 *  \brief A set of LCNode
 */
/////////////////////////////////////////////////////////////////////////
template<class PL>
struct LCNodeSet { };

// forward declaration
template<class T, unsigned D> class CrystalLattice;

/////////////////////////////////////////////////////////////////////////
// specialization with CrystalLattice<T,D>
/////////////////////////////////////////////////////////////////////////
template<class T, unsigned D>
struct LCNodeSet<CrystalLattice<T,D> > { };

/////////////////////////////////////////////////////////////////////////
// specialization with CrystalLattice<T,1>
/////////////////////////////////////////////////////////////////////////
template<class T>
struct LCNodeSet<CrystalLattice<T,1> >
{

  typedef CrystalLattice<T,1> PL_t;
  typedef typename PL_t::SingleParticlePos_t SingleParticlePos_t;
  typedef LCNode<PL_t>                       LCNode_t;

  int nbox, max_nnbox;
  T   dinv;
  std::vector<LCNode_t* > Nodes;                  //!< List of nodes

  void reset(int nat)
  {
    for(int ic=0; ic<Nodes.size(); ic++)
      Nodes[ic]->reset(nat);
  }

  inline int size() const
  {
    return Nodes.size();
  }
  //!< Returns the number of nodes

  inline int size(int ic) const
  {
    return Nodes[ic]->Next.size();
  }
  //!< Returns the number of adjacent nodes of the ic-th node

  inline int key(int ic, int loc) const
  {
    return Nodes[ic]->Next[loc];
  }
  //!< Returns the key for the loc-th adjacent node of the ic-th node

  inline int numptcl(int ic) const
  {
    return Nodes[ic]->LocalNum;
  }
  //!< Returns the number of particles of the ic-th node

  inline int ptclID(int ic, int loc) const
  {
    return Nodes[ic]->PtclID[loc];
  }
  //!< Returns the ID of the loc-th particle of the ic-th node

  inline SingleParticlePos_t shift(int ic, int loc) const
  {
    return Nodes[ic]->Del[loc];
  }//!< Returns the displ vector for the loc-th adjacent node of the ic-th node

  inline int check(int i, const SingleParticlePos_t& r)
  {
    return static_cast<int>(r[0]*dinv);
  }

  void set(const PL_t& lattice)
  {
    if(Nodes.empty())
    {
      nbox = lattice.Grid[0];
      dinv =  static_cast<T>(nbox);
      Nodes.reserve(nbox);
      for(int ic=0; ic<nbox; ic++)
      {
        Nodes.push_back(new LCNode_t);
        Nodes[ic]->init(ic,2);  //assign the cell index, and reserve space
      }
      for(int ibx=0; ibx< nbox; ibx++)
      {
        SingleParticlePos_t delta(0.0);
        for(int ib=-1; ib<=1; ib++)
        {
          if(ib == 0)
            continue;
          int jbx0 = ibx+ib;
          int jbx = lattice.BConds.apply(jbx0, 0, nbox);
          delta[0] = static_cast<T>((jbx0-jbx)/nbox);
          //if(lattice.BoxBConds[0]) {
          //  if(jbx >= nbox) { jbx -= nbox; delta[0] =  1.0e0;}
          //  if(jbx <  0)    { jbx += nbox; delta[0] = -1.0e0;}
          //}
          if(jbx < 0 || jbx >= nbox)
            continue;
          Nodes[ibx]->add(jbx,delta);
        } // 0 - dimension
      } // 0-dimension
      max_nnbox = size(0);
      for(int ic=1; ic<size(); ic++)
        max_nnbox = std::max(size(ic),max_nnbox);
    }
  }

};

/////////////////////////////////////////////////////////////////////////
// specialization with CrystalLattice<T,2>
/////////////////////////////////////////////////////////////////////////
template<class T>
struct LCNodeSet<CrystalLattice<T,2> >
{

  typedef CrystalLattice<T,2> PL_t;
  typedef typename PL_t::SingleParticlePos_t SingleParticlePos_t;
  typedef LCNode<PL_t>                       LCNode_t;

  int nbox[2];
  T   dinv[2];
  int max_nnbox;
  std::vector<LCNode_t* > Nodes;                  //!< List of nodes

  void reset(int nat)
  {
    for(int ic=0; ic<Nodes.size(); ic++)
      Nodes[ic]->reset(nat);
  }

  inline int size() const
  {
    return Nodes.size();
  }

  inline int size(int ic) const
  {
    return Nodes[ic]->Next.size();
  }

  inline int key(int ic, int loc) const
  {
    return Nodes[ic]->Next[loc];
  }

  inline int numptcl(int ic) const
  {
    return Nodes[ic]->LocalNum;
  }

  inline int ptclID(int ic, int loc) const
  {
    return Nodes[ic]->PtclID[loc];
  }

  inline SingleParticlePos_t shift(int ic, int loc) const
  {
    return Nodes[ic]->Del[loc];
  }

  inline int check(int i, const SingleParticlePos_t& r) const
  {
    return static_cast<int>(r[1]*dinv[1]) + static_cast<int>(r[0]*dinv[0])*nbox[1];
  } //!< Returns the cell-id a particle belongs to

  void set(const PL_t& lattice)
  {
    if(Nodes.empty())
    {
      nbox[0] = lattice.Grid[0];
      nbox[1] = lattice.Grid[1];
      int nctot = nbox[0]*nbox[1];
      dinv[0] = static_cast<T>(nbox[0]);
      dinv[1] = static_cast<T>(nbox[1]);
      Nodes.reserve(nctot);
      for(int ic=0; ic<nctot; ic++)
      {
        Nodes.push_back(new LCNode_t);
        Nodes[ic]->init(ic,9);  //assign the cell index, and reserve space
      }
      int iii=0;
      for(int ibx=0; ibx< nbox[0]; ibx++)
      {
        for(int iby=0; iby < nbox[1]; iby++)
        {
          SingleParticlePos_t delta(0.0);
          for(int jb=-1; jb<=1; jb++)
          {
            int jby =iby+jb;
            delta[1] = 0.0e0;
            if(lattice.BoxBConds[1])
            {
              if(jby >= nbox[1])
              {
                jby -= nbox[1];
                delta[1] =  1.0e0;
              }
              if(jby <  0)
              {
                jby += nbox[1];
                delta[1] = -1.0e0;
              }
            }
            if(jby < 0 || jby >= nbox[1])
              continue;
            for(int ib=-1; ib<=1; ib++)
            {
              if(ib == 0 && jb == 0)
                continue;
              int jbx = ibx+ib;
              delta[0] = 0.0e0;
              if(lattice.BoxBConds[0])
              {
                if(jbx >= nbox[0])
                {
                  jbx -= nbox[0];
                  delta[0] =  1.0e0;
                }
                if(jbx <  0)
                {
                  jbx += nbox[0];
                  delta[0] = -1.0e0;
                }
              }
              if(jbx < 0 || jbx >= nbox[0])
                continue;
              Nodes[iii]->add(jby + jbx*nbox[1],delta);
            } // 0 - dimension
          } // 1 - dimension
          //cout << "-----------------------------------------" << std::endl;
          iii++;
        } // 1-dimension
      } // 0-dimension
      max_nnbox = size(0);
      for(int ic=1; ic<size(); ic++)
        max_nnbox = std::max(size(ic),max_nnbox);
    }
  }
};

/////////////////////////////////////////////////////////////////////////
// specialization with CrystalLattice<T,3>
/////////////////////////////////////////////////////////////////////////
template<class T>
struct LCNodeSet<CrystalLattice<T,3> >
{

  typedef CrystalLattice<T,3> PL_t;
  typedef typename PL_t::SingleParticlePos_t SingleParticlePos_t;
  typedef LCNode<PL_t>                       LCNode_t;

  T dinv[3];
  int nbox[3];
  int max_nnbox;

  std::vector<LCNode_t* > Nodes;//!< List of nodes

  void reset(int nat)
  {
    for(int ic=0; ic<Nodes.size(); ic++)
      Nodes[ic]->reset(nat);
  }

  inline int size() const
  {
    return Nodes.size();
  }

  inline int size(int ic) const
  {
    return Nodes[ic]->Next.size();
  }

  inline int key(int ic, int loc) const
  {
    return Nodes[ic]->Next[loc];
  }

  inline int numptcl(int ic) const
  {
    return Nodes[ic]->LocalNum;
  }

  inline int ptclID(int ic, int loc) const
  {
    return Nodes[ic]->PtclID[loc];
  }

  inline SingleParticlePos_t shift(int ic, int loc) const
  {
    return Nodes[ic]->Del[loc];
  }

  inline int check(int i, const SingleParticlePos_t& r) const
  {
    return static_cast<int>(r[2]*dinv[2])
           +nbox[2]*( static_cast<int>(r[1]*dinv[1])
                      +static_cast<int>(r[0]*dinv[0])*nbox[1]);
  } //!< Returns the id of a cell a particle belongs to

  void set(const PL_t& lattice)
  {
    if(Nodes.empty())
    {
      nbox[0] = lattice.Grid[0];
      nbox[1] = lattice.Grid[1];
      nbox[2] = lattice.Grid[2];
      int nctot = nbox[0]*nbox[1]*nbox[2];
      dinv[0] = static_cast<T>(nbox[0]);
      dinv[1] = static_cast<T>(nbox[1]);
      dinv[2] = static_cast<T>(nbox[2]);
      Nodes.reserve(nctot);
      for(int ic=0; ic<nctot; ic++)
      {
        Nodes.push_back(new LCNode_t);
        Nodes[ic]->init(ic,27);  //assign the cell index, and reserve space
      }
      int iii=0;
      for(int ibx=0; ibx< nbox[0]; ibx++)
      {
        for(int iby=0; iby < nbox[1]; iby++)
        {
          for(int ibz=0; ibz < nbox[2]; ibz++)
          {
            SingleParticlePos_t delta(0.0);
            for(int kb=-1; kb<=1; kb++)
            {
              int jbz = ibz+kb;
              delta[2] = 0.0e0;
              if(lattice.BoxBConds[2])
              {
                if(jbz >= nbox[2])
                {
                  jbz -= nbox[2];
                  delta[2] =  1.0e0;
                }
                if(jbz <  0)
                {
                  jbz += nbox[2];
                  delta[2] = -1.0e0;
                }
              }
              if(jbz < 0 || jbz >= nbox[2])
                continue;
              for(int jb=-1; jb<=1; jb++)
              {
                int jby =iby+jb;
                delta[1] = 0.0e0;
                if(lattice.BoxBConds[1])
                {
                  if(jby >= nbox[1])
                  {
                    jby -= nbox[1];
                    delta[1] =  1.0e0;
                  }
                  if(jby <  0)
                  {
                    jby += nbox[1];
                    delta[1] = -1.0e0;
                  }
                }
                if(jby < 0 || jby >= nbox[1])
                  continue;
                for(int ib=-1; ib<=1; ib++)
                {
                  if(ib == 0 && jb == 0 && kb == 0)
                    continue;
                  int jbx = ibx+ib;
                  delta[0] = 0.0e0;
                  if(lattice.BoxBConds[0])
                  {
                    if(jbx >= nbox[0])
                    {
                      jbx -= nbox[0];
                      delta[0] =  1.0e0;
                    }
                    if(jbx <  0)
                    {
                      jbx += nbox[0];
                      delta[0] = -1.0e0;
                    }
                  }
                  if(jbx < 0 || jbx >= nbox[0])
                    continue;
                  Nodes[iii]->add(jbz+nbox[2]*(jby + jbx*nbox[1]),delta);
                } // 0 - dimension
              } // 1 - dimension
            } // 2 -dimension
            //cout << "-----------------------------------------" << std::endl;
            iii++;
          } // 2-dimension
        } // 1-dimension
      } // 0-dimension
//       for(int ic=0; ic<nctot;ic++) {
// 	cout << "Key = " << Nodes[ic]->ID << std::endl;
// 	for(int next=0; next<Nodes[ic]->Next.size(); next++) {
// 	  std::cout << Nodes[ic]->Next[next] <<" "
// 	       << Nodes[ic]->Del[next] <<  std::endl;
// 	}
//       }
      max_nnbox = size(0);
      for(int ic=1; ic<size(); ic++)
        max_nnbox = std::max(size(ic),max_nnbox);
    }
  }

};
#endif // OHMMS_LINKCEDCELL_H

