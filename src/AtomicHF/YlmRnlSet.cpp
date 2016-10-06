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
    
    


#ifndef OHMMS_YLMRNLSET_ONGRID_CPP_H
#define OHMMS_YLMRNLSET_ONGRID_CPP_H

/*! \fn template<class GT>
  bool YlmRnlSet<GT>::add(int n, int l, int m, int s, value_type occ)
  * \param n principal quantum number
  * \param l angular momentum
  * \param m z-component of angular momentum
  * \param s spin
  * \param occ occupation
  * \return true if succeeds
  * \brief Add a new orbital with quantum numbers \f$ (n,l,m,s) \f$ to
  the list of orbitals.  The orbitals are sorted by their restriction
  type, i.e. if the restriction type is \f$ (n,l) \f$ "spin+space" then
  all the orbitals with the same \f$ (n,l) \f$ are grouped together.
  Each orbital is assigned an id, with the possibility of several
  orbitals sharing the same id if they are restricted.
*/

template<class GT>
bool YlmRnlSet<GT>::add(int n, int l, int m, int s, value_type occ)
{
  if(Restriction == "spin+space")
  {
    NLIndex nl(n,l);
    NL_Map_t::iterator it = NL.find(nl);
    ///if the orbital is new add it to the end of the list
    if(it == NL.end())
    {
      ///assign the orbital a new id
      ID.push_back(NumUniqueOrb);
      ///the id counter is set to 1
      IDcount.push_back(1);
      ///add a new element to the map
      NL[nl] = NumUniqueOrb;
      ///increment the number of unique orbitals
      NumUniqueOrb++;
      ///now add the radial grid orbital
      psi.push_back(RadialOrbital_t(m_grid));
      ///add the quantum numbers to the list
      N.push_back(n);
      L.push_back(l);
      M.push_back(m);
      S.push_back(s);
      Occ.push_back(occ);
    }
    else
    {
      ///if an orbital of the same restriction type has already
      ///been added, add the orbital such that all the
      ///orbitals with the same restriction type are grouped
      ///together
      ///increment the id counter
      IDcount[(*it).second]++;
      ///locate the position in the array where the orbital
      ///will be added
      //	IDmap.resize(IDcount.size());
      std::vector<int> IDmap(IDcount.size());
      IDmap[0] = 0;
      int sum = 0;
      for(int i=1; i < IDmap.size(); i++)
      {
        sum += IDcount[i-1];
        IDmap[i] = sum;
      }
      ID.insert(ID.begin()+IDmap[(*it).second],(*it).second);
      psi.insert(psi.begin()+IDmap[(*it).second],RadialOrbital_t(m_grid));
      N.insert(N.begin()+IDmap[(*it).second],n);
      L.insert(L.begin()+IDmap[(*it).second],l);
      M.insert(M.begin()+IDmap[(*it).second],m);
      S.insert(S.begin()+IDmap[(*it).second],s);
      Occ.insert(Occ.begin()+IDmap[(*it).second],occ);
      IDmap.clear();
    }
  }
  if(Restriction == "spin")
  {
    NLMIndex nlm(n,l,m);
    NLM_Map_t::iterator it = NLM.find(nlm);
    if(it == NLM.end())
    {
      ID.push_back(NumUniqueOrb);
      IDcount.push_back(1);
      NLM[nlm] = NumUniqueOrb;
      NumUniqueOrb++;
      psi.push_back(RadialOrbital_t(m_grid));
      N.push_back(n);
      L.push_back(l);
      M.push_back(m);
      S.push_back(s);
      Occ.push_back(occ);
    }
    else
    {
      IDcount[(*it).second]++;
      std::vector<int> IDmap(IDcount.size());
      //	IDmap.resize(IDcount.size());
      IDmap[0] = 0;
      int sum = 0;
      for(int i=1; i < IDmap.size(); i++)
      {
        sum += IDcount[i-1];
        IDmap[i] = sum;
      }
      ID.insert(ID.begin()+IDmap[(*it).second],(*it).second);
      psi.insert(psi.begin()+IDmap[(*it).second],RadialOrbital_t(m_grid));
      N.insert(N.begin()+IDmap[(*it).second],n);
      L.insert(L.begin()+IDmap[(*it).second],l);
      M.insert(M.begin()+IDmap[(*it).second],m);
      S.insert(S.begin()+IDmap[(*it).second],s);
      Occ.insert(Occ.begin()+IDmap[(*it).second],occ);
      IDmap.clear();
    }
  }
  if(Restriction == "none")
  {
    ///add the orbital at the end of the list
    ID.push_back(NumUniqueOrb);
    IDcount.push_back(1);
    NumUniqueOrb++;
    psi.push_back(RadialOrbital_t(m_grid));
    N.push_back(n);
    L.push_back(l);
    M.push_back(m);
    S.push_back(s);
    Occ.push_back(occ);
  }
  return true;
}



template<class GT>
void YlmRnlSet<GT>::applyRestriction(int norb)
{
  static std::vector<value_type> sum;
  if(sum.empty())
  {
    sum.resize(m_grid->size());
    for(int ig=0; ig < m_grid->size(); ig++)
      sum[ig] = 0.0;
  }
  int o_start = 0;
  int o_end = 0;
  int orb = 0;
  while (orb < norb)
  {
    for(int uorb=0; uorb < NumUniqueOrb; uorb++)
    {
      for(int i=0; i < IDcount[uorb]; i++)
      {
        for(int ig=0; ig < m_grid->size(); ig++)
        {
          sum[ig] += psi[orb](ig);
        }
        orb++;
      }
      int o_end = o_start+IDcount[uorb];
      for(int o = o_start; o < o_end; o++)
      {
        for(int ig=0; ig < m_grid->size(); ig++)
        {
          psi[o](ig) = sum[ig]/IDcount[uorb];
        }
      }
      o_start = o_end;
      for(int ig=0; ig < m_grid->size(); ig++)
        sum[ig] = 0.0;
    }
  }
}
#endif
