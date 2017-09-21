//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 */
template<class FT>
class DiffTwoBodyJastrowOrbital: public DiffOrbitalBase
{
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///number of groups, e.g., for the up/down electrons
  int NumGroups;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  std::vector<FT*> F;
  ///offset for the optimizable variables
  std::vector<std::pair<int,int> > OffSet;
  Vector<RealType> dLogPsi;
  std::vector<GradVectorType*> gradLogPsi;
  std::vector<ValueVectorType*> lapLogPsi;
  std::map<std::string,FT*> J2Unique;

public:

  ///constructor
  DiffTwoBodyJastrowOrbital(ParticleSet& p):NumVars(0)
  {
    NumPtcls=p.getTotalNum();
    NumGroups=p.groups();
    F.resize(NumGroups*NumGroups,0);
  }

  ~DiffTwoBodyJastrowOrbital()
  {
    delete_iter(gradLogPsi.begin(),gradLogPsi.end());
    delete_iter(lapLogPsi.begin(),lapLogPsi.end());
  }

  void addFunc(int ia, int ib, FT* j)
  {
    // make all pair terms equal to uu initially
    //   in case some terms are not provided explicitly
    if(ia==ib)
    {
      if(ia==0)//first time, assign everything
      {
        int ij=0;
        for(int ig=0; ig<NumGroups; ++ig)
          for(int jg=0; jg<NumGroups; ++jg, ++ij)
            if(F[ij]==nullptr) F[ij]=j;
      }
      else
        F[ia*NumGroups+ib]=j;
    }
    else
    {
      if(NumPtcls==2)
      {
        // a very special case, 1 up + 1 down
        // uu/dd was prevented by the builder
        for(int ig=0; ig<NumGroups; ++ig)
          for(int jg=0; jg<NumGroups; ++jg)
            F[ig*NumGroups+jg]=j;
      }
      else
      {
        // generic case
        F[ia*NumGroups+ib]=j;
        F[ib*NumGroups+ia]=j;
      }
    }
    std::stringstream aname;
    aname<<ia<<ib;
    J2Unique[aname.str()]=j;
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it++).second->resetParameters(active);
    }
  }

  ///reset the distance table
  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->myVars.getIndex(active);
      myVars.insertFrom((*it).second->myVars);
      ++it;
    }
    myVars.getIndex(active);
    NumVars=myVars.size();

    //myVars.print(std::cout);

    if (NumVars && dLogPsi.size()==0)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars,0);
      lapLogPsi.resize(NumVars,0);
      for (int i=0; i<NumVars; ++i)
      {
        gradLogPsi[i]=new GradVectorType(NumPtcls);
        lapLogPsi[i]=new ValueVectorType(NumPtcls);
      }
      OffSet.resize(F.size());
      int varoffset=myVars.Index[0];
      for (int i=0; i<F.size(); ++i)
      {
        if(F[i] && F[i]->myVars.Index.size())
        {
          OffSet[i].first=F[i]->myVars.Index.front()-varoffset;
          OffSet[i].second=F[i]->myVars.Index.size()+OffSet[i].first;
        }
        else
        {
          OffSet[i].first=OffSet[i].second=-1;
        }
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi)
  {
    if(myVars.size()==0) return;
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(),false);
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0)
        continue;
      if (active.recompute(kk))
        recalculate=true;
      rcsingles[k]=true;
    }
    if (recalculate)
    {
      ///precomputed recalculation switch
      std::vector<bool> RecalcSwitch(F.size(),false);
      for (int i=0; i<F.size(); ++i)
      {
        if(OffSet[i].first<0)
        {
          // nothing to optimize
          RecalcSwitch[i]=false;
        }
        else
        {
          bool recalcFunc(false);
          for (int rcs=OffSet[i].first; rcs<OffSet[i].second; rcs++)
            if (rcsingles[rcs]==true) recalcFunc=true;
          RecalcSwitch[i]=recalcFunc;
        }
      }
      dLogPsi=0.0;
      for (int p=0; p<NumVars; ++p)
        (*gradLogPsi[p])=0.0;
      for (int p=0; p<NumVars; ++p)
        (*lapLogPsi[p])=0.0;
      std::vector<TinyVector<RealType,3> > derivs(NumVars);
      const DistanceTableData* d_table=P.DistTables[0];
      if(d_table->DTType == DT_SOA)
      {
        constexpr RealType cone(1);
        constexpr RealType lapfac(OHMMS_DIM-cone);
        const size_t n=d_table->size(SourceIndex);
        const size_t ng=P.groups();
        for(size_t i=1; i<n; ++i)
        {
          const size_t ig=P.GroupID[i]*ng;
          const RealType* dist=d_table->Distances[i];
          const auto& displ=d_table->Displacements[i];
          for(size_t j=0; j<i; ++j)
          {
            const size_t ptype=ig+P.GroupID[j];
            if (RecalcSwitch[ptype])
            {
              std::fill(derivs.begin(),derivs.end(),0.0);
              if (!F[ptype]->evaluateDerivatives(dist[j],derivs)) continue;
              RealType rinv(cone/dist[j]);
              PosType dr(displ[j]);
              for (int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
              {
                RealType dudr(rinv*derivs[ip][1]);
                RealType lap(derivs[ip][2]+lapfac*dudr);
                //RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
                PosType gr(dudr*dr);
                dLogPsi[p]-=derivs[ip][0];
                (*gradLogPsi[p])[i] += gr;
                (*gradLogPsi[p])[j] -= gr;
                (*lapLogPsi[p])[i] -=lap;
                (*lapLogPsi[p])[j] -=lap;
              }
            }
          }
        }
      }
      else
      {
        for (int i=0; i<d_table->size(SourceIndex); ++i)
        {
          for (int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn)
          {
            int ptype=d_table->PairID[nn];
            if (RecalcSwitch[ptype])
            {
              std::fill(derivs.begin(),derivs.end(),0.0);
              if (!F[ptype]->evaluateDerivatives(d_table->r(nn),derivs)) continue;
              int j = d_table->J[nn];
              RealType rinv(d_table->rinv(nn));
              PosType dr(d_table->dr(nn));
              for (int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
              {
                RealType dudr(rinv*derivs[ip][1]);
                RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
                PosType gr(dudr*dr);
                dLogPsi[p]-=derivs[ip][0];
                (*gradLogPsi[p])[i] += gr;
                (*gradLogPsi[p])[j] -= gr;
                (*lapLogPsi[p])[i] -=lap;
                (*lapLogPsi[p])[j] -=lap;
              }
            }
          }
        }
      }
      for (int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if (kk<0)
          continue;
        if (rcsingles[k])
        {
          dlogpsi[kk]=dLogPsi[k];
          dhpsioverpsi[kk]=-0.5*Sum(*lapLogPsi[k])-Dot(P.G,*gradLogPsi[k]);
        }
        //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
      }
    }
  }

  DiffOrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    DiffTwoBodyJastrowOrbital<FT>* j2copy=new DiffTwoBodyJastrowOrbital<FT>(tqp);
    std::map<const FT*,FT*> fcmap;
    for (int ig=0; ig<NumGroups; ++ig)
      for (int jg=ig; jg<NumGroups; ++jg)
      {
        int ij=ig*NumGroups+jg;
        if (F[ij]==0)
          continue;
        typename std::map<const FT*,FT*>::iterator fit=fcmap.find(F[ij]);
        if (fit == fcmap.end())
        {
          FT* fc=new FT(*F[ij]);
          j2copy->addFunc(ig,jg,fc);
          fcmap[F[ij]]=fc;
        }
      }
    j2copy->myVars.clear();
    j2copy->myVars.insertFrom(myVars);
    j2copy->NumVars=NumVars;
    j2copy->NumPtcls=NumPtcls;
    j2copy->NumGroups=NumGroups;
    j2copy->dLogPsi.resize(NumVars);
    j2copy->gradLogPsi.resize(NumVars,0);
    j2copy->lapLogPsi.resize(NumVars,0);
    for (int i=0; i<NumVars; ++i)
    {
      j2copy->gradLogPsi[i]=new GradVectorType(NumPtcls);
      j2copy->lapLogPsi[i]=new ValueVectorType(NumPtcls);
    }
    j2copy->OffSet=OffSet;
    return j2copy;
  }

};
}
#endif

