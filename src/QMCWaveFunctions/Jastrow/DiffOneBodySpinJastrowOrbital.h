//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DIFFERENTIAL_ONEBODYSPINJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_ONEBODYSPINJASTROW_H
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
class DiffOneBodySpinJastrowOrbital: public DiffOrbitalBase
{
  bool Spin;
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///starting index
  int VarOffset;
  ///index of the table
  int myTableIndex;
  ///reference to the ions
  const ParticleSet& CenterRef;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  Matrix<FT*> F;
  ///container for the unique Jastrow functions
  Matrix<int> Fmask;
  std::vector<int> s_offset;
  std::vector<int> t_offset;
  Vector<RealType> dLogPsi;
  std::vector<GradVectorType*> gradLogPsi;
  std::vector<ValueVectorType*> lapLogPsi;

public:

  ///constructor
  DiffOneBodySpinJastrowOrbital(const ParticleSet& centers, ParticleSet& els)
    :Spin(false),CenterRef(centers),NumVars(0),VarOffset(0)
  {
    myTableIndex=els.addTable(CenterRef,DT_SOA_PREFERRED);
    NumPtcls=els.getTotalNum();
    F.resize(CenterRef.groups(), els.groups());
    for(int i=0; i<F.size(); ++i)
      F(i)=0;
    Fmask.resize(CenterRef.groups(), els.groups());
    Fmask=-1;
    s_offset.resize(CenterRef.groups()+1,0);
    t_offset.resize(els.groups()+1,0);
    for(int s=0; s<F.rows(); ++s)
      s_offset[s+1]=centers.last(s);
    for(int t=0; t<F.cols(); ++t)
      t_offset[t+1]=els.last(t);
  }

  ~DiffOneBodySpinJastrowOrbital()
  {
    delete_iter(gradLogPsi.begin(),gradLogPsi.end());
    delete_iter(lapLogPsi.begin(),lapLogPsi.end());
    if(Spin)
    {
      for(int sg=0; sg<F.rows(); ++sg)
        for(int tg=0; tg<F.cols(); ++tg)
          if(F(sg,tg))
          {
            delete F(sg,tg);
          }
          else
          {
            for(int sg=0; sg<F.rows(); ++sg)
              if(F(sg,0))
                delete F(sg,0);
          }
    }
  }

  /** Add a radial functor for a group
   * @param source_type group index of the center species
   * @param afunc radial functor
   */
  void addFunc(int source_g, FT* afunc, int target_g=-1)
  {
    if(target_g<0)
    {
      if(Spin)
      {
        APP_ABORT("Cannot mix spin-depdent with spin-indepdentent Jastrow");
      }
      int pid=source_g*F.cols();
      for(int ig=0; ig<F.cols(); ++ig)
      {
        F(source_g,ig)=afunc;
        Fmask(source_g,ig)=pid;
      }
    }
    else
    {
      Spin=true;
      F(source_g,target_g)=afunc;
      Fmask(source_g,target_g)=source_g*F.cols()+target_g;
    }
  }


  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<F.size(); ++i)
      if(Fmask(i) == i)
        F(i)->resetParameters(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.clear();
    for(int i=0; i<F.size(); ++i)
      if(Fmask(i) == i)
      {
        F(i)->checkOutVariables(active);
        myVars.insertFrom(F(i)->myVars);
      }
    myVars.getIndex(active);
    NumVars=myVars.size();
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
    }
    //int varoffset=myVars.Index[0];
    //OffSet.resize(F.size());
    //for(int i=0; i<F.size(); ++i)
    //{
    //  if(F(i))
    //  {
    //    OffSet[i].first=F(i)->myVars.Index.front()-varoffset;
    //    OffSet[i].second=F(i)->myVars.Index.back()-varoffset+1;
    //  }
    //}
  }

  ///reset the distance table
  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi)
  {
    if (myVars.Index.size()==0)
      return;
    dLogPsi=0.0;
    for (int p=0; p<NumVars; ++p)
      (*gradLogPsi[p])=0.0;
    for (int p=0; p<NumVars; ++p)
      (*lapLogPsi[p])=0.0;
    std::vector<TinyVector<RealType,3> > derivs(NumVars);
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int varoffset=myVars.Index[0];
    for(int ig=0; ig<F.rows(); ++ig)//species
    {
      for(int iat=s_offset[ig]; iat< s_offset[ig+1]; ++iat)//
      {
        int nn=d_table->M[iat];//starting nn for the iat-th source
        for(int jg=0; jg<F.cols(); ++jg)
        {
          FT* func=F(ig,jg);
          if(func && func->myVars.is_optimizable())
          {
            int first=func->myVars.Index.front()-varoffset;
            int last=func->myVars.Index.back()-varoffset+1;
            for(int jat=t_offset[jg]; jat< t_offset[jg+1]; ++jat,++nn)
            {
              std::fill(derivs.begin(),derivs.end(),0.0);
              if (!func->evaluateDerivatives(d_table->r(nn),derivs))
                continue;
              RealType rinv(d_table->rinv(nn));
              PosType dr(d_table->dr(nn));
              for (int p=first, ip=0; p<last; ++p,++ip)
              {
                dLogPsi[p] -= derivs[ip][0];
                RealType dudr(rinv*derivs[ip][1]);
                (*gradLogPsi[p])[jat] -= dudr*dr;
                (*lapLogPsi[p])[jat] -= derivs[ip][2]+2.0*dudr;
              }
            }
          }
          else
          {
            nn+=t_offset[jg+1]-t_offset[jg];
          }
        }//j groups
      }//iat in the ig-th group
    }//ig
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0)
        continue;
      dlogpsi[kk]=dLogPsi[k];
      dhpsioverpsi[kk]=-0.5*Sum(*lapLogPsi[k])-Dot(P.G,*gradLogPsi[k]);
    }
  }

  inline void setVars(const opt_variables_type& vars)
  {
    NumVars=vars.size();
    if(NumVars==0)
      return;
    myVars=vars;
    dLogPsi.resize(NumVars);
    gradLogPsi.resize(NumVars,0);
    lapLogPsi.resize(NumVars,0);
    for (int i=0; i<NumVars; ++i)
    {
      gradLogPsi[i]=new GradVectorType(NumPtcls);
      lapLogPsi[i]=new ValueVectorType(NumPtcls);
    }
  }

  DiffOrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    DiffOneBodySpinJastrowOrbital<FT>* j1copy=new DiffOneBodySpinJastrowOrbital<FT>(CenterRef,tqp);
    if(Spin)
    {
      for(int sg=0; sg<F.rows(); ++sg)
        for(int tg=0; tg<F.cols(); ++tg)
          if(F(sg,tg))
            j1copy->addFunc(sg,new FT(*F(sg,tg)),tg);
    }
    else
    {
      for(int sg=0; sg<F.rows(); ++sg)
        if(F(sg,0))
          j1copy->addFunc(sg,new FT(*F(sg,0)),-1);
    }
    j1copy->Spin=Spin;
    j1copy->setVars(myVars);
    return j1copy;
  }
};
}
#endif

