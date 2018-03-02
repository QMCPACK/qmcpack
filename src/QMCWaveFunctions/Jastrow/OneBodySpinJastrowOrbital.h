//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GENERIC_ONEBODYJASTROWSPIN_H
#define QMCPLUSPLUS_GENERIC_ONEBODYJASTROWSPIN_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
//#include "QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/** generic implementation of spin-dependent one-body Jastrow function.
 *
 * Template parameter FT denotes the scaling function, e.g., BsplineFunctor.
 * Based on the OneBodyJastrowOrbital for the grouped particle sets.
 */
template<class FT>
class OneBodySpinJastrowOrbital: public OrbitalBase
{
  bool Spin;
  int myTableIndex;
  const ParticleSet& CenterRef;
  RealType curVal;
  RealType curLap;
  PosType curGrad;
  ParticleAttrib<RealType> U,d2U;
  ParticleAttrib<PosType> dU;
  RealType *FirstAddressOfdU;
  RealType *LastAddressOfdU;
  /** scaling functors
   *
   * F(i,j) where i= source species index and j=target group index
   */
  Matrix<FT*> F;
  /** mask for cloning and optimization
   *
   * Fmask(i,j) is set to -1 and updated by addFunc.
   */
  Matrix<int> Fmask;
  /** index offset for the source particles
   *
   * The loop over the particles is replaced by a double loop as
   * for(int sg=0;sg<F.rows();++sg)
   *   for(int s=s_offset[sg]; s<s_offset[sg+1];++s)
   */
  std::vector<int> s_offset;
  /** index offset for the target particles
   *
   * The loop over the particles is replaced by a double loop as
   * for(int tg=0;tg<F.cols();++tg)
   *   for(int t=t_offset[tg]; t<t_offset[tg+1];++t)
   */
  std::vector<int> t_offset;
public:
  typedef FT FuncType;

  ///constructor
  OneBodySpinJastrowOrbital(const ParticleSet& centers, ParticleSet& els)
    : Spin(false), CenterRef(centers), FirstAddressOfdU(0), LastAddressOfdU(0)
  {
    OrbitalName = "OneBodySpinJastrow";
    U.resize(els.getTotalNum());
    myTableIndex=els.addTable(CenterRef,DT_AOS);
    //allocate vector of proper size  and set them to 0
    F.resize(CenterRef.groups(), els.groups());
    for(int i=0; i<F.size(); ++i)
      F(i)=0;
    //initialize mask to handle cloning and variable updates
    Fmask.resize(CenterRef.groups(), els.groups());
    Fmask=-1;
    s_offset.resize(CenterRef.groups()+1,0);
    t_offset.resize(els.groups()+1,0);
    for(int s=0; s<F.rows(); ++s)
      s_offset[s+1]=centers.last(s);
    for(int t=0; t<F.cols(); ++t)
      t_offset[t+1]=els.last(t);
  }

  ~OneBodySpinJastrowOrbital()
  {
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

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
    if (dPsi)
      dPsi->resetTargetParticleSet(P);
  }

  /** add a functor
   * @param source_g index of the source species
   * @param afunc a scaling functor
   * @param target_g index of the target species
   *
   * When target_g is negative,  F(source_g,*)=afunc.
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
      app_log() << " Adding functor of type "  << source_g << " for all the target. " << std::endl;
    }
    else
    {
      Spin=true;
      F(source_g,target_g)=afunc;
      Fmask(source_g,target_g)=source_g*F.cols()+target_g;
      app_log() << " Adding functor of type "  << source_g << " for the target type " << target_g << std::endl;
    }
  }

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    for(int i=0; i<F.size(); ++i)
      if(Fmask(i)==i)
      {
        F(i)->checkInVariables(active);
        F(i)->checkInVariables(myVars);
      }
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    Optimizable=myVars.is_optimizable();
    for(int i=0; i<F.size(); ++i)
      if(Fmask(i) == i)
        F(i)->checkOutVariables(active);
    if (dPsi)
      dPsi->checkOutVariables(active);
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    if (!Optimizable)
      return;
    for(int i=0; i<F.size(); ++i)
      if(Fmask(i) == i)
        F(i)->resetParameters(active);
    //for (int i=0; i<myVars.size(); ++i)
    //{
    //  int ii=myVars.Index[i];
    //  if (ii>=0) myVars[i]= active[ii];
    //}
    if (dPsi)
      dPsi->resetParameters(active);
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
    for(int i=0; i<F.size(); ++i)
    {
      if(Fmask(i) ==i)
      {
        os << "  One-body Jastrow for the group pair " << i/F.cols() << "-" << i%F.cols() << std::endl;
        F(i)->myVars.print(os);
      }
    }
  }

  RealType evaluateLog(ParticleSet& P
                       , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    LogValue=0.0;
    U=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    RealType dudr, d2udr2;
    for(int sg=0; sg<F.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=d_table->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<F.cols(); ++tg)
        {
          FT* func=F(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2);
              LogValue -= uij;
              U[jat] += uij;
              dudr *= d_table->rinv(nn);
              G[jat] -= dudr*d_table->dr(nn);
              L[jat] -= d2udr2+2.0*dudr;
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P
                     , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
   * @param P active particle set
   * @param iat particle that has been moved.
   */
  inline ValueType ratio(ParticleSet& P, int iat)
  {
    curVal=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int tg=P.GroupID[iat];
    for(int sg=0; sg<F.rows(); ++sg)
    {
      FT* func=F(sg,tg);
      if(!func)
        continue;
      for(int s=s_offset[sg]; s<s_offset[sg+1]; ++s)
        curVal += func->evaluate(d_table->Temp[s].r1);
    }
    return std::exp(U[iat]-curVal);
  }

  inline void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    std::vector<RealType> myr(ratios.size(),U[VP.activePtcl]);
    const DistanceTableData* d_table=VP.DistTables[myTableIndex];
    int tg=VP.GroupID[VP.activePtcl];
    for(int sg=0; sg<F.rows(); ++sg)
    {
      FT* func=F(sg,tg);
      if(func)
      {
        for(int s=s_offset[sg]; s<s_offset[sg+1]; ++s)
          for (int nn=d_table->M[s],j=0; nn<d_table->M[s+1]; ++nn,++j)
            myr[j]-=func->evaluate(d_table->r(nn));
      }
    }
    for(int k=0; k<ratios.size(); ++k) ratios[k]=std::exp(myr[k]);
    //RealType x=U[VP.activePtcl];
    //for(int k=0; k<ratios.size(); ++k)
    //  ratios[k]=std::exp(x-myr[k]);
  }

  /** evaluate the ratio
   */
  inline void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
  {
    std::fill(ratios.begin(),ratios.end(),0.0);
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    for(int sg=0; sg<F.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=d_table->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<F.cols(); ++tg)
        {
          FT* func=F(sg,tg);
          if(func)
	  {
	    RealType up=func->evaluate(d_table->Temp[iat].r1);
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
              ratios[jat]+=func->evaluate(d_table->r(nn))-up;
	  }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
	  }
        }
      }  
    
    
    for(int i=0; i<ratios.size(); ++i)
      ratios[i] = std::exp(ratios[i]);
  }  
  

  inline GradType evalGrad(ParticleSet& P, int iat)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int n=d_table->size(VisitorIndex);
    int tg=P.GroupID[iat];
    curGrad = 0.0;
    RealType ur,dudr, d2udr2;
    int nn=iat;
    for(int sg=0; sg<F.rows(); ++sg)
    {
      FT* func=F(sg,tg);
      if(func)
        for(int s=s_offset[sg]; s< s_offset[sg+1]; ++s,nn+=n)
        {
          ur=func->evaluate(d_table->r(nn),dudr,d2udr2);
          dudr *= d_table->rinv(nn);
          curGrad -= dudr*d_table->dr(nn);
        }
      else
        nn += n*(s_offset[sg+1]-s_offset[sg]);
    }
    return curGrad;
  }

  inline GradType evalGradSource(ParticleSet& P, ParticleSet& source, int isrc)
  {
    APP_ABORT("NOT DONE YET");
    return GradType();
  }


  inline GradType
  evalGradSource(ParticleSet& P, ParticleSet& source, int isrc
                 , TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad
                 , TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    APP_ABORT("NOT DONE YET");
    return GradType();
  }

  inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int tg=P.GroupID[iat];//pick the target group
    curVal=0.0;
    curGrad = 0.0;
    RealType dudr, d2udr2;
    for(int sg=0; sg<F.rows(); ++sg)
    {
      FT* func=F(sg,tg);
      if(func)
        for(int s=s_offset[sg]; s< s_offset[sg+1]; ++s)
        {
          curVal += func->evaluate(d_table->Temp[s].r1,dudr,d2udr2);
          dudr *= d_table->Temp[s].rinv1;
          curGrad -= dudr*d_table->Temp[s].dr1;
        }
    }
    grad_iat += curGrad;
    return std::exp(U[iat]-curVal);
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
    U[iat] = curVal;
    dU[iat]=curGrad;
    d2U[iat]=curLap;
  }

  void evaluateLogAndStore(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& dG,
                           ParticleSet::ParticleLaplacian_t& dL)
  {
    LogValue = 0.0;
    U=0.0;
    dU=0.0;
    d2U=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    RealType uij, dudr, d2udr2;
    for(int ig=0; ig<F.rows(); ++ig)
    {
      for(int iat=s_offset[ig]; iat<s_offset[ig+1]; ++iat)
      {
        int nn=d_table->M[iat];
        for(int jg=0; jg<F.cols(); ++jg)
        {
          FT* func=F(ig,jg);
          if(func)
            for(int jat=t_offset[jg]; jat<t_offset[jg+1]; ++jat,++nn)
            {
              uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
              LogValue-=uij;
              U[jat]+=uij;
              dudr *= d_table->rinv(nn);
              dU[jat] -= dudr*d_table->dr(nn);
              d2U[jat] -= d2udr2+2.0*dudr;
              dG[jat] -= dudr*d_table->dr(nn);
              dL[jat] -= d2udr2+2.0*dudr;
            }
          else
            nn += t_offset[jg+1]-t_offset[jg];
        }
      }
    }
  }

  /** equivalent to evalaute with additional data management */
  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    d2U.resize(d_table->size(VisitorIndex));
    dU.resize(d_table->size(VisitorIndex));
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    evaluateLogAndStore(P,P.G,P.L);
    //add U, d2U and dU. Keep the order!!!
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::registerData ",buf.current());
    buf.add(U.begin(), U.end());
    buf.add(d2U.begin(), d2U.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::registerData ",buf.current());
  }

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    evaluateLogAndStore(P,P.G,P.L);
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::updateBuffer ",buf.current());
    buf.put(U.first_address(), U.last_address());
    buf.put(d2U.first_address(), d2U.last_address());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::updateBuffer ",buf.current());
    return LogValue;
  }

  /** copy the current data from a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   *
   *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::copyFromBuffer ",buf.current());
    buf.get(U.first_address(), U.last_address());
    buf.get(d2U.first_address(), d2U.last_address());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodySpinJastrow::copyFromBuffer ",buf.current());
  }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    OneBodySpinJastrowOrbital<FT>* j1copy=new OneBodySpinJastrowOrbital<FT>(CenterRef,tqp);
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
    if(dPsi)
      j1copy->dPsi =  dPsi->makeClone(tqp);
    return j1copy;
  }

  void copyFrom(const OrbitalBase& old)
  {
    //nothing to do
  }
};

}
#endif

