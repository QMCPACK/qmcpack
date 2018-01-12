//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl, MultiDiracDeterminantBase* up, MultiDiracDeterminantBase* dn):
  C2node_up(nullptr),C2node_dn(nullptr),C(nullptr),
  CSFcoeff(nullptr),DetsPerCSF(nullptr),CSFexpansion(nullptr),
  IsCloned(false),
  RatioTimer("MultiSlaterDeterminantFast::ratio"),
  RatioGradTimer("MultiSlaterDeterminantFast::ratioGrad"),
  RatioAllTimer("MultiSlaterDeterminantFast::ratio(all)"),
  Ratio1Timer("MultiSlaterDeterminantFast::detEval_ratio"),
  Ratio1GradTimer("MultiSlaterDeterminantFast::detEval_ratioGrad"),
  Ratio1AllTimer("MultiSlaterDeterminantFast::detEval_ratio(all)"),
  UpdateTimer("MultiSlaterDeterminantFast::updateBuffer"),
  EvaluateTimer("MultiSlaterDeterminantFast::evaluate"),
  AccRejTimer("MultiSlaterDeterminantFast::Accept_Reject")
{
  registerTimers();
  //Optimizable=true;
  Optimizable=true;
  OrbitalName="MultiSlaterDeterminantFast";
  usingCSF=false;
  NP = targetPtcl.getTotalNum();
  nels_up = targetPtcl.last(0)-targetPtcl.first(0);
  nels_dn = targetPtcl.last(1)-targetPtcl.first(1);
  FirstIndex_up=targetPtcl.first(0);
  FirstIndex_dn=targetPtcl.first(1);
  Dets.resize(2);
  Dets[0]=up;
  Dets[1]=dn;
  myG.resize(NP);
  myL.resize(NP);
  myG_temp.resize(NP);
  myL_temp.resize(NP);

  usingBF=false;
  BFTrans=0;
}

void MultiSlaterDeterminantFast::initialize()
{
  if(C==nullptr)
  {
    C2node_up=new std::vector<size_t>;
    C2node_dn=new std::vector<size_t>;
    C=new std::vector<RealType>;
    CSFcoeff=new std::vector<RealType>;
    DetsPerCSF=new std::vector<size_t>;
    CSFexpansion=new std::vector<RealType>;
    myVars=new opt_variables_type;
    IsCloned=false;
  }
}

OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
{
  MultiDiracDeterminantBase* up_clone = new MultiDiracDeterminantBase(*Dets[0]);
  MultiDiracDeterminantBase* dn_clone = new MultiDiracDeterminantBase(*Dets[1]);
  MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp,up_clone,dn_clone);
  if(usingBF)
  {
    BackflowTransformation *tr = BFTrans->makeClone(tqp);
    clone->setBF(tr);
  }
  clone->resetTargetParticleSet(tqp);

  //Set IsCloned so that only the main object handles the optimizable data
  clone->IsCloned=true;

  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->C=C;
  clone->myVars=myVars;

  clone->Optimizable=Optimizable;
  clone->usingCSF=usingCSF;
  clone->usingBF=usingBF;

  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
  return clone;
}

MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() 
{ 
  if(!IsCloned)
  {
    delete myVars;
    delete CSFexpansion;
    delete DetsPerCSF;
    delete CSFcoeff;
    delete C;
    delete C2node_dn;
    delete C2node_up;
  }
  //clean up determinants too!
}

void MultiSlaterDeterminantFast::resetTargetParticleSet(ParticleSet& P)
{
  if(usingBF)
  {
    BFTrans->resetTargetParticleSet(P);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(BFTrans->QP);
  }
  else
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }
}


void MultiSlaterDeterminantFast::testMSD(ParticleSet& P, int iat)
{
//     APP_ABORT("Testing disabled for safety");
  app_log() <<"Testing MSDFast. \n";
  int n = nels_up+nels_dn;
  ParticleSet::ParticleGradient_t G(n),G0(n);
  ParticleSet::ParticleLaplacian_t L(n),L0(n);
  ValueType log0;
  GradType G1;
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  /*
       app_log() <<"Testing evaluate(P,G,L). \n";
       std::cout << std::endl << std::endl;
       std::cout <<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 << std::endl;

       for(int i=0; i<n; i++) {
         std::cout <<i  <<"\n"
             <<"  x: " <<G[i][0]-G0[i][0] <<"\n"
             <<"  y: " <<G[i][1]-G0[i][1] <<"\n"
             <<"  z: " <<G[i][2]-G0[i][2] <<"\n"
             <<"  d2: " <<L(i)-L0(i) <<"\n"
             << std::endl;
       }
       std::cout << std::endl << std::endl;
       APP_ABORT("end of test 1");
  */
  Walker_t::WFBuffer_t wbuffer;
  wbuffer.clear();
  registerData(P,wbuffer);
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  PosType dr;
  dr[0] = 0.1;
  dr[1]=0.05;
  dr[2] = -0.01;
  PosType newpos(P.makeMove(iat,dr));
  app_log() <<"Testing ratio(P,dG,dL). \n";
  G=0;
  G0=0;
  L=0;
  log0 = ratioGrad(P,iat,G1);
  G0[iat]=G1;
  std::cout <<"Psi: " << log0 << std::endl;
  for(int i=0; i<n; i++)
  {
    std::cout <<i  <<"\n"
        <<"  x: " <<G[i][0]-G0[i][0] <<"  " <<G[i][0]   <<"\n"
        <<"  y: " <<G[i][1]-G0[i][1] <<"  " <<G[i][1] <<"\n"
        <<"  z: " <<G[i][2]-G0[i][2] <<"  " <<G[i][2] <<"\n"
        << std::endl;
  }
  std::cout << std::endl << std::endl;
  APP_ABORT("After MultiSlaterDeterminantFast::testMSD()");
}

/** Compute VGL of this MultiSlaterDeterminantFast
 *
 * THis is introduced to remove redundant code in 
 * - evaluate(P,G,L)
 * - evaluateLog(P,G,L,buf,fillbuffer)
 * Miguel's note: can this change over time??? I don't know yet
 */
OrbitalBase::ValueType MultiSlaterDeterminantFast::evaluate_vgl_impl(ParticleSet& P
    , ParticleSet::ParticleGradient_t& g_tmp, ParticleSet::ParticleLaplacian_t& l_tmp)
{
  const ValueVector_t& detValues_up = Dets[0]->detValues;
  const ValueVector_t& detValues_dn = Dets[1]->detValues;
  const GradMatrix_t& grads_up = Dets[0]->grads;
  const GradMatrix_t& grads_dn = Dets[1]->grads;
  const ValueMatrix_t& lapls_up = Dets[0]->lapls;
  const ValueMatrix_t& lapls_dn = Dets[1]->lapls;
  const size_t N1 = Dets[0]->FirstIndex;
  const size_t N2 = Dets[1]->FirstIndex;
  const size_t NP1 = Dets[0]->NumPtcls;
  const size_t NP2 = Dets[1]->NumPtcls;
  CONSTEXPR ValueType czero(0);
  ValueType psi=czero;
  g_tmp=czero;
  l_tmp=czero;

  const RealType* restrict cptr=C->data();
  const size_t nc=C->size();
  const size_t* restrict upC=C2node_up->data();
  const size_t* restrict dnC=C2node_dn->data();
  for(size_t i=0; i<nc; ++i)
  {
    const RealType c=cptr[i];
    const size_t up=upC[i];
    const size_t down=dnC[i];
    psi += c*detValues_up[up]*detValues_dn[down];
    const ValueType c_up=c*detValues_dn[down];
    const ValueType c_dn=c*detValues_up[up];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      g_tmp[n] += c_up*grads_up(up,k); //myG[n] += c*grads_up(up,k)*detValues_dn[down];
      l_tmp[n] += c_up*lapls_up(up,k); //myL[n] += c*lapls_up(up,k)*detValues_dn[down];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      g_tmp[n] += c_dn*grads_dn(down,k); //myG[n] += c*grads_dn(down,k)*detValues_up[up];
      l_tmp[n] += c_dn*lapls_dn(down,k); //myL[n] += c*lapls_dn(down,k)*detValues_up[up];
    }
  }
  ValueType psiinv = RealType(1)/psi;
  g_tmp *= psiinv;
  l_tmp *= psiinv;

  return psi;
}

OrbitalBase::ValueType MultiSlaterDeterminantFast::evaluate(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
 
  Dets[0]->evaluateForWalkerMove(P);
  Dets[1]->evaluateForWalkerMove(P);

  psiCurrent=evaluate_vgl_impl(P,myG,myL);

  G += myG;
  for(size_t i=0; i<L.size(); i++)
    L[i] += myL[i] - dot(myG[i],myG[i]);

  EvaluateTimer.stop();
  return psiCurrent;
}

OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = evaluate(P,G,L);
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
}

OrbitalBase::ValueType
MultiSlaterDeterminantFast::evalGrad_impl(ParticleSet& P, int iat, bool newpos, GradType& g_at)
{
  const bool upspin=(iat<FirstIndex_dn);
  const int spin0=(upspin)? 0: 1;
  const int spin1=(upspin)? 1: 0;

  if(newpos)
    Dets[spin0]->evaluateDetsAndGradsForPtclMove(P,iat);
  else
    Dets[spin0]->evaluateGrads(P,iat);

  const GradMatrix_t& grads = (newpos)? Dets[spin0]->new_grads:Dets[spin0]->grads;
  const ValueType *restrict detValues0 = (newpos)? Dets[spin0]->new_detValues.data(): Dets[spin0]->detValues.data();
  const ValueType *restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t *restrict det0=(upspin)? C2node_up->data():C2node_dn->data();
  const size_t *restrict det1=(upspin)? C2node_dn->data():C2node_up->data();
  const RealType *restrict cptr=C->data();
  const size_t nc=C->size();
  const size_t noffset=Dets[spin0]->FirstIndex;
  ValueType psi=ValueType(0);
  for(size_t i=0; i<nc; ++i)
  {
    const size_t d0=det0[i];
    //const size_t d1=det1[i];
    //psi +=  cptr[i]*detValues0[d0]        * detValues1[d1];
    //g_at += cptr[i]*grads(d0,iat-noffset) * detValues1[d1];
    const ValueType t=cptr[i]*detValues1[ det1[i] ];
    psi +=  t*detValues0[d0];
    g_at += t*grads(d0,iat-noffset);
  }
  return psi;
}

OrbitalBase::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: evalGrad not implemented. \n");
  }
  CONSTEXPR RealType cone(1);
  GradType grad_iat;
  ValueType psi=evalGrad_impl(P,iat,false,grad_iat);;
  grad_iat*= (cone/psi);
  return grad_iat;
}

OrbitalBase::ValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratioGrad not implemented. \n");
  }
  UpdateMode=ORB_PBYP_PARTIAL;

  CONSTEXPR RealType cone(1);
  GradType dummy;
  ValueType psiNew=evalGrad_impl(P,iat,true,dummy);
  grad_iat+=(cone/psiNew)*dummy;
  curRatio=psiNew/psiCurrent;
  return curRatio;
}

OrbitalBase::ValueType 
MultiSlaterDeterminantFast::ratio_impl(ParticleSet& P, int iat)
{
  const bool upspin=(iat<FirstIndex_dn);
  const int spin0=(upspin)? 0: 1;
  const int spin1=(upspin)? 1: 0;

  Dets[spin0]->evaluateDetsForPtclMove(P,iat);

  const ValueType *restrict detValues0 = Dets[spin0]->new_detValues.data(); //always new
  const ValueType *restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t *restrict det0=(upspin)? C2node_up->data():C2node_dn->data();
  const size_t *restrict det1=(upspin)? C2node_dn->data():C2node_up->data();
  const RealType *restrict cptr=C->data();
  const size_t nc=C->size();

  ValueType psi=0;
  for(size_t i=0; i<nc; ++i)
    psi += cptr[i]*detValues0[ det0[i] ]*detValues1[ det1[i] ];
  return psi;
}

// use ci_node for this routine only
OrbitalBase::ValueType MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio not implemented. \n");
  }
  UpdateMode=ORB_PBYP_RATIO;
  ValueType psiNew=ratio_impl(P,iat);
  curRatio = psiNew/psiCurrent;
  return curRatio;
}

void MultiSlaterDeterminantFast::acceptMove(ParticleSet& P, int iat)
{
// this should depend on the type of update, ratio / ratioGrad
// for now is incorrect fot ratio(P,iat,dG,dL) updates
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: acceptMove not implemented. \n");
  }
// update psiCurrent,myG_temp,myL_temp
  AccRejTimer.start();
  psiCurrent *= curRatio;
  curRatio=1.0;

  Dets[iat>=nels_up]->acceptMove(P,iat);
  //Dets[DetID[iat]]->acceptMove(P,iat);

  AccRejTimer.stop();
}

void MultiSlaterDeterminantFast::restore(int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
  AccRejTimer.start();

  Dets[iat>=nels_up]->restore(iat);
  //Dets[DetID[iat]]->restore(iat);
  curRatio=1.0;
  AccRejTimer.stop();
}

void MultiSlaterDeterminantFast::registerData(ParticleSet& P, WFBufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }

  Dets[0]->registerData(P,buf);
  Dets[1]->registerData(P,buf);
  LogValue = evaluateLog(P,P.G,P.L);

  buf.add(psiCurrent);
}

// FIX FIX FIX
OrbitalBase::RealType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{

  UpdateTimer.start();

  Dets[0]->updateBuffer(P,buf,fromscratch);
  Dets[1]->updateBuffer(P,buf,fromscratch);

  psiCurrent=evaluate_vgl_impl(P,myG,myL);
  
  P.G += myG;
  for(int i=0; i<P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i],myG[i]);

  buf.put(psiCurrent);

  UpdateTimer.stop();
  return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
}

void MultiSlaterDeterminantFast::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: copyFromBuffer not implemented. \n");
  }
  Dets[0]->copyFromBuffer(P,buf);
  Dets[1]->copyFromBuffer(P,buf);

  buf.get(psiCurrent);
}


void MultiSlaterDeterminantFast::checkInVariables(opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
  {
    if(myVars->size())
      active.insertFrom(*myVars);
    else
      Optimizable=false;
  }
}

void MultiSlaterDeterminantFast::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
    myVars->getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
  {
    if(usingCSF)
    {
      RealType *restrict CSFcoeff_p=CSFcoeff->data();
      for(int i=0; i<CSFcoeff->size()-1; i++)
      {
        int loc=myVars->where(i);
        if(loc>=0)
        {
          CSFcoeff_p[i+1]= (*myVars)[i]=active[loc];
        }
      }
      int cnt=0;
      RealType *restrict C_p=C->data();
      const RealType *restrict CSFexpansion_p=CSFexpansion->data();
      for(int i=0; i<DetsPerCSF->size(); i++)
      {
        for(int k=0; k<(*DetsPerCSF)[i]; k++)
        {
          C_p[cnt] = CSFcoeff_p[i]*CSFexpansion_p[cnt];
          cnt++;
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
    else
    {
      RealType *restrict C_p=C->data();
      for(int i=0; i<C->size()-1; i++)
      {
        int loc=myVars->where(i);
        if(loc>=0)
        {
          C_p[i+1]=(*myVars)[i]=active[loc];
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
  }
}
void MultiSlaterDeterminantFast::reportStatus(std::ostream& os)
{
}


void MultiSlaterDeterminantFast::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
  bool recalculate(false);
  for (int k=0; k<myVars->size(); ++k)
  {
    int kk=myVars->where(k);
    if (kk<0)
      continue;
    if (optvars.recompute(kk))
      recalculate=true;
  }
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
  {
    if(usingCSF)
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      const size_t N1 = Dets[0]->FirstIndex;
      const size_t N2 = Dets[1]->FirstIndex;
      const size_t NP1 = Dets[0]->NumPtcls;
      const size_t NP2 = Dets[1]->NumPtcls;
// myG,myL should already be calculated
      const size_t n = P.getTotalNum();

      ValueType psiinv = (RealType)1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }

      const RealType *restrict C_p=C->data();
      for(size_t i=0; i<C->size(); i++)
      {
        size_t upC = (*C2node_up)[i];
        size_t dnC = (*C2node_dn)[i];
        ValueType tmp1 = C_p[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C_p[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(size_t k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(size_t k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(size_t i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
//       for(int i=0; i<C.size(); i++){
      num=CSFcoeff->size()-1;
      int cnt=0;
//        this one is not optable
      cnt+=(*DetsPerCSF)[0];
      int ip(1);
      for(int i=0; i<num; i++,ip++)
      {
        int kk=myVars->where(i);
        if (kk<0)
        {
          cnt+=(*DetsPerCSF)[ip];
          continue;
        }
        ValueType cdet=0.0,q0=0.0,v1=0.0,v2=0.0;
        const RealType *restrict CSFexpansion_p=CSFexpansion->data();
        for(int k=0; k<(*DetsPerCSF)[ip]; k++)
        {
          size_t upC = (*C2node_up)[cnt];
          size_t dnC = (*C2node_dn)[cnt];
          ValueType tmp1=CSFexpansion_p[cnt]*detValues_dn[dnC]*psiinv;
          ValueType tmp2=CSFexpansion_p[cnt]*detValues_up[upC]*psiinv;
          cdet+=CSFexpansion_p[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
          q0 += (tmp1*laplSum_up[upC] + tmp2*laplSum_dn[dnC]);
          for(size_t l=0,j=N1; l<NP1; l++,j++)
            v1 += tmp1*static_cast<ValueType>(dot(P.G[j],grads_up(upC,l))-dot(myG_temp[j],grads_up(upC,l)));
          for(size_t l=0,j=N2; l<NP2; l++,j++)
            v2 += tmp2*static_cast<ValueType>(dot(P.G[j],grads_dn(dnC,l))-dot(myG_temp[j],grads_dn(dnC,l)));
          cnt++;
        }
        convert(cdet,dlogpsi[kk]);
        ValueType dhpsi =  (RealType)-0.5*(q0-cdet*lapl_sum)
                           -cdet*gg-v1-v2;
        //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
        //                         -cdet*lapl_sum)
        //                   -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    }
    else
      //usingCSF
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;
      int n = P.getTotalNum();
      ValueType psiinv = (RealType)1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(size_t k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      const RealType *restrict C_p=C->data();
      for(size_t i=0; i<C->size(); i++)
      {
        size_t upC = (*C2node_up)[i];
        size_t dnC = (*C2node_dn)[i];
        ValueType tmp1 = C_p[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C_p[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(size_t k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(size_t k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(size_t i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
      for(size_t i=1; i<C->size(); i++)
      {
        int kk=myVars->where(i-1);
        if (kk<0)
          continue;
        const size_t upC = (*C2node_up)[i];
        const size_t dnC = (*C2node_dn)[i];
        ValueType cdet=detValues_up[upC]*detValues_dn[dnC]*psiinv;
        ValueType tmp1=detValues_dn[dnC]*psiinv;
        ValueType tmp2=detValues_up[upC]*psiinv;
        convert(cdet,dlogpsi[kk]);
        ValueType v1=0.0,v2=0.0;
        for(size_t k=0,j=N1; k<NP1; k++,j++)
          v1 += (dot(P.G[j],grads_up(upC,k))-dot(myG_temp[j],grads_up(upC,k)) );
        for(size_t k=0,j=N2; k<NP2; k++,j++)
          v2 += (dot(P.G[j],grads_dn(dnC,k))-dot(myG_temp[j],grads_dn(dnC,k)));
        ValueType dhpsi =  (RealType)-0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
                                 -cdet*lapl_sum)
                           -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    } // usingCSF
  }
}

void MultiSlaterDeterminantFast::registerTimers()
{
  RatioTimer.reset();
  RatioGradTimer.reset();
  RatioAllTimer.reset();
  Ratio1Timer.reset();
  Ratio1GradTimer.reset();
  Ratio1AllTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&RatioGradTimer);
  TimerManager.addTimer (&RatioAllTimer);
  TimerManager.addTimer (&Ratio1Timer);
  TimerManager.addTimer (&Ratio1GradTimer);
  TimerManager.addTimer (&Ratio1AllTimer);
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&EvaluateTimer);
  TimerManager.addTimer (&AccRejTimer);
}


}
