//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h>
#include <LongRange/StructFact.h>
#include <config/stdlib/math.h>

namespace qmcplusplus
{

LRTwoBodyJastrow::LRTwoBodyJastrow(ParticleSet& p ):
  NumPtcls(0), NumSpecies(0), skRef(0)
{
  Optimizable=false;
  NumSpecies=p.groups();
  skRef=p.SK;
  if(skRef)
  {
    p.SK->DoUpdate=true;
    CellVolume = p.Lattice.Volume;
    Rs =std::pow(3.0/4.0/M_PI*p.Lattice.Volume/static_cast<RealType>(p.getTotalNum()),1.0/3.0);
    NormConstant=4.0*M_PI*Rs*NumPtcls*(NumPtcls-1)*0.5;
    NumPtcls=p.getTotalNum();
    NumKpts=skRef->KLists.numk;
    MaxKshell=skRef->KLists.kshell.size();
    resize();
  }
}


void LRTwoBodyJastrow::resize()
{
  //set the  maximum k-shell
#if defined(USE_REAL_STRUCT_FACTOR)
  Rhok_r.resize(NumKpts);
  Rhok_i.resize(NumKpts);
  rokbyF_r.resize(NumPtcls,NumKpts);
  rokbyF_i.resize(NumPtcls,NumKpts);
#else
  Rhok.resize(NumKpts);
  rokbyF.resize(NumPtcls,NumKpts);
#endif
  U.resize(NumPtcls);
  dU.resize(NumPtcls);
  d2U.resize(NumPtcls);
  FirstAddressOfdU=&(dU[0][0]);
  LastAddressOfdU = FirstAddressOfdU+NumPtcls*DIM;
  offU.resize(NumPtcls);
  offdU.resize(NumPtcls);
  offd2U.resize(NumPtcls);
}


void LRTwoBodyJastrow::checkInVariables(opt_variables_type& active)
{
}

void LRTwoBodyJastrow::checkOutVariables(const opt_variables_type& active)
{
}

void LRTwoBodyJastrow::resetParameters(const opt_variables_type& active)
{
}

void LRTwoBodyJastrow::reportStatus(std::ostream& os)
{
}

void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P)
{
  skRef=P.SK;
}

LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::evaluateLog(ParticleSet& P,
                              ParticleSet::ParticleGradient_t& G,
                              ParticleSet::ParticleLaplacian_t& L)
{
  RealType sum(0.0);
#if defined(USE_REAL_STRUCT_FACTOR)
  std::copy(P.SK->rhok_r[0],P.SK->rhok_r[0]+MaxK,Rhok_r.data());
  std::copy(P.SK->rhok_i[0],P.SK->rhok_i[0]+MaxK,Rhok_i.data());
  
  for(int spec1=1; spec1<NumSpecies; spec1++)
  {
    accumulate_elements(P.SK->rhok_r[spec1],P.SK->rhok_r[spec1]+MaxK,Rhok_r.data());
    accumulate_elements(P.SK->rhok_i[spec1],P.SK->rhok_i[spec1]+MaxK,Rhok_i.data());
  }
#else
  //memcopy if necessary but this is not so critcal
  std::copy(P.SK->rhok[0],P.SK->rhok[0]+MaxK,Rhok.data());
  for(int spec1=1; spec1<NumSpecies; spec1++)
    accumulate_elements(P.SK->rhok[spec1],P.SK->rhok[spec1]+MaxK,Rhok.data());
#endif
  const KContainer::VContainer_t& Kcart(P.SK->KLists.kpts_cart);
  for(int iat=0; iat<NumPtcls; iat++)
  {
    RealType res(0.0),l(0.0);
    PosType g;
    #if defined(USE_REAL_STRUCT_FACTOR)
    const RealType* restrict eikr_r_ptr(P.SK->eikr_r[iat]);
    const RealType* restrict eikr_i_ptr(P.SK->eikr_i[iat]);
    const RealType* restrict rhok_r_ptr(Rhok_r.data());
    const RealType* restrict rhok_i_ptr(Rhok_i.data());
    #else
    const ComplexType* restrict eikr_ptr(P.SK->eikr[iat]);
    const ComplexType* restrict rhok_ptr(Rhok.data());
    #endif
    int ki=0;
    for(int ks=0; ks<MaxKshell; ks++)
    {
      RealType res_k(0.0),l_k(0.0);
      PosType g_k;
      for(; ki<Kshell[ks+1]; ki++)
      {
        #if defined(USE_REAL_STRUCT_FACTOR)
        RealType rr=((*eikr_r_ptr)*(*rhok_r_ptr)+(*eikr_i_ptr)*(*rhok_i_ptr));
        RealType ii=((*eikr_r_ptr)*(*rhok_i_ptr)-(*eikr_i_ptr)*(*rhok_r_ptr));
        eikr_r_ptr++;
        eikr_i_ptr++;
        rhok_r_ptr++;
        rhok_i_ptr++;
        #else
        RealType rr=((*eikr_ptr).real()*(*rhok_ptr).real()+(*eikr_ptr).imag()*(*rhok_ptr).imag());
        RealType ii=((*eikr_ptr).real()*(*rhok_ptr).imag()-(*eikr_ptr).imag()*(*rhok_ptr).real());
        eikr_ptr++;
        rhok_ptr++;
        #endif
        res_k +=  rr;
        l_k += (1.0-rr);
        g_k += ii*Kcart[ki];
      }
      res += Fk_symm[ks]*res_k;
      l += FkbyKK[ks]*l_k;
      g += Fk_symm[ks]*g_k;
    }
    sum+=(U[iat]=res);
    G[iat]+=(dU[iat]=g);
    L[iat]+=(d2U[iat]=l);
  }
  return 0.5*sum;
}


/* evaluate the ratio with P.R[iat]
 *
 */
LRTwoBodyJastrow::ValueType
LRTwoBodyJastrow::ratio(ParticleSet& P, int iat)
{
  //restore, if called should do nothing
  NeedToRestore=false;
#if defined(USE_REAL_STRUCT_FACTOR)
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  
  const RealType* restrict eikr_r_ptr(P.SK->eikr_r[iat]);
  const RealType* restrict eikr_i_ptr(P.SK->eikr_i[iat]);
  const RealType* restrict rhok_r_ptr(Rhok_r.data());
  const RealType* restrict rhok_i_ptr(Rhok_i.data());

  const PosType &pos(P.activePos);
  curVal=0.0;
  int ki=0;
  for(int ks=0; ks<MaxKshell; ks++)
  {
    RealType dd=0.0,s,c;
    for(; ki<Kshell[ks+1]; ki++)
    {
      
      sincos(dot(kpts[ki],pos),&s,&c);
      dd += c*(c+(*rhok_r_ptr)-(*eikr_r_ptr))
            + s*(s+(*rhok_i_ptr)-(*eikr_i_ptr));

      eikr_r_ptr++;
      eikr_i_ptr++;
      rhok_r_ptr++;
      rhok_i_ptr++;
    }
    curVal += Fk_symm[ks]*dd;
  }
#else
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  const ComplexType* restrict eikr_ptr(P.SK->eikr[iat]);
  const ComplexType* restrict rhok_ptr(Rhok.data());
  const PosType &pos(P.activePos);
  curVal=0.0;
  int ki=0;
  for(int ks=0; ks<MaxKshell; ks++)
  {
    RealType dd=0.0,s,c;
    for(; ki<Kshell[ks+1]; ki++,eikr_ptr++,rhok_ptr++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      dd += c*(c+(*rhok_ptr).real()-(*eikr_ptr).real())
            + s*(s+(*rhok_ptr).imag()-(*eikr_ptr).imag());
    }
    curVal += Fk_symm[ks]*dd;
  }
#endif
  return std::exp(curVal-U[iat]);
}


LRTwoBodyJastrow::ValueType
LRTwoBodyJastrow::ratioGrad(ParticleSet& P, int iat, GradType & g)
{
  NeedToRestore=true;
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  {
    const PosType &pos(P.activePos);
    RealType c,s;
  #if defined(USE_REAL_STRUCT_FACTOR)
    RealType* restrict eikr1_r(eikr_new_r.data());
    RealType* restrict eikr1_i(eikr_new_i.data());
    RealType* restrict deikr_r(delta_eikr_r.data());
    RealType* restrict deikr_i(delta_eikr_i.data());
    const RealType* restrict eikr0_r(eikr_r[iat]);
    const RealType* restrict eikr0_i(eikr_i[iat]);
    for(int ki=0; ki<MaxK; ki++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      (*eikr1_r)=c;
      (*eikr1_i)=s;
      (*deikr_r++)=(*eikr1_r++)-(*eikr0_r++);
      (*deikr_i++)=(*eikr1_i++)-(*eikr0_i++);
    }
  #else  
    ComplexType* restrict eikr1(eikr_new.data());
    ComplexType* restrict deikr(delta_eikr.data());
    const ComplexType* restrict eikr0(eikr[iat]);
    for(int ki=0; ki<MaxK; ki++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      (*eikr1)=ComplexType(c,s);
      (*deikr++)=(*eikr1++)-(*eikr0++);
    }
  #endif
  }
  curVal=0.0;
  curGrad=0.0;
#if defined(USE_REAL_STRUCT_FACTOR)
  //new Rhok: restored by rejectMove
  Rhok_r += delta_eikr_r;
  Rhok_i += delta_eikr_i;
  const RealType* restrict rhok_ptr_r(Rhok_r.data());
  const RealType* restrict rhok_ptr_i(Rhok_i.data());
  const RealType* restrict eikr1_r(eikr_new_r.data());
  const RealType* restrict eikr1_i(eikr_new_i.data());
#else
  //new Rhok: restored by rejectMove
  Rhok += delta_eikr;
  const ComplexType* restrict rhok_ptr(Rhok.data());
  const ComplexType* restrict eikr1(eikr_new.data());
#endif
  for(int ks=0,ki=0; ks<MaxKshell; ks++)
  {
    RealType v(0.0);
    PosType g;
    for(; ki<Kshell[ks+1]; ki++)
    {
    #if defined(USE_REAL_STRUCT_FACTOR)
      RealType rr=((*eikr1_r)*(*rhok_ptr_r)+(*eikr1_i)*(*rhok_ptr_i));
      RealType ii=((*eikr1_r)*(*rhok_ptr_i)-(*eikr1_i)*(*rhok_ptr_r));
      eikr1_r++;
      eikr1_i++;
      rhok_ptr_r++;
      rhok_ptr_i++;
    #else
      RealType rr=((*eikr1).real()*(*rhok_ptr).real()+(*eikr1).imag()*(*rhok_ptr).imag());
      RealType ii=((*eikr1).real()*(*rhok_ptr).imag()-(*eikr1).imag()*(*rhok_ptr).real());
      eikr1++;
      rhok_ptr++;
    #endif
      v +=  rr;
      g += ii*kpts[ki];
    }
    curVal += Fk_symm[ks]*v;
    curGrad += Fk_symm[ks]*g;
  }
  
  g+=curGrad;
  return std::exp(curVal-U[iat]);

}

LRTwoBodyJastrow::GradType
LRTwoBodyJastrow::evalGrad(ParticleSet& P, int iat)
{
  NeedToRestore=true;
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  {
    const PosType &pos(P.R[iat]);
    RealType c,s;
  #if defined(USE_REAL_STRUCT_FACTOR)
    RealType* restrict eikr1_r(eikr_new_r.data());
    RealType* restrict eikr1_i(eikr_new_i.data());
    RealType* restrict deikr_r(delta_eikr_r.data());
    RealType* restrict deikr_i(delta_eikr_i.data());
    const RealType* restrict eikr0_r(eikr_r[iat]);
    const RealType* restrict eikr0_i(eikr_i[iat]);
    for(int ki=0; ki<MaxK; ki++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      (*eikr1_r)=c;
      (*eikr1_i)=s;
      (*deikr_r++)=(*eikr1_r++)-(*eikr0_r++);
      (*deikr_i++)=(*eikr1_i++)-(*eikr0_i++);
    }
  #else  
    ComplexType* restrict eikr1(eikr_new.data());
    ComplexType* restrict deikr(delta_eikr.data());
    const ComplexType* restrict eikr0(eikr[iat]);
    for(int ki=0; ki<MaxK; ki++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      (*eikr1)=ComplexType(c,s);
      (*deikr++)=(*eikr1++)-(*eikr0++);
    }
  #endif
  }
  curGrad=0.0;
#if defined(USE_REAL_STRUCT_FACTOR)
  //new Rhok: restored by rejectMove
  Rhok_r += delta_eikr_r;
  Rhok_i += delta_eikr_i;
  const RealType* restrict rhok_ptr_r(Rhok_r.data());
  const RealType* restrict rhok_ptr_i(Rhok_i.data());
  const RealType* restrict eikr1_r(eikr_new_r.data());
  const RealType* restrict eikr1_i(eikr_new_i.data());
#else
  //new Rhok: restored by rejectMove
  Rhok += delta_eikr;
  const ComplexType* restrict rhok_ptr(Rhok.data());
  const ComplexType* restrict eikr1(eikr_new.data());
#endif
  for(int ks=0,ki=0; ks<MaxKshell; ks++)
  {
    PosType g;
    for(; ki<Kshell[ks+1]; ki++)
    {
    #if defined(USE_REAL_STRUCT_FACTOR)
      RealType rr=((*eikr1_r)*(*rhok_ptr_r)+(*eikr1_i)*(*rhok_ptr_i));
      RealType ii=((*eikr1_r)*(*rhok_ptr_i)-(*eikr1_i)*(*rhok_ptr_r));
      eikr1_r++;
      eikr1_i++;
      rhok_ptr_r++;
      rhok_ptr_i++;
    #else
      RealType rr=((*eikr1).real()*(*rhok_ptr).real()+(*eikr1).imag()*(*rhok_ptr).imag());
      RealType ii=((*eikr1).real()*(*rhok_ptr).imag()-(*eikr1).imag()*(*rhok_ptr).real());
      eikr1++;
      rhok_ptr++;
    #endif
      g += ii*kpts[ki];
    }
    curGrad += Fk_symm[ks]*g;
  }
  
  return curGrad;

}

void LRTwoBodyJastrow::restore(int iat)
{
  //substract the addition in logRatio
  if (NeedToRestore)
  {
  #if defined(USE_REAL_STRUCT_FACTOR)
    Rhok_r -= delta_eikr_r;
    Rhok_i -= delta_eikr_i;
  #else
    Rhok -= delta_eikr;
  #endif
  }
}

void LRTwoBodyJastrow::acceptMove(ParticleSet& P, int iat)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  std::copy(eikr_new_r.data(),eikr_new_r.data()+MaxK,eikr_r[iat]);
  std::copy(eikr_new_i.data(),eikr_new_i.data()+MaxK,eikr_i[iat]);
#else
  std::copy(eikr_new.data(),eikr_new.data()+MaxK,eikr[iat]);
#endif
  U += offU;
  dU += offdU;
  d2U += offd2U;
}

void
LRTwoBodyJastrow::registerData(ParticleSet& P, WFBufferType& buf)
{
  LogValue=evaluateLog(P,P.G,P.L);
#if defined(USE_REAL_STRUCT_FACTOR)
  eikr_r.resize(NumPtcls,MaxK);
  eikr_i.resize(NumPtcls,MaxK);
  eikr_new_r.resize(MaxK);
  eikr_new_i.resize(MaxK);
  delta_eikr_r.resize(MaxK);
  delta_eikr_i.resize(MaxK);

  for(int iat=0; iat<NumPtcls; iat++)
  {
    std::copy(P.SK->eikr_r[iat],P.SK->eikr_r[iat]+MaxK,eikr_r[iat]);
    std::copy(P.SK->eikr_i[iat],P.SK->eikr_i[iat]+MaxK,eikr_i[iat]);
  }
  buf.add(Rhok_r.first_address(), Rhok_r.last_address());
  buf.add(Rhok_i.first_address(), Rhok_i.last_address());
#else
  eikr.resize(NumPtcls,MaxK);
  eikr_new.resize(MaxK);
  delta_eikr.resize(MaxK);

  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
  buf.add(Rhok.first_address(), Rhok.last_address());
#endif

  buf.add(U.first_address(), U.last_address());
  buf.add(d2U.first_address(), d2U.last_address());
  buf.add(FirstAddressOfdU,LastAddressOfdU);
}

LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf,
                               bool fromscratch)
{
  LogValue=evaluateLog(P,P.G,P.L);
#if defined(USE_REAL_STRUCT_FACTOR)
  for(int iat=0; iat<NumPtcls; iat++)
  {
    std::copy(P.SK->eikr_r[iat],P.SK->eikr_r[iat]+MaxK,eikr_r[iat]);
    std::copy(P.SK->eikr_i[iat],P.SK->eikr_i[iat]+MaxK,eikr_i[iat]);
  }
  buf.put(Rhok_r.first_address(), Rhok_r.last_address());
  buf.put(Rhok_i.first_address(), Rhok_i.last_address());
#else
  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
  buf.put(Rhok.first_address(), Rhok.last_address());
#endif
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}

void LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  buf.get(Rhok_r.first_address(), Rhok_r.last_address());
  buf.get(Rhok_i.first_address(), Rhok_i.last_address());
#else
  buf.get(Rhok.first_address(), Rhok.last_address());
#endif
  buf.get(U.first_address(), U.last_address());
  buf.get(d2U.first_address(), d2U.last_address());
  buf.get(FirstAddressOfdU,LastAddressOfdU);
#if defined(USE_REAL_STRUCT_FACTOR)
  for(int iat=0; iat<NumPtcls; iat++)
  {  
    std::copy(P.SK->eikr_r[iat],P.SK->eikr_r[iat]+MaxK,eikr_r[iat]);
    std::copy(P.SK->eikr_i[iat],P.SK->eikr_i[iat]+MaxK,eikr_i[iat]);
  } 
#else
  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
#endif
}


bool LRTwoBodyJastrow::put(xmlNodePtr cur)
{
  if(skRef == 0)
  {
    app_error() << "  LRTowBodyJastrow should not be used for non periodic systems." << std::endl;
    return false;
  }
  return true;
}

void LRTwoBodyJastrow::resetByHandler(HandlerType* handler)
{
  MaxKshell=handler->MaxKshell;
  Fk_symm.resize(MaxKshell);
  FkbyKK.resize(MaxKshell);
  Fk_0.resize(handler->Fk.size());
  Fk_0 = -1.0 * handler->Fk;
  Fk.resize(Fk_0.size());
  Fk=Fk_0;
  int ki=0;
  int ish=0;
  while(ish<MaxKshell && ki<NumKpts)
  {
    Fk_symm[ish]=Fk[ki];
    FkbyKK[ish]=Fk_symm[ish]*skRef->KLists.ksq[ki];
    for(; ki<skRef->KLists.kshell[ish+1]; ki++)
      Fk[ki]=Fk_symm[ish];
    ++ish;
  }
  MaxK=skRef->KLists.kshell[MaxKshell];
  Kshell.resize(MaxKshell+1);
  std::copy(skRef->KLists.kshell.begin(),skRef->KLists.kshell.begin()+MaxKshell+1, Kshell.begin());

#if defined(USE_REAL_STRUCT_FACTOR)
  if (Rhok_r.size()!=MaxK)
    Rhok_r.resize(MaxK);
  if (Rhok_i.size()!=MaxK)
    Rhok_i.resize(MaxK);
#else
  if (Rhok.size()!=MaxK)
    Rhok.resize(MaxK);
#endif
}

OrbitalBasePtr LRTwoBodyJastrow::makeClone(ParticleSet& tqp) const
{
  LRTwoBodyJastrow* myclone=new LRTwoBodyJastrow(*this);
  myclone->skRef=tqp.SK;
  return myclone;
}
}

