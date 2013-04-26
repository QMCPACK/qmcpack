//////////////////////////////////////////////////////////////////
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
// -*- C++ -*-
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
  Rhok.resize(NumKpts);
  rokbyF.resize(NumPtcls,NumKpts);
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
//     if(handler) resetByHandler();
}

void LRTwoBodyJastrow::reportStatus(ostream& os)
{
}

void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P)
{
  // update handler as well, should there also be a reset?
  skRef=P.SK;
  //Do not update handler and breakup
  //if(handler)
  //{
  //  handler->initBreakup(P);
  //  resetInternals();
  //}
}

LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::evaluateLog(ParticleSet& P,
                              ParticleSet::ParticleGradient_t& G,
                              ParticleSet::ParticleLaplacian_t& L)
{
  RealType sum(0.0);
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("LRTwoBodyJastrow::evaluateLog ");
#else
  //memcopy if necessary but this is not so critcal
  std::copy(P.SK->rhok[0],P.SK->rhok[0]+MaxK,Rhok.data());
  for(int spec1=1; spec1<NumSpecies; spec1++)
    accumulate_elements(P.SK->rhok[spec1],P.SK->rhok[spec1]+MaxK,Rhok.data());
  //Rhok=0.0;
  //for(int spec1=0; spec1<NumSpecies; spec1++)
  //{
  //  const ComplexType* restrict rhok(P.SK->rhok[spec1]);
  //  for(int ki=0; ki<MaxK; ki++) Rhok[ki] += rhok[ki];
  //}
  const KContainer::VContainer_t& Kcart(P.SK->KLists.kpts_cart);
  for(int iat=0; iat<NumPtcls; iat++)
  {
    RealType res(0.0),l(0.0);
    PosType g;
    const ComplexType* restrict eikr_ptr(P.SK->eikr[iat]);
    const ComplexType* restrict rhok_ptr(Rhok.data());
    int ki=0;
    for(int ks=0; ks<MaxKshell; ks++)
    {
      RealType res_k(0.0),l_k(0.0);
      PosType g_k;
      for(; ki<Kshell[ks+1]; ki++,eikr_ptr++,rhok_ptr++)
      {
        RealType rr=((*eikr_ptr).real()*(*rhok_ptr).real()+(*eikr_ptr).imag()*(*rhok_ptr).imag());
        RealType ii=((*eikr_ptr).real()*(*rhok_ptr).imag()-(*eikr_ptr).imag()*(*rhok_ptr).real());
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
#endif
  return 0.5*sum;
//      const KContainer::SContainer_t& ksq(P.SK->KLists.ksq);
//      ValueType sum(0.0);
//      for(int iat=0; iat<NumPtcls; iat++) {
//        ValueType res(0.0),l(0.0);
//        GradType g;
//        const ComplexType* restrict eikr(P.SK->eikr[iat]);
//        for(int ki=0; ki<MaxK; ki++) {
//          ComplexType skp((Fk[ki]*conj(eikr[ki])*Rhok[ki]));
//#if defined(QMC_COMPLEX)
//          res +=  skp;
//          l += ksq[ki]*(Fk[ki]-skp);
//          g += ComplexType(skp.imag(),-skp.real())*kpts[ki];
//#else
//          res +=  skp.real();
//          g += Kcart[ki]*skp.imag();
//          l += ksq[ki]*(Fk[ki]-skp.real());
//#endif
//        }
//        sum+=(U[iat]=res);
//        G[iat]+=(dU[iat]=g);
//        L[iat]+=(d2U[iat]=l);
//      }
//
//      return sum*0.5;
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
  APP_ABORT("LRTwoBodyJastrow::ratio(ParticleSet& P, int iat)");
#else
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  const ComplexType* restrict eikr_ptr(P.SK->eikr[iat]);
  const ComplexType* restrict rhok_ptr(Rhok.data());
  PosType pos(P.R[iat]);
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
//      const Vector<ComplexType>& eikr1(P.SK->eikr_new);
//      const Vector<ComplexType>& del_eikr(P.SK->delta_eikr);
//      //Rhok += del_eikr;
//      for(int ki=0; ki<MaxK; ki++) {
//        //ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
//        ComplexType skp((Fk[ki]*conj(eikr1[ki])*(Rhok[ki]+del_eikr[ki])));
//#if defined(QMC_COMPLEX)
//        curVal +=  skp;
//#else
//        curVal +=  skp.real();
//#endif
//      }
#endif
  return std::exp(curVal-U[iat]);
}


LRTwoBodyJastrow::ValueType
LRTwoBodyJastrow::logRatio(ParticleSet& P, int iat,
                           ParticleSet::ParticleGradient_t& dG,
                           ParticleSet::ParticleLaplacian_t& dL)
{
  NeedToRestore=true;
  const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
  {
    ComplexType* restrict eikr1(eikr_new.data());
    ComplexType* restrict deikr(delta_eikr.data());
    const ComplexType* restrict eikr0(eikr[iat]);
    PosType pos(P.R[iat]);
    RealType c,s;
    for(int ki=0; ki<MaxK; ki++)
    {
      sincos(dot(kpts[ki],pos),&s,&c);
      (*eikr1)=ComplexType(c,s);
      (*deikr++)=(*eikr1++)-(*eikr0++);
    }
  }
  //new Rhok: restored by rejectMove
  Rhok += delta_eikr;
  curVal=0.0;
  curLap=0.0;
  curGrad=0.0;
  const ComplexType* restrict rhok_ptr(Rhok.data());
  const ComplexType* restrict eikr1(eikr_new.data());
  for(int ks=0,ki=0; ks<MaxKshell; ks++)
  {
    RealType v(0.0),l(0.0);
    PosType g;
    for(; ki<Kshell[ks+1]; ki++,eikr1++,rhok_ptr++)
    {
      RealType rr=((*eikr1).real()*(*rhok_ptr).real()+(*eikr1).imag()*(*rhok_ptr).imag());
      RealType ii=((*eikr1).real()*(*rhok_ptr).imag()-(*eikr1).imag()*(*rhok_ptr).real());
      v +=  rr;
      l += 1.0-rr;
      g += ii*kpts[ki];
    }
    curVal += Fk_symm[ks]*v;
    curGrad += Fk_symm[ks]*g;
    curLap += FkbyKK[ks]*l;
  }
  for(int jat=0; jat<NumPtcls; jat++)
  {
    if(jat == iat)
      continue;
    const ComplexType* restrict eikri(delta_eikr.data());
    const ComplexType* restrict eikrj(eikr[jat]);
    RealType v(0.0),l(0.0);
    PosType g;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      RealType v_k(0.0),l_k(0.0);
      PosType g_k;
      for(; ki<Kshell[ks+1]; ki++,eikri++,eikrj++)
      {
        RealType rr=(*eikrj).real()*(*eikri).real()+(*eikrj).imag()*(*eikri).imag();
        RealType ii=(*eikrj).real()*(*eikri).imag()-(*eikrj).imag()*(*eikri).real();
        v_k += rr;
        l_k -= rr;
        g_k += ii*kpts[ki];
      }
      v += Fk_symm[ks]*v_k;
      g += Fk_symm[ks]*g_k;
      l += FkbyKK[ks]*l_k;
    }
    offU[jat]=v;
    offdU[jat]=g;
    offd2U[jat]=l;
    dG[jat] += g;
    dL[jat] += l;
  }
//      for(int jat=0;jat<NumPtcls; jat++) {
//        if(iat==jat) {
//          for(int ki=0; ki<MaxK; ki++) {
//            //ComplexType rhok_new(Rhok[ki]+del_eikr[ki]);
//            //ComplexType skp((Fk[ki]*conj(eikr1[ki])*rhok_new));
//#if defined(QMC_COMPLEX)
//            ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
//            curVal +=  skp;
//            curGrad += ComplexType(skp.imag(),-skp.real())*kpts[ki];
//            curLap += ksq[ki]*(Fk[ki]-skp);
//#else
//            RealType skp_r=Fk[ki]*(eikr1[ki].real()*Rhok[ki].real()+eikr1[ki].imag()*Rhok[ki].imag());
//            RealType skp_i=Fk[ki]*(eikr1[ki].real()*Rhok[ki].imag()-eikr1[ki].imag()*Rhok[ki].real());
//            curVal +=  skp_r;
//            curLap += ksq[ki]*(Fk[ki]-skp_r);
//            curGrad += kpts[ki]*skp_i;
//            //curVal +=  skp.real();
//            //curLap += ksq[ki]*(Fk[ki]-skp.real());
//            //curGrad += skp.imag()*kpts[ki];
//#endif
//          }
//        } else {
//          const ComplexType* restrict eikrj(P.SK->eikr[jat]);
//          GradType g;
//          ValueType l(0.0), v(0.0);
//          for(int ki=0; ki<MaxK; ki++) {
//#if defined(QMC_COMPLEX)
//            ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
//            GradType dg(skp.imag()*kpts[ki]);
//            ValueType dl(skp.real()*ksq[ki]);
//            v += skp.real();
//            g +=dg;
//            l -= dl;
//#else
//            ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
//            //GradType dg(skp.imag()*kpts[ki]);
//            //ValueType dl(skp.real()*ksq[ki]);
//            v += skp.real();
//            l -= skp.real()*ksq[ki];
//            g += skp.imag()*kpts[ki];
//#endif
//            //dG[jat] += Fk[ki]*skp.imag()*kpts[ki];
//            //dL[jat] -= Fk[ki]*skp.real()*ksq[ki];
//          }
//          offU[jat]=v;
//          offdU[jat]=g;
//          offd2U[jat]=l;
//          dG[jat] += g;
//          dL[jat] += l;
//        }
//      }
//
  dG[iat] += offdU[iat] = curGrad-dU[iat];
  dL[iat] += offd2U[iat] = curLap-d2U[iat];
  return offU[iat] = curVal-U[iat];
}

void LRTwoBodyJastrow::restore(int iat)
{
  //substract the addition in logRatio
  if(NeedToRestore)
    Rhok -= delta_eikr;
}

void LRTwoBodyJastrow::acceptMove(ParticleSet& P, int iat)
{
  std::copy(eikr_new.data(),eikr_new.data()+MaxK,eikr[iat]);
  U += offU;
  dU += offdU;
  d2U += offd2U;
}

void LRTwoBodyJastrow::update(ParticleSet& P,
                              ParticleSet::ParticleGradient_t& dG,
                              ParticleSet::ParticleLaplacian_t& dL,
                              int iat)
{
  app_error() << "LRTwoBodyJastrow::update is INCOMPLETE " << endl;
}


LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  LogValue=evaluateLog(P,P.G,P.L);
  eikr.resize(NumPtcls,MaxK);
  eikr_new.resize(MaxK);
  delta_eikr.resize(MaxK);
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("LRTwoBodyJastrow::registerData");
#else
  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
#endif
  buf.add(Rhok.first_address(), Rhok.last_address());
  buf.add(U.first_address(), U.last_address());
  buf.add(d2U.first_address(), d2U.last_address());
  buf.add(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}

LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
                               bool fromscratch)
{
  LogValue=evaluateLog(P,P.G,P.L);
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("LRTwoBodyJastrow::updateBuffer");
#else
  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
#endif
  buf.put(Rhok.first_address(), Rhok.last_address());
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}

void LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.get(Rhok.first_address(), Rhok.last_address());
  buf.get(U.first_address(), U.last_address());
  buf.get(d2U.first_address(), d2U.last_address());
  buf.get(FirstAddressOfdU,LastAddressOfdU);
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("LRTwoBodyJastrow::copyFromBuffer");
#else
  for(int iat=0; iat<NumPtcls; iat++)
    std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
#endif
}

LRTwoBodyJastrow::RealType
LRTwoBodyJastrow::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.put(Rhok.first_address(), Rhok.last_address());
  buf.put(U.first_address(), U.last_address());
  buf.put(d2U.first_address(), d2U.last_address());
  buf.put(FirstAddressOfdU,LastAddressOfdU);
  return LogValue;
}


bool LRTwoBodyJastrow::put(xmlNodePtr cur)
{
  if(skRef == 0)
  {
    app_error() << "  LRTowBodyJastrow should not be used for non periodic systems." << endl;
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
    Fk_symm[ish]=Fk[ki];//-1.0*handler->Fk[ki];
    FkbyKK[ish]=Fk_symm[ish]*skRef->KLists.ksq[ki];
    for(; ki<skRef->KLists.kshell[ish+1]; ki++)
      Fk[ki]=Fk_symm[ish];
    ++ish;
  }
  MaxK=skRef->KLists.kshell[MaxKshell];
  Kshell.resize(MaxKshell+1);
  std::copy(skRef->KLists.kshell.begin(),skRef->KLists.kshell.begin()+MaxKshell+1, Kshell.begin());
  if (Rhok.size()!=MaxK)
    Rhok.resize(MaxK);
}

OrbitalBasePtr LRTwoBodyJastrow::makeClone(ParticleSet& tqp) const
{
  LRTwoBodyJastrow* myclone=new LRTwoBodyJastrow(*this);
  myclone->skRef=tqp.SK;
  return myclone;
}
}

