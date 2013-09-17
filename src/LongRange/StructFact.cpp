#include <LongRange/StructFact.h>
#include <config/stdlib/math.h>
#include <Numerics/e2iphi.h>
#include <simd/vmath.hpp>
#include <Numerics/OhmmsBlas.h>
#include <qmc_common.h>

namespace qmcplusplus
{

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& P, RealType kc):
  DoUpdate(false),SuperCellEnum(SUPERCELL_BULK)
{
  UpdateNewCell(P,kc);

  if(qmc_common.use_ewald && P.LRBox.SuperCellEnum == SUPERCELL_SLAB)
  {
    app_log() << "  Setting StructFact::SuperCellEnum=SUPERCELL_SLAB " << endl;
    SuperCellEnum=SUPERCELL_SLAB;
  }
}

//Destructor
StructFact::~StructFact() { }

void
StructFact::UpdateNewCell(ParticleSet& P, RealType kc)
{
  //Generate the lists of k-vectors
  KLists.UpdateKLists(P.LRBox,kc);
  //resize any arrary
  resize(P.getSpeciesSet().size(),P.getTotalNum(),KLists.numk);
  //Compute the entire Rhok
  FillRhok(P);
}

//void StructFact::resize()
//{
//  SpeciesSet& tspecies(PtclRef.getSpeciesSet());
//  phiV.resize(KLists.numk);
//  phiM.resize(PtclRef.getTotalNum(),KLists.numk);
//#if defined(USE_REAL_STRUCT_FACTOR)
//  rhok_r.resize(tspecies.TotalNum,KLists.numk);
//  rhok_i.resize(tspecies.TotalNum,KLists.numk);
//  eikr_r.resize(PtclRef.getTotalNum(),KLists.numk);
//  eikr_i.resize(PtclRef.getTotalNum(),KLists.numk);
//  eikr_r_temp.resize(KLists.numk);
//  eikr_i_temp.resize(KLists.numk);
//#else
//  rhok.resize(tspecies.TotalNum,KLists.numk);
//  eikr.resize(PtclRef.getTotalNum(),KLists.numk);
//  eikr_temp.resize(KLists.numk);
//#endif
//  //int maxdim = std::max(KLists.mmax[0],std::max(KLists.mmax[1],KLists.mmax[2]));
//  int maxdim=KLists.mmax[DIM];
//  C.resize(DIM,2*maxdim+1);
//}

void StructFact::resize(int ns, int nptcl, int nkpts)
{
  phiV.resize(nkpts);
  phiM.resize(nptcl,nkpts);
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r.resize(ns,nkpts);
  rhok_i.resize(ns,nkpts);
  eikr_r.resize(nptcl,nkpts);
  eikr_i.resize(nptcl,nkpts);
  eikr_r_temp.resize(nkpts);
  eikr_i_temp.resize(nkpts);
#else
  rhok.resize(ns,nkpts);
  eikr.resize(nptcl,nkpts);
  eikr_temp.resize(nkpts);
#endif
  //int maxdim = std::max(KLists.mmax[0],std::max(KLists.mmax[1],KLists.mmax[2]));
  int maxdim=KLists.mmax[DIM];
  C.resize(DIM,2*maxdim+1);
}


//void
//StructFact::Update1Part(const PosType& rold,const PosType& rnew,int iat,int GroupID) {
//  UpdateRhok(rold,rnew,iat,GroupID);
//}

void
StructFact::UpdateAllPart(ParticleSet& P)
{
  //if(!DoUpdate) FillRhok();
  FillRhok(P);
}

///** Experimental functions to support real storage for structure factor
// */
//namespace simd
//{
//  template<typename T>
//    inline void add(int n, const T* restrict in, T* restrict out)
//    {
//      for(int i=0; i<n; ++i) out[i]+=in[i];
//    }
//
//  template<typename T>
//    inline void get_phase(int n, const T* restrict kpts, const T* restrict xyz, T* restrict phi)
//    {
//      T x=xyz[0]; T y=xyz[1]; T z=xyz[2];
//      for(int i=0; i<n; ++i)
//        phi[i]=x*kpts[i*3]+y*kpts[i*3+1]+z*kpts[i*3+2];
//    }
//
//  template<typename AT, typename BT, typename CT>
//  inline void get_phase(const AT& kpts, const BT& pos, CT& phase)
//  {
//    const char transa = 'T';
//    const char transb = 'N';
//    const double zone(1.0);
//    const double zero(0.0);
//    dgemm(transa, transb
//        , phase.cols(), phase.rows(), 3
//        , zone
//        , &(kpts[0][0]), 3
//        , &(pos[0][0]), 3
//        , zero, phase.data(), phase.rows());
//  }
//}

/** evaluate rok per species, eikr  per particle
 */
void
StructFact::FillRhok(ParticleSet& P)
{
  int npart = P.getTotalNum();
#if defined(QMC_SK_USE_RECURSIVE)
  rhok=0.0;
  for(int i=0; i<npart; i++)
  {
    //operate with a reduced positon
    PosType tau_red=P.LRBox.toUnit(P.R[i]);
    for(int idim=0; idim<DIM; idim++)
    {
      RealType phi=TWOPI*tau_red[idim];
      ComplexType ctemp(std::cos(phi),std::sin(phi));
      C(idim,KLists.mmax[idim])=1.0;
      for(int n=1; n<=KLists.mmax[idim]; n++)
      {
        C(idim,KLists.mmax[idim]+n) = ctemp*C(idim,KLists.mmax[idim]+n-1);
        C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
      }
    }
    ComplexType* restrict eikr_ref=eikr[i];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      eikr_ref[ki]=C(0,KLists.kpts[ki][0]+KLists.mmax[0]);
      for(idim=1; idim<DIM; id++)
        eikr_ref[ki] *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    }
    accumulate_elements(eikr_ref,eikr_ref+KLists.numk,rhok[P.GroupID[i]]);
    //valid version only with orthorohmbic cell, generalized to any cell above
    //  for(int idim=0; idim<3; idim++){
    //    complex<double> Ctemp;
    //    //start the recursion with the 111 vector.
    //    double phi = (P.R[i])[idim] * k111[idim];
    //    Ctemp = complex<double>(std::cos(phi), std::sin(phi));
    //    C(idim,KLists.mmax[idim]) = 1.0; // K=0 term
    //    //Recursively generate all Cs.
    //    for(int n=1; n<=KLists.mmax[idim]; n++){
    //      C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
    //      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
    //    }
    //  }
    //Now add the contribution to Rhok for this particle
    //for(int ki=0; ki<KLists.numk; ki++){
    //  eikr(i,ki) = ComplexType(1.0,0.0); //Initialize
    //  for(int idim=0; idim<3; idim++)
    //    eikr(i,ki) *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    //  rhok(P.GroupID[i],ki) += eikr(i,ki);
    //}
  } //End particle loop
#else
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r=0.0;
  rhok_i=0.0;
  //algorithmA
  const int nk=KLists.numk;
  for(int i=0; i<npart; ++i)
  {
    //defined in this file
    simd::get_phase(nk,&(KLists.kpts_cart[0][0]), P.R[i].data(), phiV.data());
    //get_phase simply encapsulate this
    //PosType pos(P.R[i]);
    //for(int ki=0; ki<KLists.numk; ki++)
    //  phiV[ki]=dot(KLists.kpts_cart[ki],pos);
    eval_e2iphi(nk, phiV.data(), eikr_r[i], eikr_i[i]);
    simd::add(nk,eikr_r[i],rhok_r[P.GroupID[i]]);
    simd::add(nk,eikr_i[i],rhok_i[P.GroupID[i]]);
  }
  //use dgemm: vtune shows algorithmA is better
  //simd::get_phase(KLists.kpts_cart,P.R,phiM);
  //eval_e2iphi(phiM.size(), phiM.data(), eikr_r.data(), eikr_i.data());
  //for(int i=0; i<npart; ++i)
  //  simd::add(nk,eikr_r[i],rhok_r[P.GroupID[i]]);
  //for(int i=0; i<npart; ++i)
  //  simd::add(nk,eikr_i[i],rhok_i[P.GroupID[i]]);
#else
  rhok=0.0;
  for(int i=0; i<npart; i++)
  {
    PosType pos(P.R[i]);
    RealType s,c;//get sin and cos
    ComplexType* restrict eikr_ref=eikr[i];
    ComplexType* restrict rhok_ref=rhok[P.GroupID[i]];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      //RealType phi(dot(KLists.kpts_cart[ki],pos));
      //eikr_ref[ki] = ComplexType(std::cos(phi),std::sin(phi));
      sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
      eikr_ref[ki]=ComplexType(c,s);
      rhok_ref[ki]+= eikr_ref[ki];
    }
  }
#endif
#endif
}


void
StructFact::UpdateRhok(const PosType& rold,const PosType& rnew,int iat,int GroupID)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("WHO IS USING UpdateRhok");
#else
//  TinyVector<double,DIM> k111; //k=1*b1 + 1*b2 + 1*b3
//  //Convert to Cartesian
//  for(int idim=0; idim<DIM; idim++)
//  {
//    k111[idim] = 0.0;
//    for(int idir=0; idir<DIM; idir++)
//    {
//      k111[idim] += P.LRBox.b(idir)[idim];
//    }
//    k111[idim] *= TWOPI;
//  }
//  //Prepare for subtracting old position
//  for(unsigned int idim=0; idim<DIM; idim++)
//  {
//    complex<double> Ctemp;
//    //start the recursion with the 111 vector.
//    double phi = rold[idim] * k111[idim];
//    Ctemp = complex<double>(std::cos(phi), std::sin(phi));
//    C(idim,KLists.mmax[idim]) = 1.0; // K=0
//    //Recursively generate all Cs.
//    for(int n=1; n<=KLists.mmax[idim]; n++)
//    {
//      C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
//      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
//    }
//  }
//  //Subtract old position
//  for(int ki=0; ki<KLists.numk; ki++)
//  {
//    complex<double> temp = 1.0;
//    for(int idim=0; idim<DIM; idim++)
//      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
//    rhok(GroupID,ki) -= temp;
//  }
//  //Prepare for adding new position
//  for(unsigned int idim=0; idim<DIM; idim++)
//  {
//    complex<double> Ctemp;
//    //start the recursion with the 111 vector.
//    double phi = rnew[idim] * k111[idim];
//    Ctemp = complex<double>(std::cos(phi), std::sin(phi));
//    C(idim,KLists.mmax[idim]) = 1.0; // K=0
//    //Recursively generate all Cs.
//    for(int n=1; n<=KLists.mmax[idim]; n++)
//    {
//      C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
//      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
//    }
//  }
//  //Add new position
//  for(int ki=0; ki<KLists.numk; ki++)
//  {
//    complex<double> temp = 1.0;
//    for(int idim=0; idim<DIM; idim++)
//      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
//    rhok(GroupID,ki) += temp;
//    eikr(iat,ki) = temp;
//  }
#endif
}

void StructFact::makeMove(int active, const PosType& pos)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  for(int ki=0; ki<KLists.numk; ki++)
    phiV[ki]=dot(KLists.kpts_cart[ki],pos);
  eval_e2iphi(KLists.numk, phiV.data(), eikr_r_temp.data(), eikr_i_temp.data());
#else
  RealType s,c;//get sin and cos
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
    eikr_temp[ki]=ComplexType(c,s);
  }
#endif
}

void StructFact::acceptMove(int active, int gid)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("NOT DONE WITH StructFact::acceptMove");
#else
  //cout << "StructFact::acceptMove " << active << endl;
  //APP_ABORT("StructFact::acceptMove should not be used yet");
  ComplexType* restrict eikr_ptr=eikr[active];
  ComplexType* restrict rhok_ptr(rhok[gid]);
  //ComplexType* restrict rhok_ptr(rhok[P.GroupID[active]]);
  //const ComplexType* restrict t(eikr_ref.data());
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    //(*rho_ptr++) += (*t)-(*eikr_ptr);
    //*eikr_ptr++ = *t++;
    rhok_ptr[ki] += (eikr_temp[ki]-eikr_ptr[ki]);
    eikr_ptr[ki]=eikr_temp[ki];
  }
#endif
}

void StructFact::rejectMove(int active, int gid)
{
  //APP_ABORT("StructFact::rejectMove should not be used yet");
  //do nothing
}
}
