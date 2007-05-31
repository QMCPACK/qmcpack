#include "LongRange/StructFact.h"
using namespace qmcplusplus;

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& ref, RealType kc): PtclRef(ref), KLists(ref.Lattice) 
{
  //Update Rhok with new "Lattice" information.
  UpdateNewCell(kc);
}

/** Copy Constructor
 *
 * Lattices are forced to match in PtclRef initialization.
 * The KLists constructor doesn't generate the lists. It merely sets the lattice reference.
 * "=" is defined for all data members that need to be copied.
 */
StructFact::StructFact(const StructFact &ref): PtclRef(ref.PtclRef), KLists(ref.PtclRef.Lattice) 
{
  KLists = ref.KLists; //= checks for same cutoff and returns with no cost if equal.
  resize();
  rhok = ref.rhok;
}

//Destructor
StructFact::~StructFact() { }

//Overload the assignment operator
StructFact& 
StructFact::operator=(const StructFact &ref) {
  if(this != &ref){
    //Copy data from ref.
    //Check that PtclRefs are the same and k-shells
    //If PtclRef match then Lattice match. All KLists should then 
    //match if cutoffs do.
    if(&PtclRef != &ref.PtclRef){
      LOGMSG("ERROR: tried to copy SK with different PtclRef");
      return *this;
    }
    // "=" is defined for all data members that need to be copied.
    KLists = ref.KLists; //= checks for same cutoff and returns with no cost if equal.
    resize();
    rhok = ref.rhok;
  }
  return *this; //Allows assignment chaining 
}

//Public Methods:
// UpdateNewCell - recompute Rhok if lattice changed
// Update1Part - update Rhok if 1 particle moved
// UpdateAllPart - update Rhok if all particles moved
void
StructFact::UpdateNewCell(RealType kc) {
  //Generate the lists of k-vectors
  KLists.UpdateKLists(PtclRef.Lattice,kc);
  //resize any arrary
  resize();
  //Compute the entire Rhok 
  FillRhok();
}

void StructFact::resize()
{
  SpeciesSet& tspecies(PtclRef.getSpeciesSet());
  rhok.resize(tspecies.TotalNum,KLists.numk);
  eikr.resize(PtclRef.getTotalNum(),KLists.numk);
  int maxdim = std::max(KLists.mmax[0],std::max(KLists.mmax[1],KLists.mmax[2]));
  C.resize(3,2*maxdim+1);
}

//void 
//StructFact::Update1Part(const PosType& rold,const PosType& rnew,int iat,int GroupID) {
//  UpdateRhok(rold,rnew,iat,GroupID);
//}

void 
StructFact::UpdateAllPart() {
  FillRhok();
}


/** evaluate rok per species, eikr  per particle
 */
void 
StructFact::FillRhok() {
  rhok=0.0;
  int npart = PtclRef.getTotalNum();
#if defined(QMC_SK_USE_RECURSIVE)
  for(int i=0; i<npart; i++)
  {
    //operate with a reduced positon
    PosType tau_red=PtclRef.Lattice.toUnit(PtclRef.R[i]);
    for(int idim=0; idim<3; idim++)
    {
      RealType phi=TWOPI*tau_red[idim];
      ComplexType ctemp(std::cos(phi),std::sin(phi));
      C(idim,KLists.mmax[idim])=1.0;
      for(int n=1; n<=KLists.mmax[idim]; n++){
        C(idim,KLists.mmax[idim]+n) = ctemp*C(idim,KLists.mmax[idim]+n-1);
        C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
      }
    }
    ComplexType* restrict eikr_ref=eikr[i];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      eikr_ref[ki]=C(0,KLists.kpts[ki][0]+KLists.mmax[0])
        *C(1,KLists.kpts[ki][1]+KLists.mmax[1])
        *C(2,KLists.kpts[ki][2]+KLists.mmax[2]);
    }
    accumulate_elements(eikr_ref,eikr_ref+KLists.numk,rhok[PtclRef.GroupID[i]]);
    //valid version only with orthorohmbic cell, generalized to any cell above
    //  for(int idim=0; idim<3; idim++){
    //    complex<double> Ctemp;
    //    //start the recursion with the 111 vector.
    //    double phi = (PtclRef.R[i])[idim] * k111[idim];
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
    //  rhok(PtclRef.GroupID[i],ki) += eikr(i,ki);
    //}
  } //End particle loop
#else
  for(int i=0; i<npart; i++)
  {
    PosType pos(PtclRef.R[i]);
    ComplexType* restrict eikr_ref=eikr[i];
    ComplexType* restrict rhok_ref=rhok[PtclRef.GroupID[i]];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      RealType phi(dot(KLists.kpts_cart[ki],pos));
//#if defined(HAVE_SINCOS)
//      sincos(phi,&(eikr_ref[ki].imag()),&(eikr_ref[ki].real()));
//#else
//      eikr_ref[ki] = ComplexType(std::cos(phi),std::sin(phi));
//#endif
      eikr_ref[ki] = ComplexType(std::cos(phi),std::sin(phi));
      rhok_ref[ki]+= eikr_ref[ki];
    }
  }
#endif
}


void 
StructFact::UpdateRhok(const PosType& rold,const PosType& rnew,int iat,int GroupID){

  cout << "##### StructFact::UpdateRhok(const PosType& rold,const PosType& rnew USED" << endl;
  TinyVector<double,3> k111; //k=1*b1 + 1*b2 + 1*b3
  //Convert to Cartesian
  
  for(int idim=0; idim<3; idim++){
    k111[idim] = 0.0;
    for(int idir=0; idir<3; idir++){
      k111[idim] += PtclRef.Lattice.b(idir)[idim];
    }
    k111[idim] *= TWOPI;
  }

  //Prepare for subtracting old position
  for(unsigned int idim=0; idim<3; idim++){
    complex<double> Ctemp;
    //start the recursion with the 111 vector.
    double phi = rold[idim] * k111[idim];
    Ctemp = complex<double>(std::cos(phi), std::sin(phi));
    C(idim,KLists.mmax[idim]) = 1.0; // K=0
    //Recursively generate all Cs.
    for(int n=1; n<=KLists.mmax[idim]; n++){
      C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
    }
  }

  //Subtract old position
  for(int ki=0; ki<KLists.numk; ki++){
    complex<double> temp = 1.0;
    for(int idim=0; idim<3; idim++)
      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(GroupID,ki) -= temp;
  }

  //Prepare for adding new position
  for(unsigned int idim=0; idim<3; idim++){
    complex<double> Ctemp;
    //start the recursion with the 111 vector.
    double phi = rnew[idim] * k111[idim];
    Ctemp = complex<double>(std::cos(phi), std::sin(phi));
    C(idim,KLists.mmax[idim]) = 1.0; // K=0
    //Recursively generate all Cs.
    for(int n=1; n<=KLists.mmax[idim]; n++){
      C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
    }
  }

  //Add new position
  for(int ki=0; ki<KLists.numk; ki++){
    complex<double> temp = 1.0;
    for(int idim=0; idim<3; idim++)
      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(GroupID,ki) += temp;
    eikr(iat,ki) = temp;
  }
}

//void StructFact::makeMove(int iat, const PosType& pos) {
//  cout << "Nobody should call this! " << endl;
//  const ComplexType* restrict eikr0(eikr[iat]);
//#if defined(QMC_SK_USE_RECURSIVE)
//  PosType tau_red=PtclRef.Lattice.toUnit(pos);
//  for(int idim=0; idim<3; idim++)
//  {
//    RealType phi=TWOPI*tau_red[idim];
//    ComplexType ctemp(std::cos(phi),std::sin(phi));
//    C(idim,KLists.mmax[idim])=1.0;
//    for(int n=1; n<=KLists.mmax[idim]; n++){
//      C(idim,KLists.mmax[idim]+n) = ctemp*C(idim,KLists.mmax[idim]+n-1);
//      C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
//    }
//  }
//  for(int ki=0; ki<KLists.numk; ki++)
//  {
//    eikr_new[ki]=C(0,KLists.kpts[ki][0]+KLists.mmax[0])
//      *C(1,KLists.kpts[ki][1]+KLists.mmax[1])
//      *C(2,KLists.kpts[ki][2]+KLists.mmax[2]);
//    delta_eikr[ki]=eikr_new[ki]-eikr0[ki];
//  }
//#else
//  for(int ki=0; ki<KLists.numk; ki++){
//    RealType kdotr(dot(KLists.kpts_cart[ki],pos));
//    eikr_new[ki]=ComplexType(std::cos(kdotr),std::sin(kdotr));
//    delta_eikr[ki]=eikr_new[ki]-eikr0[ki];
//  }
//#endif
//}
//
//void StructFact::acceptMove(int iat) {
//  std::copy(eikr_new.begin(),eikr_new.end(),eikr[iat]);
//  ComplexType* restrict rhok_ptr(rhok[PtclRef.GroupID[iat]]);
//  for(int ki=0; ki<KLists.numk; ki++){
//    rhok_ptr[ki]+= delta_eikr[ki];
//  }
//}
//
//void StructFact::rejectMove(int iat) {
//  //do nothing
//}
