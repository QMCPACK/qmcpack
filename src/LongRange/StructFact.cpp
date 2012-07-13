#include <LongRange/StructFact.h>
#include <config/stdlib/math.h>
namespace qmcplusplus
{

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& ref, RealType kc): 
  DoUpdate(false),PtclRef(ref), KLists(ref.Lattice) 
{
  //Update Rhok with new "Lattice" information.
  UpdateNewCell(kc);
}

///** Copy Constructor
// *
// * Lattices are forced to match in PtclRef initialization.
// * The KLists constructor doesn't generate the lists. It merely sets the lattice reference.
// * "=" is defined for all data members that need to be copied.
// */
//StructFact::StructFact(const StructFact &ref): 
// DoUpdate(ref.DoUpdate),PtclRef(ref.PtclRef), KLists(ref.PtclRef.Lattice) 
//{
//  KLists = ref.KLists; //= checks for same cutoff and returns with no cost if equal.
//  resize();
//  rhok = ref.rhok;
//}

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
  eikr_temp.resize(KLists.numk);
  //int maxdim = std::max(KLists.mmax[0],std::max(KLists.mmax[1],KLists.mmax[2]));
  int maxdim=KLists.mmax[DIM];
  C.resize(DIM,2*maxdim+1);
}

//void 
//StructFact::Update1Part(const PosType& rold,const PosType& rnew,int iat,int GroupID) {
//  UpdateRhok(rold,rnew,iat,GroupID);
//}

void 
StructFact::UpdateAllPart() {
  //if(!DoUpdate) FillRhok();
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
    for(int idim=0; idim<DIM; idim++)
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
      eikr_ref[ki]=C(0,KLists.kpts[ki][0]+KLists.mmax[0]);
      for(idim=1;idim<DIM; id++)
        eikr_ref[ki] *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
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
    RealType s,c;//get sin and cos
    ComplexType* restrict eikr_ref=eikr[i];
    ComplexType* restrict rhok_ref=rhok[PtclRef.GroupID[i]];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      //RealType phi(dot(KLists.kpts_cart[ki],pos));
      //eikr_ref[ki] = ComplexType(std::cos(phi),std::sin(phi));
      sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
      eikr_ref[ki]=ComplexType(c,s);
      rhok_ref[ki]+= eikr_ref[ki];
    }
  }

  //RealType s,c;//get sin and cos
  //for(int i=0; i<npart; i++)
  //{
  //  const PosType& pos(PtclRef.R[i]);
  //  ComplexType* restrict eikr_ref=eikr[i];
  //  for(int ki=0; ki<KLists.numk; ki++)
  //  {
  //    sincos( dot(KLists.kpts_cart[ki],pos) , &s, &c);
  //    eikr_ref[ki]=ComplexType(c,s);
  //  }
  //  accumulate_elements(eikr[i],eikr[i]+KLists.numk,rhok[PtclRef.GroupID[i]]);
  //}

  //RealType KdotP[KLists.numk];
  //for(int ig=0; ig<PtclRef.groups(); ++ig)
  //{
  //  for(int i=PtclRef.first(ig); i<PtclRef.last(ig); ++i)
  //  {
  //    for(int ki=0; ki<KLists.numk; ki++)
  //      KdotP[ki]=dot(KLists.kpts_cart[ki],PtclRef.R[i]);
  //    eval_e2iphi(KLists.numk,KdotP,eikr[i]);
  //    accumulate_elements(eikr[i],eikr[i]+KLists.numk,rhok[ig]);
  //  }

#endif
}


void 
StructFact::UpdateRhok(const PosType& rold,const PosType& rnew,int iat,int GroupID){

  TinyVector<double,DIM> k111; //k=1*b1 + 1*b2 + 1*b3
  //Convert to Cartesian
  
  for(int idim=0; idim<DIM; idim++){
    k111[idim] = 0.0;
    for(int idir=0; idir<DIM; idir++){
      k111[idim] += PtclRef.Lattice.b(idir)[idim];
    }
    k111[idim] *= TWOPI;
  }

  //Prepare for subtracting old position
  for(unsigned int idim=0; idim<DIM; idim++){
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
    for(int idim=0; idim<DIM; idim++)
      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(GroupID,ki) -= temp;
  }

  //Prepare for adding new position
  for(unsigned int idim=0; idim<DIM; idim++){
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
    for(int idim=0; idim<DIM; idim++)
      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(GroupID,ki) += temp;
    eikr(iat,ki) = temp;
  }
}

void StructFact::makeMove(int active, const PosType& pos) 
{
  //APP_ABORT("StructFact::makeMove should not be used yet");
  //cout << "StructFact::makeMove " << active << " " << pos << endl;
  RealType s,c;//get sin and cos
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
    eikr_temp[ki]=ComplexType(c,s);
  }
}

void StructFact::acceptMove(int active) 
{
  //cout << "StructFact::acceptMove " << active << endl;
  //APP_ABORT("StructFact::acceptMove should not be used yet");
  ComplexType* restrict eikr_ptr=eikr[active];
  ComplexType* restrict rhok_ptr(rhok[PtclRef.GroupID[active]]);
  //const ComplexType* restrict t(eikr_ref.data());
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    //(*rho_ptr++) += (*t)-(*eikr_ptr);
    //*eikr_ptr++ = *t++;
    rhok_ptr[ki] += (eikr_temp[ki]-eikr_ptr[ki]);
    eikr_ptr[ki]=eikr_temp[ki];
  }
}

void StructFact::rejectMove(int active) {
  //APP_ABORT("StructFact::rejectMove should not be used yet");
  //do nothing
}
}
