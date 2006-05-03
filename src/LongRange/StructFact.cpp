#include "LongRange/StructFact.h"

using namespace qmcplusplus;

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& ref, RealType kc): PtclRef(ref), KLists(ref.Lattice) {
  //Update Rhok with new "Lattice" information.
  UpdateNewCell(kc);
}

//Copy Constructor
StructFact::StructFact(const StructFact &ref): PtclRef(ref.PtclRef), KLists(ref.PtclRef.Lattice) {
  // Lattices are forced to match in PtclRef initialization.
  // The KLists constructor doesn't generate the lists. It merely sets the lattice reference.
  // "=" is defined for all data members that need to be copied.
  KLists = ref.KLists; //= checks for same cutoff and returns with no cost if equal.
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
  //Compute the entire Rhok 
  FillRhok();
}

void 
StructFact::Update1Part(Position_t rold,Position_t rnew,int iat,int GroupID) {
  UpdateRhok(rold,rnew,iat,GroupID);
}

void 
StructFact::UpdateAllPart() {
  FillRhok();
}

//Private Methods
// FillRhok
// UpdateRhok
void 
StructFact::FillRhok() {
  SpeciesSet& tspecies(PtclRef.getSpeciesSet());
	
  //Evaluate "Rho_k" using fast method
  //Ken's breakup doc., section 5.1.
  //This is the structure-factor of the ion coordinates.
  //rhok.resize(KLists.numk,tspecies.TotalNum);
  rhok.resize(tspecies.TotalNum,KLists.numk);
  eikr.resize(PtclRef.getTotalNum(),KLists.numk);

  //Zero out rhok
  for(int ki=0; ki<KLists.numk; ki++) 
    for(int t=0; t<tspecies.TotalNum; t++)
      rhok(t,ki) = complex<RealType>(0.0,0.0);
     
  TinyVector<double,3> k111; //k=1*b1 + 1*b2 + 1*b3
  //Convert to Cartesian  
  for(int idim=0; idim<3; idim++){
    k111[idim] = 0.0;
    for(int idir=0; idir<3; idir++){
      k111[idim] += PtclRef.Lattice.b(idir)[idim];
    }
    k111[idim] *= TWOPI;
  }
  
  
  //Allocate the 'C' arrays.
  //C needs to be allocated from [0..2][-max..max]
  //Actually, map -max..max to 0..2max
  //where max is the largest of KLists.mmax[3]
  int maxdim = max(KLists.mmax[0],max(KLists.mmax[1],KLists.mmax[2]));
  Matrix<complex<double> > C;
  C.resize(3,2*maxdim+1);
  
  int npart = PtclRef.getTotalNum();
  for(int i=0; i<npart; i++){
    for(int idim=0; idim<3; idim++){
      complex<double> Ctemp;
      //start the recursion with the 111 vector.
      double phi = (PtclRef.R[i])[idim] * k111[idim];
      Ctemp = complex<double>(cos(phi), sin(phi));
      C(idim,KLists.mmax[idim]) = 1.0; // K=0 term
      //Recursively generate all Cs.
      for(int n=1; n<=KLists.mmax[idim]; n++){
	C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
	C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
      }
    }
    
    //Now add the contribution to Rhok for this particle
    for(int ki=0; ki<KLists.numk; ki++){
      eikr(i,ki) = complex<RealType>(1.0,0.0); //Initialize 
      for(int idim=0; idim<3; idim++)
	eikr(i,ki) *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
      rhok(PtclRef.GroupID[i],ki) += eikr(i,ki);
    }
  } //End particle loop
}


void 
StructFact::UpdateRhok(Position_t rold,Position_t rnew,int iat,int GroupID){
  TinyVector<double,3> k111; //k=1*b1 + 1*b2 + 1*b3
  //Convert to Cartesian
  
  for(int idim=0; idim<3; idim++){
    k111[idim] = 0.0;
    for(int idir=0; idir<3; idir++){
      k111[idim] += PtclRef.Lattice.b(idir)[idim];
    }
    k111[idim] *= TWOPI;
  }

  //Allocate the 'C' arrays.
  //C needs to be allocated from [0..2][-max..max]
  //Actually, map -max..max to 0..2max
  Matrix<complex<double> > C;
  C.resize(3,2*KLists.mmax[3]+1);

  //Prepare for subtracting old position
  for(unsigned int idim=0; idim<3; idim++){
    complex<double> Ctemp;
    //start the recursion with the 111 vector.
    double phi = rold[idim] * k111[idim];
    Ctemp = complex<double>(cos(phi), sin(phi));
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
    Ctemp = complex<double>(cos(phi), sin(phi));
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
