#include "LongRange/StructFact.h"

using namespace ohmmsqmc;

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& ref, RealType kc): PtclRef(ref), 
						       KLists(ref.Lattice) {
  //Update Rhok with new "Lattice" information.
  UpdateNewCell(kc);
}

//Destructor
StructFact::~StructFact() { }

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
StructFact::Update1Part(Position_t rold,Position_t rnew,int GroupID) {
  UpdateRhok(rold,rnew,GroupID);
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
  //Currently only for 1 species!!! Extend this.
  rhok.resize(KLists.numk,tspecies.TotalNum);

  //Zero out rhok
  for(int ki=0; ki<KLists.numk; ki++)
    for(int t=0; t<tspecies.TotalNum; t++)
      rhok(ki,t) = complex<RealType>(0.0,0.0);
  
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
  
  int npart = PtclRef.getTotalNum();
  for(int i=0; i<npart; i++){
    for(int idim=0; idim<3; idim++){
      complex<double> Ctemp;
      //start the recursion with the 111 vector.
      double phi = (PtclRef.R[i])[idim] * k111[idim];
      Ctemp = complex<double>(cos(phi), sin(phi));
      C(idim,KLists.mmax[idim]) = 1.0; // K=0
      //Recursively generate all Cs.
      for(int n=1; n<=KLists.mmax[idim]; n++){
	C(idim,KLists.mmax[idim]+n) = Ctemp*C(idim,KLists.mmax[idim]+n-1);
	C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
      }
    }
    
    //Now add the contribution to Rhok for this particle
    for(int ki=0; ki<KLists.numk; ki++){
      complex<double> temp = 1.0;
      for(int idim=0; idim<3; idim++)
	temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
      rhok(ki,PtclRef.GroupID[i]) += temp;
    }
  } //End particle loop
}


void 
StructFact::UpdateRhok(Position_t rold,Position_t rnew,int GroupID){
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
      temp *= -C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(ki,GroupID) += temp;
  }

  //Prepare for subtracting new position
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

  //Subtract old position
  for(int ki=0; ki<KLists.numk; ki++){
    complex<double> temp = 1.0;
    for(int idim=0; idim<3; idim++)
      temp *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    rhok(ki,GroupID) += temp;
  }
}
