#include "Message/Communicate.h"
#include "LongRange/KContainer.h"

using namespace qmcplusplus;

//Constructor
KContainer::KContainer(ParticleLayout_t& ref): Lattice(ref) { }

//Destructor
KContainer::~KContainer() { }

//Public Methods
// UpdateKLists - call for new k or when lattice changed.
void
KContainer::UpdateKLists(ParticleLayout_t& ref, RealType kc) {
  kcutoff = kc;
  kcut2 = kc*kc;
  Lattice = ref;

  if(kcutoff < 0.0){
    LOGMSG("KContainer initialised with cutoff " << kcutoff);
    OHMMS::Controller->abort();
  }

  FindApproxMMax();
  BuildKLists();
}

//Private Methods:
// FindApproxMMax - compute approximate parallelpiped that surrounds kcutoff
// BuildKLists - Correct mmax and fill lists of k-vectors.
void 
KContainer::FindApproxMMax() {
  //Estimate the size of the parallelpiped that encompases a sphere of kcutoff.
  //mmax is stored as integer translations of the reciprocal cell vectors.
  //Does not depend on orthorhombic cells.
  Matrix<RealType> mmat;
  mmat.resize(3,3);
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++){
      mmat[i][j] = 0.0;
      for(int k=0;k<3;k++)
	mmat[i][j] = mmat[i][j] + Lattice.b(k)[i]* Lattice.b(j)[k];
    }
	  
  TinyVector<RealType,3> x,temp;
  RealType tempr;
  for(int idim=0;idim<3;idim++){
    int i = ((idim-1)%3);
    int j = (idim%3);
    int k = ((idim+1)%3);
    
    x[i] = 1.0;
    x[j] = (mmat[j][k]*mmat[k][i] - mmat[k][k]*mmat[i][j]);
    x[j]/= (mmat[j][j]*mmat[k][k] - mmat[j][k]*mmat[j][k]);
    x[k] = -(mmat[k][i] + mmat[j][k]*x[j])/mmat[k][k];
    
    for(i=0;i<3;i++){
      temp[i] = 0.0;
	for(j=0;j<3;j++)
	  temp[i] += mmat[i][j]*x[j];
    }
    tempr = dot(x,temp);

    mmax[idim] = static_cast<int>(sqrt(4.0*kcut2)) + 1;  
  }
}

void 
KContainer::BuildKLists() {
  TinyVector<int,4> TempActualMax;
  TinyVector<int,3> kvec;
  TinyVector<RealType,3> kvec_cart;
  RealType modk2;

  kpts.clear();
  kpts_cart.clear();
  
  for(int i=0; i<4; i++)
    TempActualMax[i] = 0;
  
  //Loop over guesses for valid k-points.
  for(int i=-mmax[0]; i<=mmax[0]; i++){
    kvec[0] = i;
    for(int j=-mmax[1]; j<=mmax[1]; j++){
      kvec[1] = j;
      for(int k=-mmax[2]; k<=mmax[2]; k++){
	kvec[2] = k;
	      
	//Do not include k=0 in evaluations.
	if(i==0 && j==0 && k==0)continue;

	//Convert kvec to Cartesian
	/*
      	for(int idim=0; idim<3; idim++){
	  kvec_cart[idim] = 0.0;
	  for(int idir=0; idir<3; idir++){
	    kvec_cart[idim]+=kvec[idir]*Lattice.b(idir)[idim];
	  }
	  kvec_cart[idim]*=TWOPI;
	}
	*/
	kvec_cart = Lattice.k_cart(kvec);
	

	//Find modk
	modk2 = dot(kvec_cart,kvec_cart);

	if(modk2>kcut2)continue; //Inside cutoff?

	//This k-point should be added to the list
	kpts.push_back(kvec);
	kpts_cart.push_back(kvec_cart);
	
	//Update record of the allowed maximum translation.
	for(int idim=0; idim<3; idim++)
	  if(abs(kvec[idim]) > TempActualMax[idim])
	    TempActualMax[idim] = abs(kvec[idim]);
	
      }
    }
  }
  
  //Finished searching k-points. Copy list of maximum translations.
  mmax[3] = 0;
  for(int idim=0; idim<3; idim++) {
    mmax[idim] = TempActualMax[idim];
    if(mmax[idim] > mmax[3])
      mmax[3] = mmax[idim];
  }

  //Update a record of the number of k vectors
  numk = kpts.size();

  //Now fill the array that returns the index of -k when given the index of k.
  minusk.resize(numk);
  for(int ki=0; ki<numk; ki++) {
    int kj=0;
    do {
      if(kpts[ki][0] == -kpts[kj][0] 
	 &&kpts[ki][1] == -kpts[kj][1] 
	 &&kpts[ki][2] == -kpts[kj][2])
	break;
      kj++;
    } while(kj<numk);
    if(kj==numk){
      LOGMSG("Error finding -K in list");
      OHMMS::Controller->abort();
    }
    minusk[ki] = kj;
  }
}
