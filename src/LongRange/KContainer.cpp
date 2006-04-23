#include "Message/Communicate.h"
#include "LongRange/KContainer.h"
#include <map>

using namespace qmcplusplus;

//Constructor
KContainer::KContainer(ParticleLayout_t& ref): Lattice(ref),kcutoff(0.0) { }

//Destructor
KContainer::~KContainer() { }

//Overloaded assignment operator
KContainer&
KContainer::operator=(const KContainer& ref) {
  //Lattices should be equal. 
  if(&Lattice != &ref.Lattice){
    LOGMSG("ERROR: tried to copy KContainer with different lattices");
    OHMMS::Controller->abort();
  }
  //Now, if kcutoffs are the same then we can be sure that the lists are identical.
  //otherwise the STL containers must have contents copied.
  if(this!=&ref && kcutoff!=ref.kcutoff){
    //All components have a valid '=' defined
    kcutoff = ref.kcutoff;
    kcut2 = ref.kcut2;
    mmax = ref.mmax;
    kpts = ref.kpts;
    kpts_cart = ref.kpts_cart;
    minusk = ref.minusk;
    numk = ref.numk;
  }

  return *this;
}

//Public Methods
// UpdateKLists - call for new k or when lattice changed.
void
KContainer::UpdateKLists(ParticleLayout_t& ref, RealType kc) {
  kcutoff = kc;
  kcut2 = kc*kc;
  Lattice = ref;

  if(kcutoff <= 0.0){
    LOGMSG("KContainer initialised with cutoff " << kcutoff);
    OHMMS::Controller->abort();
  }

  FindApproxMMax();
  BuildKLists();
}

// UpdateKLists - call for new k or when lattice changed.
void
KContainer::UpdateKLists(RealType kc) {
  kcutoff = kc;
  kcut2 = kc*kc;

  if(kcutoff <= 0.0){
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
  //Does not require an orthorhombic cell. 
  /* Old method.
  //2pi is not included in Lattice.b
  Matrix<RealType> mmat;
  mmat.resize(3,3);
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++){
      mmat[i][j] = 0.0;
      for(int k=0;k<3;k++)
	mmat[i][j] = mmat[i][j] + 4.0*M_PI*M_PI*Lattice.b(k)[i]*Lattice.b(j)[k];
    }

  TinyVector<RealType,3> x,temp;
  RealType tempr;
  for(int idim=0;idim<3;idim++){
    int i = ((idim)%3);
    int j = ((idim+1)%3);
    int k = ((idim+2)%3);
    
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
    mmax[idim] = static_cast<int>(sqrt(4.0*kcut2/tempr)) + 1;  
  }
  */
  // see rmm, Electronic Structure, p. 85 for details
  for (int i = 0; i < 3; i++) 
    mmax[i] = static_cast<int>(floor(sqrt(dot(Lattice.a(i),Lattice.a(i))) * kcutoff / (2 * M_PI))) + 1;
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

  // Create a map from the hash value for each k vector to the index
  std::map<int, int> hashToIndex;
  for (int ki=0; ki<numk; ki++) {
    hashToIndex[GetHashOfVec(kpts[ki], numk)] = ki;
  }

  // Use the map to find the index of -k from the index of k
  for(int ki=0; ki<numk; ki++) {
    minusk[ki] = hashToIndex[ GetHashOfVec(-1 * kpts[ki], numk) ];
  }
}
