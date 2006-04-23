// -*- C++ -*-
#ifndef QMCPLUSPLUS_PLANEWAVEBASIS_H
#define QMCPLUSPLUS_PLANEWAVEBASIS_H

#include "Configuration.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include <complex>

namespace qmcplusplus {

  class PlaneWaveBasis: public QMCTraits {
  private:
    typedef ParticleSet::ParticleLayout_t ParticleLayout_t;    
    ///the number of particles
    int NumPtcls;

    //The PlaneWave data
    double ecut;
    TinyVector<double,3> twist;
    vector<double> modkplusg;
    vector<TinyVector<int,3> > gvecs; //Reduced coordinates
    vector<TinyVector<double,3> > kplusgvecs_cart; //Cartesian.
    //Need to store the maximum translation in each dimension to use recursive PW generation.
    TinyVector<int,3> maxg;
    int maxmaxg;

    //Storage for basis functions evaluated at given particle positions.
#if defined(QMC_COMPLEX)
    ///matrix to store values \f$ Y[i,j] = \exp(iG_jr_i) \f$
    Matrix<complex<double> > Y;
    ///matrix to store gradients \f$ dY[i,j] = {\bf \nabla}_i Y[i,j] \f$
    Matrix<TinyVector<complex<double>,3> >  dY;
    ///matrix to store laplacians \f$ d2Y[i,j] = \nabla^2_i Y[i,j] \f$
    Matrix<complex<double> > d2Y;
#else
    //Real wavefunctions here. Now the basis states are cos(Gr) or sin(Gr), not exp(iGr) 
    //We need a way of switching between them for G -> -G, otherwise the
    //determinant will have multiple rows that are equal (to within a constant) of
    //others, giving a zero determinant. For this, we build a vector (negative) which 
    //stores whether a vector is "+" or "-" (with some criterion, to be defined). We
    //the switch from cos() to sin() based on the value of this input.
    vector<int> negative;
    ///matrix to store values \f$ Y[i,j] = \cos(Gr) \mathrm{or} \sin(Gr) \f$
    Matrix<double> Y;
    Matrix<TinyVector<double,3> >  dY;
    Matrix<double> d2Y;
#endif
  public:
    /* inputmap is used for a memory efficient way of 
    ** importing the basis-set and coefficients when the desired energy cutoff may be
    ** lower than that represented by all data in the wavefunction input file.
    ** The steps taken are:
    **  * Read all basis data.
    **  * Create map. inputmap[i] = j; j is correct PW index, i is input coef index.
    **    For basis elements outside cutoff, inputmap[i] = gvecs.size(); 
    **  * Coefficients are in same order as PWs in inputfile => simply file into
    **    storage matrix using the map as the input. All excess coefficients are 
    **    put into [gvecs.size()] and not used. i.e. coefs need to be allocated 1 higher.
    ** Such an approach is not needed for Gamma-point only calculations because the
    ** basis is spherically ordered. However, when a twist-angle is used, the "sphere"
    ** of allowed planewaves is shifted.
    */
    vector<int> inputmap;

    ///total number of basis functions
    int NumPlaneWaves;

    ///constructor
    PlaneWaveBasis(TinyVector<double,3> twistangle, int nptcl): NumPtcls(nptcl), twist(twistangle) {
      NumPlaneWaves=0;
      //Initialize containers for maximum plane-wave translations (reduced coordinates)
      for(int i=0; i<3; i++)
	maxg[i] = 0;
    }

    /**
       @param nptcl number of particles
       @brief resize the containers for data 
       for nptcl particles and TotalBasis basis functions
    */
    inline void resize(int nptcl) {
      NumPtcls = nptcl;
      Y.resize(nptcl,NumPlaneWaves);
      dY.resize(nptcl,NumPlaneWaves);
      d2Y.resize(nptcl,NumPlaneWaves);
    }

    ///Read basisset from hdf5 file. Apply ecut. Resize internal storage.
    inline void
    readbasis(hid_t h5basisgroup,double ecutoff,ParticleLayout_t &Lattice) {
      //First, read the total number of planewaves:
      int idata;
      HDFAttribIO<int> hdfint(idata);
      hdfint.read(h5basisgroup,"num_planewaves");
      NumPlaneWaves = idata;

      //Resize the storage for the G-Vectors in reduced coordinates:
      gvecs.resize(NumPlaneWaves);

      //Read ALL of the planewaves (even those > ecut):
      //      HDFAttribIO<int*>hdfvtv(&gvecs[0][0],NumPlaneWaves*3);
      HDFAttribIO<std::vector<TinyVector<int,3> > >hdfvtv(gvecs);
      hdfvtv.read(h5basisgroup,"planewaves");

      //Now remove elements outside Ecut. At the same time, fill k+G and |k+G| lists.
      //Also keep track of the original index ordering (using indexmap[]) so that
      //orbital coefficients can be ordered and trimmed for ecut in the same way.
      ecut = ecutoff;
      trimforecut(Lattice);

      //Resize storage for the basis-elemtents evaluated at the 
      //particle coordinates, since basis has been trimmed.
      resize(NumPtcls);

      //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
      for(int ig=0; ig<NumPlaneWaves; ig++)
        for(int i=0; i<3; i++)
           if(abs(gvecs[ig][i]) > maxg[i]) 
	 maxg[i] = abs(gvecs[ig][i]);
      maxmaxg = max(maxg[0],max(maxg[1],maxg[2]));

      LOGMSG("\n\tBasisset energy cutoff = " << ecut);
      LOGMSG("\tNumber of planewaves = " << NumPlaneWaves<<"\n");

      //Check that we actually kept some elements within ecut.
      if(NumPlaneWaves < 1){
	LOGMSG("No planewaves exist within ecut (="<<ecut<<")");
	OHMMS::Controller->abort();
      }
    }


    ///Remove basis elements if kinetic energy > ecut.
    //Keep and indexmap so we know how to match coefficients on read.
    void trimforecut(ParticleLayout_t &Lattice) { 
      TinyVector<double,3> tempvec;
      double mod2, mod;
      double kcutoff = std::sqrt(2.0*ecut);

      //resize inputmap
      NumPlaneWaves = gvecs.size();
      inputmap.resize(NumPlaneWaves);

      //ig is the loop index to access the member of gvecs for testing.
      //newig is the index showing where ig exists in the new (truncated) basis.
      //oldig is the index showing where ig came from...differs from ig after gvecs 
      // has at least one element truncated.
      for(int ig=0, newig=0, oldig=0; ig<NumPlaneWaves; ig++,oldig++) {
	//Check size of this g-vector
	tempvec = Lattice.k_cart(gvecs[ig]+twist);
	mod2 = dot(tempvec,tempvec);
	mod = std::sqrt(mod2);
	if(mod<=kcutoff){
	  //Keep this element
	  kplusgvecs_cart.push_back(tempvec);
	  modkplusg.push_back(mod);
	  //Remember which position in the HDF5 file this came from...for coefficients
	  inputmap[oldig] = newig;
	  newig++;
#if !defined(QMC_COMPLEX)
	  //Build the negative vector. See comment at declaration (above) for details.
	  if(gvecs[ig][0] < 0)
	    negative.push_back(0);
	  else if(gvecs[ig][0] > 0)
	    negative.push_back(1);
	  else { //gx == 0, test gy
	    if(gvecs[ig][1] < 0)
	      negative.push_back(0);
	    else if(gvecs[ig][1] > 0)
	      negative.push_back(1);
	    else { //gx == gy == 0; test gz. If gz==0 also, take negative=1 (arbitrary)
	      if(gvecs[ig][2] < 0)
		negative.push_back(0);
	      else
		negative.push_back(1);
	    }
	  }
#endif
	} else {
	  //Remove this element. Remember to set ig back by one element so 
	  //removal doesn't lead to a skipping
	  inputmap[oldig] = -1; //Temporary value...need to know final NumPlaneWaves.
	  gvecs.erase(gvecs.begin()+ig,gvecs.begin()+ig+1);
	  ig--; NumPlaneWaves--;
	}
      }

      //Finalize the basis. Fix temporary values of inputmap.
      for(int ig=0; ig<inputmap.size(); ig++) 
        if(inputmap[ig] == -1)
	  inputmap[ig] = NumPlaneWaves; //For dumping coefficients of PWs>ecut
    }

    ///Evaluate all planewaves for current particle coordinates
    inline void 
    evaluate(const ParticleSet& P) {
      if(P.getTotalNum() != NumPtcls) {
	cout << "PlaneWaveBasis Error: Storage not allocated correctly" << endl;
	OHMMS::Controller->abort();
      }
      //Evaluate the plane-waves at current particle coordinates using a fast
      //recursion algorithm. Only Y is evaluated in this routine.
      //These can be "dotted" with coefficients later to complete
      //orbital evaluations.
      TinyVector<double,3> G111; // Cartesian of twist+G for 1,1,1 (reduced coordinates)
      for(int idim=0; idim<3; idim++)
	G111[idim] = 1.0 + twist[idim]; //Reduced
      G111 = P.Lattice.k_cart(G111); //Cartesian

      //Allocate the 'C' temporary arrays. Fill with information to decouple dimensions.
      Matrix<complex<double> > C;
      complex<double> pw;
      C.resize(3,2*maxmaxg+1);
      for(int ip=0; ip<NumPtcls; ip++){
	//Precompute a small number of complex factors (PWs along b1,b2.b3 lines)
	//using a fast recursion algorithm
	for(int idim=0; idim<3; idim++){
	  complex<double> Ctemp;
	  //start the recursion with the 111 vector.
	  double phi = (P.R[ip])[idim] * G111[idim];
	  Ctemp = complex<double>(std::cos(phi), std::sin(phi));
	  C(idim,maxg[idim]) = 1.0; // G=0. Index is shifted: [0..2*max]. Zero at 'max'.
	  //Recursively generate all Cs for each dimension independently.
	  for(int n=1; n<=maxg[idim]; n++){
	    C(idim,maxg[idim]+n) = Ctemp*C(idim,maxg[idim]+n-1);
	    C(idim,maxg[idim]-n) = conj(C(idim,maxg[idim]+n));
	  }
	}


	//Evaluate the planewaves and derivatives.
	for(int ig=0; ig<NumPlaneWaves; ig++) {
	  pw = 1.0;
	  for(int idim=0; idim<3; idim++)
	    pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
	  Y(ip,ig) = pw;
#else
	  Y(ip,ig) = negative[ig]*pw.real() + (1-negative[ig])*pw.imag();
#endif
	}
      }
      /* Old, slow way
      double gr;
      for(int ig=0; ig<NumPlaneWaves; ig++) {
	for(int ip=0; ip<NumPtcls; ip++) {
	  //G.r...
	  gr = dot(kplusgvecs_cart[ig],P.R[ip]);
#if defined(QMC_COMPLEX)
	  Y(ip,ig) = complex<double>(std::cos(gr),std::sin(gr));
#else
	  Y(ip,ig) = negative[ig]*std::cos(gr) + (1-negative[ig])*std::sin(gr);
#endif
	}
      }
      */
    }

    ///Evaluate all planewaves and derivatives for current coordinates
    inline void 
    evaluateAll(const ParticleSet& P) {
      if(P.getTotalNum() != NumPtcls) {
	LOGMSG("PlaneWaveBasis Error: Storage not allocated correctly");
	OHMMS::Controller->abort();
      }

      //Evaluate the plane-waves at current particle coordinates using a fast
      //recursion algorithm. Order of Y,dY and d2Y is kept correct.
      //These can be "dotted" with coefficients later to complete
      //orbital evaluations.
      TinyVector<double,3> G111; // Cartesian of twist+G for 1,1,1 (reduced coordinates)
      for(int idim=0; idim<3; idim++)
	G111[idim] = 1.0 + twist[idim]; //Reduced
      G111 = P.Lattice.k_cart(G111); //Cartesian

      //Allocate the 'C' temporary arrays. Fill with information to decouple dimensions.
      Matrix<complex<double> > C;
      complex<double> pw;
      C.resize(3,2*maxmaxg+1);
      for(int ip=0; ip<NumPtcls; ip++){
	//Precompute a small number of complex factors (PWs along b1,b2.b3 lines)
	//using a fast recursion algorithm
	for(int idim=0; idim<3; idim++){
	  complex<double> Ctemp;
	  //start the recursion with the 111 vector.
	  double phi = (P.R[ip])[idim] * G111[idim];
	  Ctemp = complex<double>(std::cos(phi), std::sin(phi));
	  C(idim,maxg[idim]) = 1.0; // G=0. Index is shifted: [0..2*max]. Zero at 'max'.
	  //Recursively generate all Cs for each dimension independently.
	  for(int n=1; n<=maxg[idim]; n++){
	    C(idim,maxg[idim]+n) = Ctemp*C(idim,maxg[idim]+n-1);
	    C(idim,maxg[idim]-n) = conj(C(idim,maxg[idim]+n));
	  }
	}


	//Evaluate the planewaves and derivatives.
	for(int ig=0; ig<NumPlaneWaves; ig++) {
	  pw = 1.0;
	  for(int idim=0; idim<3; idim++)
	    pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
	  Y(ip,ig) = pw;
	  for(int idim=0; idim<3; idim++)
	    dY(ip,ig)[idim] = Y(ip,ig)*complex<double>(0.0,kplusgvecs_cart[ig][idim]);
#else
	  Y(ip,ig) = negative[ig]*pw.real() + (1-negative[ig])*pw.imag();
	  for(int idim=0; idim<3; idim++) {
	    dY(ip,ig)[idim] = (1-negative[ig])*pw.real();
	    dY(ip,ig)[idim]-= negative[ig]*pw.imag();
	    dY(ip,ig)[idim]*= kplusgvecs_cart[ig][idim];
	  }

#endif
	  d2Y(ip,ig) = -modkplusg[ig]*modkplusg[ig]*Y(ip,ig);
	}
      }

      /*
      //DEBUGGING CODE:
      //Test the accuracy of the recursively-generated PWs vs. the explicit (slow) evaluation
      double gr;
#if defined(QMC_COMPLEX)
      complex<double> testY,testd2Y;
      TinyVector<complex<double>,3> testdY;
#else
      double testY,testd2Y;
      TinyVector<double,3> testdY;
#endif
      for(int ig=0; ig<NumPlaneWaves; ig++)
	for(int ip=0; ip<NumPtcls; ip++) {
	  gr = dot(kplusgvecs_cart[ig],P.R[ip]);
#if defined(QMC_COMPLEX)
	  testY = complex<double>(std::cos(gr),std::sin(gr));
	  for(int idim=0; idim<3; idim++)
	    testdY[idim] = Y(ip,ig) * complex<double>(0.0,kplusgvecs_cart[ig][idim]);
#else
	  testY = negative[ig]*std::cos(gr) + (1-negative[ig])*std::sin(gr);
	  for(int idim=0; idim<3; idim++)
	    testdY[idim] = kplusgvecs_cart[ig][idim]*((1-negative[ig])*std::cos(gr)-negative[ig]*std::sin(gr));
#endif
	  testd2Y = -modkplusg[ig]*modkplusg[ig]*Y(ip,ig);
	  
	  testY -= Y(ip,ig);
	  for(int idim=0; idim<3; idim++)
	    testdY[idim] -= dY(ip,ig)[idim];
	  //	  testd2Y -= d2Y(ip,ig);

	  if(abs(testY) > 1.e-06) {
	    cout << "Planewave evaluation error:" << endl;
	    cout << "ip = " << ip << endl;
	    cout << "ig = " << ig << endl;
	    cout << "G.r = " << gr << endl;
	    cout << "cos(Gr) = " << std::cos(gr) << endl;
	    cout << "sin(Gr) = " << std::sin(gr) << endl;
	    cout << "Y(ip,ig) = " << Y(ip,ig) << endl;
	    abort();
	  }
	  for(int idim=0; idim<3; idim++)
	    if(abs(testdY[idim]) > 1.e-06) {
	      cout << "Planewave evaluation error:" << endl;
	      cout << "ip = " << ip << endl;
	      cout << "ig = " << ig << endl;
	      cout << "idim = " << idim << endl;
	      cout << "G.r = " << gr << endl;
	      cout << "cos(Gr) = " << std::cos(gr) << endl;
	      cout << "sin(Gr) = " << std::sin(gr) << endl;
	      cout << "dY(ip,ig) = " << dY(ip,ig)[idim] << endl;
	      abort();
	    }
	  if(abs(testd2Y-d2Y(ip,ig)) > 1.e-06) { 
	    cout << "Planewave evaluation error:" << endl;
	    cout << "ip = " << ip << endl;
	    cout << "ig = " << ig << endl;
	    cout << "G.r = " << gr << endl;
	    cout << "cos(Gr) = " << std::cos(gr) << endl;
	    cout << "sin(Gr) = " << std::sin(gr) << endl;
	    cout << "d2Y(ip,ig) fast = " << d2Y(ip,ig) << endl;
	    cout << "d2Y(ip,ig) direct = " << testd2Y << endl;
	    abort();
	  }

	}
      exit(0);
      //END DEBUGGING CODE.
      */
    }

#if defined(QMC_COMPLEX)
    ///row i of matrix Y
    inline const complex<double>* restrict y(int i){ return &Y(i,0);}
    ///row i of vector matrix dY
    inline const TinyVector<complex<double>,3>* restrict dy(int i){ return &dY(i,0);}
    ///row i of matrix d2Y
    inline const complex<double>* restrict d2y(int i){ return &d2Y(i,0);}
#else
    ///row i of matrix Y
    inline const double* restrict y(int i){ return &Y(i,0);}
    ///row i of vector matrix dY
    inline const TinyVector<double,3>* restrict dy(int i){ return &dY(i,0);}
    ///row i of matrix d2Y
    inline const double* restrict d2y(int i){ return &d2Y(i,0);}
#endif
  };

}
#endif
