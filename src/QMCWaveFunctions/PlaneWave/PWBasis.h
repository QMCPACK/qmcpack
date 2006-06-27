// -*- C++ -*-
#ifndef QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H
#define QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include <complex>

namespace qmcplusplus {

  /** Plane-wave basis set
   *
   * Rewrite of PlaneWaveBasis to utilize blas II or III
   */
  class PWBasis: public QMCTraits {
  private:
    typedef ParticleSet::ParticleLayout_t ParticleLayout_t;    
    //The PlaneWave data - keep all of these strictly private to prevent inconsistencies.
    RealType ecut;
    TinyVector<RealType,3> twist,twist_cart; //Twist angle in reduced and Cartesian.
    vector<RealType> minusModKplusG2;
    vector<TinyVector<int,3> > gvecs; //Reduced coordinates
    vector<PosType> kplusgvecs_cart; //Cartesian.
    //Need to store the maximum translation in each dimension to use recursive PW generation.
    TinyVector<int,3> maxg;
    int maxmaxg;
    Matrix<ValueType> C;
#if !defined(QMC_COMPLEX)
    //Real wavefunctions here. Now the basis states are cos(Gr) or sin(Gr), not exp(iGr) 
    //We need a way of switching between them for G -> -G, otherwise the
    //determinant will have multiple rows that are equal (to within a constant factor) 
    //of others, giving a zero determinant. For this, we build a vector (negative) which 
    //stores whether a vector is "+" or "-" (with some criterion, to be defined). We
    //the switch from cos() to sin() based on the value of this input.
    vector<int> negative;
#endif
  public:
    //enumeration for the value, laplacian, gradients and size
    enum {PW_VALUE, PW_LAP, PW_GRADX, PW_GRADY, PW_GRADZ, PW_MAXINDEX};

    Matrix<ValueType> Z;
    Vector<ValueType> Zv;
    /* inputmap is used for a memory efficient way of 
     *
     * importing the basis-set and coefficients when the desired energy cutoff may be
     * lower than that represented by all data in the wavefunction input file.
     * The steps taken are:
     *  - Read all basis data.
     *  - Create map. inputmap[i] = j; j is correct PW index, i is input coef index.
     *    For basis elements outside cutoff, inputmap[i] = gvecs.size(); 
     *  - Coefficients are in same order as PWs in inputfile => simply file into
     *    storage matrix using the map as the input. All excess coefficients are 
     *    put into [gvecs.size()] and not used. i.e. coefs need to be allocated 1 higher.
     * Such an approach is not needed for Gamma-point only calculations because the
     * basis is spherically ordered. However, when a twist-angle is used, the "sphere"
     * of allowed planewaves is shifted.
     */
    vector<int> inputmap;

    ///total number of basis functions
    int NumPlaneWaves;

    ///constructor
    PWBasis(const TinyVector<RealType,3>& twistangle): twist(twistangle) {
      NumPlaneWaves=0;
      maxg=0.0;
    }

    ///Read basisset from hdf5 file. Apply ecut. Resize internal storage.
    inline void
    readbasis(hid_t h5basisgroup,RealType ecutoff,int& nh5gvecs, 
              ParticleLayout_t &Lattice) {
      //First, read the total number of planewaves:
      int idata;
      HDFAttribIO<int> hdfint(idata);
      hdfint.read(h5basisgroup,"num_planewaves");
      NumPlaneWaves = idata;
      nh5gvecs = idata;

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

      //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
      for(int ig=0; ig<NumPlaneWaves; ig++)
        for(int i=0; i<3; i++)
           if(abs(gvecs[ig][i]) > maxg[i]) 
             maxg[i] = abs(gvecs[ig][i]);
      maxmaxg = max(maxg[0],max(maxg[1],maxg[2]));

      //changes the order????
      C.resize(3,2*maxmaxg+1);

      LOGMSG("\n\tBasisset energy cutoff = " << ecut);
      LOGMSG("\tNumber of planewaves = " << NumPlaneWaves<<"\n");

      //Check that we actually kept some elements within ecut.
      if(NumPlaneWaves < 1){
        LOGMSG("No planewaves exist within ecut (="<<ecut<<")");
        OHMMS::Controller->abort();
      }

      Z.resize(NumPlaneWaves,2+DIM);
      Zv.resize(NumPlaneWaves);
    }


    ///Remove basis elements if kinetic energy > ecut.
    //Keep and indexmap so we know how to match coefficients on read.
    void trimforecut(ParticleLayout_t &Lattice) { 
      TinyVector<RealType,3> tempvec;
      RealType mod2, mod;
      RealType kcutoff = std::sqrt(2.0*ecut);

      //resize inputmap
      NumPlaneWaves = gvecs.size();
      inputmap.resize(NumPlaneWaves);

      //Convert the twist angle to Cartesian coordinates.
      twist_cart = Lattice.k_cart(twist);

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
          minusModKplusG2.push_back(-mod2);
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

    void BuildRecursionCoefs(const ParticleSet &P, int iat) {
      //Fill the recursion coefficients matrix.
      TinyVector<RealType,3> G111; // Cartesian of twist+G for 1,1,1 (reduced coordinates)
      for(int idim=0; idim<3; idim++)
        G111[idim] = 1.0 + twist[idim]; //Reduced
      G111 = P.Lattice.k_cart(G111); //Cartesian

      //Precompute a small number of complex factors (PWs along b1,b2,b3 lines)
      //using a fast recursion algorithm
      for(int idim=0; idim<3; idim++){
        //start the recursion with the 111 vector.
        RealType phi = (P.R[iat])[idim] * G111[idim];
        register ComplexType Ctemp(std::cos(phi), std::sin(phi));
        register int ng=maxg[idim];
        ComplexType* restrict cp_ptr=C[idim]+ng;
        ComplexType* restrict cn_ptr=C[idim]+ng-1;
        *cp_ptr=1.0;
        //add INTEL vectorization
        for(int n=1; n<=ng; n++,cn_ptr--){
          ComplexType t(Ctemp*(*cp_ptr++));
          *cp_ptr = t;
          *cn_ptr = conj(t);
        }

        //C(idim,maxg[idim]) = 1.0; // G=0. Index is shifted: [0..2*max]. Zero at 'max'.
        ////Recursively generate all Cs for each dimension independently.
        //for(int n=1; n<=maxg[idim]; n++){
        //  C(idim,maxg[idim]+n) = Ctemp*C(idim,maxg[idim]+n-1);
        //  C(idim,maxg[idim]-n) = conj(C(idim,maxg[idim]+n));
        //}
      }
    }

    inline void 
    evaluate(const ParticleSet& P, int iat) {
      BuildRecursionCoefs(P,iat);
      RealType twistdotr = dot(twist_cart,P.R[iat]);
      ComplexType pw0(std::cos(twistdotr),std::sin(twistdotr));
      //Evaluate the planewaves for particle iat.
      for(int ig=0; ig<NumPlaneWaves; ig++) {
        //PW is initialized as exp(i*twist.r) so that the final basis evaluations
        //are for (twist+G).r
        ComplexType pw(pw0); //std::cos(twistdotr),std::sin(twistdotr));
        for(int idim=0; idim<3; idim++)
          pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
        Zv[ig]=pw; 
#else
        Zv[ig]= negative[ig]*pw.real() + (1-negative[ig])*pw.imag();
#endif
      }
    }

    /** Evaluate all planewaves and derivatives for the iat-th particle
     *
     * The basis functions are evaluated for particles iat: first <= iat < last
     * Evaluate the plane-waves at current particle coordinates using a fast
     * recursion algorithm. Order of Y,dY and d2Y is kept correct.
     * These can be "dotted" with coefficients later to complete orbital evaluations.
     */
    inline void 
    evaluateAll(const ParticleSet& P, int iat) {
      BuildRecursionCoefs(P,iat);
      RealType twistdotr = dot(twist_cart,P.R[iat]);
      ComplexType pw0(std::cos(twistdotr),std::sin(twistdotr));
      //Evaluate the planewaves and derivatives.
      ComplexType* restrict zptr=Z.data();
      for(int ig=0; ig<NumPlaneWaves; ig++,zptr+=5) {
        //PW is initialized as exp(i*twist.r) so that the final basis evaluations
        //are for (twist+G).r
        ComplexType pw(pw0);
        // THE INDEX ORDER OF C DOESN'T LOOK TOO GOOD: this could be fixed
        for(int idim=0; idim<3; idim++)
          pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
        zptr[0]= pw;
        zptr[1]= minusModKplusG2[ig]*pw;
        zptr[2]= pw*ComplexType(0.0,kplusgvecs_cart[ig][0]);
        zptr[3]= pw*ComplexType(0.0,kplusgvecs_cart[ig][1]);
        zptr[4]= pw*ComplexType(0.0,kplusgvecs_cart[ig][2]);
#else
        ERROR_MSG("DO NOT USE THIS UNTIL TESTED")
#endif
      }
    }
  };
}
#endif
