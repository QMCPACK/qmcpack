//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_PLANEWAVEBASIS_H
#define QMCPLUSPLUS_PLANEWAVEBASIS_H

#include "Configuration.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include <complex>

namespace qmcplusplus
{

class PlaneWaveBasis: public QMCTraits
{
private:
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  ///the number of particles
  int NumPtcls;

  //The PlaneWave data - keep all of these strictly private to prevent inconsistencies.
  RealType ecut;
  TinyVector<RealType,3> twist,twist_cart; //Twist angle in reduced and Cartesian.
  std::vector<RealType> modkplusg;
  std::vector<TinyVector<int,3> > gvecs; //Reduced coordinates
  std::vector<TinyVector<RealType,3> > kplusgvecs_cart; //Cartesian.
  //Need to store the maximum translation in each dimension to use recursive PW generation.
  TinyVector<int,3> maxg;
  int maxmaxg;

  //Storage for basis functions evaluated at given particle positions.
  //valuetype handles changing from real to complex in Configuration.h
  ///matrix to store values \f$ Y[i,j] = \exp(iG_jr_i) \f$
  Matrix<ValueType> Y;
  ///matrix to store gradients \f$ dY[i,j] = {\bf \nabla}_i Y[i,j] \f$
  Matrix<TinyVector<ValueType,3> > dY;
  ///matrix to store laplacians \f$ d2Y[i,j] = \nabla^2_i Y[i,j] \f$
  Matrix<ValueType> d2Y;
#if !defined(QMC_COMPLEX)
  //Real wavefunctions here. Now the basis states are cos(Gr) or sin(Gr), not exp(iGr)
  //We need a way of switching between them for G -> -G, otherwise the
  //determinant will have multiple rows that are equal (to within a constant factor)
  //of others, giving a zero determinant. For this, we build a vector (negative) which
  //stores whether a vector is "+" or "-" (with some criterion, to be defined). We
  //the switch from cos() to sin() based on the value of this input.
  std::vector<int> negative;
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
  std::vector<int> inputmap;

  ///total number of basis functions
  int NumPlaneWaves;

  ///constructor
  PlaneWaveBasis(TinyVector<RealType,3> twistangle, int nptcl): NumPtcls(nptcl), twist(twistangle)
  {
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
  inline void resize(int nptcl)
  {
    NumPtcls = nptcl;
    Y.resize(nptcl,NumPlaneWaves);
    dY.resize(nptcl,NumPlaneWaves);
    d2Y.resize(nptcl,NumPlaneWaves);
  }

  ///Read basisset from hdf5 file. Apply ecut. Resize internal storage.
  inline void
  readbasis(hid_t h5basisgroup,RealType ecutoff,int& nh5gvecs,
            ParticleLayout_t &Lattice)
  {
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
    //Resize storage for the basis-elemtents evaluated at the
    //particle coordinates, since basis has been trimmed.
    resize(NumPtcls);
    //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
    for(int ig=0; ig<NumPlaneWaves; ig++)
      for(int i=0; i<3; i++)
        if(std::abs(gvecs[ig][i]) > maxg[i])
          maxg[i] = std::abs(gvecs[ig][i]);
    maxmaxg = std::max(maxg[0],max(maxg[1],maxg[2]));
    LOGMSG("\n\tBasisset energy cutoff = " << ecut);
    LOGMSG("\tNumber of planewaves = " << NumPlaneWaves<<"\n");
    //Check that we actually kept some elements within ecut.
    if(NumPlaneWaves < 1)
    {
      LOGMSG("No planewaves exist within ecut (="<<ecut<<")");
      OHMMS::Controller->abort();
    }
  }


  ///Remove basis elements if kinetic energy > ecut.
  //Keep and indexmap so we know how to match coefficients on read.
  void trimforecut(ParticleLayout_t &Lattice)
  {
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
    for(int ig=0, newig=0, oldig=0; ig<NumPlaneWaves; ig++,oldig++)
    {
      //Check size of this g-vector
      tempvec = Lattice.k_cart(gvecs[ig]+twist);
      mod2 = dot(tempvec,tempvec);
      mod = std::sqrt(mod2);
      if(mod<=kcutoff)
      {
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
        else
          if(gvecs[ig][0] > 0)
            negative.push_back(1);
          else
            //gx == 0, test gy
          {
            if(gvecs[ig][1] < 0)
              negative.push_back(0);
            else
              if(gvecs[ig][1] > 0)
                negative.push_back(1);
              else
                //gx == gy == 0; test gz. If gz==0 also, take negative=1 (arbitrary)
              {
                if(gvecs[ig][2] < 0)
                  negative.push_back(0);
                else
                  negative.push_back(1);
              }
          }
#endif
      }
      else
      {
        //Remove this element. Remember to set ig back by one element so
        //removal doesn't lead to a skipping
        inputmap[oldig] = -1; //Temporary value...need to know final NumPlaneWaves.
        gvecs.erase(gvecs.begin()+ig,gvecs.begin()+ig+1);
        ig--;
        NumPlaneWaves--;
      }
    }
    //Finalize the basis. Fix temporary values of inputmap.
    for(int ig=0; ig<inputmap.size(); ig++)
      if(inputmap[ig] == -1)
        inputmap[ig] = NumPlaneWaves; //For dumping coefficients of PWs>ecut
  }

  void BuildRecursionCoefs(Matrix<std::complex<RealType> > &C, const ParticleSet &P,
                           int iat)
  {
    ////Fill the recursion coefficients matrix.
    //for(int idim=0; idim<3; idim++)
    //  G111[idim] = 1.0 + twist[idim]; //Reduced
    //G111 = P.Lattice.k_cart(G111); //Cartesian
    TinyVector<RealType,3> G111(1.0,1.0,1.0);
    G111 = P.Lattice.k_cart(G111);
    //Precompute a small number of complex factors (PWs along b1,b2.b3 lines)
    //using a fast recursion algorithm
    for(int idim=0; idim<3; idim++)
    {
      std::complex<RealType> Ctemp;
      //start the recursion with the 111 vector.
      RealType phi = (P.R[iat])[idim] * G111[idim];
      Ctemp = std::complex<RealType>(std::cos(phi), std::sin(phi));
      C(idim,maxg[idim]) = 1.0; // G=0. Index is shifted: [0..2*max]. Zero at 'max'.
      //Recursively generate all Cs for each dimension independently.
      for(int n=1; n<=maxg[idim]; n++)
      {
        C(idim,maxg[idim]+n) = Ctemp*C(idim,maxg[idim]+n-1);
        C(idim,maxg[idim]-n) = conj(C(idim,maxg[idim]+n));
      }
    }
  }



  ///Evaluate all planewaves for current particle coordinates
  //This function does not evaluate first or second derivatives -> use evaluateAll.
  //The basis functions are evaluated for particles iat: first <= iat < last
  inline void
  evaluate(const ParticleSet& P, int first, int last)
  {
    if(P.getTotalNum() != NumPtcls)
    {
      std::cout << "PlaneWaveBasis Error: Storage not allocated correctly" << std::endl;
      OHMMS::Controller->abort();
    }
    //Evaluate the plane-waves at current particle coordinates using a fast
    //recursion algorithm. Only Y is evaluated in this routine.
    //These can be "dotted" with coefficients later to complete
    //orbital evaluations.
    //Allocate the 'C' temporary arrays. Fill with information to decouple dimensions.
    Matrix<std::complex<RealType> > C;
    RealType twistdotr;
    std::complex<RealType> pw;
    C.resize(3,2*maxmaxg+1);
    for(int iat=first; iat<last; iat++)
    {
      BuildRecursionCoefs(C,P,iat);
      twistdotr = dot(twist_cart,P.R[iat]);
      //Evaluate the planewaves for particle iat.
      for(int ig=0; ig<NumPlaneWaves; ig++)
      {
        //PW is initialized as exp(i*twist.r) so that the final basis evaluations
        //are for (twist+G).r
        pw = std::complex<RealType>(std::cos(twistdotr),std::sin(twistdotr));
        for(int idim=0; idim<3; idim++)
          pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
        Y(iat,ig) = pw;
#else
        Y(iat,ig) = negative[ig]*pw.real() + (1-negative[ig])*pw.imag();
#endif
      }
    }
  }



  ///Evaluate all planewaves and derivatives for current coordinates
  //The basis functions are evaluated for particles iat: first <= iat < last
  inline void
  evaluateAll(const ParticleSet& P, int first, int last)
  {
    if(P.getTotalNum() != NumPtcls)
    {
      LOGMSG("PlaneWaveBasis Error: Storage not allocated correctly");
      OHMMS::Controller->abort();
    }
    //Evaluate the plane-waves at current particle coordinates using a fast
    //recursion algorithm. Order of Y,dY and d2Y is kept correct.
    //These can be "dotted" with coefficients later to complete
    //orbital evaluations.
    //Allocate the 'C' temporary arrays. Fill with information to decouple dimensions.
    Matrix<std::complex<RealType> > C;
    RealType twistdotr;
    std::complex<RealType> pw;
    C.resize(3,2*maxmaxg+1);
    for(int iat=first; iat<last; iat++)
    {
      BuildRecursionCoefs(C,P,iat);
      twistdotr = dot(twist_cart,P.R[iat]);
      //Evaluate the planewaves and derivatives.
      for(int ig=0; ig<NumPlaneWaves; ig++)
      {
        //PW is initialized as exp(i*twist.r) so that the final basis evaluations
        //are for (twist+G).r
        pw = std::complex<RealType>(std::cos(twistdotr),std::sin(twistdotr));
        for(int idim=0; idim<3; idim++)
          pw *= C(idim,gvecs[ig][idim]+maxg[idim]);
#if defined(QMC_COMPLEX)
        Y(iat,ig) = pw;
        for(int idim=0; idim<3; idim++)
          dY(iat,ig)[idim] = Y(iat,ig)*complex<RealType>(0.0,kplusgvecs_cart[ig][idim]);
#else
        Y(iat,ig) = negative[ig]*pw.real() + (1-negative[ig])*pw.imag();
        for(int idim=0; idim<3; idim++)
        {
          dY(iat,ig)[idim] = (1-negative[ig])*pw.real();
          dY(iat,ig)[idim]-= negative[ig]*pw.imag();
          dY(iat,ig)[idim]*= kplusgvecs_cart[ig][idim];
        }
#endif
        d2Y(iat,ig) = -modkplusg[ig]*modkplusg[ig]*Y(iat,ig);
      }
    }
  }

  ///row i of matrix Y
  inline const ValueType* __restrict__ y(int i)
  {
    return &Y(i,0);
  }
  ///row i of vector matrix dY
  inline const TinyVector<ValueType,3>* __restrict__ dy(int i)
  {
    return &dY(i,0);
  }
  ///row i of matrix d2Y
  inline const ValueType* __restrict__ d2y(int i)
  {
    return &d2Y(i,0);
  }
};
}
#endif
