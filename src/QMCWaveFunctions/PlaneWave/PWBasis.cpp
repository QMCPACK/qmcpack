//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file PWBasis.cpp
 * @brief Definition of member functions of Plane-wave basis set
 */
#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

  int PWBasis::readbasis(hid_t h5basisgroup, RealType ecutoff, ParticleLayout_t &lat,
      const string& pwname, bool resizeContainer) 
  {
    ///make a local copy
    Lattice=lat;

    ecut = ecutoff;

    HDFAttribIO<std::vector<GIndex_t> >hdfvtv(gvecs);
    hdfvtv.read(h5basisgroup,pwname.c_str()); //"planewaves");
    NumPlaneWaves=gvecs.size();
    //Now remove elements outside Ecut. At the same time, fill k+G and |k+G| lists.
    //Also keep track of the original index ordering (using indexmap[]) so that
    //orbital coefficients can be ordered and trimmed for ecut in the same way.

    //support older parser
    if(resizeContainer) reset();

    //std::copy(gvecs.begin(),gvecs.end(),ostream_iterator<GIndex_t>(cout,"\n"));
    return NumPlaneWaves;
  }

  void PWBasis::setTwistAngle(const PosType& tang)
  {
    PosType dang=twist-tang;
    bool sameTwist = dot(dang,dang)<numeric_limits<RealType>::epsilon();
    if(maxmaxg && sameTwist) return;
    twist=tang;
    reset();
  }

  void PWBasis::reset() 
  {
    trimforecut();
    //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
    for(int ig=0; ig<NumPlaneWaves; ig++)
      for(int i=0; i<3; i++)
        if(abs(gvecs[ig][i]) > maxg[i]) 
          maxg[i] = abs(gvecs[ig][i]);
    maxmaxg = max(maxg[0],max(maxg[1],maxg[2]));

    //changes the order????
    C.resize(3,2*maxmaxg+1);
    //logC.resize(3,2*maxmaxg+1);

    Z.resize(NumPlaneWaves,2+DIM);
    Zv.resize(NumPlaneWaves);
  }

  /** Remove basis elements if kinetic energy > ecut.
   *
   * Keep and indexmap so we know how to match coefficients on read.
   */
  void PWBasis::trimforecut() { 
    //Convert the twist angle to Cartesian coordinates.
    twist_cart = Lattice.k_cart(twist);

    //resize inputmap
    NumPlaneWaves = gvecs.size();
    inputmap.resize(NumPlaneWaves);

    app_log() << "  PWBasis::trimforecut NumPlaneWaves (before) =" << NumPlaneWaves << endl;

    //make a copy of input to gvecCopy
    vector<GIndex_t> gvecCopy(gvecs);
    gvecs.clear();
    gvecs.reserve(gvecCopy.size());

    RealType kcutoff2 = 2.0*ecut; //std::sqrt(2.0*ecut);
    int ngIn=NumPlaneWaves;
    for(int ig=0, newig=0; ig<ngIn; ig++) {
      //Check size of this g-vector
      PosType tempvec = Lattice.k_cart(gvecCopy[ig]+twist);
      RealType mod2 = dot(tempvec,tempvec);
      if(mod2<=kcutoff2){ //Keep this element
        gvecs.push_back(gvecCopy[ig]);
        kplusgvecs_cart.push_back(tempvec);
        minusModKplusG2.push_back(-mod2);
        //Remember which position in the HDF5 file this came from...for coefficients
        inputmap[ig] = newig++;
#if !defined(QMC_COMPLEX)
        //Build the negative vector. See comment at declaration (above) for details.
        if(gvecCopy[ig][0] < 0)
          negative.push_back(0);
        else if(gvecCopy[ig][0] > 0)
          negative.push_back(1);
        else { //gx == 0, test gy
          if(gvecCopy[ig][1] < 0)
            negative.push_back(0);
          else if(gvecCopy[ig][1] > 0)
            negative.push_back(1);
          else { //gx == gy == 0; test gz. If gz==0 also, take negative=1 (arbitrary)
            if(gvecCopy[ig][2] < 0)
              negative.push_back(0);
            else
              negative.push_back(1);
          }
        }
#endif
      } else {
        inputmap[ig] = -1; //Temporary value...need to know final NumPlaneWaves.
        NumPlaneWaves--;
      }
    }
    //Finalize the basis. Fix temporary values of inputmap.
    //for(int ig=0; ig<inputmap.size(); ig++) 
    //  if(inputmap[ig] == -1)
    //    inputmap[ig] = NumPlaneWaves; //For dumping coefficients of PWs>ecut
    app_log() << "                       NumPlaneWaves (after)  =" <<NumPlaneWaves << endl;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
