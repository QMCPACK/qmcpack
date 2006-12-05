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

  ///Remove basis elements if kinetic energy > ecut.
  //Keep and indexmap so we know how to match coefficients on read.
  void PWBasis::trimforecut() { 
    PosType tempvec;
    RealType mod2, mod;
    RealType kcutoff = std::sqrt(2.0*ecut);


    //resize inputmap
    NumPlaneWaves = gvecs.size();
    inputmap.resize(NumPlaneWaves);

    app_log() << "  PWBasis::trimforecut NumPlaneWaves (before) =" << NumPlaneWaves << endl;

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

    app_log() << "                       NumPlaneWaves (after)  =" <<NumPlaneWaves << endl;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
