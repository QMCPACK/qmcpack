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

namespace qmcplusplus
{

int PWBasis::readbasis(hid_t h5basisgroup, RealType ecutoff, ParticleLayout_t &lat,
                       const string& pwname,
                       const string& pwmultname,
                       bool resizeContainer)
{
  ///make a local copy
  Lattice=lat;
  ecut = ecutoff;
  if(pwmultname[0] != '0')
  {
    app_log() << "  PWBasis::" << pwmultname << " is found " << endl;
    HDFAttribIO<std::vector<GIndex_t> >hdfvtv(gvecs);
    hdfvtv.read(h5basisgroup,pwmultname.c_str());
  }
  if(pwname[0] != '0')
  {
    app_log() << "  PWBasis::" << pwname << " is found " << endl;
    HDFAttribIO<vector<PosType> > hdfg(kplusgvecs_cart);
    hdfg.read(h5basisgroup,pwname.c_str()); //"planewaves");
  }
  //hdfvtv.read(h5basisgroup,pwname.c_str()); //"planewaves");
  NumPlaneWaves=std::max(gvecs.size(),kplusgvecs_cart.size());
  if(NumPlaneWaves ==0)
  {
    app_error() << "  PWBasis::readbasis Basis is missing. Abort " << endl;
    abort();//FIX_ABORT
  }
  if(kplusgvecs_cart.empty())
  {
    kplusgvecs_cart.resize(NumPlaneWaves);
    for(int i=0; i<NumPlaneWaves; i++)
      kplusgvecs_cart[i]=Lattice.k_cart(gvecs[i]);
  }
  //app_log() << "  Gx Gy Gz " << endl;
  //for(int i=0; i<kplusgvecs_cart.size(); i++)
  //{
  //  app_log() << kplusgvecs_cart[i] << endl;
  //}
  //Now remove elements outside Ecut. At the same time, fill k+G and |k+G| lists.
  //Also keep track of the original index ordering (using indexmap[]) so that
  //orbital coefficients can be ordered and trimmed for ecut in the same way.
  //support older parser
  if(resizeContainer)
    reset();
  //std::copy(gvecs.begin(),gvecs.end(),ostream_iterator<GIndex_t>(cout,"\n"));
  return NumPlaneWaves;
}

void PWBasis::setTwistAngle(const PosType& tang)
{
  PosType dang=twist-tang;
  bool sameTwist = dot(dang,dang)<numeric_limits<RealType>::epsilon();
  if(maxmaxg && sameTwist)
    return;
  twist=tang;
  reset();
}

void PWBasis::reset()
{
  trimforecut();
  //logC.resize(3,2*maxmaxg+1);
  Z.resize(NumPlaneWaves,2+DIM);
  Zv.resize(NumPlaneWaves);
}

/** Remove basis elements if kinetic energy > ecut.
 *
 * Keep and indexmap so we know how to match coefficients on read.
 */
void PWBasis::trimforecut()
{
  //Convert the twist angle to Cartesian coordinates.
  twist_cart = Lattice.k_cart(twist);
  inputmap.resize(NumPlaneWaves);
  app_log() << "  PWBasis::TwistAngle (unit) =" << twist << endl;
  app_log() << "  PWBasis::TwistAngle (cart) =" << twist_cart << endl;
  app_log() << "  PWBasis::trimforecut NumPlaneWaves (before) =" << NumPlaneWaves << endl;
  vector<GIndex_t> gvecCopy(gvecs);
  vector<PosType> gcartCopy(kplusgvecs_cart);
  gvecs.clear();
  kplusgvecs_cart.clear();
  minusModKplusG2.reserve(NumPlaneWaves);
  RealType kcutoff2 = 2.0*ecut; //std::sqrt(2.0*ecut);
  int ngIn=NumPlaneWaves;
  for(int ig=0, newig=0; ig<ngIn; ig++)
  {
    //PosType tempvec = Lattice.k_cart(gvecCopy[ig]+twist);
    PosType tempvec = gcartCopy[ig]+twist_cart;
    RealType mod2 = dot(tempvec,tempvec);
    if(mod2<=kcutoff2)
    {
      gvecs.push_back(gvecCopy[ig]);
      kplusgvecs_cart.push_back(tempvec);
      minusModKplusG2.push_back(-mod2);
      //Remember which position in the HDF5 file this came from...for coefficients
      inputmap[ig] = newig++;
    }
    else
    {
      inputmap[ig] = -1; //Temporary value...need to know final NumPlaneWaves.
      NumPlaneWaves--;
    }
  }
#if defined(PWBASIS_USE_RECURSIVE)
  //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
  for(int ig=0; ig<NumPlaneWaves; ig++)
    for(int i=0; i<OHMMS_DIM; i++)
      if(abs(gvecs[ig][i]) > maxg[i])
        maxg[i] = abs(gvecs[ig][i]);
  gvecs_shifted.resize(NumPlaneWaves);
  for(int ig=0; ig<NumPlaneWaves; ig++)
    gvecs_shifted[ig]=gvecs[ig]+maxg;
  maxmaxg = std::max(maxg[0],std::max(maxg[1],maxg[2]));
  //changes the order???? ok
  C.resize(3,2*maxmaxg+2);
#else
  maxmaxg=1;
#endif
//    //make a copy of input to gvecCopy
////    for(int ig=0, newig=0; ig<ngIn; ig++) {
//      //Check size of this g-vector
//      PosType tempvec = Lattice.k_cart(gvecCopy[ig]+twist);
//      RealType mod2 = dot(tempvec,tempvec);
//      if(mod2<=kcutoff2){ //Keep this element
//        gvecs.push_back(gvecCopy[ig]);
//        kplusgvecs_cart.push_back(tempvec);
//        minusModKplusG2.push_back(-mod2);
//        //Remember which position in the HDF5 file this came from...for coefficients
//        inputmap[ig] = newig++;
////#if !defined(QMC_COMPLEX)
////        //Build the negative vector. See comment at declaration (above) for details.
////        if(gvecCopy[ig][0] < 0)
////          negative.push_back(0);
////        else if(gvecCopy[ig][0] > 0)
////          negative.push_back(1);
////        else { //gx == 0, test gy
////          if(gvecCopy[ig][1] < 0)
////            negative.push_back(0);
////          else if(gvecCopy[ig][1] > 0)
////            negative.push_back(1);
////          else { //gx == gy == 0; test gz. If gz==0 also, take negative=1 (arbitrary)
////            if(gvecCopy[ig][2] < 0)
////              negative.push_back(0);
////            else
////              negative.push_back(1);
////          }
////        }
////#endif
//      } else {
//        inputmap[ig] = -1; //Temporary value...need to know final NumPlaneWaves.
//        NumPlaneWaves--;
//      }
//    }
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
