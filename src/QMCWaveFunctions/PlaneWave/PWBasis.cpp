//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file PWBasis.cpp
 * @brief Definition of member functions of Plane-wave basis set
 */
#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"

namespace qmcplusplus
{

int PWBasis::readbasis(hid_t h5basisgroup, RealType ecutoff, ParticleLayout_t &lat,
                       const std::string& pwname,
                       const std::string& pwmultname,
                       bool resizeContainer)
{
  ///make a local copy
  Lattice=lat;
  ecut = ecutoff;
  app_log() << "  PWBasis::" << pwmultname << " is found " << std::endl;
  HDFAttribIO<std::vector<GIndex_t> >hdfvtv(gvecs);
  hdfvtv.read(h5basisgroup,"/electrons/kpoint_0/gvectors");
  NumPlaneWaves=std::max(gvecs.size(),kplusgvecs_cart.size());
  if(NumPlaneWaves ==0)
  {
    app_error() << "  PWBasis::readbasis Basis is missing. Abort " << std::endl;
    abort();//FIX_ABORT
  }
  if(kplusgvecs_cart.empty())
  {
    kplusgvecs_cart.resize(NumPlaneWaves);
    for(int i=0; i<NumPlaneWaves; i++)
      kplusgvecs_cart[i]=Lattice.k_cart(gvecs[i]);
  }
  //app_log() << "  Gx Gy Gz " << std::endl;
  //for(int i=0; i<kplusgvecs_cart.size(); i++)
  //{
  //  app_log() << kplusgvecs_cart[i] << std::endl;
  //}
  //Now remove elements outside Ecut. At the same time, fill k+G and |k+G| lists.
  //Also keep track of the original index ordering (using indexmap[]) so that
  //orbital coefficients can be ordered and trimmed for ecut in the same way.
  //support older parser
  if(resizeContainer)
    reset();
  //std::copy(gvecs.begin(),gvecs.end(),std::ostream_iterator<GIndex_t>(std::cout,"\n"));
  return NumPlaneWaves;
}

void PWBasis::setTwistAngle(const PosType& tang)
{
  PosType dang=twist-tang;
  bool sameTwist = dot(dang,dang)<std::numeric_limits<RealType>::epsilon();
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
  phi.resize(NumPlaneWaves);
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
  app_log() << "  PWBasis::TwistAngle (unit) =" << twist << std::endl;
  app_log() << "  PWBasis::TwistAngle (cart) =" << twist_cart << std::endl;
  app_log() << "  PWBasis::trimforecut NumPlaneWaves (before) =" << NumPlaneWaves << std::endl;
  std::vector<GIndex_t> gvecCopy(gvecs);
  std::vector<PosType> gcartCopy(kplusgvecs_cart);
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

   // Keep all the g-vectors
   // The cutoff energy is not stored in the HDF file now.
   // Is truncating the gvectors to a spherical shell necessary?
   if (true)
   {
      gvecs.push_back(gvecCopy[ig]);
      kplusgvecs_cart.push_back(tempvec);
      minusModKplusG2.push_back(-mod2);
      //Remember which position in the HDF5 file this came from...for coefficients
      inputmap[ig] = newig++;
    }
#if 0
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
#endif
  }
#if defined(PWBASIS_USE_RECURSIVE)
  //Store the maximum number of translations, within ecut, of any reciprocal cell vector.
  for(int ig=0; ig<NumPlaneWaves; ig++)
    for(int i=0; i<OHMMS_DIM; i++)
      if(std::abs(gvecs[ig][i]) > maxg[i])
        maxg[i] = std::abs(gvecs[ig][i]);
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
  app_log() << "                       NumPlaneWaves (after)  =" <<NumPlaneWaves << std::endl;
}
}
