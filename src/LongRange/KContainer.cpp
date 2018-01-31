//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Message/Communicate.h"
#include "LongRange/KContainer.h"
#include <qmc_common.h>
#include <map>

namespace qmcplusplus
{

void
KContainer::UpdateKLists(ParticleLayout_t& lattice, RealType kc, bool useSphere)
{
  kcutoff = kc;
  kcut2 = kc*kc;
  if(kcutoff <= 0.0)
  {
    APP_ABORT("  Illegal cutoff for KContainer");
  }
  FindApproxMMax(lattice);
  BuildKLists(lattice,useSphere);

  app_log() << "  KContainer initialised with cutoff " << kcutoff << std::endl;
  app_log() << "   # of K-shell  = " << kshell.size() << std::endl;
  app_log() << "   # of K points = " << kpts.size() << std::endl;
  app_log() << std::endl;
}

void
KContainer::FindApproxMMax(ParticleLayout_t& lattice)
{
  //Estimate the size of the parallelpiped that encompases a sphere of kcutoff.
  //mmax is stored as integer translations of the reciprocal cell vectors.
  //Does not require an orthorhombic cell.
  /* Old method.
  //2pi is not included in lattice.b
  Matrix<RealType> mmat;
  mmat.resize(3,3);
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++){
      mmat[i][j] = 0.0;
      for(int k=0;k<3;k++)
  mmat[i][j] = mmat[i][j] + 4.0*M_PI*M_PI*lattice.b(k)[i]*lattice.b(j)[k];
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
  for (int i = 0; i < DIM; i++)
    mmax[i] = static_cast<int>(std::floor(std::sqrt(dot(lattice.a(i),lattice.a(i))) * kcutoff / (2 * M_PI))) + 1;


//overwrite the non-periodic directon to be zero
if(qmc_common.use_ewald)
{
  app_log() << "  Using Ewald sum for the slab " << std::endl;
#if OHMMS_DIM==3
  if(lattice.SuperCellEnum == SUPERCELL_SLAB) mmax[2]=0;
//  if(lattice.SuperCellEnum == SUPERCELL_WIRE) mmax[1]=mmax[2]=0;
//#elif OHMMS_DIM==2
//  if(lattice.SuperCellEnum == SUPERCELL_WIRE)
//    mmax[1]=0;
#endif
}

  mmax[DIM]=mmax[0];
  for(int i=1; i<DIM; ++i) mmax[DIM]=std::max(mmax[i],mmax[DIM]);
}

void
KContainer::BuildKLists(ParticleLayout_t& lattice, bool useSphere)
{
  TinyVector<int,DIM+1> TempActualMax;
  TinyVector<int,DIM> kvec;
  TinyVector<RealType,DIM> kvec_cart;
  RealType modk2;
  std::vector<TinyVector<int,DIM> > kpts_tmp;
  VContainer_t kpts_cart_tmp;
  SContainer_t ksq_tmp;
  // reserve the space for memory efficiency
#if OHMMS_DIM ==3
  int numGuess=(2*mmax[0]+1)*(2*mmax[1]+1)*(2*mmax[2]+1);
  if(useSphere)
  {
    //Loop over guesses for valid k-points.
    for(int i=-mmax[0]; i<=mmax[0]; i++)
    {
      kvec[0] = i;
      for(int j=-mmax[1]; j<=mmax[1]; j++)
      {
        kvec[1] = j;
        for(int k=-mmax[2]; k<=mmax[2]; k++)
        {
          kvec[2] = k;
          //Do not include k=0 in evaluations.
          if(i==0 && j==0 && k==0)
            continue;
          //Convert kvec to Cartesian
          kvec_cart = lattice.k_cart(kvec);
          //Find modk
          modk2 = dot(kvec_cart,kvec_cart);
          if(modk2>kcut2)
            continue; //Inside cutoff?
          //This k-point should be added to the list
          kpts_tmp.push_back(kvec);
          kpts_cart_tmp.push_back(kvec_cart);
          ksq_tmp.push_back(modk2);
          //Update record of the allowed maximum translation.
          for(int idim=0; idim<3; idim++)
            if(std::abs(kvec[idim]) > TempActualMax[idim])
              TempActualMax[idim] = std::abs(kvec[idim]);
        }
      }
    }
  }
  else
  {
    // Loop over all k-points in the parallelpiped and add them to kcontainer
    // note layout is for interfacing with fft, so for each dimension, the
    // positive indexes come first then the negative indexes backwards
    // e.g.    0, 1, .... mmax, -mmax+1, -mmax+2, ... -1
    const int idimsize = mmax[0]*2;
    const int jdimsize = mmax[1]*2;
    const int kdimsize = mmax[2]*2;
    for (int i = 0; i < idimsize; i++)
    {
      kvec[0] = i;
      if (kvec[0] > mmax[0])
        kvec[0] -= idimsize;
      for (int j = 0; j < jdimsize; j++)
      {
        kvec[1] = j;
        if (kvec[1] > mmax[1])
          kvec[1] -= jdimsize;
        for (int k = 0; k < kdimsize; k++)
        {
          kvec[2] = k;
          if (kvec[2] > mmax[2])
            kvec[2] -= kdimsize;
          // get cartesian location and modk2
          kvec_cart = lattice.k_cart(kvec);
          modk2 = dot(kvec_cart, kvec_cart);
          // add k-point to lists
          kpts_tmp.push_back(kvec);
          kpts_cart_tmp.push_back(kvec_cart);
          ksq_tmp.push_back(modk2);
        }
      }
    }
    // set allowed maximum translation
    TempActualMax[0] = mmax[0];
    TempActualMax[1] = mmax[1];
    TempActualMax[2] = mmax[2];
  }
#elif OHMMS_DIM == 2
  int numGuess=(2*mmax[0]+1)*(2*mmax[1]+1);
  if(useSphere)
  {
    //Loop over guesses for valid k-points.
    for(int i=-mmax[0]; i<=mmax[0]; i++)
    {
      kvec[0] = i;
      for(int j=-mmax[1]; j<=mmax[1]; j++)
      {
        kvec[1] = j;
        //Do not include k=0 in evaluations.
        if(i==0 && j==0)
          continue;
        //Convert kvec to Cartesian
        kvec_cart = lattice.k_cart(kvec);
        //Find modk
        modk2 = dot(kvec_cart,kvec_cart);
        if(modk2>kcut2)
          continue; //Inside cutoff?
        //This k-point should be added to the list
        kpts_tmp.push_back(kvec);
        kpts_cart_tmp.push_back(kvec_cart);
        ksq_tmp.push_back(modk2);
        //Update record of the allowed maximum translation.
        for(int idim=0; idim<3; idim++)
          if(std::abs(kvec[idim]) > TempActualMax[idim])
            TempActualMax[idim] = std::abs(kvec[idim]);
      }
    }
  }
  else
  {
    // Loop over all k-points in the parallelpiped and add them to kcontainer
    // note layout is for interfacing with fft, so for each dimension, the
    // positive indexes come first then the negative indexes backwards
    // e.g.    0, 1, .... mmax, -mmax+1, -mmax+2, ... -1
    const int idimsize = mmax[0]*2;
    const int jdimsize = mmax[1]*2;
    for (int i = 0; i < idimsize; i++)
    {
      kvec[0] = i;
      if (kvec[0] > mmax[0])
        kvec[0] -= idimsize;
      for (int j = 0; j < jdimsize; j++)
      {
        kvec[1] = j;
        if (kvec[1] > mmax[1])
          kvec[1] -= jdimsize;
        // get cartesian location and modk2
        kvec_cart = lattice.k_cart(kvec);
        modk2 = dot(kvec_cart, kvec_cart);
        // add k-point to lists
        kpts_tmp.push_back(kvec);
        kpts_cart_tmp.push_back(kvec_cart);
        ksq_tmp.push_back(modk2);
      }
    }
    // set allowed maximum translation
    TempActualMax[0] = mmax[0];
    TempActualMax[1] = mmax[1];
  }
//#elif OHMMS_DIM == 1
//add one-dimension
#else
#error "OHMMS_DIM != 2 || OHMMS_DIM != 3"
#endif
  //Update a record of the number of k vectors
  numk = kpts_tmp.size();
  std::map<int,std::vector<int>*>  kpts_sorted;
//create the map: use simple integer with resolution of 0.001 in ksq
  for(int ik=0; ik<numk; ik++)
  {
    int k_ind=static_cast<int>(ksq_tmp[ik]*1000);
    std::map<int,std::vector<int>*>::iterator it(kpts_sorted.find(k_ind));
    if(it == kpts_sorted.end())
    {
      std::vector<int>* newSet=new std::vector<int>;
      kpts_sorted[k_ind]=newSet;
      newSet->push_back(ik);
    }
    else
    {
      (*it).second->push_back(ik);
    }
  }
  std::map<int,std::vector<int>*>::iterator it(kpts_sorted.begin());
  kpts.resize(numk);
  kpts_cart.resize(numk);
  ksq.resize(numk);
  kshell.resize(kpts_sorted.size()+1,0);
  int ok=0, ish=0;
  while(it != kpts_sorted.end())
  {
    std::vector<int>::iterator vit((*it).second->begin());
    while(vit != (*it).second->end())
    {
      int ik=(*vit);
      kpts[ok]=kpts_tmp[ik];
      kpts_cart[ok]=kpts_cart_tmp[ik];
      ksq[ok]=ksq_tmp[ik];
      ++vit;
      ++ok;
    }
    kshell[ish+1]=kshell[ish]+(*it).second->size();
    ++it;
    ++ish;
  }
  it=kpts_sorted.begin();
  std::map<int,std::vector<int>*>::iterator e_it(kpts_sorted.end());
  while(it != e_it)
  {
    delete it->second;
    it++;
  }
  //Finished searching k-points. Copy list of maximum translations.
  mmax[DIM] = 0;
  for(int idim=0; idim<DIM; idim++)
  {
    mmax[idim] = TempActualMax[idim];
    mmax[DIM]=std::max(mmax[idim],mmax[DIM]);
    //if(mmax[idim] > mmax[DIM]) mmax[DIM] = mmax[idim];
  }
  //Now fill the array that returns the index of -k when given the index of k.
  minusk.resize(numk);
  // Create a map from the hash value for each k vector to the index
  std::map<int, int> hashToIndex;
  for (int ki=0; ki<numk; ki++)
  {
    hashToIndex[GetHashOfVec(kpts[ki], numk)] = ki;
  }
  // Use the map to find the index of -k from the index of k
  for(int ki=0; ki<numk; ki++)
  {
    minusk[ki] = hashToIndex[ GetHashOfVec(-1 * kpts[ki], numk) ];
  }
}

}
