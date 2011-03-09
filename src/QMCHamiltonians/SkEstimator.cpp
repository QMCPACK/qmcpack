//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <QMCHamiltonians/SkEstimator.h>
#include <LongRange/StructFact.h>
#include <Utilities/IteratorUtility.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus 
{

  SkEstimator::SkEstimator(ParticleSet& source)
  {
    sourcePtcl = &source;
    UpdateMode.set(COLLECTABLE,1);
    NumSpecies=source.getSpeciesSet().getTotalNum();
    NumK=source.SK->KLists.numk;
    OneOverN=1.0/static_cast<RealType>(source.getTotalNum());
    Kshell=source.SK->KLists.kshell;
    MaxKshell=Kshell.size()-1;
    RhokTot.resize(NumK);
    values.resize(NumK);

    Kmag.resize(MaxKshell);
    OneOverDnk.resize(MaxKshell);
    for(int ks=0, k=0; ks<MaxKshell; ks++)
    {
      Kmag[ks]=std::sqrt(source.SK->KLists.ksq[Kshell[ks]]);
      OneOverDnk[ks]=1.0/static_cast<RealType>(Kshell[ks+1]-Kshell[ks]);
    }
    hdf5_out=true;
  }

  void SkEstimator::resetTargetParticleSet(ParticleSet& P)
  {
    sourcePtcl = &P;
  }

  SkEstimator::Return_t SkEstimator::evaluate(ParticleSet& P)
  {
    //sum over species
    std::copy(P.SK->rhok[0],P.SK->rhok[0]+NumK,RhokTot.begin());
    for(int i=1; i<NumSpecies; ++i)
      accumulate_elements(P.SK->rhok[i],P.SK->rhok[i]+NumK,RhokTot.begin());

    if(hdf5_out)
    {
      Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
      for(int i=myIndex;iit != iit_end;++iit,++i)
        P.Collectables[i]+=OneOverN*((*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag());
    }
    else
    {
      Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
      for(int i=0;iit != iit_end;++iit,++i)
        values[i]=OneOverN*((*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag());
    }
    
    return 0.0;
  }

  void SkEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    if(hdf5_out)
    {
      myIndex=collectables.size();
      vector<RealType> tmp(NumK);
      collectables.add(tmp.begin(),tmp.end());
    }
    else
    {
      myIndex=plist.size();
      for (int i=0;i<NumK;i++)
      {
        std::stringstream sstr;
        sstr << "sk_" <<i;
        int id=plist.add(sstr.str());
      }
    }
  }
  
    void SkEstimator::addObservables(PropertySetType& plist )
  { 
    myIndex=plist.size();
    for (int i=0;i<NumK;i++)
    {
      std::stringstream sstr;
      sstr << "sk_" <<i;
      int id=plist.add(sstr.str());
    } 
  }
  
    void SkEstimator::setObservables(PropertySetType& plist)
  {
    if (!hdf5_out)
      std::copy(values.begin(),values.end(),plist.begin()+myIndex);
  }

  void SkEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    if (!hdf5_out)
      std::copy(values.begin(),values.end(),plist.begin()+myIndex+offset);
  }
  

  void SkEstimator::registerCollectables(vector<observable_helper*>& h5desc
      , hid_t gid) const
  {
    if (hdf5_out)
    {
      vector<int> ndim(1,NumK);
      observable_helper* h5o=new observable_helper(myName);
      h5o->set_dimensions(ndim,myIndex);
      h5o->open(gid);
      h5desc.push_back(h5o);
      
      hsize_t kdims[2];
      kdims[0] = NumK;
      kdims[1] = OHMMS_DIM;
      string kpath = myName + "/kpoints";
      hid_t k_space = H5Screate_simple(2,kdims, NULL);
      hid_t k_set   = H5Dcreate (gid, kpath.c_str(), H5T_NATIVE_DOUBLE, k_space, H5P_DEFAULT);
      hid_t mem_space = H5Screate_simple (2, kdims, NULL);
      double *ptr = &(sourcePtcl->SK->KLists.kpts_cart[0][0]);
      herr_t ret = H5Dwrite(k_set, H5T_NATIVE_DOUBLE, mem_space, k_space, H5P_DEFAULT, ptr);
      H5Dclose (k_set);
      H5Sclose (mem_space);
      H5Sclose (k_space);
      H5Fflush(gid, H5F_SCOPE_GLOBAL);
    }
  }

  bool SkEstimator::put(xmlNodePtr cur)
  {
    OhmmsAttributeSet pAttrib;
    string hdf5_flag="no";
    pAttrib.add(hdf5_flag,"hdf5");
    pAttrib.put(cur);

    if (hdf5_flag=="yes") 
      hdf5_out=true;
    else 
      hdf5_out=false;
    
    return true;
  }

  bool SkEstimator::get(std::ostream& os) const
  {
    return true;
  }

  QMCHamiltonianBase* SkEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    SkEstimator* myclone = new SkEstimator(*this);
    myclone->hdf5_out=hdf5_out;
    myclone->myIndex=myIndex;
    return myclone;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
