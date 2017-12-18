//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
#if defined(USE_REAL_STRUCT_FACTOR)
  RhokTot_r.resize(NumK);
  RhokTot_i.resize(NumK);
#else
  RhokTot.resize(NumK);
#endif
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
#if defined(USE_REAL_STRUCT_FACTOR)
  //sum over species
  std::copy(P.SK->rhok_r[0],P.SK->rhok_r[0]+NumK,RhokTot_r.begin());
  std::copy(P.SK->rhok_i[0],P.SK->rhok_i[0]+NumK,RhokTot_i.begin());
  for(int i=1; i<NumSpecies; ++i)
    accumulate_elements(P.SK->rhok_r[i],P.SK->rhok_r[i]+NumK,RhokTot_r.begin());
  for(int i=1; i<NumSpecies; ++i)
    accumulate_elements(P.SK->rhok_i[i],P.SK->rhok_i[i]+NumK,RhokTot_i.begin());
  if(hdf5_out)
  {
    Vector<RealType>::const_iterator iit_r(RhokTot_r.begin()),iit_r_end(RhokTot_r.end());
    Vector<RealType>::const_iterator iit_i(RhokTot_i.begin()),iit_i_end(RhokTot_i.end());
    for(int i=myIndex; iit_r != iit_r_end; ++iit_r,++iit_i,++i)
      P.Collectables[i]+=OneOverN*((*iit_r)*(*iit_r)+(*iit_i)*(*iit_i));
  }
  else
  {
    Vector<RealType>::const_iterator iit_r(RhokTot_r.begin()),iit_r_end(RhokTot_r.end());
    Vector<RealType>::const_iterator iit_i(RhokTot_i.begin()),iit_i_end(RhokTot_i.end());
    for(int i=0; iit_r != iit_r_end; ++iit_r,++iit_i,++i)
      values[i]=OneOverN*((*iit_r)*(*iit_r)+(*iit_i)*(*iit_i));
  }
#else
  //sum over species
  std::copy(P.SK->rhok[0],P.SK->rhok[0]+NumK,RhokTot.begin());
  for(int i=1; i<NumSpecies; ++i)
    accumulate_elements(P.SK->rhok[i],P.SK->rhok[i]+NumK,RhokTot.begin());
  if(hdf5_out)
  {
    Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
    for(int i=myIndex; iit != iit_end; ++iit,++i)
      P.Collectables[i]+=OneOverN*((*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag());
  }
  else
  {
    Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
    for(int i=0; iit != iit_end; ++iit,++i)
      values[i]=OneOverN*((*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag());
  }
#endif
  return 0.0;
}

void SkEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if(hdf5_out)
  {
    myIndex=collectables.size();
    std::vector<RealType> tmp(NumK);
    collectables.add(tmp.begin(),tmp.end());
  }
  else
  {
    myIndex=plist.size();
    for (int i=0; i<NumK; i++)
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
  for (int i=0; i<NumK; i++)
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


void SkEstimator::registerCollectables(std::vector<observable_helper*>& h5desc
                                       , hid_t gid) const
{
  if (hdf5_out)
  {
    std::vector<int> ndim(1,NumK);
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,myIndex);
    h5o->open(gid);
    h5desc.push_back(h5o);
    hsize_t kdims[2];
    kdims[0] = NumK;
    kdims[1] = OHMMS_DIM;
    std::string kpath = myName + "/kpoints";
    hid_t k_space = H5Screate_simple(2,kdims, NULL);
    hid_t k_set   = H5Dcreate (gid, kpath.c_str(), H5T_NATIVE_DOUBLE, k_space, H5P_DEFAULT);
    hid_t mem_space = H5Screate_simple (2, kdims, NULL);
    RealType *ptr = &(sourcePtcl->SK->KLists.kpts_cart[0][0]);
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
  std::string hdf5_flag="no";
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

