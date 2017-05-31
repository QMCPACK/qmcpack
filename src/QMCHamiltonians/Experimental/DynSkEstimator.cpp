//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCHamiltonians/DynSkEstimator.h>
#include <LongRange/StructFact.h>
#include <Utilities/IteratorUtility.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

DynSkEstimator::DynSkEstimator(ParticleSet& source)
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
  hdf5_out=false;
}

void DynSkEstimator::resetTargetParticleSet(ParticleSet& P)
{
  sourcePtcl = &P;
}

DynSkEstimator::Return_t DynSkEstimator::calculate(ParticleSet& P)
{
  for(int j(0); j<NumT-MinT; j++)
  {
//       j steps ago, all pindex synched
    int oindex=tWalker->PHindex[pindx]-1-MinT;
    if(oindex<0)
      oindex+=NumT;
    cindx = oindex-j;
    if(cindx<0)
      cindx+=NumT;
    for(int i=0; i<NumK; ++i)
    {
      RealType then_real=tWalker->PropertyHistory[pindx+2*i][cindx];
      RealType now_real=tWalker->PropertyHistory[pindx+2*i][oindex];
      RealType then_imag=tWalker->PropertyHistory[pindx+2*i+1][cindx];
      RealType now_imag=tWalker->PropertyHistory[pindx+2*i+1][oindex];
      values[i+NumK*j]=OneOverN*(now_real*then_real+now_imag*then_imag );
    }
  }
  return 0.0;
}


DynSkEstimator::Return_t DynSkEstimator::evaluate(ParticleSet& P)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT(" DynSkEstimator::evaluate");
#else
  //sum over species for t=0 (now)
  copy(P.SK->rhok[0],P.SK->rhok[0]+NumK,RhokTot.begin());
  for(int i=1; i<NumSpecies; ++i)
    accumulate_elements(P.SK->rhok[i],P.SK->rhok[i]+NumK,RhokTot.begin());
  for(int i=0; i<NumK; ++i)
  {
    tWalker->addPropertyHistoryPoint(pindx+2*i  ,  RhokTot[i].real());
    tWalker->addPropertyHistoryPoint(pindx+2*i+1,  RhokTot[i].imag());
  }
#endif
  calculate(P);
  return 0.0;
}

void DynSkEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
//     if(hdf5_out)
//     {
//       myIndex=collectables.size();
//       std::vector<RealType> tmp(NumK);
//       collectables.add(tmp.begin(),tmp.end());
//     }
//     else
//     {
  myIndex=plist.size();
  for (int i=MinT; i<NumT; i++)
  {
    for(int j(0); j<NumK; j++)
    {
      std::stringstream sstr;
      sstr << "dsk_" <<j<<"_"<<i;
      int id=plist.add(sstr.str());
    }
  }
//     }
}

void DynSkEstimator::addObservables(PropertySetType& plist )
{
  myIndex=plist.size();
  for (int i=MinT; i<NumT; i++)
  {
    for(int j(0); j<NumK; j++)
    {
      std::stringstream sstr;
      sstr << "dsk_" <<j<<"_"<<i;
      int id=plist.add(sstr.str());
    }
  }
}

void DynSkEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
    copy(values.begin(),values.end(),plist.begin()+myIndex);
}

void DynSkEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  if (!hdf5_out)
    copy(values.begin(),values.end(),plist.begin()+myIndex+offset);
}


void DynSkEstimator::registerCollectables(std::vector<observable_helper*>& h5desc
    , hid_t gid) const
{
//     if (hdf5_out)
//     {
//       std::vector<int> ndim(1,NumK);
//       observable_helper* h5o=new observable_helper(myName);
//       h5o->set_dimensions(ndim,myIndex);
//       h5o->open(gid);
//       h5desc.push_back(h5o);
//
//       hsize_t kdims[2];
//       kdims[0] = NumK;
//       kdims[1] = OHMMS_DIM;
//       std::string kpath = myName + "/kpoints";
//       hid_t k_space = H5Screate_simple(2,kdims, NULL);
//       hid_t k_set   = H5Dcreate (gid, kpath.c_str(), H5T_NATIVE_DOUBLE, k_space, H5P_DEFAULT);
//       hid_t mem_space = H5Screate_simple (2, kdims, NULL);
//       double *ptr = &(sourcePtcl->SK->KLists.kpts_cart[0][0]);
//       herr_t ret = H5Dwrite(k_set, H5T_NATIVE_DOUBLE, mem_space, k_space, H5P_DEFAULT, ptr);
//       H5Dclose (k_set);
//       H5Sclose (mem_space);
//       H5Sclose (k_space);
//       H5Fflush(gid, H5F_SCOPE_GLOBAL);
//     }
}

bool DynSkEstimator::putSpecial(xmlNodePtr cur, ParticleSet& P)
{
  std::string tagName("dSK");
  NumT=10;
  MinT=0;
  OhmmsAttributeSet Tattrib;
  Tattrib.add(tagName,"name");
  Tattrib.add(NumT,"max");
  Tattrib.add(MinT,"min");
  Tattrib.put(cur);
  app_log()<<"   "<<NumT<<" dynamic SK values will be calculated"<< std::endl;
//     property history is real, so we must add two (complex)
  pindx = P.addPropertyHistory(NumT);
  for(int j(1); j<NumK*2; j++)
    P.addPropertyHistory(NumT);
// we need NumK*NumT values for all observables
  values.resize(NumK*(NumT-MinT));
  cindx=0;
  return true;
}

bool DynSkEstimator::get(std::ostream& os) const
{
  return true;
}

QMCHamiltonianBase* DynSkEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  DynSkEstimator* myclone = new DynSkEstimator(*this);
//     myclone->hdf5_out=hdf5_out;
//     myclone->myIndex=myIndex;
  return myclone;
}
}

