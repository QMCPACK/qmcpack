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
#include <QMCHamiltonians/SkAllEstimator.h>
#include <LongRange/StructFact.h>
#include <Utilities/IteratorUtility.h>
#include <OhmmsData/AttributeSet.h>
#include "Configuration.h"
#include <fstream>
#include <sstream>
#include <string>
namespace qmcplusplus
{

SkAllEstimator::SkAllEstimator(ParticleSet& source, ParticleSet& target)
{
 app_log()<<"SkAllEstimator Constructor\n"; 
 elns= &target;
  ions= &source;
  NumeSpecies=elns->getSpeciesSet().getTotalNum();
  NumIonSpecies=ions->getSpeciesSet().getTotalNum();
  UpdateMode.set(COLLECTABLE,1);
  
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
  //for values, we are including e-e structure factor, and e-Ion.  So a total of NumIonSpecies+1 structure factors.  
  //+2 for the real and imaginary parts of rho_k^e
  values.resize(3*NumK);
  Kmag.resize(MaxKshell);
  OneOverDnk.resize(MaxKshell);
  for(int ks=0, k=0; ks<MaxKshell; ks++)
  {
    Kmag[ks]=std::sqrt(source.SK->KLists.ksq[Kshell[ks]]);
    OneOverDnk[ks]=1.0/static_cast<RealType>(Kshell[ks+1]-Kshell[ks]);
  }
  hdf5_out=false;

}

void SkAllEstimator::resetTargetParticleSet(ParticleSet& P)
{
  elns = &P;
}


void SkAllEstimator::evaluateIonIon()
{
  std::stringstream ss;
  ss << " kx ky kz";

  std::ofstream skfile;
  std::stringstream filebuffer;
  skfile.open("rhok_IonIon.dat");
  for(int i=0; i<NumIonSpecies; i++)
    ss<<" rho_"<<i<<"_r"<<" rho_"<<i<<"_i";
  
  filebuffer<<ss.str()<<std::endl;
 
  for (int k=0; k<NumK; k++)
  {
    PosType kvec= ions->SK->KLists.kpts_cart[k];
   
    filebuffer<<kvec;
    for (int i=0; i<NumIonSpecies; i++) 
    {
      double rho_i(0.0);
      double rho_r(0.0);
       
#if defined(USE_REAL_STRUCT_FACTOR)
      rho_r=ions->SK->rhok_r[i][k];
      rho_i=ions->SK->rhok_i[i][k];
#else
      rho_r=ions->SK->rhok[i][k].real();
      rho_i=ions->SK->rhok[i][k].imag();
#endif
      filebuffer<<" "<<rho_r<<" "<<rho_i;
    } 
   
    filebuffer<<std::endl;
  }
  
  skfile<<filebuffer.str();

  skfile.close();	
}

SkAllEstimator::Return_t SkAllEstimator::evaluate(ParticleSet& P)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  //sum over species
  std::copy(P.SK->rhok_r[0],P.SK->rhok_r[0]+NumK,RhokTot_r.begin());
  std::copy(P.SK->rhok_i[0],P.SK->rhok_i[0]+NumK,RhokTot_i.begin());
  for(int i=1; i<NumeSpecies; ++i)
    accumulate_elements(P.SK->rhok_r[i],P.SK->rhok_r[i]+NumK,RhokTot_r.begin());
  for(int i=1; i<NumeSpecies; ++i)
    accumulate_elements(P.SK->rhok_i[i],P.SK->rhok_i[i]+NumK,RhokTot_i.begin());
    
    for(int k=0; k<NumK; k++)
      values[k]=RhokTot_r[k]*RhokTot_r[k]+RhokTot_i[k]*RhokTot_i[k];

//    for (int ionSpec=0; ionSpec<NumIonSpecies; ionSpec++)
//    {
//   	 for(int k=0; k<NumK; k++)
//   	 {
//		RealType rhok_A_r(ions->SK->rhok_r[ionSpec][k]), rhok_A_i(ions->SK->rhok_i[ionSpec][k]);
//		values[(ionSpec+1)*NumK+k]=RhokTot_r[k]*rhok_A_r+RhokTot_i[k]*rhok_A_i;
//   	 }	
//    } 
    
    for(int k=0,count=0; k<NumK; k++)
    {
	values[NumK+2*k]=RhokTot_r[k];
	values[NumK+2*k+1]=RhokTot_i[k];
    }
#else
  //sum over species
  std::copy(P.SK->rhok[0],P.SK->rhok[0]+NumK,RhokTot.begin());
  for(int i=1; i<NumeSpecies; ++i)
    accumulate_elements(P.SK->rhok[i],P.SK->rhok[i]+NumK,RhokTot.begin());
    Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
    for(int k=0; k<NumK; k++)
      values[k]=rhok[k].real()*rhok[k].real()+rhok[k].imag()*rhok[k].imag();

//    for (int ionSpec=0; ionSpec<NumIonSpecies; ionSpec++)
//    {
//   	 for(int k=0; k<NumK; k++)
//   	 {
//		RealType rhok_A_r(ions->SK->rhok[ionSpec][k].real()), rhok_A_i(ions->SK->rhok[ionSpec][k].imag());
//		values[(ionSpec+1)*NumK+k]=rhok[k].real()*rhok_A_r+rho_k[k].imag()*rhok_A_i;
//  	 }	
//    } 
    for(int k=0,count=0; k<NumK; k++)
    {
	value[NumK+count]=rhok[k].real();
	count++;
	value[NumK+count]=rhok[k].imag();
	count++;
    }
#endif
  return 0.0;
}

void SkAllEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
//  if(hdf5_out)
//  {
//    myIndex=collectables.size();
//    vector<RealType> tmp(NumK);
//    collectables.add(tmp.begin(),tmp.end());
//  }
//  else
//  {
    
    myIndex=plist.size();
//First the electron structure factor
    for (int i=0; i<NumK; i++)
    {
      std::stringstream sstr;
      sstr << "rhok_e_e_" <<i;
      int id=plist.add(sstr.str());
    }
//Now the e-Ion structure factors.  IonIon are dumped to file, and not evaluated.  
//    for (int ionSpec=0; ionSpec<NumIonSpecies; ionSpec++)
//    { 
//       for (int i=0; i<NumK; i++)
//       {
//         std::stringstream sstr;
//         sstr << "rhok_e_" <<ionSpec<<"_"<<i;
//         int id=plist.add(sstr.str());
//       }
//    }

    for(int k=0; k<NumK; k++)
    {
         std::stringstream sstr1,sstr2;
         sstr1 << "rhok_e_r_"<<k;
	 sstr2 << "rhok_e_i_"<<k;
         int id=plist.add(sstr1.str());
	 int id2=plist.add(sstr2.str());

    }

//  }
}

void SkAllEstimator::addObservables(PropertySetType& plist )
{
    myIndex=plist.size();
//First the electron structure factor
    for (int i=0; i<NumK; i++)
    {
      std::stringstream sstr;
      sstr << "rhok_e_e_" <<i;
      int id=plist.add(sstr.str());
    }
//Now the e-Ion structure factors.  IonIon are dumped to file, and not evaluated.  
//    for (int ionSpec=0; ionSpec<NumIonSpecies; ionSpec++)
//    { 
//       for (int i=0; i<NumK; i++)
//       {
//         std::stringstream sstr;
//         sstr << "rhok_e_" <<ionSpec<<"_"<<i;
//         int id=plist.add(sstr.str());
//       }
//    }
    for(int k=0; k<NumK; k++)
    {
         std::stringstream sstr1,sstr2;
         sstr1 << "rhok_e_r_"<<k;
	 sstr2 << "rhok_e_i_"<<k;
         int id=plist.add(sstr1.str());
	 int id2=plist.add(sstr2.str());

    }

}

void SkAllEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
    std::copy(values.begin(),values.end(),plist.begin()+myIndex);
}

void SkAllEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  if (!hdf5_out)
    std::copy(values.begin(),values.end(),plist.begin()+myIndex+offset);
}


void SkAllEstimator::registerCollectables(std::vector<observable_helper*>& h5desc
                                       , hid_t gid) const
{
/*  if (hdf5_out)
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
  } */
}

bool SkAllEstimator::put(xmlNodePtr cur)
{
  OhmmsAttributeSet pAttrib;
  std::string hdf5_flag="no";
  std::string write_ionion_flag="no";

  pAttrib.add(hdf5_flag,"hdf5");
  pAttrib.add(hdf5_flag,"hdf5");
  pAttrib.add(write_ionion_flag,"writeionion");

  pAttrib.put(cur);
  if (hdf5_flag=="yes")
    hdf5_out=true;
  else
    hdf5_out=false;

  if (write_ionion_flag=="Yes" || write_ionion_flag=="yes")
  {
	app_log()<<"SkAll:  evaluateIonIon()\n";
	 evaluateIonIon();
  }

  return true;
}

bool SkAllEstimator::get(std::ostream& os) const
{
  return true;
}

QMCHamiltonianBase* SkAllEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  SkAllEstimator* myclone = new SkAllEstimator(*this);
  myclone->hdf5_out=hdf5_out;
  myclone->myIndex=myIndex;
  return myclone;
}
}

