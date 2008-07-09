//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/MultiChain.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
namespace qmcplusplus {

/** specialization for PooledData<double>
 */
template<>
struct HDFAttribIO<PooledData<double> >: public HDFAttribIOBase {

  typedef PooledData<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }
 
  inline void write(hid_t grp, const char* name) {
    hsize_t dim=ref.size();
    int rank = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void overwrite(hid_t grp, const char* name) 
  {
    hid_t dataset =  H5Dopen(grp, name);
    hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
  
};

  MultiChain::MultiChain(Bead* abead,int len, int direction, int npsi): 
GrowthDirection(direction), nPsi(npsi), h_config(-1)
{
  //always add number of beads
  Middle = len/2;
  Last = len-1;
  for(int i=0; i<len; i++) Beads.push_back(new Bead(*abead));
  GlobalAction.resize(npsi);   GlobalAction=0.0;
  UmbrellaWeight.resize(npsi); UmbrellaWeight=1.0;
  GlobalSignWgt.resize(npsi);  GlobalSignWgt=0;
  RefSign.resize(npsi); RefSign=0;
  
  Age=0;
  
}

MultiChain::~MultiChain() 
{
  delete_iter(Beads.begin(),Beads.end());
}

void MultiChain::copyFromBuffer(Buffer_t& buf) 
{
  int n(Beads.size());
  buf.get(n);
  buf.get(GrowthDirection);
  buf.get(Middle);
  buf.get(Last);
  buf.get(nPsi);
  buf.get(Age);
  buf.get(GlobalWgt);
  buf.get(GlobalAction.begin(),GlobalAction.end());
  buf.get(UmbrellaWeight.begin(),UmbrellaWeight.end());
  for(int i=0; i<GlobalSignWgt.size(); i++) buf.get(GlobalSignWgt[i]);
  for(int i=0; i<RefSign.size(); i++) buf.get(RefSign[i]);
  //buf.get(GlobalSignWgt.begin(),GlobalSignWgt.end());
  //buf.get(RefSign.begin(),RefSign.end());
}

void MultiChain::copyToBuffer(Buffer_t& buf) {
  double n= static_cast<double>(Beads.size());
  buf.add(n);
  buf.add(GrowthDirection);
  buf.add(Middle);
  buf.add(Last);
  buf.add(nPsi);
  buf.add(Age);
  buf.add(GlobalWgt);
  buf.add(GlobalAction.begin(),GlobalAction.end());
  buf.add(UmbrellaWeight.begin(),UmbrellaWeight.end());
  buf.add(GlobalSignWgt.begin(),GlobalSignWgt.end());
  buf.add(RefSign.begin(),RefSign.end());
}

/** 
 * - MultiChain
 *   -- Version
 *   -- NumberOfBeads
 *   -- state: number of beads, grotwh direction etc, see MultiChain::copyToBuffer 
 *   -- bead0000
 *      -- state:  everything is stored in a buffer using PooledData<t>
 *      R, Drift, Multiple Gradients,
 *      properties, weights etc, see Bead::copyToBuffer
 */
void MultiChain::open(const string& aroot) {
  hid_t h_file=-1;
  string h5file=aroot+".config.h5";
  h_file  =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  h_config = H5Gcreate(h_file,"MultiChain",0);
  int m_version=1;
  HDFAttribIO<int> v(m_version);
  v.write(h_config,"Version");

  int nc=Beads.size();
  HDFAttribIO<int> c(nc);
  c.write(h_config,"NumberOfBeads");

  //typedef Bead::PosType PosType;
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;

  Buffer_t chain_buffer,bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);
  nc=bead_buffer.size();
  HDFAttribIO<int> c2(nc);
  c2.write(h_config,"BufferSize");

  chain_buffer.rewind();
  copyToBuffer(chain_buffer);

  HDFAttribIO<Buffer_t> mcout(chain_buffer);
  mcout.write(h_config,"state");

  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end) {
    sprintf(GrpName,"bead%04d",ibead);
    hid_t bead_id = H5Gcreate(h_config,GrpName,0);

    bead_buffer.rewind();
    Bead& bead(**bead_it);
    bead.copyToBuffer(bead_buffer);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.write(bead_id,"state");

    H5Gclose(bead_id);
    ++bead_it;
    ++ibead;
  }
  if(h_file>-1) H5Fclose(h_file);
}

void MultiChain::record() {
  //h_config = H5Gopen(h_file,"MultiChain");
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;

  Buffer_t chain_buffer,bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);

  chain_buffer.rewind();
  copyToBuffer(chain_buffer);

  HDFAttribIO<Buffer_t> mcout(chain_buffer);
  mcout.overwrite(h_config,"state");

  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end) {
    sprintf(GrpName,"bead%04d",ibead);
    hid_t bead_id = H5Gopen(h_config,GrpName);
    bead_buffer.rewind();
    Bead& bead(**bead_it);
    bead.copyToBuffer(bead_buffer);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.overwrite(bead_id,"state");
    H5Gclose(bead_id);
    ++bead_it;
    ++ibead;
  }
}

bool MultiChain::read(const string& aroot){

  string h5file = aroot;
  string ext=getExtension(h5file);
  if(ext != "h5") { //if the filename does not h5 extension, add the extension
    h5file.append(".config.h5");
  }

  hid_t  h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  bool success = read(h_file);
  H5Fclose(h_file);
  return success;
}

void MultiChain::close() 
{
  H5Gclose(h_config);
}

bool MultiChain::read(hid_t grp){

  hid_t hgrp = H5Gopen(grp,"MultiChain");

  int m_version=1;
  HDFAttribIO<int> v(m_version);
  v.read(hgrp,"Version");

  int nc(0);
  HDFAttribIO<int> c(nc);
  c.read(hgrp,"NumberOfBeads");
  if(nc != Beads.size()) {
    WARNMSG("The number of chains is different. Previous = "  << nc << " Current = " << Beads.size())
  }

  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;

  Buffer_t chain_buffer,bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);

  HDFAttribIO<int> c2(nc);
  c2.read(hgrp,"BufferSize");
  if(nc != bead_buffer.size()) {
    ERRORMSG("The buffer size is different. Ignore restart data")
        H5Gclose(hgrp);
    return false;
  }

  chain_buffer.rewind();
  copyToBuffer(chain_buffer);

  HDFAttribIO<Buffer_t> mcin(chain_buffer);
  mcin.read(hgrp,"state");
  chain_buffer.rewind();
  copyFromBuffer(chain_buffer);

  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end) {
    sprintf(GrpName,"bead%04d",ibead);
    hid_t bead_id = H5Gopen(hgrp,GrpName);

    Bead& bead(**bead_it);

    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.read(bead_id,"state");

    bead_buffer.rewind();
    bead.copyFromBuffer(bead_buffer);

    H5Gclose(bead_id);
    ++bead_it;
    ++ibead;
  }

  H5Gclose(hgrp);

  return true;
}
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
