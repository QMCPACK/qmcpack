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
using namespace qmcplusplus;

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

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
  
};

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
void MultiChain::write(hid_t grp) {

  hid_t h_config = H5Gcreate(grp,"MultiChain",0);

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

  H5Gclose(h_config);
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

bool MultiChain::read(hid_t grp){

  hid_t h_config = H5Gopen(grp,"MultiChain");

  int m_version=1;
  HDFAttribIO<int> v(m_version);
  v.read(h_config,"Version");

  int nc(0);
  HDFAttribIO<int> c(nc);
  c.read(h_config,"NumberOfBeads");
  if(nc != Beads.size()) {
    WARNMSG("The number of chains is different. Previous = "  << nc << " Current = " << Beads.size())
  }

  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;

  Buffer_t chain_buffer,bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);

  HDFAttribIO<int> c2(nc);
  c2.read(h_config,"BufferSize");
  if(nc != bead_buffer.size()) {
    ERRORMSG("The buffer size is different. Ignore restart data")
    H5Gclose(h_config);
    return false;
  }

  chain_buffer.rewind();
  copyToBuffer(chain_buffer);

  HDFAttribIO<Buffer_t> mcin(chain_buffer);
  mcin.read(h_config,"state");
  chain_buffer.rewind();
  copyFromBuffer(chain_buffer);

  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end) {
    sprintf(GrpName,"bead%04d",ibead);
    hid_t bead_id = H5Gopen(h_config,GrpName);

    Bead& bead(**bead_it);

    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.read(bead_id,"state");

    bead_buffer.rewind();
    bead.copyFromBuffer(bead_buffer);

    H5Gclose(bead_id);
    ++bead_it;
    ++ibead;
  }

  H5Gclose(h_config);

  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
