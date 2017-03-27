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
    
    


#include "QMCDrivers/MultiChainIO.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
using namespace qmcplusplus;

/** specialization for PooledData<double>
 */
template<>
struct HDFAttribIO<PooledData<double> >: public HDFAttribIOBase
{

  typedef PooledData<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim=ref.size();
    int rank = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }

};
/** Create the HDF5 file "aroot.config.h5" for output.
 *@param aroot the root file name
 *@param append
 *
 * @if append == true
 * the configuration is appended
 * @else
 * the configuration is overwritten
 * The constructor
 * - opens a hdf5 file
 * - create the main groups /multichain and /random_state
 */
HDFMultiChainOutput::HDFMultiChainOutput(const std::string& aroot, int count) :
  Counter(count) , m_version(1)
{
  h5file = aroot;
  h5file.append(".config.h5");
}

/** Destructor closes the HDF5 file and main group. */

HDFMultiChainOutput::~HDFMultiChainOutput()
{
}


/** Write the set of walker configurations to the HDF5 file.
 *@param W set of walker configurations
 *
 * hdf5 file conttains
 * - random_state
 * - counter
 * - multichain
 *   -- state: number of beads, grotwh direction etc, see MultiChain::copyToBuffer
 *   -- bead0000
 *      -- state: properties, weights etc, see Bead::copyToBuffer
 *      -- R
 *      -- Drift
 *      -- Multiple Gradients
 */
bool HDFMultiChainOutput::get(MultiChain& W)
{
  hid_t h_file = H5Fcreate(h5file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  HDFAttribIO<int> v(m_version);
  v.write(h_file,"version");
  hid_t h_random = H5Gcreate(h_file,"random_state",0);
  Random.write(h_random,false);
  H5Gclose(h_random);
  HDFAttribIO<int> c(Counter);
  c.write(h_file,"counter");
  hid_t h_config = H5Gcreate(h_file,"multichain",0);
  //typedef Bead::PosType PosType;
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;
  Buffer_t chain_buffer,bead_buffer;
  chain_buffer.rewind();
  W.copyToBuffer(chain_buffer);
  (*W.begin())->registerData(bead_buffer);
  HDFAttribIO<Buffer_t> mcout(chain_buffer);
  mcout.write(h_config,"state");
  std::deque<Bead*>::iterator bead_it(W.begin());
  std::deque<Bead*>::iterator bead_end(W.end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end)
  {
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
  H5Fclose(h_file);
  return true;
}

/** Open the HDF5 file "aroot.config.h5" for reading.
 *@param aroot the root file name
 *
 *@note The main group is "/config_collection"
 */
HDFMultiChainInput::HDFMultiChainInput(const std::string& aroot)
{
  std::string h5file = aroot;
  std::string ext=getExtension(h5file);
  if(ext != "h5")
    //if the filename does not h5 extension, add the extension
  {
    h5file.append(".config.h5");
  }
}

/** Destructor closes the HDF5 file and main group. */
HDFMultiChainInput::~HDFMultiChainInput() { }

bool
HDFMultiChainInput::put(MultiChain& W)
{
  hid_t h_file = H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  HDFAttribIO<int> v(m_version);
  v.read(h_file,"version");
  hid_t h_random = H5Gopen(h_file,"random_state");
  Random.read(h_random);
  H5Gclose(h_random);
  HDFAttribIO<int> c(Counter);
  c.read(h_file,"counter");
  hid_t h_config = H5Gopen(h_file,"multichain");
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;
  Buffer_t chain_buffer,bead_buffer;
  chain_buffer.rewind();
  W.copyToBuffer(chain_buffer);
  (*W.begin())->registerData(bead_buffer);
  HDFAttribIO<Buffer_t> mcin(chain_buffer);
  mcin.read(h_config,"state");
  std::deque<Bead*>::iterator bead_it(W.begin());
  std::deque<Bead*>::iterator bead_end(W.end());
  //create the group and increment counter
  char GrpName[128];
  int ibead=0;
  while(bead_it != bead_end)
  {
    sprintf(GrpName,"bead%04d",ibead);
    hid_t bead_id = H5Gopen(h_config,GrpName);
    bead_buffer.rewind();
    Bead& bead(**bead_it);
    bead.copyToBuffer(bead_buffer);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.read(bead_id,"state");
    H5Gclose(bead_id);
    ++bead_it;
    ++ibead;
  }
  H5Gclose(h_config);
  H5Fclose(h_file);
  //  hid_t h_random = H5Gopen(h_file,"random_state");
  //  if(h_random>-1) {
  //    LOGMSG("Reading the state of the random number generator from the configuration file")
  //    Random.read(h_random);
  //    H5Gclose(h_random);
  //  }
  return true;
}

