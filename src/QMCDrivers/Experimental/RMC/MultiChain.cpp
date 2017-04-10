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
    
    


#include "QMCDrivers/MultiChain.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
namespace qmcplusplus
{

void Bead::registerData(Buffer_t& buf)
{
  buf.clear();
  buf.rewind();
  buf.add(get_first_address(R),get_last_address(R));
  buf.add(get_first_address(Drift),get_last_address(Drift));
  for(int i=0; i<Gradients.size(); ++i)
    buf.add(get_first_address(*Gradients[i]),get_last_address(*Gradients[i]));
  for(int i=0; i<Laplacians.size(); ++i)
    buf.add(Laplacians[i]->first_address(),Laplacians[i]->last_address());
  //buf.add(Laplacians[i]->begin(),Laplacians[i]->end());
  for(int i=0; i<DriftVectors.size(); ++i)
    buf.add(get_first_address(*DriftVectors[i]),get_last_address(*DriftVectors[i]));
  buf.add(BeadSignWgt.begin(),BeadSignWgt.end());
  buf.add(TransProb[0]);
  buf.add(TransProb[1]);
  buf.add(Action.begin(),Action.end());
  buf.add(Properties.begin(),Properties.end());
  buf.add(deltaRSquared.begin(),deltaRSquared.end());
}

void Bead::copyFromBuffer(Buffer_t& buf)
{
  buf.rewind();
  buf.get(get_first_address(R),get_last_address(R));
  buf.get(get_first_address(Drift),get_last_address(Drift));
  for(int i=0; i<Gradients.size(); ++i)
    buf.get(get_first_address(*Gradients[i]),get_last_address(*Gradients[i]));
  for(int i=0; i<Laplacians.size(); ++i)
    buf.get(Laplacians[i]->begin(),Laplacians[i]->end());
  for(int i=0; i<DriftVectors.size(); ++i)
    buf.get(get_first_address(*DriftVectors[i]),get_last_address(*DriftVectors[i]));
  buf.get(BeadSignWgt.begin(),BeadSignWgt.end());
  buf.get(TransProb[0]);
  buf.get(TransProb[1]);
  buf.get(Action.begin(),Action.end());
  buf.get(Properties.begin(),Properties.end());
  buf.get(deltaRSquared.begin(),deltaRSquared.end());
}

void Bead::copyToBuffer(Buffer_t& buf)
{
  buf.rewind();
  buf.put(get_first_address(R),get_last_address(R));
  buf.put(get_first_address(Drift),get_last_address(Drift));
  for(int i=0; i<Gradients.size(); ++i)
    buf.put(get_first_address(*Gradients[i]),get_last_address(*Gradients[i]));
  for(int i=0; i<Laplacians.size(); ++i)
    buf.add(Laplacians[i]->first_address(),Laplacians[i]->last_address());
//      buf.put(Laplacians[i]->begin(),Laplacians[i]->end());
  for(int i=0; i<DriftVectors.size(); ++i)
    buf.put(get_first_address(*DriftVectors[i]),get_last_address(*DriftVectors[i]));
  buf.put(BeadSignWgt.begin(),BeadSignWgt.end());
  buf.put(TransProb[0]);
  buf.put(TransProb[1]);
  buf.put(Action.begin(),Action.end());
  buf.put(Properties.begin(),Properties.end());
  buf.put(deltaRSquared.begin(),deltaRSquared.end());
}

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

  inline void overwrite(hid_t grp, const char* name)
  {
    hid_t dataset =  H5Dopen(grp, name);
    hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
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
  for(int i=0; i<len; i++)
    Beads.push_back(new Bead(*abead));
  GlobalAction.resize(npsi);
  GlobalAction=0.0;
  UmbrellaWeight.resize(npsi);
  UmbrellaWeight=1.0;
  GlobalSignWgt.resize(npsi);
  GlobalSignWgt=0;
  RefSign.resize(npsi);
  RefSign=0;
  Age=0;
}

MultiChain::~MultiChain()
{
  delete_iter(Beads.begin(),Beads.end());
}

void MultiChain::copyFromBuffer(Buffer_t& buf)
{
  buf.rewind();
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
  //for(int i=0; i<GlobalSignWgt.size(); i++) buf.get(GlobalSignWgt[i]);
  //for(int i=0; i<RefSign.size(); i++) buf.get(RefSign[i]);
  buf.get(GlobalSignWgt.begin(),GlobalSignWgt.end());
  buf.get(RefSign.begin(),RefSign.end());
}

void MultiChain::copyToBuffer(Buffer_t& buf)
{
  buf.clear();
  buf.rewind();
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
 *   -- BeadBufferSize
 *   -- state: number of beads, grotwh direction etc, see MultiChain::copyToBuffer
 *   -- beada_#
 *      -- state:  everything is stored in a buffer using PooledData<t>
 *      R, Drift, Multiple Gradients,
 *      properties, weights etc, see Bead::copyToBuffer
 */
void MultiChain::open(const std::string& aroot)
{
  hid_t h_file=-1;
  std::string h5file=aroot+".config.h5";
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
  //this buffer is used only here and nothing to do with other buffer.
  //reused by all the beads in the chain
  Buffer_t bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);
  BeadBufferSize=bead_buffer.size();
  app_log() << "  MultiChain::Bead Buffer size = " << BeadBufferSize << std::endl;
  HDFAttribIO<int> c2(BeadBufferSize);
  c2.write(h_config,"BeadBufferSize");
  Buffer_t chain_buffer;
  copyToBuffer(chain_buffer);
  MyBufferSize=chain_buffer.size();
  HDFAttribIO<Buffer_t> mcout(chain_buffer);
  mcout.write(h_config,"state");
  app_log() << "  MultiChain state Buffer size = " << MyBufferSize << std::endl;
  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  //create the group and increment counter
  char GrpName[128];
  for(int ibead=0; bead_it!= bead_end; ++ibead, ++bead_it)
  {
    (*bead_it)->copyToBuffer(bead_buffer);
    sprintf(GrpName,"bead_%d",ibead);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.write(h_config,GrpName);
  }
  if(h_file>-1)
    H5Fclose(h_file);
}

void MultiChain::record()
{
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;
  Buffer_t chain_buffer;
  chain_buffer.reserve(MyBufferSize);
  this->copyToBuffer(chain_buffer);
  Buffer_t bead_buffer;
  bead_buffer.resize(BeadBufferSize);
  HDFAttribIO<Buffer_t> mcout(chain_buffer);
  mcout.overwrite(h_config,"state");
  int nc=Beads.size();
  HDFAttribIO<int> c(nc);
  c.write(h_config,"NumberOfBeads");
  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  char GrpName[128];
  for(int ibead=0; bead_it!=bead_end; ++ibead,++bead_it)
  {
    (*bead_it)->copyToBuffer(bead_buffer);
    sprintf(GrpName,"bead_%d",ibead);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.overwrite(h_config,GrpName);
  }
}

bool MultiChain::read(hid_t grp)
{
  hid_t hgrp = H5Gopen(grp,"MultiChain");
  int m_version=1;
  HDFAttribIO<int> v(m_version);
  v.read(hgrp,"Version");
  int nc(0);
  HDFAttribIO<int> c(nc);
  c.read(hgrp,"NumberOfBeads");
  if(nc != Beads.size())
  {
    app_warning() << "MultiChain::read stopped due to the difference in the number of beads"<< std::endl;
    return false;
  }
  typedef Bead::RealType PosType;
  typedef Bead::Buffer_t Buffer_t;
  Buffer_t bead_buffer;
  (*(this->begin()))->registerData(bead_buffer);
  BeadBufferSize=bead_buffer.size();
  HDFAttribIO<int> c2(nc);
  c2.read(hgrp,"BeadBufferSize");
  if(nc != bead_buffer.size())
  {
    APP_ABORT("MultiChain::read due to the difference in the buffer bead size");
  }
  Buffer_t chain_buffer;
  copyToBuffer(chain_buffer);//this is just to measure the size
  HDFAttribIO<Buffer_t> mcin(chain_buffer);
  mcin.read(hgrp,"state");
  copyFromBuffer(chain_buffer);
  std::deque<Bead*>::iterator bead_it(this->begin());
  std::deque<Bead*>::iterator bead_end(this->end());
  char GrpName[128];
  for(int ibead=0; bead_it != bead_end; ++ibead, ++bead_it)
  {
    sprintf(GrpName,"bead_%d",ibead);
    HDFAttribIO<Buffer_t> bout(bead_buffer);
    bout.read(h_config, GrpName);
    (*bead_it)->copyFromBuffer(bead_buffer);
  }
  H5Gclose(hgrp);
  return true;
}
bool MultiChain::read(const std::string& aroot)
{
  std::string h5file = aroot;
  std::string ext=getExtension(h5file);
  if(ext != "h5")
    //if the filename does not h5 extension, add the extension
  {
    h5file.append(".config.h5");
  }
  hid_t  h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  bool success = read(h_file);
  H5Fclose(h_file);
  return success;
}

void MultiChain::close()
{
  if(h_config != -1)
  {
    H5Gclose(h_config);
    h_config=-1;
  }
}

}

