//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim and Jordan Vincent
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
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include <blitz/array.h>
using namespace ohmmsqmc;

/**
 *@param aroot the root file name
 *@brief Create the HDF5 file "aroot.config.h5" for output. 
 *@note The main group is "/config_collection"
 */

HDFWalkerOutput::HDFWalkerOutput(const string& aroot, bool append):
  Counter(0), AppendMode(append) {

  string h5file = aroot;
  h5file.append(".config.h5");
  h_file = H5Fcreate(h5file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  h_config = H5Gcreate(h_file,"config_collection",0);
}

/** Destructor closes the HDF5 file and main group. */

HDFWalkerOutput::~HDFWalkerOutput() {

  if(AppendMode)  H5Gclose(h_config);
  H5Fclose(h_file);
}


/**
 *@param W set of walker configurations
 *@brief Write the set of walker configurations to the
 HDF5 file.  
*/

bool HDFWalkerOutput::get(MCWalkerConfiguration& W) {

  typedef MCWalkerConfiguration::PosType PosType;
  typedef MCWalkerConfiguration::RealType RealType;
  typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;

  typedef blitz::Array<PosType,2> Array_t;
  typedef Vector<RealType>        Vector_t;

  //OverWrite -> create a new group, close it later
  //if(!AppendMode)  {
  //  cout << "Create a new group " << endl;
  // h_config = H5Gcreate(h_file,"config_collection",0);
  //}

  PropertyContainer_t Properties;
  //2D array of PosTypes (x,y,z) indexed by (walker,particle)
  Array_t Pos;
  Vector_t sample_1, sample_2;

  Pos.resize(W.getActiveWalkers(),W.R.size());
  sample_1.resize(W.getActiveWalkers());
  sample_2.resize(W.getActiveWalkers());

  //store walkers in a temporary array
  int nw = 0; int item=0;
  for (MCWalkerConfiguration::iterator it = W.begin(); 
       it != W.end(); ++it, ++nw) {
    sample_1(nw) = (*it)->Properties(PSISQ);
    sample_2(nw) = (*it)->Properties(LOCALPOTENTIAL);
    for(int np=0; np < W.getParticleNum(); ++np)
      Pos(nw,np) = (*it)->R(np);    
  }
  //create the group and increment counter
  char GrpName[128];
  sprintf(GrpName,"config%04d",Counter++);
  hid_t group_id = H5Gcreate(h_config,GrpName,0);

  //write the dataset
  HDFAttribIO<Array_t> Pos_out(Pos);
  Pos_out.write(group_id,"coord");
  HDFAttribIO<Vector_t> sample_out(sample_1);
  sample_out.write(group_id,"psisq");
  HDFAttribIO<Vector_t> sample_out2(sample_2);
  sample_out2.write(group_id,"localpotential");

  H5Gclose(group_id);

  //closing h_config if overwriting
  //if(!AppendMode)  H5Gclose(h_config);
  //XMLReport("Printing " << W.getActiveWalkers() << " Walkers to file")
  
  return true;
}

/**
 *@param aroot the root file name
 *@brief Open the HDF5 file "aroot.config.h5" for reading. 
 *@note The main group is "/config_collection"
 */

HDFWalkerInput::HDFWalkerInput(const string& aroot):
  Counter(0), NumSets(0) {
  string h5file = aroot;
  h5file.append(".config.h5");
  h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  h_config = H5Gopen(h_file,"config_collection");
  H5Gget_num_objs(h_config,&NumSets);
  if(!NumSets) {
    ERRORMSG("File does not contain walkers!")
   }
}

/** Destructor closes the HDF5 file and main group. */

HDFWalkerInput::~HDFWalkerInput() {
  H5Gclose(h_config);
  H5Fclose(h_file);
}

/**
 *@param W set of walker configurations
 *@brief Write the set of walker configurations to the
 HDF5 file.  
*/

bool HDFWalkerInput::put(MCWalkerConfiguration& W){

  if(Counter >= NumSets) return false;

  int ic = Counter++;
  return put(W,ic);
}

/**
 *@param W set of walker configurations
 *@param ic  
 *@brief Write the set of walker configurations to the
 HDF5 file.  
*/

bool  
HDFWalkerInput::put(MCWalkerConfiguration& W, int ic){

  int selected = ic;
  if(ic<0) {
    XMLReport("Will use the last set from " << NumSets << " of configurations.")
    selected = NumSets-1;
  }

  typedef MCWalkerConfiguration::PosType PosType;
  typedef MCWalkerConfiguration::RealType RealType;
  typedef MCWalkerConfiguration::PropertyContainer_t ProtertyContainer_t;

  typedef blitz::Array<PosType,2> Array_t;
  typedef Vector<RealType>        Vector_t;

  int nwt = 0;
  int npt = 0;
  //2D array of PosTypes (x,y,z) indexed by (walker,particle)
  Array_t Pos_temp;
  Vector_t psisq_in, localene_in;

  //open the group
  char GrpName[128];
  sprintf(GrpName,"config%04d",selected);
  hid_t group_id = H5Gopen(h_config,GrpName);
    
  HDFAttribIO<Array_t> Pos_in(Pos_temp);
  HDFAttribIO<Vector_t> sample1(psisq_in), sample2(localene_in);
  //read the dataset
  Pos_in.read(group_id,"coord");
  sample1.read(group_id,"psisq");
  sample2.read(group_id,"localpotential");
  //close the group
  H5Gclose(group_id);

  /*check to see if the number of walkers and particles is  consistent with W */
  int nptcl = Pos_temp.extent(1);
  nwt = Pos_temp.extent(0);
  if(nwt != W.getActiveWalkers() || nptcl != W.getParticleNum()) {
    W.resize(nwt,nptcl); 
  }

  //assign configurations to W
  int iw=0;
  for(MCWalkerConfiguration::iterator it = W.begin(); it != W.end(); 
      ++it, iw++) {
    (*it)->Properties(PSISQ) = psisq_in[iw];
    for(int np=0; np < W.getParticleNum(); ++np){
      (*it)->R(np) = Pos_temp(iw,np);
    }
  }

  if(localene_in.size()>0) {
    iw=0;
    for(MCWalkerConfiguration::iterator it = W.begin(); it != W.end(); 
	++it, iw++) 
      (*it)->Properties(LOCALPOTENTIAL) = localene_in[iw];
  }
  //XMLReport("Read in " << W.getActiveWalkers() << " Walkers from file" << GrpName)
  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
