//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_BOOSTRANDOM_H
#define OHMMS_BOOSTRANDOM_H

#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif

#include <ctime>        
#include <sstream>        
#include <boost/random.hpp>
#include "Message/Communicate.h"
#include "OhmmsData/HDFAttribIO.h"

template<class T>
class BoostRandom {

public:

  typedef T Return_t;
  typedef boost::minstd_rand      base_generator_type;
  typedef boost::mt19937          generator_type;
  typedef boost::variate_generator<generator_type,boost::uniform_real<T> > uniform_generator_type;

  BoostRandom(): baseSeed(0), state_array(0), unit_dist(0,1), generator(0), uni(0) { }

  BoostRandom(int iseed): 
    baseSeed(iseed), unit_dist(0,1), generator(0), uni(0)
  { }
 
  ~BoostRandom() {
    if(uni) delete uni;
    if(generator) delete generator;
    if(state_array) delete [] state_array;
  }

  /** initialize the generator
   * @param i thread index
   * @param nstr number of threads
   * @param iseed_in input seed
   *
   * Initialize generator with the seed. 
   */
  void init(int i, int nstr, int iseed_in) {
    baseSeed=iseed_in;
    if(iseed_in<=0) {
      baseSeed=static_cast<uint32_t>(std::time(0))%16081+(i+1)*nstr+i;
    } 

    //use the constructor with the seed
    //generator->seed(baseSeed) is not working!!!
    if(generator == 0) {
      generator = new generator_type(baseSeed);
      uni = new uniform_generator_type(*generator,unit_dist);
      std::stringstream a;
      a << uni->engine();
      state_size=a.str().size();
      state_array=new char[state_size+128];
      state_name="mt19937";
    }

  }

  inline uniform_generator_type& getGenerator() { return *uni;}

  //!< return [0,1)
  inline T getRandom() { return (*uni)(); }

  inline Return_t operator()() { return (*uni)();} 
  inline int irand() { return 1;}

  inline void bivariate(Return_t& g1, Return_t &g2) {
    Return_t v1, v2, r;
    do {
    v1 = 2.0e0*((*uni)()) - 1.0e0;
    v2 = 2.0e0*((*uni)()) - 1.0e0;
    r = v1*v1+v2*v2;
    } while(r > 1.0e0);
    Return_t fac = sqrt(-2.0e0*log(r)/r);
    g1 = v1*fac;
    g2 = v2*fac;
  }

  inline void read(std::istream& rin) {
    rin >> uni->engine();
  }

  inline void write(std::ostream& rout) const {
    rout << uni->engine();
  }

  /** read the state from a hdf5 group
   * @param gid group id
   */
  inline void read(hid_t gid) {
    hid_t dataset = H5Dopen(gid,state_name.c_str());
    hid_t ret= H5Dread(dataset,H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,state_array);
    H5Dclose(dataset);

    std::istringstream a(state_array);
    a >> uni->engine();
  }

  /** write the state to a hdf5 group
   * @param gid group id
   * @param overwrite if true, open the dataset
   */
  inline void write(hid_t gid, bool overwrite) {
    std::ostringstream a;
    a << uni->engine();
    if(overwrite) {
      hid_t dataset = H5Dopen(gid,state_name.c_str());
      hid_t ret = H5Dwrite(dataset,H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,a.str().c_str());
      H5Dclose(dataset);
    } else {
      hid_t dataspace = H5Screate_simple(1,&state_size,NULL);
      hid_t dataset = H5Dcreate(gid,state_name.c_str(),H5T_NATIVE_CHAR,dataspace,H5P_DEFAULT);
      hid_t ret = H5Dwrite(dataset,H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,a.str().c_str());
      H5Dclose(dataset);
      H5Sclose(dataspace);
    }
  }

private:
  uint32_t baseSeed;
  char* state_array;
  std::string state_name;
  hsize_t state_size;
  base_generator_type    base_generator;
  boost::uniform_real<T> unit_dist;
  generator_type         *generator;
  uniform_generator_type *uni;
};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
