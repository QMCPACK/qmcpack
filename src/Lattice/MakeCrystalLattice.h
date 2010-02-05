//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
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
/**@file MakeCrystalLattice.h
 *@brief Functors to create a lattice with command-line options.
 *
 *The arguments are stored in std::vector<string>. 
 *
 */

#include <cstdlib>

namespace APPNAMESPACE {
/** dummy template class to be specialized */
template<class CL>
struct makelattice { };


/** Specialization of makelattice<CL> for CrystalLattice<T,D>
 *
 * Does nothing but enables specialization for D-dimensional lattice.
 */
template<class T, unsigned D>
struct makelattice<CrystalLattice<T,D> > {
  inline static void apply(CrystalLattice<T,D>& , std::vector<string>& argv) {
  }
};

/** Specialization of makelattice<CL> for CrystalLattice<T,1>*/
template<class T>
struct makelattice<CrystalLattice<T,1> > {

  inline static 
  void 
  apply(CrystalLattice<T,1>& lat, std::vector<string>& argv) {
    int i=0; 
    int argc = argv.size();
    while(i<argc) {
      if(argv[i] == "a0") {
        lat.R(0,0) = std::atof(argv[++i].c_str());
      }
      i++;
    }
    lat.reset();
  }
};

/** Specialization of makelattice<CL> for CrystalLattice<T,2>*/
template<class T>
struct makelattice<CrystalLattice<T,2> > {

  inline static 
  void 
  apply(CrystalLattice<T,2>& lat, std::vector<string>& argv) {
    T a0 = 1.0e0;
    int i=0; 
    int argc = argv.size();
    while(i<argc) {
      if(argv[i] == "cubic") {
        a0 = std::atof(argv[++i].c_str());
        lat.R.diagonal(1.0);
      } else if(argv[i] == "orthorombic") { 
        lat.R = 0.0e0;
        lat.R(0,0) = std::atof(argv[++i].c_str());
        lat.R(1,1) = std::atof(argv[++i].c_str());
      } else if(argv[i] == "general") {
        lat.R = 0.0e0;
        lat.R(0,0) = std::atof(argv[++i].c_str());
        lat.R(0,1) = std::atof(argv[++i].c_str());
        lat.R(1,0) = std::atof(argv[++i].c_str());
        lat.R(1,1) = std::atof(argv[++i].c_str());
      } 
      i++;
    }
    lat.R *= a0;
    lat.reset();
  }
};

/** Specialization of makelattice<CL> for CrystalLattice<T,3>*/
template<class T>
struct makelattice<CrystalLattice<T,3> > {

  /*! \fn makelattic<CrystalLattice<T,3> >
   *  ::apply(CrystalLattice<T,3>& lattice, vector<string>& argv)
   *  \param lattice an CrystalLattice to be set
   *  \param argv   input parameters
   *  \note Keywords to set a speical 3D primitive cell.
   *  \li \p lattice \p cubic \p a
   *  \li \p lattice \p fcc \p a
   *  \li \p lattice \p bcc \p a
   *  \li \p lattice \p hcp \p a \p [c/a]
   */
  inline static 
  void 
  apply(CrystalLattice<T,3>& lat, std::vector<string>& argv) {
    T a0 = 1.0e0;
    int i=0; 
    int argc = argv.size();
    while(i<argc) {
      if(argv[i] == "cubic") {
        a0 = std::atof(argv[++i].c_str());
        lat.R.diagonal(1.0);
      } else if(argv[i] == "orthorombic") { 
        lat.R = 0.0e0;
        lat.R(0,0) = std::atof(argv[++i].c_str());
        lat.R(1,1) = std::atof(argv[++i].c_str());
        lat.R(2,2) = std::atof(argv[++i].c_str());
      } else if(argv[i] == "fcc") {
        a0 = std::atof(argv[++i].c_str());
        lat.R(0,0) = 0.0; lat.R(0,1) = 0.5; lat.R(0,2) = 0.5;
        lat.R(1,0) = 0.5; lat.R(1,1) = 0.0; lat.R(1,2) = 0.5;
        lat.R(2,0) = 0.5; lat.R(2,1) = 0.5; lat.R(2,2) = 0.0;
      } else if(argv[i] == "bcc") {
        a0 = std::atof(argv[++i].c_str());
        lat.R(0,0) = -0.5; lat.R(0,1) =  0.5; lat.R(0,2) =  0.5;
        lat.R(1,0) =  0.5; lat.R(1,1) = -0.5; lat.R(1,2) =  0.5;
        lat.R(2,0) =  0.5; lat.R(2,1) =  0.5; lat.R(2,2) = -0.5;
      } else if(argv[i] == "hcp") {
        a0 = std::atof(argv[++i].c_str());
        double covera = std::sqrt(8.0/3.0);
        if(argc-i > 1) covera = std::atof(argv[++i].c_str());
        lat.R(0,0) = 0.5*a0; lat.R(0,1) = -std::sqrt(3.0)*0.5*a0; lat.R(0,2) = 0.0;
        lat.R(1,0) = 0.5*a0; lat.R(1,1) =  std::sqrt(3.0)*0.5*a0; lat.R(1,2) = 0.0;
        lat.R(2,0) = 0.0;    lat.R(2,1) =   0.0;             lat.R(2,2) = covera*a0;
        a0 = 1.0e0;
      } 
      i++;
    }
    lat.R *= a0;
    lat.reset();
  }
};
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
