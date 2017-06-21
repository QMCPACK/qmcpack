//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file MakeCrystalLattice.h
 *@brief Functors to create a lattice with command-line options.
 *
 *The arguments are stored in std::vector<std::string>.
 *
 */

#include <cstdlib>

namespace qmcplusplus
{
/** dummy template class to be specialized */
template<class CL>
struct makelattice { };


/** Specialization of makelattice<CL> for CrystalLattice<T,D>
 *
 * Does nothing but enables specialization for D-dimensional lattice.
 */
template<class T, unsigned D>
struct makelattice<CrystalLattice<T,D> >
{
  inline static void apply(CrystalLattice<T,D>& , std::vector<std::string>& argv)
  {
  }
};

/** Specialization of makelattice<CL> for CrystalLattice<T,1>*/
template<class T>
struct makelattice<CrystalLattice<T,1> >
{

  inline static
  void
  apply(CrystalLattice<T,1>& lat, std::vector<std::string>& argv)
  {
    int i=0;
    int argc = argv.size();
    while(i<argc)
    {
      if(argv[i] == "a0")
      {
        lat.R(0,0) = std::atof(argv[++i].c_str());
      }
      i++;
    }
    lat.reset();
  }
};

/** Specialization of makelattice<CL> for CrystalLattice<T,2>*/
template<class T>
struct makelattice<CrystalLattice<T,2> >
{

  inline static
  void
  apply(CrystalLattice<T,2>& lat, std::vector<std::string>& argv)
  {
    T a0 = 1.0e0;
    int i=0;
    int argc = argv.size();
    while(i<argc)
    {
      if(argv[i] == "cubic")
      {
        a0 = std::atof(argv[++i].c_str());
        lat.R.diagonal(1.0);
      }
      else
        if(argv[i] == "orthorombic")
        {
          lat.R = 0.0e0;
          lat.R(0,0) = std::atof(argv[++i].c_str());
          lat.R(1,1) = std::atof(argv[++i].c_str());
        }
        else
          if(argv[i] == "general")
          {
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
struct makelattice<CrystalLattice<T,3> >
{

  /*! \fn makelattic<CrystalLattice<T,3> >
   *  ::apply(CrystalLattice<T,3>& lattice, std::vector<std::string>& argv)
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
  apply(CrystalLattice<T,3>& lat, std::vector<std::string>& argv)
  {
    T a0 = 1.0e0;
    int i=0;
    int argc = argv.size();
    while(i<argc)
    {
      if(argv[i] == "cubic")
      {
        a0 = std::atof(argv[++i].c_str());
        lat.R.diagonal(1.0);
      }
      else
        if(argv[i] == "orthorombic")
        {
          lat.R = 0.0e0;
          lat.R(0,0) = std::atof(argv[++i].c_str());
          lat.R(1,1) = std::atof(argv[++i].c_str());
          lat.R(2,2) = std::atof(argv[++i].c_str());
        }
        else
          if(argv[i] == "fcc")
          {
            a0 = std::atof(argv[++i].c_str());
            lat.R(0,0) = 0.0;
            lat.R(0,1) = 0.5;
            lat.R(0,2) = 0.5;
            lat.R(1,0) = 0.5;
            lat.R(1,1) = 0.0;
            lat.R(1,2) = 0.5;
            lat.R(2,0) = 0.5;
            lat.R(2,1) = 0.5;
            lat.R(2,2) = 0.0;
          }
          else
            if(argv[i] == "bcc")
            {
              a0 = std::atof(argv[++i].c_str());
              lat.R(0,0) = -0.5;
              lat.R(0,1) =  0.5;
              lat.R(0,2) =  0.5;
              lat.R(1,0) =  0.5;
              lat.R(1,1) = -0.5;
              lat.R(1,2) =  0.5;
              lat.R(2,0) =  0.5;
              lat.R(2,1) =  0.5;
              lat.R(2,2) = -0.5;
            }
            else
              if(argv[i] == "hcp")
              {
                a0 = std::atof(argv[++i].c_str());
                double covera = std::sqrt(8.0/3.0);
                if(argc-i > 1)
                  covera = std::atof(argv[++i].c_str());
                lat.R(0,0) = 0.5*a0;
                lat.R(0,1) = -std::sqrt(3.0)*0.5*a0;
                lat.R(0,2) = 0.0;
                lat.R(1,0) = 0.5*a0;
                lat.R(1,1) =  std::sqrt(3.0)*0.5*a0;
                lat.R(1,2) = 0.0;
                lat.R(2,0) = 0.0;
                lat.R(2,1) =   0.0;
                lat.R(2,2) = covera*a0;
                a0 = 1.0e0;
              }
      i++;
    }
    lat.R *= a0;
    lat.reset();
  }
};
}
