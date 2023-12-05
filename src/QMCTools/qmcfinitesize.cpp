//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "Message/Communicate.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "Particle/ParticleSetPool.h"

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/QMCFiniteSize.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/SkParserASCII.h"
#include "QMCTools/QMCFiniteSize/SkParserScalarDat.h"
#include "QMCTools/QMCFiniteSize/SkParserHDF5.h"

#include "Numerics/OneDimGridBase.h"

//Purpose of this routine is to compute the finite size effects
//for a given simulation cell from post-processed QMC Data.  For
//the potential, this is done by splining the structure factor
//and performing the integral given in Holzmann et al., PRB 035126 (2016)
//Assuming you have provided a long-ranged jastrow in the input xml, this will also
//calculate the kinetic energy correction from Holzmann et al.
//
//Input:  Cell geometry.  Ion positions, cell geometries, etc, are taken from main.xml file.
//                        Code recognizes ESHDF5 declarations in main.xml.
//                        Will use kspace jastrow to compute kinetic correction if available.
//        S(k):  This is the electron-electron fluctuation structure factor.
//
//Returns: (E(N=infty)-E(N)) for the given simulation cell.

using namespace qmcplusplus;
using RealType = QMCTraits::RealType;
using PosType  = QMCTraits::PosType;
using Grid_t   = SkParserBase::Grid_t;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc, argv);
  Random.init(-1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right, std::ios::adjustfield);
  std::cout.precision(12);

  // Suppress HDF5 warning and error messages.
  qmcplusplus::hdf_error_suppression hide_hdf_errors;

  std::unique_ptr<SkParserBase> skparser(nullptr);

  bool show_usage = false;
  bool show_warn  = false;
  std::string skname;
  std::string skfile;
  /* For a successful execution of the code, atleast 3 arguments will need to be
   * provided along with the executable. Therefore, print usage information if
   * argc is less than 4.
   */
  if (argc < 4)
  {
    std::cout << "Insufficient number of arguments.\n\n";
    show_usage = true;
  }
  else
  {
    int iargc      = 2;
    bool skf_found = false;
    while (iargc + 1 < argc)
    {
      std::string a(argv[iargc]);
      std::string anxt(argv[iargc + 1]);
      if (a == "--ascii")
      {
        skparser = std::make_unique<SkParserASCII>();
        skfile   = anxt;
        if (skf_found)
          show_warn = true;
        skf_found = true;
      }
      else if (a == "--scalardat")
      {
        skparser = std::make_unique<SkParserScalarDat>();
        skfile   = anxt;
        if (skf_found)
          show_warn = true;
        skf_found = true;
      }
      else if (a == "--hdf5")
      {
        skparser = std::make_unique<SkParserHDF5>();
        skfile   = anxt;
        if (skf_found)
          show_warn = true;
        skf_found = true;
      }
      else if (a == "--skname")
      {
        skname = anxt;
      }
      else
      {
        std::cout << "Unrecognized flag '" << a << "'.\n\n";
        show_usage = true;
      }
      iargc += 2;
    }
  }

  if (show_usage)
  {
    std::cout << "Usage:  qmcfinitesize [main.xml] --[skformat] [SK_FILE] --[optional_arg] [option]\n";
    std::cout << "  [main.xml]\n";
    std::cout << "    input file to qmcpack corresponding that calculated S(k) data. ";
    std::cout << "  [skformat]\n";
    std::cout << "    ascii:      S(k) given in kx ky kz sk sk_err format.  Header necessary.\n";
    std::cout << "    scalardat:  File containing skall elements with energy.pl output format.\n";
    std::cout << "    hdf5:       stat.h5 file containing skall data.\n";
    std::cout << "  [SK_FILE]\n";
    std::cout << "    filename containing the S(k) data\n";
    std::cout << "  [optional_args]\n";
    std::cout << "    --skname:     in the stat.h5, the S(k) group name. Set to \"name\" defined in qmcpack input file "
                 "where SkAll estimator was specified\n";
    std::cout << "  [option]\n";
    std::cout << "    option corresonding to the optional_arg\n";
    std::cout << "---------------------------\n";
    std::cout << "Examples:\n";
    std::cout << "  qmcfinitesize qmc.in.xml --hdf5 qmc.g000.s000.stat.h5 --skname SkAll\n";
    std::cout << "  qmcfinitesize qmc.in.xml --ascii processed_sk.dat\n";
    return 0;
  }

  if (show_warn)
    std::cout << "WARNING:  multiple skformats were provided. All but the last will be ignored.\n\n";

  if (!skname.empty())
    skparser->setName(skname);
  skparser->parse(skfile);

  if (skparser == NULL)
  {
    APP_ABORT("qmcfinitesize:  skparser failed to initialize");
  }

  QMCFiniteSize qmcfs(skparser.get());
  qmcfs.parse(std::string(argv[1]));
  qmcfs.validateXML();
  qmcfs.execute();

  // Jobs done. Clean up.
  OHMMS::Controller->finalize();
  return 0;
}
