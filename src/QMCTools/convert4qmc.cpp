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


#include "QMCTools/GaussianFCHKParser.h"
#include "QMCTools/GamesAsciiParser.h"
#include "QMCTools/LCAOHDFParser.h"
#include "QMCTools/DiracParser.h"
#include "QMCTools/RMGParser.h"
#include "Message/Communicate.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "Platforms/Host/OutputManager.h"
#include <sstream>

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller->initialize(env);
#endif
  if (argc < 2)
  {
    std::cout << "Usage: convert [-gaussian|-gamess|-orbitals|-dirac|-rmg] filename " << std::endl;
    std::cout << "[-nojastrow -hdf5 -prefix title -addCusp -production -NbImages NimageX NimageY NimageZ]" << std::endl;
    std::cout << "[-psi_tag psi0 -ion_tag ion0 -gridtype log|log0|linear -first ri -last rf]" << std::endl;
    std::cout << "[-size npts -multidet multidet.h5 -ci file.out -threshold cimin -TargetState state_number "
                 "-NaturalOrbitals NumToRead -optDetCoeffs]"
              << std::endl;
    std::cout << "Defaults : -gridtype log -first 1e-6 -last 100 -size 1001 -ci required -threshold 0.01 -TargetState "
                 "0 -prefix sample"
              << std::endl;
    std::cout << "When the input format is missing, the  extension of filename is used to determine the format "
              << std::endl;
    std::cout << " *.Fchk -> gaussian; *.out -> gamess; *.h5 -> HDF5" << std::endl;
    return 0;
  }
  else
  {
    try
    {
      if (OHMMS::Controller->rank() != 0)
      {
        outputManager.shutOff();
      }
      Random.init(0, 1, -1);
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout.setf(std::ios::right, std::ios::adjustfield);
      std::cout.precision(12);
      QMCGaussianParserBase::init();
      std::unique_ptr<QMCGaussianParserBase> parser;
      int iargc = 0;
      std::string in_file(argv[1]);


      std::string punch_file;
      std::string psi_tag("psi0");
      std::string ion_tag("ion0");
      std::string jastrow("j");
      std::string prefix;

      int TargetState = 0;
      bool addJastrow = true;
      bool usehdf5    = false;
      bool h5         = false;
      bool useprefix  = false;
      bool debug      = false;
      bool prod       = false;
      bool ci = false, zeroCI = false, orderByExcitation = false, addCusp = false, multidet = false,
           optDetCoeffs = false;
      double thres      = 1e-20;
      int readNO        = 0; // if > 0, read Natural Orbitals from gamess output
      int readGuess     = 0; // if > 0, read Initial Guess from gamess output
      std::vector<int> Image;
      while (iargc < argc)
      {
        std::string a(argv[iargc]);
        if (a == "-gaussian")
        {
          parser  = std::make_unique<GaussianFCHKParser>(argc, argv);
          in_file = argv[++iargc];
        }
        else if (a == "-gamess")
        {
          parser  = std::make_unique<GamesAsciiParser>(argc, argv);
          in_file = argv[++iargc];
        }
        else if (a == "-dirac")
        {
          parser  = std::make_unique<DiracParser>(argc, argv);
          in_file = argv[++iargc];
          usehdf5 = true;
        }
        else if (a == "-orbitals")
        {
          parser  = std::make_unique<LCAOHDFParser>(argc, argv);
          h5      = true;
          in_file = argv[++iargc];
        }
        else if (a == "-rmg")
        {
          parser  = std::make_unique<RMGParser>(argc, argv);
          h5      = true;
          in_file = argv[++iargc];
        }
        else if (a == "-hdf5")
        {
          usehdf5 = true;
        }
        else if (a == "-psi_tag")
        {
          psi_tag = argv[++iargc];
        }
        else if (a == "-production")
        {
          prod = true;
        }
        else if (a == "-ion_tag")
        {
          ion_tag = argv[++iargc];
        }
        else if (a == "-prefix")
        {
          prefix    = argv[++iargc];
          useprefix = true;
        }
        else if (a == "-ci")
        {
          ci         = true;
          punch_file = argv[++iargc];
        }
        else if (a == "-multidet")
        {
          multidet   = true;
          punch_file = argv[++iargc];
        }
        else if (a == "-NbImages")
        {
          int temp;
          temp = atoi(argv[++iargc]);
          temp += 1 - temp % 2;
          Image.push_back(temp);
          temp = atoi(argv[++iargc]);
          temp += 1 - temp % 2;
          Image.push_back(temp);
          temp = atoi(argv[++iargc]);
          temp += 1 - temp % 2;
          Image.push_back(temp);
        }
        else if (a == "-addCusp")
        {
          addCusp = true;
        }
        else if (a == "-threshold")
        {
          thres = atof(argv[++iargc]);
        }
        else if (a == "-optDetCoeffs")
        {
          optDetCoeffs = true;
        }
        else if (a == "-TargetState")
        {
          TargetState = atoi(argv[++iargc]);
        }
        else if (a == "-NaturalOrbitals")
        {
          readNO = atoi(argv[++iargc]);
        }
        else if (a == "-readInitialGuess")
        {
          readGuess = atoi(argv[++iargc]);
        }
        else if (a == "-zeroCI")
        {
          zeroCI = true;
        }
        else if (a == "-orderByExcitation")
        {
          orderByExcitation = true;
        }
        else if (a == "-cutoff")
        {
          orderByExcitation = true;
        }
        else if (a == "-debug")
        {
          debug = true;
        }
        else if (a == "-nojastrow")
        {
          addJastrow = false;
          jastrow    = "noj";
        }
        ++iargc;
      }
      if (readNO > 0 && readGuess > 0)
      {
        std::cerr << "Can only use one of: -NaturalOrbitals or -readInitialGuess. \n";
        abort();
      }
      //Failed to create a parser. Try with the extension
      std::string ext = getExtension(in_file);
      if (parser == 0)
      {
        if (ext == "Fchk")
        {
          WARNMSG("Creating GaussianFCHKParser")
          parser = std::make_unique<GaussianFCHKParser>(argc, argv);
        }
        else if (ext == "h5")
        {
          WARNMSG("Creating LCAOHDFParser")
          parser = std::make_unique<LCAOHDFParser>(argc, argv);
          h5      = true;
        }
        else if (ext == "out")
        {
          WARNMSG("Creating GamesAsciiParser")
          parser = std::make_unique<GamesAsciiParser>(argc, argv);
        }
        else
        {
          std::cerr << "Unknown extension: " << ext << std::endl;
          exit(1);
        }
      }
      if (useprefix != true)
      {
        prefix = in_file;
        std::string delimiter;
        if (ext == "h5")
          delimiter = ".h5";
        else
          delimiter = ".out";
        int pos = 0;
        std::string token;
        pos   = prefix.find(delimiter);
        token = prefix.substr(0, pos);
        prefix.erase(0, pos + delimiter.length());
        prefix = token;
      }
      std::cout << "Using " << prefix << " to name output files" << std::endl;

      parser->Title       = prefix;
      parser->debug       = debug;
      parser->DoCusp      = addCusp;
      parser->UseHDF5     = usehdf5;
      parser->singledetH5 = h5;
      if (h5)
      {
        parser->UseHDF5 = false;
        parser->h5file  = in_file;
      }
      if (usehdf5)
        parser->h5file = parser->Title + ".orbs.h5";
      parser->IonSystem.setName(ion_tag);
      if (debug)
      {
        parser->UseHDF5 = false;
        parser->h5file  = "";
      }
      parser->multideterminant = false;
      if (ci)
        parser->multideterminant = ci;
      if (multidet)
      {
        parser->multideterminant = multidet;
        parser->multidetH5       = multidet;
      }
      parser->multih5file       = punch_file;
      parser->production        = prod;
      parser->ci_threshold      = thres;
      parser->optDetCoeffs      = optDetCoeffs;
      parser->target_state      = TargetState;
      parser->readNO            = readNO;
      parser->orderByExcitation = orderByExcitation;
      parser->zeroCI            = zeroCI;
      parser->readGuess         = readGuess;
      parser->outputFile        = punch_file;
      parser->Image             = Image;
      parser->parse(in_file);
      if (prod)
      {
        parser->addJastrow = addJastrow;
        parser->WFS_name   = jastrow;
        if (parser->PBC)
        {
          std::cout << "Generating Inputs for Supertwist  with coordinates:" << parser->STwist_Coord[0] << "  "
                    << parser->STwist_Coord[1] << "  " << parser->STwist_Coord[2] << std::endl;
          parser->dumpPBC(psi_tag, ion_tag);
        }
        else
          parser->dump(psi_tag, ion_tag);
        parser->dumpStdInputProd(psi_tag, ion_tag);
      }
      else
      {
        parser->addJastrow = addJastrow;
        parser->WFS_name   = jastrow;
        if (parser->PBC)
        {
          std::cout << "Generating Inputs for Supertwist  with coordinates:" << parser->STwist_Coord[0] << "  "
                    << parser->STwist_Coord[1] << "  " << parser->STwist_Coord[2] << std::endl;
          parser->dumpPBC(psi_tag, ion_tag);
        }
        else
          parser->dump(psi_tag, ion_tag);
        parser->dumpStdInput(psi_tag, ion_tag);
      }
    }
    catch (const std::exception& e)
    {
      app_error() << e.what() << std::endl;
      APP_ABORT("Unhandled Exception");
    }
  }
  OHMMS::Controller->finalize();
  return 0;
}
