//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_DATASPORKANALYSIS_H
#define QMCPLUSPLUS_DATASPORKANALYSIS_H
#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct DataSporkAnalysis
{

  po::variables_map vm;

  int output_precision;
  int FirstIndex;
  int LastIndex;
  std::string merged_file;;
  //keywords
  std::string include_dir;
  std::string input_file;
  std::string output_file;
  std::string xslt_file;

  std::string message_txt;
  std::string generator_file;
  std::string observable;
  std::string collectable;

  std::map<std::string, int> Observables;
  std::map<std::string, int> Collectables;

  DataSporkAnalysis():include_dir("include-dir"),input_file("input-file"),
    observable("observable"),
    collectable("collectable") {}

  int getOptions(int ac, char* av[]);

  void execute();
  void printOptions(std::ostream& os);
};
#endif
