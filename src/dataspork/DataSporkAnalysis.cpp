// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//Modify for DataSporkAnalysis
#include "DataSporkAnalysis.h"
#include "ScalarDataSetManager.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/tokenizer.hpp>
#include <Platforms/sysutil.h>
//#include <boost/date_time/posix_time/posix_time.hpp>

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}


int DataSporkAnalysis::getOptions(int ac, char* av[])
{
  try
  {
    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("xsl",
     po::value<std::string>(&xslt_file)->default_value("no"),
     "xslt file to generate html using xsltproc")
    ("merge",
     po::value<std::string>(&merged_file)->default_value("no"),
     "merge files for one output")
    ("output",
     po::value<std::string>(&output_file)->default_value("generic.xml"),
     "write to a file file")
    ("first",
     po::value<int>(&FirstIndex)->default_value(0),
     "First index of the data")
    ("last",
     po::value<int>(&LastIndex)->default_value(-1),
     "Last index of the data. Default uses the last valid index.")
    ("precision",
     po::value<int>(&output_precision)->default_value(6),
     "precision of the output")
    ("qmc-input",
     po::value<std::string>(&generator_file)->default_value("missing"),
     "The main driver for qmc run")
    ("message,m",
     po::value<std::string>(&message_txt)->default_value("dataspork analysis"),
     "The main driver for qmc run")
    ;
    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
    ("observable,o",
     po::value< std::vector<std::string> >()->composing(),
     "observables")
    ("collectable,c",
     po::value< std::vector<std::string> >()->composing(),
     "collectables")
    ;
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", po::value< std::vector<std::string> >(), "input file")
    ;
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);
    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);
    po::options_description visible("Allowed options");
    visible.add(generic).add(config);
    po::positional_options_description p;
    p.add("input-file", -1);
    store(po::command_line_parser(ac, av).
          options(cmdline_options).positional(p).run(), vm);
    std::ifstream ifs("dataspork.cfg");
    if(ifs.fail())
    {
      std::stringstream qin;
      qin << "observable = LocalEnergy LocalPotential Kinetic ElecElec Coulomb\n";
      qin << "collectable =  Variance Weight NumOfWalkers TrialEnergy BlockCPU AcceptRatio\n";
      store(parse_config_file(qin, config_file_options), vm);
    }
    else
    {
      //use the configuraton file
      store(parse_config_file(ifs, config_file_options), vm);
    }
    notify(vm);
    if (vm.count("help"))
    {
      std::cout << "Usage: datasporkpp input-file+ [options]\n\n";
      std::cout << visible << "\n";
      return 0;
    }
    if (vm.count("version"))
    {
      std::cout << "datasporkpp version 1.0\n";
      return 0;
    }
    //if (vm.count("include-dir"))
    //{
    //    std::cout << "Include directories are: "
    //         << vm["include-dir"].as< std::vector<std::string> >() << "\n";
    //}
    if (vm.count(observable))
    {
      const std::vector<std::string>& vm_o(vm[observable].as< std::vector<std::string> >());
      std::vector<std::string>::const_iterator it(vm_o.begin()), it_end(vm_o.end());
      int curObservable=Observables.size();
      while(it != it_end)
      {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep(" ");
        tokenizer tok(*it, sep);
        for(tokenizer::iterator tit=tok.begin(); tit!=tok.end(); tit++)
        {
          Observables[*tit]=curObservable++;
        }
        ++it;
      }
    }
    if (vm.count(collectable))
    {
      const std::vector<std::string>& vm_o(vm[collectable].as< std::vector<std::string> >());
      std::vector<std::string>::const_iterator it(vm_o.begin()), it_end(vm_o.end());
      int curCollectable=Collectables.size();
      while(it != it_end)
      {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep(" ");
        tokenizer tok(*it, sep);
        for(tokenizer::iterator tit=tok.begin(); tit!=tok.end(); tit++)
        {
          Collectables[*tit]=curCollectable++;
        }
        ++it;
      }
    }
  }
  catch(exception& e)
  {
    std::cout << e.what() << "\n";
    return 1;
  }
  return 0;
}

void DataSporkAnalysis::printOptions(std::ostream& os)
{
  os << "    <options>\n";
  if (vm.count("input-file"))
  {
    const std::vector<std::string>& ins(vm["input-file"].as< std::vector<std::string> >());
    os << "    <input_file>\n" << ins << "\n</input_file>\n";
  }
  if (vm.count(observable))
  {
    const std::vector<std::string>& ins(vm["observable"].as< std::vector<std::string> >());
    os << "    <observable>\n" << ins << "\n</observable>\n";
  }
  os << "  </options>" << std::endl;
}

void DataSporkAnalysis::execute()
{
  ScalarDataSetManager fh;
  fh.registerObservables(Observables);
  fh.registerCollectables(Collectables);
  std::ofstream out(output_file.c_str());
  out.setf(std::ios::scientific, std::ios::floatfield);
  out.precision(output_precision);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(output_precision);
  if (vm.count("input-file"))
  {
    if(merged_file != "no")
    {
      fh.addDataSet(vm["input-file"].as< std::vector<std::string> >(), merged_file, FirstIndex, LastIndex);
    }
    else
    {
      const std::vector<std::string>& ins(vm["input-file"].as< std::vector<std::string> >());
      for(int i=0; i<ins.size(); i++)
      {
        fh.addDataSet(ins[i], FirstIndex, LastIndex);
      }
    }
  }
  /* create the local_time to report*/
  //using namespace boost::posix_time;
  //ptime now=second_clock::local_time();
  //time_facet* facet(new time_facet("%Y-%m-%d %T"));
  //out.imbue(std::locale(out.getloc(), facet));
  std::string now=getDateAndTime("%Y-%m-%d %T");
  out << "<?xml version=\"1.0\"?>\n"
      << "<?xml-stylesheet type=\"text/xsl\" href=\"dataspork.xsl\"?>\n"
      << "<dataspork>" << std::endl;
  out
      << "  <generated by=\"datasporkpp\" version=\"1.0\">\n"
      << "    <when>" << now << "</when>\n";
  printOptions(out);
  out << "  </generated>\n";
  out << "  <original_input>"<<generator_file<< "\n</original_input>\n";
  out << "  <comment>"<<message_txt << "\n</comment>\n";
  fh.write(out);
  out <<"</dataspork>" << std::endl;
  out.close();
  if(xslt_file != "no")
  {
    std::cout << "Generating " << output_file << ".html using " << xslt_file << std::endl;
    char cmds[128];
    sprintf(cmds,"xsltproc -o %s.html %s %s",
            output_file.c_str(),xslt_file.c_str(),output_file.c_str());
    system(cmds);
  }
}
