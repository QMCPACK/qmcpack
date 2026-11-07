//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "CommandLineOptions.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "config.h"
#include "Message/UniformCommunicateError.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/qmc_common.h"
#include "qmcpack_version.h"

namespace qmcplusplus
{

void CommandLineOptions::parse(int argc, char** argv)
{
  *this = CommandLineOptions{};

  auto is_xml = [](const std::string& str) {
    return str.size() >= 4 && str.compare(str.size() - 4, 4, ".xml") == 0;
  };

  bool save_wfs = false;
  bool vacuum   = false;
  std::vector<std::string> unrecognized_options;

  std::vector<std::string> fgroup_names_cmd;
  std::vector<std::string> fgroup_names_txt;

  int i = 1;
  while (i < argc)
  {
    std::string c(argv[i]);
    if (c[0] == '-')
    {
      if (c == "-debug" || c == "--debug")
        debug_ = true;
      else if (c == "-enable-timers" || c == "--enable-timers")
      {
        enable_timers_ = true;
      }
      else if (c.rfind("-enable-timers=", 0) == 0 || c.rfind("--enable-timers=", 0) == 0)
      {
        enable_timers_ = true;
        timer_level_   = c.substr(c.find('=') + 1);
      }
      else if (c == "-verbosity" || c == "--verbosity" || c.rfind("-verbosity=", 0) == 0 ||
               c.rfind("--verbosity=", 0) == 0)
      {
        const auto pos = c.find("=");
        std::string verbose_level;
        if (pos != std::string::npos)
          verbose_level = c.substr(pos + 1);
        else if (i + 1 < argc && argv[i + 1][0] != '-')
          verbose_level = argv[++i];
        else
          throw UniformCommunicateError("Missing argument for --verbosity");

        if (verbose_level == "low")
          verbosity_ = Verbosity::LOW;
        else if (verbose_level == "high")
          verbosity_ = Verbosity::HIGH;
        else if (verbose_level == "debug")
          verbosity_ = Verbosity::DEBUG;
        else
          throw UniformCommunicateError("Unknown verbosity level: " + verbose_level);
      }
      else if (c == "--vnode-ranks" || c == "-vnode-ranks" || c.rfind("--vnode-ranks=", 0) == 0 ||
               c.rfind("-vnode-ranks=", 0) == 0)
      {
        const auto pos = c.find("=");
        try
        {
          if (pos != std::string::npos)
            vnode_ranks_ = std::stoi(c.substr(pos + 1));
          else if (i + 1 < argc)
            vnode_ranks_ = std::stoi(argv[++i]);
          else
            throw UniformCommunicateError("Missing argument for --vnode-ranks");
        }
        catch (const std::invalid_argument&)
        {
          throw UniformCommunicateError("Invalid argument for --vnode-ranks");
        }
        catch (const std::out_of_range&)
        {
          throw UniformCommunicateError("Argument out of range for --vnode-ranks");
        }
        if (*vnode_ranks_ <= 0)
          throw UniformCommunicateError("--vnode-ranks must be a positive integer.");
      }
      else if (c == "--dryrun")
      {
        dryrun_ = true;
      }
      else if (c == "--noprint")
      {
        noprint_ = true;
      }
      else if (c == "--save_wfs" || c.rfind("--save_wfs=", 0) == 0)
      {
        save_wfs = true;
      }
      else if (c == "--help" || c == "-h")
      {
        help_ = true;
      }
      else if (c == "--version")
      {
        version_ = true;
      }
      else if (c == "--vacuum" || c.rfind("--vacuum=", 0) == 0)
      {
        vacuum = true;
      }
      else
      {
        unrecognized_options.push_back(c);
      }
    }
    else
    {
      if (is_xml(c))
        fgroup_names_cmd.push_back(argv[i]);
      else
      {
        std::ifstream fin(argv[i], std::ifstream::in);
        bool valid = !fin.fail();
        while (valid)
        {
          std::vector<std::string> words;
          getwords(words, fin);
          if (words.size())
          {
            if (is_xml(words[0]))
              fgroup_names_txt.push_back(words[0]);
          }
          else
            valid = false;
        }
      }
    }
    ++i;
  }

  inputs_ = std::move(fgroup_names_txt);
  inputs_.insert(inputs_.end(), fgroup_names_cmd.begin(), fgroup_names_cmd.end());

  if (!should_normal_stop())
  {
    std::string error_msg;
    for (const auto& opt : unrecognized_options)
      error_msg += "ERROR: command line option " + opt + " is unrecognized.\n";
    if (save_wfs)
      error_msg +=
          "ERROR: command line option --save_wfs has been removed. Use save_coefs input tag as described in the manual.\n";
    if (vacuum)
      error_msg +=
          "ERROR: command line option --vacuum has been removed. Use vacuum input tag as described in the manual.\n";
    if (inputs_.empty())
      error_msg += "No valid input file is given.\n"
                   "Usage: qmcpack [options] <input-files.xml>\n"
                   "Ensemble runs may be initialized by specifying either multiple .xml files or text files "
                   "containing lists of .xml input files.\n";

    if (!error_msg.empty())
    {
      if (error_msg.back() == '\n')
        error_msg.pop_back();
      throw UniformCommunicateError(error_msg);
    }
  }
}

void CommandLineOptions::print_version_and_usage(std::ostream& os) const
{
  os << std::endl
     << "QMCPACK version " << QMCPACK_VERSION_MAJOR << "." << QMCPACK_VERSION_MINOR << "."
     << QMCPACK_VERSION_PATCH << " built on " << __DATE__ << std::endl;
  qmc_common.print_git_info_if_present(os);
  os << std::endl
     << "Usage: qmcpack [options] <input-files.xml> [input2.xml ...]\n"
     << "       qmcpack [options] <input-list.txt>\n"
     << "\n"
     << "Options:\n"
     << "  --help, -h               Print this help message and exit.\n"
     << "  --version                Print version information and exit.\n"
     << "  --dryrun                 Validate input files without performing the simulation.\n"
     << "                           All QMC and loop sections will be skipped.\n"
     << "  --noprint                Do not write extra pseudopotential or Jastrow .dat files.\n"
     << "  -debug, --debug          Enable detailed debug output reporting.\n"
     << "  --enable-timers[=level]  Set timer granularity threshold (none, coarse, medium, fine).\n"
     << "                           Default is 'coarse'.\n"
     << "  --verbosity=level        Set output verbosity (low, high, debug).\n"
     << "                           Default is 'low'.\n"
     << "  --vnode-ranks=N          Number of MPI ranks per virtual node for ensemble runs.\n"
     << "                           Must be a positive integer and a factor of the node rank count.\n"
     << std::endl;
}

} // namespace qmcplusplus
