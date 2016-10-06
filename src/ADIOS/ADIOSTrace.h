//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Norbert Podhorszki, pnorbert@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Norbert Podhorszki, pnorbert@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef ADIOS_TRACE_H
#define ADIOS_TRACE_H

#ifdef HAVE_ADIOS
#include <string>
#include <fstream>
#include <sstream>
#include <Configuration.h>
//#include <OhmmsData/OhmmsElementBase.h>
#include <OhmmsData/AttributeSet.h>
#include <stdint.h>
#include <strings.h>
#include <vector>
#include <map>

#include "adios.h"
//#include "adios_read.h"


namespace ADIOS
{

/** ADIOS Trace helper class to define and write the tracing data with ADIOS.
 */
class Trace
{
public:

    Trace ( std::string group_name, MPI_Comm comm, xmlNodePtr adios_options);
    ~Trace();
    int define_var( std::string path, int ndim, int *dims, std::string type);

    void open ( std::string filename);
    void set_group_size (int nrows);
    void write ( std::string varname, void *data);
    void close();

private:

    int64_t  group;         // a pointer to ADIOS' internal group
    uint64_t scalars_size;  // size of all adios auxiliary variables to be written in bytes
    uint64_t rowsize;       // size of one row of all trace variables in bytes
    std::string   groupname;     // name of the adios group defining the group of variables
    int64_t  f;             // file descriptor
    MPI_Comm t_comm;        // the participants in the parallel trace writing

    // options from QMC XML (under <traces><adios_options>...</>)
    std::string   method_name;   // ADIOS output method
    std::string   method_args;   // ADIOS output method's arguments
    std::map<std::string,std::string> transforms;
    void process_options (xmlNodePtr adios_options);

};

} //namespace

#endif // HAVE_ADIOS
#endif // ADIOS_TRACE_H


















