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

    Trace (string group_name, MPI_Comm comm, xmlNodePtr adios_options);
    ~Trace();
    int define_var(string path, int ndim, int *dims, string type);

    void open (string filename);
    void set_group_size (int nrows);
    void write (string varname, void *data);
    void close();

private:

    int64_t  group;         // a pointer to ADIOS' internal group
    uint64_t scalars_size;  // size of all adios auxiliary variables to be written in bytes
    uint64_t rowsize;       // size of one row of all trace variables in bytes
    string   groupname;     // name of the adios group defining the group of variables
    int64_t  f;             // file descriptor
    MPI_Comm t_comm;        // the participants in the parallel trace writing

    // options from QMC XML (under <traces><adios_options>...</>)
    string   method_name;   // ADIOS output method
    string   method_args;   // ADIOS output method's arguments
    std::map<std::string,std::string> transforms;
    void process_options (xmlNodePtr adios_options);

};

} //namespace

#endif // HAVE_ADIOS
#endif // ADIOS_TRACE_H


















