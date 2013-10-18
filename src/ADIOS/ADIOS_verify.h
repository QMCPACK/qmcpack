
#ifndef ADIOS_ADIOS_VERIFY_H
#define ADIOS_ADIOS_VERIFY_H

#if defined(ADIOS_VERIFY) && defined(HAVE_ADIOS)

#include <string>
#include <Configuration.h>
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/BranchIO.h"
#include <stdint.h>
#include "Utilities/RandomGenerator.h"

namespace IO_VERIFY
{
typedef qmcplusplus::SimpleFixedNodeBranch::RealType RealType;
typedef qmcplusplus::RandomGenerator_t::uint_type uint_type;
void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, unsigned long origin);
void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, RealType* origin);
void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, int* origin);
void adios_checkpoint_verify_intarray_variables(ADIOS_FILE* fp, const char* name, int* origin);
void adios_checkpoint_verify_random_variables(ADIOS_FILE* fp, const char* name, uint_type* origin);
void adios_checkpoint_verify_local_variables(ADIOS_FILE* fp, const char* name, OHMMS_PRECISION* origin);
void adios_trace_verify_local_variables(ADIOS_FILE* fp, const char* name, double* origin);
};

#endif
#endif

