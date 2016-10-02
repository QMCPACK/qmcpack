//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory 
//
// File created by: Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


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

