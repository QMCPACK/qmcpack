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
    
    


adios_groupsize = 8 \
                  + 4 * (4) \
                  + 8 * (4) \
                  + 8 * (4) \
                  + 8 * (4) \
                  + 4 * (8) \
                  + strlen(iparem_def) \
                  + 4 * (16) \
                  + strlen(vparem_def);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "branchmode", &branchmode);
adios_write (adios_handle, "energy", energy);
adios_write (adios_handle, "r2accepted", r2accepted);
adios_write (adios_handle, "r2proposed", r2proposed);
adios_write (adios_handle, "variance", variance);
adios_write (adios_handle, "iparam", iparam);
adios_write (adios_handle, "iparem_def", iparem_def);
adios_write (adios_handle, "vparam", vparam);
adios_write (adios_handle, "vparem_def", vparem_def);
