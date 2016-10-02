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
    
    


adios_groupsize = 4 \
                  + 4 \
                  + 4 \
                  + 8 * (walker_num) * (particle_num) * (walker_dim_num);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "walker_num", &walker_num);
adios_write (adios_handle, "particle_num", &particle_num);
adios_write (adios_handle, "walker_dim_num", &walker_dim_num);
adios_write (adios_handle, "walkers", walkers);
