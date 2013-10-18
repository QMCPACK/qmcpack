adios_groupsize = 4 \
                  + 4 \
                  + 4 \
                  + 8 * (walker_num) * (particle_num) * (walker_dim_num);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "walker_num", &walker_num);
adios_write (adios_handle, "particle_num", &particle_num);
adios_write (adios_handle, "walker_dim_num", &walker_dim_num);
adios_write (adios_handle, "walkers", walkers);
