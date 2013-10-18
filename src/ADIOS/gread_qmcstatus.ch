s = adios_selection_writeblock (rank);
adios_schedule_read (fp, s, "energy", 1, 1, energy);
adios_schedule_read (fp, s, "r2accepted", 1, 1, r2accepted);
adios_schedule_read (fp, s, "r2proposed", 1, 1, r2proposed);
adios_schedule_read (fp, s, "variance", 1, 1, variance);
adios_schedule_read (fp, s, "iparam", 1, 1, iparam);
adios_schedule_read (fp, s, "vparam", 1, 1, vparam);
adios_perform_reads (fp, 1);
adios_selection_delete (s);
