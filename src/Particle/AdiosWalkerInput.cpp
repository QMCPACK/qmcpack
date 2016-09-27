//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifdef HAVE_ADIOS
#include "OhmmsData/AttributeSet.h"
#include "Configuration.h"
#include <cstdlib>
#include "Particle/AdiosWalkerInput.h"


namespace qmcplusplus
{

void check_adios_error()
{
  if (adios_errno < 0)
  {
    app_error() << adios_errmsg() << std::endl;
  }
}


AdiosWalkerInput::AdiosWalkerInput(MCWalkerConfiguration& w, Communicate* c):
  targetW(w),myComm(c)
{
  FileRoot = "";
}

AdiosWalkerInput::~AdiosWalkerInput()
{
}

string AdiosWalkerInput::getFileRoot()
{
  return FileRoot;
}

bool AdiosWalkerInput::put(std::vector<xmlNodePtr>& wset)
{
  //Get our mpi proc rank
  int rank = myComm->rank();
  //Check if their are any wsetwalker tags in the xml file
  if (wset.empty())
    return false;
  int walker_win = 0;
  //This will accumalate the walkers we will be reading from the file
  std::string bpFileName;
  //Number of walkers per process
  std::vector<int> nw(myComm->size(), 0);
  //iterate through the mcwalkerset tags
  for (int i = 0; i < wset.size(); i++)
  {
    //We have multiple tags corresponding to mcWalkerset
    xmlNodePtr curr = wset[i];
    //set up a parse of the xml node curr
    std::string froot, file, href;
    int nprocs;
    OhmmsAttributeSet pAttrib;
    //full file name could either be href or file
    pAttrib.add(href, "href");
    pAttrib.add(file, "file");
    //<fileroot>.config.bp
    pAttrib.add(froot, "fileroot");
    //Get the attributes
    pAttrib.put(curr);
    //Check to see if any file attributes were set
    if (froot.empty() && href.empty() && file.empty())
      return false;
    //We default to using froot if all of them were supplied
    //followed by file and then href
    if (!froot.empty())
    {
      bpFileName = froot + ".config.bp";
      app_log() << "Using froot: ignoring href and file tags"
                << std::endl;
    }
    else if (!file.empty())
    {
      bpFileName = file;
      app_log() << "Using file: ignoring href" << std::endl;
    }
    else if (!href.empty())
    {
      bpFileName = href;
      app_log() << "Using href tag" << std::endl;
    }
    else
      app_error() << "No file associated tag in mcwalkerset tag" << std::endl;
    app_log() << "Reading walker configurations from: " << bpFileName << std::endl;
    //Open the bp file
    ADIOS_FILE* adios_file_handle = adios_read_open_file(bpFileName.c_str(),
                                    ADIOS_READ_METHOD_BP,
                                    myComm->getMPI());
    //Did the bp file open successfully
    check_adios_error();
    //Inquire about the number of proccess
    ADIOS_VARINFO* var_info = adios_inq_var(adios_file_handle, "walkers");
    nprocs = *var_info->nblocks;
    app_log() << "Number of procs that wrote " << nprocs << std::endl;
    adios_free_varinfo(var_info);
    //read in the data
    read(nprocs, adios_file_handle, walker_win, nw);
  }
  FileRoot = bpFileName;
  //Set mcwalker data
  setMCWalker( nw);
  return true;
}

void AdiosWalkerInput::setMCWalker(std::vector<int> nw)
{
  //Now we generate the walker offsets
  int np = myComm->size();
  std::vector<int> nwoff(myComm->size() + 1, 0);
  for(int ip=0; ip<np; ++ip)
    nwoff[ip+1]=nwoff[ip]+nw[ip];
  app_log() << "Number of walkers " << nwoff[np] << std::endl;
  targetW.setGlobalNumWalkers(nwoff[np]);
  targetW.setWalkerOffsets(nwoff);
}

void AdiosWalkerInput::read(int nprocs,
                            ADIOS_FILE* file_handle,
                            int& walker_win,
                            std::vector<int>& nw)
{
  //iterate over the number of blocks in the adios file
  std::vector<int> walker_num(nprocs, 0);
  int total_walker_num = 0;
  ADIOS_SELECTION* sel;
  for (int i = 0; i < nprocs; i++)
  {
    //Find out how many walkers are in each writeblock from the last proc
    sel = adios_selection_writeblock(i);
    adios_schedule_read(file_handle, sel, "walker_num", 0, 1, &(walker_num[i]));
  }
  //Force reads to finish
  adios_perform_reads(file_handle, 1);
  check_adios_error();
  for (int i = 0; i < nprocs; i++)
    total_walker_num += walker_num[i];
  for (int j = 0; j < walker_num.size(); j++)
    app_log() << walker_num[j] << ", ";
  app_log() << std::endl;
  //The number of proccess to add 1 extra walker too
  int walker_wrap = total_walker_num % myComm->size();
  //Buffer for each block in the adios
  R_t block_buff((total_walker_num / nprocs + 1) * targetW.getParticleNum());
  //Calculate how many walkers each proccess should read and from which files
  // this proccess should read
  int current_adios_block = 0;
  int current_offset = 0;
  int next_offset;
  int next_adios_block;
  for (int i = 0; i < myComm->size(); i++)
  {
    int walkers_to_read;
    if ( (walker_win <= i && i < walker_win + walker_wrap) ||
         (i < (walker_win + walker_wrap) % myComm->size() &&
          (walker_win + walker_wrap) % myComm->size() <= walker_win))
      walkers_to_read = total_walker_num / myComm->size() + 1;
    else
      walkers_to_read = total_walker_num / myComm->size();
    app_log() << "walkers_to_read " << walkers_to_read << " proc=" << i << std::endl;
    //Keep track of how many walkers each proccess has
    nw[i] += walkers_to_read;
    while (walkers_to_read != 0)
    {
      for (int j = 0; j < walker_num.size(); j++)
        app_log() << walker_num[j] << ", ";
      app_log() << std::endl;
      int read_size;
      if (walker_num[current_adios_block] > walkers_to_read)
      {
        read_size = walkers_to_read;
        next_offset = read_size + current_offset;
        walker_num[current_adios_block] -= read_size;
        next_adios_block = current_adios_block;
      }
      else
      {
        read_size = walker_num[current_adios_block];
        next_offset = 0;
        walker_num[current_adios_block] = 0;
        next_adios_block = current_adios_block + 1;
      }
      walkers_to_read -= read_size;
      app_log() << "Read Size=" << read_size << std::endl;
      app_log() << "Offset=" << current_offset << std::endl;
      app_log() << "Adios Block=" << current_adios_block << std::endl;
      app_log() << "Next Block=" << next_adios_block << std::endl;
      app_log() << " Next Offset=" << next_offset << std::endl;
      if (i == myComm->rank())
      {
        //Select a region of the walker buffer
        sel = adios_selection_writeblock(current_adios_block);
        adios_schedule_read(file_handle, sel, "walkers", 0, 1, &(block_buff[0]));
        adios_selection_delete(sel);
        adios_perform_reads(file_handle, 1);
        append_walkers(block_buff, read_size, current_offset);
      }
      current_adios_block = next_adios_block;
      current_offset = next_offset;
    }
  }
  check_adios_error();
  //Move the window
  walker_win = (walker_win + walker_wrap) % myComm->size();
  //Append the walkers we read from the buffer to the walkers list
}

void AdiosWalkerInput::append_walkers(R_t& walker_buff, int read_size,
                                      int current_offset)
{
  int offset = current_offset * targetW.getParticleNum();
  int cw = targetW.getActiveWalkers();
  targetW.createWalkers(read_size);
  int pn = targetW.getParticleNum();
  //Now lets append all those walkers to the walker list
  for (int w = cw; w < (cw + read_size); w++) //through walkers
    for ( int i = 0; i < pn; i++)   //through R
    {
      targetW.WalkerList[w]->R[i] = walker_buff[offset + w * pn + i];
    }
}

}
#endif
