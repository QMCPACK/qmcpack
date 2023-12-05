//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_WALKERIO_HPP
#define AFQMC_WALKERIO_HPP

#include <cassert>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include "type_traits/container_traits_multi.h"

#include "Configuration.h"
#include "AFQMC/config.h"

namespace qmcplusplus
{
namespace afqmc
{
template<class WalkerSet, typename = typename std::enable_if<(WalkerSet::contiguous_walker)>::type>
bool dumpSamplesHDF5(WalkerSet& wset, hdf_archive& dump, int nW_to_file)
{
  return true;
  APP_ABORT("Finish \n");
  return true;
  /*
  if(nW_to_file==0) return true;
  if(head) { 

    int nW = numWalkers();

    std::vector<int> from(nproc_heads);
    MPI_Allgather(&nW,1,MPI_INT,from.data(),1,MPI_INT,TG.TG_heads().impl_);

    int nWtot = std::accumulate(from.begin(),from.end(),int(0));
    int w0 = std::accumulate(from.begin(),from.begin()+rank_heads,int(0));

    if(nW_to_file<0) nW_to_file=nWtot;
    nW_to_file = std::min(nW_to_file,nWtot);

    // careful here, avoid sending extra information (e.g. B mats for back propg)
    int wlk_nterms = (9 + nrow*ncol);
    int wlk_sz = wlk_nterms*sizeof(ComplexType);
    int nwlk_per_block = std::min(std::max(1,hdf_block_size/wlk_sz),nW_to_file);
    int nblks = nW_to_file/nwlk_per_block + ((nW_to_file%nwlk_per_block>0)?1:0);
    std::vector<int> wlk_per_blk;
  
    if(rank_heads==0) {

      // check that restart data doesnot exist 
      counts.resize(nproc_heads);
      displ.resize(nproc_heads);
      wlk_per_blk.reserve(nblks);

      commBuff.reserve(wlk_nterms*nwlk_per_block);

      std::vector<int> Idata(6);
      Idata[0]=nW_to_file;
      Idata[1]=nblks;
      Idata[2]=wlk_nterms;
      Idata[3]=wlk_sz;
      Idata[4]=nrow;
      Idata[5]=ncol;

      dump.push("Walkers");
      dump.write(Idata,"dims");
      
    }

    // to extract a single walker from walker storage 
    MPI_Datatype rtype, stype;
    {
      MPI_Datatype stype_;
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &stype_);
      MPI_Type_commit(&stype_);
      MPI_Aint loc0, loc1;
      MPI_Get_address(walkers.values(), &loc0);
      MPI_Get_address(walkers.values()+walker_size, &loc1);
      MPI_Aint dis = loc1-loc0;
      MPI_Type_create_resized(stype_,0,dis,&stype); 
      MPI_Type_commit(&stype);
      MPI_Type_free(&stype_);

      // to deposit a single walker on contiguous memory
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &rtype);
      MPI_Type_commit(&rtype);
    }

    int nsent=0;
    // ready to send walkers to head in blocks
    for(int i=0, ndone=0; i<nblks; i++, ndone+=nwlk_per_block) {

      int nwlk_tot = std::min(nwlk_per_block,nW_to_file-ndone);  
      int nw_to_send=0; 
      if( w0+nsent >= ndone && w0+nsent < ndone + nwlk_tot) 
        nw_to_send = std::min(nW-nsent,(ndone+nwlk_tot)-(w0+nsent)); 
 
      if(rank_heads==0) {

        for(int p=0, nt=0; p<nproc_heads; p++) {
   
          int n_ = 0;
          int nn = nt + nW;   
          if( ndone+nwlk_tot > nt && ndone < nt+nW ) {
            if(ndone <= nt)
              n_ = std::min(nW,(ndone+nwlk_tot)-nt);    
            else    
              n_ = std::min(nt+nW-ndone,nwlk_tot);    
          }    

          counts[p]=n_;
          nt+=from[p];

        }    
        displ[0]=0;
        for(int p=1, nt=0; p<nproc_heads; p++) {
          nt += counts[p-1];
          displ[p]=nt;
        }

        commBuff.resize(wlk_nterms*nwlk_tot);
      }  

      // assumes current storage structure
      MPI_Gatherv( getSM(nsent), nw_to_send, stype,
                    commBuff.data(), counts.data(), displ.data(), rtype, 
                    0, TG.TG_heads().impl_);
      nsent += nw_to_send;
      
      if(rank_heads==0) {  
        dump.write(commBuff,std::string("walkers_")+std::to_string(i));
        wlk_per_blk.push_back(nwlk_tot);
      }  
  
    }

    if(rank_heads==0) {
      dump.write(wlk_per_blk,"wlk_per_blk");
      dump.pop();
    }

    MPI_Type_free(&rtype);
    MPI_Type_free(&stype);

  }

  TG.global_barrier();  

  return true;
*/
}

template<class WalkerSet, typename = typename std::enable_if<(WalkerSet::contiguous_walker)>::type>
bool restartFromHDF5(WalkerSet& wset,
                     int nW_per_tg,
                     std::string hdf_read_restart,
                     hdf_archive& read,
                     bool set_to_target)
{
  TaskGroup_& TG = wset.getTG();

  std::vector<int> Idata(7);
  if (read.is_parallel())
  {
    app_error() << " Error: hdf_archive can't be parallel in restartFromHDF5().\n";
    APP_ABORT("");
  }
  if (TG.TG_local().root())
  {
    if (!read.open(hdf_read_restart, H5F_ACC_RDONLY))
      APP_ABORT("Problem opening restart file");

    std::string path = "/Walkers/WalkerSet";

    if (!read.is_group(path))
    {
      app_error()
          << " ERROR: H5Group  could not find /Walkers/WalkerSet group in file. No restart data for walkers. \n";
      return false;
    }

    read.push("Walkers");
    read.push("WalkerSet");

    if (!read.readEntry(Idata, "dims"))
      return false;
  }
  TG.TG_local().broadcast_n(Idata.begin(), Idata.size());

  auto walker_type = wset.getWalkerType();

  int nWtot      = Idata[0];
  int wlk_nterms = Idata[2];
  int NMO        = Idata[4];
  int NAEA       = Idata[5];
  int NAEB       = Idata[6];
  if (wlk_nterms != wset.walkerSizeIO())
  {
    app_error() << " Inconsistent walker restart file: IO size, NMO, NAEA, NAEB, WalkerType:" << wset.walkerSizeIO()
                << " " << NMO << " " << NAEA << " " << NAEB << " " << wset.getWalkerType() << std::endl;
    APP_ABORT("");
  }

  // walker range belonging to this TG
  int nW0, nWN;
  if (set_to_target)
  {
    if (nWtot < nW_per_tg * TG.TG_heads().size())
      APP_ABORT(" Error: Not enough walkers in restart file.\n");
    nW0 = nW_per_tg * TG.TG_heads().rank();
    nWN = nW0 + nW_per_tg;
  }
  else
  {
    if (nWtot % TG.TG_heads().size() != 0)
      APP_ABORT(" Error: Number of walkers in restart file must be divisible by number of task groups.\n");
    nW0 = (nWtot / TG.TG_heads().size()) * TG.TG_heads().rank();
    nWN = nW0 + nWtot / TG.TG_heads().size();
  }
  int nw_local = nWN - nW0;
  { // to limit scope
    boost::multi::array<ComplexType, 2> PsiA, PsiB;
    int NMO2 = ((walker_type == NONCOLLINEAR) ? 2 * NMO : NMO);
    if (TG.TG_local().root())
    {
      PsiA.reextent({NMO2, NAEA});
      if (wset.getWalkerType() == COLLINEAR)
        PsiB.reextent({NMO, NAEB});
    }
    // PsiA/B only meaningful at root
    wset.resize(nw_local, PsiA, PsiB);
  }

  // only head of WalkerSet reads
  if (TG.TG_local().root())
  {
    std::vector<int> wlk_per_blk;
    read.read(wlk_per_blk, "wlk_per_blk");

    boost::multi::array<ComplexType, 2> Data;

    // loop through blocks and read when necessary
    int ni = 0, nread = 0, bi = 0;
    while (nread < nw_local)
    {
      if (ni + wlk_per_blk[bi] > nW0)
      {
        // determine block of walkers to read
        int w0  = std::max(0, nW0 - ni);
        int nw_ = std::min(ni + wlk_per_blk[bi], nWN) - std::max(ni, nW0);
        Data.reextent({nw_, wlk_nterms});
        hyperslab_proxy<boost::multi::array_ref<ComplexType, 2>, 2> hslab(Data,
                                                                          std::array<size_t, 2>{static_cast<size_t>(wlk_per_blk[bi]),
                                                                                             static_cast<size_t>(wlk_nterms)},
                                                                          std::array<size_t, 2>{static_cast<size_t>(nw_), static_cast<size_t>(wlk_nterms)},
                                                                          std::array<size_t, 2>{static_cast<size_t>(w0), 0});
        read.read(hslab, std::string("walkers_") + std::to_string(bi));
        for (int n = 0; n < nw_; n++, nread++)
          wset.copyFromIO(Data[n], nread);
      }
      ni += wlk_per_blk[bi++];
    }
  }
  TG.global_barrier();

  read.close();
  return true;
}

template<class WalkerSet, typename = typename std::enable_if<(WalkerSet::contiguous_walker)>::type>
bool dumpToHDF5(WalkerSet& wset, hdf_archive& dump)
{
  TaskGroup_& TG = wset.getTG();

  if (TG.TG_local().root())
  {
    int nW = wset.size();

    auto nw_per_tg = TG.TG_heads().all_gather_value(nW);

    int nWtot = std::accumulate(nw_per_tg.begin(), nw_per_tg.end(), int(0));
    int w0    = std::accumulate(nw_per_tg.begin(), nw_per_tg.begin() + TG.TG_heads().rank(), int(0));

    auto walker_type = wset.getWalkerType();

    // careful here, avoid sending extra information (e.g. B mats for back propg)
    int wlk_nterms = wset.walkerSizeIO();
    int wlk_sz     = wlk_nterms * sizeof(ComplexType);

#if defined(__ENABLE_PHDF5__)

    // parallel I/O
    APP_ABORT("Restarting with parallel HDF5 not implemented yet.\n");
    int nwlk_per_block = std::min(std::max(1, WALKER_HDF_BLOCK_SIZE / wlk_sz), nW);
    int nblks          = (nW - 1) / nwlk_per_block + 1;
    auto nblks_per_tg  = TG.TG_heads().all_gather_value(nblks);
    int nblkTot        = std::accumulate(nblks_per_tg.begin(), nblks_per_tg.end(), int(0));
    int blk0           = std::accumulate(nblks_per_tg.begin(), nblks_per_tg.begin() + TG.TG_heads().rank(), int(0));
    std::vector<int> wlk_per_blk(nblks);

    //    if(TG.TG_heads().root()) {

    // check that restart data doesnot exist
    std::string path = "/Walkers/WalkerSet";
    if (dump.is_group(path))
    {
      app_error() << " ERROR: H5Group /Walkers/WalkerSet already exists in restart file. This is a bug and should not "
                     "happen. Contact a developer.\n";
      return false;
    }

    int NMO, NAEA, NAEB = 0;
    { // to limit the scope
      auto w = wset[0];
      NMO    = (*w.SlaterMatrix(Alpha)).size(0);
      NAEA   = (*w.SlaterMatrix(Alpha)).size(1);
      if (walker_type == COLLINEAR)
        NAEB = (*w.SlaterMatrix(Beta)).size(1) if (walker_type == NONCOLLINEAR) NMO /= 2;
    }

    std::vector<int> Idata(7);
    Idata[0] = nWtot;
    Idata[1] = nblkTot;
    Idata[2] = wlk_nterms;
    Idata[3] = wlk_sz;
    Idata[4] = NMO;
    Idata[5] = NAEA;
    Idata[6] = NAEB;

    dump.push("Walkers");
    dump.push("WalkerSet");
    dump.write(Idata, "dims");

    //    }
    APP_ABORT("FINISH.\n");

    // loop through blocks and use double hyperslabs

#else

    // communicate to root
    int nwlk_per_block = std::min(std::max(1, WALKER_HDF_BLOCK_SIZE / wlk_sz), nWtot);
    int nblks          = (nWtot - 1) / nwlk_per_block + 1;
    std::vector<int> wlk_per_blk;

    boost::multi::array<ComplexType, 2> RecvBuff;
    boost::multi::array<int, 1> counts, displ;

    if (TG.TG_heads().root())
    {
      // check that restart data doesnot exist
      std::string path = "/Walkers/WalkerSet";
      if (dump.is_group(path))
      {
        app_error() << " ERROR: H5Group /Walkers/WalkerSet already exists in restart file. This is a bug and should "
                       "not happen. Contact a developer.\n";
        return false;
      }

      counts.reextent({TG.TG_heads().size()});
      displ.reextent({TG.TG_heads().size()});
      wlk_per_blk.reserve(nblks);

      int NMO, NAEA, NAEB = 0;
      { // to limit the scope
        auto w = wset[0];
        NMO    = std::get<0>((*w.SlaterMatrix(Alpha)).sizes());
        NAEA   = std::get<1>((*w.SlaterMatrix(Alpha)).sizes());
        if (walker_type == COLLINEAR)
          NAEB = std::get<1>((*w.SlaterMatrix(Beta)).sizes());
        if (walker_type == NONCOLLINEAR)
          NMO /= 2;
      }

      std::vector<int> Idata(7);
      Idata[0] = nWtot;
      Idata[1] = nblks;
      Idata[2] = wlk_nterms;
      Idata[3] = wlk_sz;
      Idata[4] = NMO;
      Idata[5] = NAEA;
      Idata[6] = NAEB;

      dump.push("Walkers");
      dump.push("WalkerSet");
      dump.write(Idata, "dims");
    }

    int nsent = 0;
    // ready to send walkers to head in blocks
    for (int i = 0, ndone = 0; i < nblks; i++, ndone += nwlk_per_block)
    {
      boost::multi::array<ComplexType, 2> SendBuff;
      int nwlk_tot   = std::min(nwlk_per_block, nWtot - ndone);
      int nw_to_send = 0;
      if (w0 + nsent >= ndone && w0 + nsent < ndone + nwlk_tot)
        nw_to_send = std::min(nW - nsent, (ndone + nwlk_tot) - (w0 + nsent));

      if (TG.TG_heads().root())
      {
        for (int p = 0, nt = 0; p < TG.TG_heads().size(); p++)
        {
          int n_ = 0;
          if (ndone + nwlk_tot > nt && ndone < nt + nW)
          {
            if (ndone <= nt)
              n_ = std::min(nW, (ndone + nwlk_tot) - nt);
            else
              n_ = std::min(nt + nW - ndone, nwlk_tot);
          }

          counts[p] = n_ * wlk_nterms;
          nt += nw_per_tg[p];
        }
        displ[0] = 0;
        for (int p = 1, nt = 0; p < TG.TG_heads().size(); p++)
        {
          nt += counts[p - 1];
          displ[p] = nt;
        }

        RecvBuff.reextent({nwlk_tot, wlk_nterms});
      }

      if (nw_to_send > 0)
      {
        SendBuff.reextent({nw_to_send, wlk_nterms});
        for (int p = 0; p < nw_to_send; p++)
        {
          wset.copyToIO(SendBuff[p], nsent + p);
        }
      }

      TG.TG_heads().gatherv_n(SendBuff.origin(), SendBuff.num_elements(), RecvBuff.origin(), counts.data(),
                              displ.data(), 0);
      nsent += nw_to_send;

      if (TG.TG_heads().root())
      {
        dump.write(RecvBuff, std::string("walkers_") + std::to_string(i));
        wlk_per_blk.push_back(nwlk_tot);
      }

      // not sure if necessary, but avoids avalanche of messages on head node
      TG.TG_heads().barrier();
    }

    if (TG.TG_heads().root())
    {
      dump.write(wlk_per_blk, "wlk_per_blk");
      dump.pop();
      dump.pop();
    }
#endif
  }

  TG.global_barrier();
  return true;
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
