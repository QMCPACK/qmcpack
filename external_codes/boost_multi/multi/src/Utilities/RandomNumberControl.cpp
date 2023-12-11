//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <Configuration.h>
#include "Concurrency/OpenMP.h"
#include "OhmmsData/AttributeSet.h"
#include "RandomNumberControl.h"
#include "Utilities/Timer.h"
#include "hdf/HDFVersion.h"
#include "hdf/hdf_archive.h"
#include "mpi/collectives.h"
#include "Utilities/SimpleParser.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
///initialize the static data members
PrimeNumberSet<RandomBase<QMCTraits::FullPrecRealType>::uint_type> RandomNumberControl::PrimeNumbers;
UPtrVector<RandomBase<QMCTraits::FullPrecRealType>> RandomNumberControl::Children;
RandomBase<QMCTraits::FullPrecRealType>::uint_type RandomNumberControl::Offset = 11u;

/// constructors and destructors
RandomNumberControl::RandomNumberControl(const char* aname)
    : OhmmsElementBase(aname), NeverBeenInitialized(true), myCur(NULL) //, Offset(5)
{}

/// generic output
bool RandomNumberControl::get(std::ostream& os) const
{
  if (omp_get_max_threads() > 1)
  {
    for (int ip = 0; ip < omp_get_max_threads(); ip++)
    {
      Children[ip]->write(os);
      os << std::endl;
    }
  }
  else
  {
    Random.write(os);
  }
  return true;
}

/// generic input
bool RandomNumberControl::put(std::istream& is) { return true; }

/// reset the generator
void RandomNumberControl::reset() { make_seeds(); }

/// reset the generator
void RandomNumberControl::make_seeds()
{
  int pid         = OHMMS::Controller->rank();
  int nprocs      = OHMMS::Controller->size();
  uint_type iseed = static_cast<uint_type>(std::time(0)) % 1024;
  mpi::bcast(*OHMMS::Controller, iseed);
  //OHMMS::Controller->bcast(iseed);//broadcast the seed
  Offset = iseed;
  std::vector<uint_type> mySeeds;
  RandomNumberControl::PrimeNumbers.get(Offset, nprocs * (omp_get_max_threads() + 2), mySeeds);
  Random.init(mySeeds[pid]);
  //change children as well
  make_children();
}

void RandomNumberControl::make_children()
{
  int nthreads = omp_get_max_threads();
  int n        = nthreads - Children.size();
  while (n)
  {
    Children.push_back(std::make_unique<RandomGenerator>());
    n--;
  }
  int rank       = OHMMS::Controller->rank();
  int nprocs     = OHMMS::Controller->size();
  int baseoffset = Offset + nprocs + nthreads * rank;
  std::vector<uint_type> myprimes;
  PrimeNumbers.get(baseoffset, nthreads, myprimes);
  for (int ip = 0; ip < nthreads; ip++)
    Children[ip]->init(myprimes[ip]);
}

xmlNodePtr RandomNumberControl::initialize(xmlXPathContextPtr acontext)
{
  OhmmsXPathObject rg_request("//random", acontext);
  put(rg_request[0]);
  return myCur;
}

void RandomNumberControl::test()
{
  /* Add random number generator tester
  */
  int nthreads = omp_get_max_threads();
  std::vector<double> avg(nthreads), avg2(nthreads);
#pragma omp parallel for
  for (int ip = 0; ip < nthreads; ++ip)
  {
    const int n = 1000000;
    double sum = 0.0, sum2 = 0.0;
    auto& myrand(*Children[ip]);
    for (int i = 0; i < n; ++i)
    {
      double r = myrand();
      sum += r;
      sum2 += r * r;
    }
    avg[ip]  = sum / static_cast<double>(n);
    avg2[ip] = sum2 / static_cast<double>(n);
  }
  std::vector<double> avg_tot(nthreads * OHMMS::Controller->size()), avg2_tot(nthreads * OHMMS::Controller->size());
  mpi::gather(*OHMMS::Controller, avg, avg_tot);
  mpi::gather(*OHMMS::Controller, avg2, avg2_tot);
  double avg_g  = 0.0;
  double avg2_g = 0.0;
  for (int i = 0, ii = 0; i < OHMMS::Controller->size(); ++i)
  {
    for (int ip = 0; ip < nthreads; ++ip, ++ii)
    {
      app_log() << "RNGTest " << std::setw(4) << i << std::setw(4) << ip << std::setw(20) << avg_tot[ii]
                << std::setw(20) << avg2_tot[ii] - avg_tot[ii] * avg_tot[ii] << std::endl;
      avg_g += avg_tot[ii];
      avg2_g += avg2_tot[ii];
    }
  }
  avg_g /= static_cast<double>(nthreads * OHMMS::Controller->size());
  avg2_g /= static_cast<double>(nthreads * OHMMS::Controller->size());
  app_log() << "RNGTest " << std::setw(4) << OHMMS::Controller->size() << std::setw(4) << nthreads << std::setw(20)
            << avg_g << std::setw(20) << avg2_g - avg_g * avg_g << std::endl;
  app_log().flush();
}

bool RandomNumberControl::put(xmlNodePtr cur)
{
  if (NeverBeenInitialized)
  {
    bool init_mpi = true;
    int offset_in = -1; // default is to generate by Wall-clock
    if (cur != NULL)
    {
      std::string pname("yes");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(pname, "parallel");
      oAttrib.add(offset_in, "seed");
      oAttrib.put(cur);
      if (pname == "0" || pname == "false" || pname == "no")
        init_mpi = false;
    }
    int nprocs = 1;
    int pid    = 0;
    if (init_mpi)
    {
      pid    = OHMMS::Controller->rank();
      nprocs = OHMMS::Controller->size();
    }

    app_summary() << std::endl;
    app_summary() << " Random Number" << std::endl;
    app_summary() << " -------------" << std::endl;
    if (offset_in < 0)
    {
      offset_in = static_cast<int>(static_cast<uint_type>(std::time(0)) % 1024);
      app_summary() << "  Offset for the random number seeds based on time: " << offset_in << std::endl;
      mpi::bcast(*OHMMS::Controller, offset_in);
    }
    else
    {
      offset_in %= 1024;
      app_summary() << "  Offset for the random number seeds from input file (mod 1024): " << offset_in << std::endl;
    }
    app_summary() << std::endl;
    Offset = offset_in;
    std::vector<uint_type> mySeeds;
    //allocate twice of what is required
    PrimeNumbers.get(Offset, nprocs * (omp_get_max_threads() + 2), mySeeds);
    Random.init(mySeeds[pid]);
    app_log() << "  Range of prime numbers to use as seeds over processors and threads = " << mySeeds[0] << "-"
              << mySeeds[nprocs * omp_get_max_threads()] << std::endl;
    app_log() << std::endl;

    make_children();
    NeverBeenInitialized = false;
  }
  else
    reset();
  return true;
}

/*New functions past this point*/
//switch between read functions
void RandomNumberControl::read(const std::string& fname, Communicate* comm)
{
  std::string h5name = fname + ".random.h5";
  hdf_archive hin(comm, true); //attempt to read in parallel
  hin.open(h5name, H5F_ACC_RDONLY);
  if (hin.is_parallel())
    read_parallel(hin, comm);
  else
    read_rank_0(hin, comm);
}

void RandomNumberControl::write(const std::string& fname, Communicate* comm)
{
  write(convertUPtrToRefVector(Children), fname, comm);
}

//switch between write functions
void RandomNumberControl::write(const RefVector<RandomBase<FullPrecRealType>>& rng,
                                const std::string& fname,
                                Communicate* comm)
{
  std::string h5name = fname + ".random.h5";
  hdf_archive hout(comm, true); //attempt to write in parallel
  hout.create(h5name);
  if (hout.is_parallel())
    write_parallel(rng, hout, comm);
  else
    write_rank_0(rng, hout, comm);
}

//Parallel read
void RandomNumberControl::read_parallel(hdf_archive& hin, Communicate* comm)
{
  // cast integer to size_t
  const size_t nthreads  = static_cast<size_t>(omp_get_max_threads());
  const size_t comm_size = static_cast<size_t>(comm->size());
  const size_t comm_rank = static_cast<size_t>(comm->rank());

  std::vector<uint_type> vt, mt;
  TinyVector<int, 3> shape_now(comm->size(), nthreads, Random.state_size()); //cur configuration
  TinyVector<int, 3> shape_hdf5(3, 0);                                       //configuration when file was written

  //grab shape and Random.state_size() used to create hdf5 file
  hin.push(hdf::main_state);
  hin.read(shape_hdf5, "nprocs_nthreads_statesize");

  //if hdf5 file's shape and the current shape don't match, abort read
  if (shape_hdf5[0] != shape_now[0] || shape_hdf5[1] != shape_now[1] || shape_hdf5[2] != shape_now[2])
  {
    app_log() << "Mismatched random number generators."
              << "\n  Number of procs in streams : old=" << shape_hdf5[0] << " new= " << shape_now[0]
              << "\n  Number of threads in streams : old=" << shape_hdf5[1] << " new= " << shape_now[1]
              << "\n  State size per stream : old=" << shape_hdf5[2] << " new= " << shape_now[2]
              << "\n  Using the random streams generated at the initialization.\n";
    return;
  }
  app_log() << "  Restart from the random number streams from the previous configuration.\n";

  vt.resize(nthreads * Random.state_size()); //buffer for children[ip]
  mt.resize(Random.state_size());            //buffer for single thread Random object of random nums

  std::array<size_t, 2> shape{comm_size * nthreads, Random.state_size()}; //global dims of children dataset
  std::array<size_t, 2> counts{nthreads, Random.state_size()};            //local dimensions of dataset
  std::array<size_t, 2> offsets{comm_rank * nthreads, 0};                 //offsets for each process to read in

  hin.push("random"); //group that holds children[ip] random nums
  hyperslab_proxy<std::vector<uint_type>, 2> slab(vt, shape, counts, offsets);
  hin.read(slab, Random.EngineName);

  hin.pop();
  hin.push("random_master"); //group that holds Random_th random nums
  shape[0]   = comm->size(); //reset shape, counts and offset for non-multiple threads
  counts[0]  = 1;
  offsets[0] = comm->rank();
  hyperslab_proxy<std::vector<uint_type>, 2> slab2(mt, shape, counts, offsets);
  hin.read(slab2, Random.EngineName);
  hin.close();

  std::vector<uint_type>::iterator vt_it(vt.begin());
  for (int ip = 0; ip < nthreads; ip++, vt_it += shape[1])
  {
    std::vector<uint_type> c(vt_it, vt_it + shape[1]);
    Children[ip]->load(c); //load random nums back to program from buffer
  }
  Random.load(mt); //load random nums back to prog from buffer
}

//Parallel write
void RandomNumberControl::write_parallel(const RefVector<RandomBase<FullPrecRealType>>& rng,
                                         hdf_archive& hout,
                                         Communicate* comm)
{
  // cast integer to size_t
  const size_t nthreads  = static_cast<size_t>(omp_get_max_threads());
  const size_t comm_size = static_cast<size_t>(comm->size());
  const size_t comm_rank = static_cast<size_t>(comm->rank());

  std::vector<uint_type> vt, mt;
  TinyVector<int, 3> shape_hdf5(comm->size(), nthreads, Random.state_size()); //configuration at write time
  vt.reserve(nthreads * Random.state_size()); //buffer for random numbers from children[ip] of each thread
  mt.reserve(Random.state_size());            //buffer for random numbers from single Random object

  std::vector<uint_type> c;
  for (int ip = 0; ip < nthreads; ++ip)
  {
    rng[ip].get().save(c);
    vt.insert(vt.end(), c.begin(), c.end()); //get nums from each thread into buffer
  }
  Random.save(mt); //get nums for single random object (no threads)

  std::array<size_t, 2> shape{comm_size * nthreads, Random.state_size()}; //global dimensions
  std::array<size_t, 2> counts{nthreads, Random.state_size()};            //local dimensions
  std::array<size_t, 2> offsets{comm_rank * nthreads, 0};                 //offset for the file write

  hout.push(hdf::main_state);
  hout.write(shape_hdf5, "nprocs_nthreads_statesize"); //save the shape of the data at write

  hout.push("random"); //group for children[ip]
  hyperslab_proxy<std::vector<uint_type>, 2> slab(vt, shape, counts, offsets);
  hout.write(slab, Random.EngineName); //write to hdf5file
  hout.pop();

  shape[0]   = comm->size(); //adjust shape, counts, offset for just one thread
  counts[0]  = 1;
  offsets[0] = comm->rank();
  hout.push("random_master"); //group for random object without threads
  hyperslab_proxy<std::vector<uint_type>, 2> slab2(mt, shape, counts, offsets);
  hout.write(slab2, Random.EngineName); //write data to hdf5 file
  hout.close();
}

//Scatter read
void RandomNumberControl::read_rank_0(hdf_archive& hin, Communicate* comm)
{
  // cast integer to size_t
  const size_t nthreads  = static_cast<size_t>(omp_get_max_threads());
  const size_t comm_size = static_cast<size_t>(comm->size());
  const size_t comm_rank = static_cast<size_t>(comm->rank());

  std::vector<uint_type> vt, vt_tot, mt, mt_tot;
  TinyVector<size_t, 3> shape_now(comm_size, nthreads, Random.state_size()); //current configuration
  TinyVector<size_t, 3> shape_hdf5;                                          //configuration when hdf5 file was written
  std::array<size_t, 2> shape{comm_size * nthreads, Random.state_size()};    //dimensions of children dataset

  //grab configuration of threads/procs and Random.state_size() in hdf5 file
  if (comm->rank() == 0)
  {
    hin.push(hdf::main_state);
    hin.read(shape_hdf5, "nprocs_nthreads_statesize");
  }

  mpi::bcast(*comm, shape_hdf5);

  //if hdf5 file's configuration and current configuration don't match, abort read
  if (shape_hdf5[0] != shape_now[0] || shape_hdf5[1] != shape_now[1] || shape_hdf5[2] != shape_now[2])
  {
    app_log() << "Mismatched random number generators."
              << "\n  Number of procs in streams : old=" << shape_hdf5[0] << " new= " << shape_now[0]
              << "\n  Number of threads in streams : old=" << shape_hdf5[1] << " new= " << shape_now[1]
              << "\n  State size per stream : old=" << shape_hdf5[2] << " new= " << shape_now[2]
              << "\n  Using the random streams generated at the initialization.\n";
    return;
  }
  app_log() << "  Restart from the random number streams from the previous configuration.\n";

  vt.resize(nthreads * Random.state_size()); //buffer for random nums in children of each thread
  mt.resize(Random.state_size());            //buffer for random numbers from single Random object

  if (comm->rank() == 0)
  {
    hin.push("random"); //group for children[ip] (Random.object for each thread)
    vt_tot.resize(nthreads * Random.state_size() * comm->size());
    hin.readSlabReshaped(vt_tot, shape, Random.EngineName);
    hin.pop();

    shape[0] = comm->size(); //reset shape to one thread per process
    mt_tot.resize(Random.state_size() * comm->size());
    hin.push("random_master"); //group for single Random object
    hin.readSlabReshaped(mt_tot, shape, Random.EngineName);
    hin.close();
  }

  if (comm->size() > 1)
  {
    mpi::scatter(*comm, vt_tot, vt); //divide big buffer into on for each proc
    mpi::scatter(*comm, mt_tot, mt);
  }
  else
  {
    copy(vt_tot.begin(), vt_tot.end(), vt.begin());
    copy(mt_tot.begin(), mt_tot.end(), mt.begin());
  }

  std::vector<uint_type>::iterator vt_it(vt.begin());
  for (int i = 0; i < nthreads; i++, vt_it += shape[1])
  {
    std::vector<uint_type> c(vt_it, vt_it + shape[1]);
    Children[i]->load(c); //read seeds for each thread from buffer back into object
  }
  Random.load(mt); //read seeds back into object
}

//scatter write
void RandomNumberControl::write_rank_0(const RefVector<RandomBase<FullPrecRealType>>& rng,
                                       hdf_archive& hout,
                                       Communicate* comm)
{
  // cast integer to size_t
  const size_t nthreads  = static_cast<size_t>(omp_get_max_threads());
  const size_t comm_size = static_cast<size_t>(comm->size());
  const size_t comm_rank = static_cast<size_t>(comm->rank());

  std::vector<uint_type> vt, vt_tot, mt, mt_tot;
  std::array<size_t, 2> shape{comm_size * nthreads, Random.state_size()};     //dimensions of children dataset
  TinyVector<size_t, 3> shape_hdf5(comm_size, nthreads, Random.state_size()); //configuration at write time
  vt.reserve(nthreads * Random.state_size()); //buffer for children[ip] (Random object of seeds for each thread)
  mt.reserve(Random.state_size()); //buffer for single Random object of seeds, one per proc regardless of thread num

  for (int i = 0; i < nthreads; ++i)
  {
    std::vector<uint_type> c;
    rng[i].get().save(c);
    vt.insert(vt.end(), c.begin(), c.end()); //copy children[nthreads] seeds to buffer
  }
  Random.save(mt); //copy random_th seeds to buffer

  if (comm->size() > 1)
  {
    vt_tot.resize(vt.size() * comm->size());
    mt_tot.resize(mt.size() * comm->size());
    mpi::gather(*comm, vt, vt_tot); //gather into one big buffer for master write
    mpi::gather(*comm, mt, mt_tot);
  }
  else
  {
    vt_tot = vt;
    mt_tot = mt;
  }

  if (comm->rank() == 0)
  {
    hout.push(hdf::main_state);
    hout.write(shape_hdf5, "nprocs_nthreads_statesize"); //configuration at write time to file

    hout.push("random"); //group for children[ip]
    hout.writeSlabReshaped(vt_tot, shape, Random.EngineName);
    hout.pop();

    shape[0] = comm->size();    //reset dims for single thread use
    hout.push("random_master"); //group for random_th object
    hout.writeSlabReshaped(mt_tot, shape, Random.EngineName);
    hout.close();
  }
}
} // namespace qmcplusplus
