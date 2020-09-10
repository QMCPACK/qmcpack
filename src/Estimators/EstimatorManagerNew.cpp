//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <functional>

#include "Estimators/EstimatorManagerNew.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommUtilities.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "Estimators/RMCLocalEnergyEstimator.h"
#include "Estimators/CollectablesEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/WalkerProperties.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include "Estimators/CSEnergyEstimator.h"
//leave it for serialization debug
//#define DEBUG_ESTIMATOR_ARCHIVE

namespace qmcplusplus
{
/** enumeration for EstimatorManagerBase.Options
 */
enum
{
  COLLECT = 0,
  MANAGE,
  RECORD,
  POSTIRECV,
  APPEND
};

//initialize the name of the primary estimator
EstimatorManagerNew::EstimatorManagerNew(Communicate* c)
    : MainEstimatorName("LocalEnergy"),
      RecordCount(0),
      h_file(-1),
      Archive(0),
      DebugArchive(0),
      my_comm_(0),
      MainEstimator(0),
      Collectables(0),
      max4ascii(8),
      FieldWidth(20)
{
  setCommunicator(c);
}

EstimatorManagerNew::EstimatorManagerNew(EstimatorManagerNew& em)
    : MainEstimatorName(em.MainEstimatorName),
      Options(em.Options),
      RecordCount(0),
      h_file(-1),
      Archive(0),
      DebugArchive(0),
      my_comm_(0),
      MainEstimator(0),
      Collectables(0),
      EstimatorMap(em.EstimatorMap),
      max4ascii(em.max4ascii),
      FieldWidth(20)
{
  //inherit communicator
  setCommunicator(em.my_comm_);

  // Here Estimators are ScalarEstimatorNew
  for (int i = 0; i < em.Estimators.size(); i++)
    Estimators.push_back(em.Estimators[i]->clone());
  MainEstimator = Estimators[EstimatorMap[MainEstimatorName]];
  if (em.Collectables)
    Collectables = em.Collectables->clone();
}

EstimatorManagerNew::~EstimatorManagerNew()
{
  delete_iter(Estimators.begin(), Estimators.end());
  delete_iter(h5desc.begin(), h5desc.end());
  if (Collectables)
    delete Collectables;
}

void EstimatorManagerNew::setCommunicator(Communicate* c)
{
  // I think this is actually checking if this is the "Main Estimator"
  if (my_comm_ && my_comm_ == c)
    return;
  my_comm_ = c ? c : OHMMS::Controller;
  //set the default options
  // This is a flag to tell manager if there is more than one rank
  // running walkers, its discovered by smelly query of my_comm_.
  // New code should not make use of these useless options
  Options.set(COLLECT, my_comm_->size() > 1);
  Options.set(MANAGE, my_comm_->rank() == 0);
  if (RemoteData.empty())
  {
    RemoteData.push_back(UPtr<FPRBuffer>(new FPRBuffer));
    RemoteData.push_back(UPtr<FPRBuffer>(new FPRBuffer));
  }
}

/** reset names of the properties
 *
 * \todo this should be in in constructor object shouldn't be reused
 * Warning this is different from some "resets" in the code, it does not clear the object
 * 
 * The number of estimators and their order can vary from the previous state.
 * reinitialized properties before setting up a new BlockAverage data list.
 *
 */
void EstimatorManagerNew::reset()
{
  //It's essential that this one is first, other code assumes weightInd always == 0
  weightInd = BlockProperties.add("BlockWeight");
  assert(weightInd == 0);
  cpuInd    = BlockProperties.add("BlockCPU");
  acceptRatioInd = BlockProperties.add("BlockAcceptRatio");
  BlockAverages.clear(); //cleaup the records
  for (int i = 0; i < Estimators.size(); i++)
    Estimators[i]->add2Record(BlockAverages);
  max4ascii += BlockAverages.size();
  if (Collectables)
    Collectables->add2Record(BlockAverages);
}

void EstimatorManagerNew::addHeader(std::ostream& o)
{
  o.setf(std::ios::scientific, std::ios::floatfield);
  o.setf(std::ios::left, std::ios::adjustfield);
  o.precision(10);
  for (int i = 0; i < BlockAverages.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockAverages.Names[i].size() + 2);
  for (int i = 0; i < BlockProperties.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockProperties.Names[i].size() + 2);
  int maxobjs = std::min(BlockAverages.size(), max4ascii);
  o << "#   index    ";
  for (int i = 0; i < maxobjs; i++)
    o << std::setw(FieldWidth) << BlockAverages.Names[i];
  for (int i = 0; i < BlockProperties.size(); i++)
    o << std::setw(FieldWidth) << BlockProperties.Names[i];
  o << std::endl;
  o.setf(std::ios::right, std::ios::adjustfield);
}

/// \todo clean up this method its a mess
void EstimatorManagerNew::start(int blocks, bool record)
{
  for (int i = 0; i < Estimators.size(); i++)
    Estimators[i]->setNumberOfBlocks(blocks);
  reset();
  RecordCount = 0;
  energyAccumulator.clear();
  varAccumulator.clear();
  int nc = (Collectables) ? Collectables->size() : 0;
  BlockAverages.setValues(0.0);
  // \todo Collectables should just have its own data structures not change the EMBS layout.
  AverageCache.resize(BlockAverages.size() + nc);
  SquaredAverageCache.resize(BlockAverages.size() + nc);
  PropertyCache.resize(BlockProperties.size());
  //count the buffer size for message
  BufferSize  = 2 * AverageCache.size() + PropertyCache.size();
  int sources = 2;
  //allocate buffer for data collection
  if (RemoteData.empty())
    for (int i = 0; i < sources; ++i)
      RemoteData.push_back(UPtr<FPRBuffer>(new FPRBuffer(BufferSize)));
  else
    for (int i = 0; i < RemoteData.size(); ++i)
      RemoteData[i]->resize(BufferSize);
#if defined(DEBUG_ESTIMATOR_ARCHIVE)
  if (record && DebugArchive == 0)
  {
    char fname[128];
    sprintf(fname, "%s.p%03d.scalar.dat", my_comm_->getName().c_str(), my_comm_->rank());
    DebugArchive = new std::ofstream(fname);
    addHeader(*DebugArchive);
  }
#endif
  //set Options[RECORD] to enable/disable output
  Options.set(RECORD, record && Options[MANAGE]);
  if (Options[RECORD])
  {
    if (Archive)
      delete Archive;
    std::string fname(my_comm_->getName());
    fname.append(".scalar.dat");
    Archive = new std::ofstream(fname.c_str());
    addHeader(*Archive);
    if (h5desc.size())
    {
      delete_iter(h5desc.begin(), h5desc.end());
      h5desc.clear();
    }
    fname  = my_comm_->getName() + ".stat.h5";
    h_file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for (int i = 0; i < Estimators.size(); i++)
      Estimators[i]->registerObservables(h5desc, h_file);
    if (Collectables)
      Collectables->registerObservables(h5desc, h_file);
  }
}

void EstimatorManagerNew::startBlock(int steps)
{
  MyTimer.restart();
  BlockWeight = 0.0;
}

void EstimatorManagerNew::stopBlock(unsigned long accept, unsigned long reject, RealType block_weight, double cpu_block_time)
{
  //take block averages and update properties per block
  PropertyCache[weightInd] = block_weight;
  PropertyCache[cpuInd]    = cpu_block_time;
  makeBlockAverages(accept, reject);
}


/** Called at end of block in Unified Driver
 *
 */
QMCTraits::FullPrecRealType EstimatorManagerNew::collectScalarEstimators(
    const RefVector<ScalarEstimatorBase>& estimators)
{
  using ScalarType = ScalarEstimatorBase::RealType;

  AverageCache        = 0.0;
  SquaredAverageCache = 0.0;

  // One scalar estimator can be accumulating many scalar values
  int num_scalars = estimators[0].get().size();
  if (AverageCache.size() != num_scalars)
    throw std::runtime_error(
        "EstimatorManagerNew and Crowd ScalarManagers do not agree on number of scalars being estimated");
  Vector<ScalarType> averages_work(num_scalars, 0.0);
  Vector<ScalarType> sq_averages_work(num_scalars, 0.0);

  auto accumulateVectorsInPlace = [](auto& vec_a, const auto& vec_b) {
    for (int i = 0; i < vec_a.size(); ++i)
      vec_a[i] += vec_b[i];
  };

  RealType tot_weight = 0.0;
  for (int i = 0; i < estimators.size(); ++i)
  {
    RealType weight = estimators[i].get().takeBlockSumsGetWeight(averages_work.begin(), sq_averages_work.begin());
    tot_weight += weight;
    accumulateVectorsInPlace(AverageCache, averages_work);
    accumulateVectorsInPlace(SquaredAverageCache, sq_averages_work);
  }
  return tot_weight;
}

// blocks don't close frequently enough that we should be sweating the mpi transfers at all.
// all this Cache stuff is premature optimization because someone wanted to be very fancy
void EstimatorManagerNew::makeBlockAverages(unsigned long accepts, unsigned long rejects)
{
  //there is only one EstimatormanagerNew per rank in the unified driver.
  //copy cached data to RemoteData[0]
  //we should not handle RemoteData elsewhere.
  std::vector<unsigned long> accepts_and_rejects(my_comm_->size() * 2, 0);
  accepts_and_rejects[my_comm_->rank()] = accepts;
  accepts_and_rejects[my_comm_->size() + my_comm_->rank()] = rejects;
  my_comm_->allreduce(accepts_and_rejects);
  unsigned long total_block_accept = std::accumulate(accepts_and_rejects.begin(), accepts_and_rejects.begin() + my_comm_->size(), 0);
  unsigned long total_block_reject = std::accumulate(accepts_and_rejects.begin() + my_comm_->size(), accepts_and_rejects.begin() + my_comm_->size() * 2, 0);

  int n1 = AverageCache.size();
  int n2 = n1 + AverageCache.size();
  int n3 = n2 + PropertyCache.size();

  // This is a hack but it needs to be the correct size

  std::vector<double> send_buffer(n3, 0.0);
  std::vector<double> recv_buffer(n3, 0.0);
  {
    auto cur = send_buffer.begin();
    copy(AverageCache.begin(), AverageCache.end(), cur);
    copy(SquaredAverageCache.begin(), SquaredAverageCache.end(), cur + n1);
    copy(PropertyCache.begin(), PropertyCache.end(), cur + n2);
  }

  // This is necessary to use mpi3's C++ style reduce
#ifdef HAVE_MPI
  my_comm_->comm.reduce_n(send_buffer.begin(), send_buffer.size(), recv_buffer.begin(), std::plus<>{}, 0);
#else
  recv_buffer = send_buffer;
#endif
  if (my_comm_->rank() == 0)
  {
    auto cur = recv_buffer.begin();
    copy(cur, cur + n1, AverageCache.begin());
    copy(cur + n1, cur + n2, SquaredAverageCache.begin());
    copy(cur + n2, cur + n3, PropertyCache.begin());
    RealType invTotWgt = 1.0 / PropertyCache[weightInd];
    AverageCache *= invTotWgt;
    SquaredAverageCache *= invTotWgt;
    //do not weight weightInd i.e. its index 0!
    for (int i = 1; i < PropertyCache.size(); i++)
      PropertyCache[i] *= invTotWgt;
  }

  // now we put the correct accept ratio in
  PropertyCache[acceptRatioInd] = static_cast<FullPrecRealType>(total_block_accept) / static_cast<FullPrecRealType>(total_block_accept + total_block_reject);
  
  //add the block average to summarize
  energyAccumulator(AverageCache[0]);
  varAccumulator(SquaredAverageCache[0]);
  if (Archive)
  {
    *Archive << std::setw(10) << RecordCount;
    int maxobjs = std::min(BlockAverages.size(), max4ascii);
    for (int j = 0; j < maxobjs; j++)
      *Archive << std::setw(FieldWidth) << AverageCache[j];
    for (int j = 0; j < PropertyCache.size(); j++)
      *Archive << std::setw(FieldWidth) << PropertyCache[j];
    *Archive << std::endl;
    for (int o = 0; o < h5desc.size(); ++o)
      h5desc[o]->write(AverageCache.data(), SquaredAverageCache.data());
    H5Fflush(h_file, H5F_SCOPE_LOCAL);
  }
  RecordCount++;
}


void EstimatorManagerNew::getApproximateEnergyVariance(RealType& e, RealType& var)
{
  if (Options[COLLECT]) //need to broadcast the value
  {
    RealType tmp[3];
    tmp[0] = energyAccumulator.count();
    tmp[1] = energyAccumulator.result();
    tmp[2] = varAccumulator.result();
    myComm->bcast(tmp, 3);
    e   = tmp[1] / tmp[0];
    var = tmp[2] / tmp[0] - e * e;
  }
  else
  {
    e   = energyAccumulator.mean();
    var = varAccumulator.mean() - e * e;
  }
}

EstimatorManagerNew::EstimatorType* EstimatorManagerNew::getEstimator(const std::string& a)
{
  std::map<std::string, int>::iterator it = EstimatorMap.find(a);
  if (it == EstimatorMap.end())
    return 0;
  else
    return Estimators[(*it).second];
}

bool EstimatorManagerNew::put(QMCHamiltonian& H, xmlNodePtr cur)
{
  std::vector<std::string> extra;
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "estimator")
    {
      std::string est_name(MainEstimatorName);
      std::string use_hdf5("yes");
      OhmmsAttributeSet hAttrib;
      hAttrib.add(est_name, "name");
      hAttrib.add(use_hdf5, "hdf5");
      hAttrib.put(cur);
      if ((est_name == MainEstimatorName) || (est_name == "elocal"))
      {
        max4ascii = H.sizeOfObservables() + 3;
        add(new LocalEnergyEstimator(H, use_hdf5 == "yes"), MainEstimatorName);
      }
      else if (est_name == "RMC")
      {
        int nobs(20);
        OhmmsAttributeSet hAttrib;
        hAttrib.add(nobs, "nobs");
        hAttrib.put(cur);
        max4ascii = nobs * H.sizeOfObservables() + 3;
        add(new RMCLocalEnergyEstimator(H, nobs), MainEstimatorName);
      }
      else if (est_name == "CSLocalEnergy")
      {
        OhmmsAttributeSet hAttrib;
        int nPsi = 1;
        hAttrib.add(nPsi, "nPsi");
        hAttrib.put(cur);
        add(new CSEnergyEstimator(H, nPsi), MainEstimatorName);
        app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
      }
      else
        extra.push_back(est_name);
    }
    cur = cur->next;
  }
  if (Estimators.empty())
  {
    app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
    max4ascii = H.sizeOfObservables() + 3;
    add(new LocalEnergyEstimator(H, true), MainEstimatorName);
  }
  //Collectables is special and should not be added to Estimators
  if (Collectables == 0 && H.sizeOfCollectables())
  {
    app_log() << "  Using CollectablesEstimator for collectables, e.g. sk, gofr, density " << std::endl;
    Collectables = new CollectablesEstimator(H);
  }
  return true;
}

int EstimatorManagerNew::add(EstimatorType* newestimator, const std::string& aname)
{
  std::map<std::string, int>::iterator it = EstimatorMap.find(aname);
  int n                                   = Estimators.size();
  if (it == EstimatorMap.end())
  {
    Estimators.push_back(newestimator);
    EstimatorMap[aname] = n;
  }
  else
  {
    n = (*it).second;
    app_log() << "  EstimatorManagerNew::add replace " << aname << " estimator." << std::endl;
    delete Estimators[n];
    Estimators[n] = newestimator;
  }
  //check the name and set the MainEstimator
  if (aname == MainEstimatorName)
    MainEstimator = newestimator;
  return n;
}

int EstimatorManagerNew::addObservable(const char* aname)
{
  int mine = BlockAverages.add(aname);
  int add  = TotalAverages.add(aname);
  if (mine < Block2Total.size())
    Block2Total[mine] = add;
  else
    Block2Total.push_back(add);
  return mine;
}

void EstimatorManagerNew::getData(int i, std::vector<RealType>& values)
{
  int entries = TotalAveragesData.rows();
  values.resize(entries);
  for (int a = 0; a < entries; a++)
    values[a] = TotalAveragesData(a, Block2Total[i]);
}

} // namespace qmcplusplus
