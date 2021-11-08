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
#include <numeric>
#include <algorithm>

#include "EstimatorManagerNew.h"
#include "SpinDensityNew.h"
#include "MomentumDistribution.h"
#include "OneBodyDensityMatrices.h"
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
#include "hdf/hdf_archive.h"
#include "OhmmsData/AttributeSet.h"
#include "Estimators/CSEnergyEstimator.h"

//leave it for serialization debug
//#define DEBUG_ESTIMATOR_ARCHIVE

namespace qmcplusplus
{
//initialize the name of the primary estimator
EstimatorManagerNew::EstimatorManagerNew(Communicate* c)
    : MainEstimatorName("LocalEnergy"), RecordCount(0), my_comm_(c), Collectables(0), max4ascii(8), FieldWidth(20)
{}

EstimatorManagerNew::~EstimatorManagerNew()
{
  if (Collectables)
    delete Collectables;
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
  cpuInd         = BlockProperties.add("BlockCPU");
  acceptRatioInd = BlockProperties.add("AcceptRatio");
  BlockAverages.clear(); //cleaup the records
  for (int i = 0; i < Estimators.size(); i++)
    Estimators[i]->add2Record(BlockAverages);
  max4ascii += BlockAverages.size();
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
void EstimatorManagerNew::startDriverRun()
{
  reset();
  RecordCount = 0;
  energyAccumulator.clear();
  varAccumulator.clear();
  int nc = (Collectables) ? Collectables->size() : 0;
  BlockAverages.setValues(0.0);
  // \todo Collectables should just have its own data structures not change the EMBS layout.
  AverageCache.resize(BlockAverages.size() + nc);
  PropertyCache.resize(BlockProperties.size());
#if defined(DEBUG_ESTIMATOR_ARCHIVE)
  if (!DebugArchive)
  {
    char fname[128];
    sprintf(fname, "%s.p%03d.scalar.dat", my_comm_->getName().c_str(), my_comm_->rank());
    DebugArchive = std::make_unique<std::ofstream>(fname);
    addHeader(*DebugArchive);
  }
#endif
  if (my_comm_->rank() == 0)
  {
    std::string fname(my_comm_->getName());
    fname.append(".scalar.dat");
    Archive = std::make_unique<std::ofstream>(fname.c_str());
    addHeader(*Archive);
    if (h5desc.size())
    {
      h5desc.clear();
    }
    fname  = my_comm_->getName() + ".stat.h5";
    h_file = std::make_unique<hdf_archive>();
    h_file->create(fname);
    for (int i = 0; i < Estimators.size(); i++)
      Estimators[i]->registerObservables(h5desc, h_file->getFileID());
    for (auto& uope : operator_ests_)
      uope->registerOperatorEstimator(h_file->getFileID());
  }
}

void EstimatorManagerNew::stopDriverRun() { h_file.reset(); }

void EstimatorManagerNew::startBlock(int steps) { block_timer_.restart(); }

void EstimatorManagerNew::stopBlock(unsigned long accept, unsigned long reject, RealType block_weight)
{
  /* Need a redesign of how accept, reject and block_weight are handled from driver to this manager.
   * DMC needs to add non-local move counters.
   * also need to add num_samples which differs from block_weight
   */
  //take block averages and update properties per block
  PropertyCache[weightInd] = block_weight;
  makeBlockAverages(accept, reject);
  reduceOperatorEstimators();
  writeOperatorEstimators();
  zeroOperatorEstimators();
  // intentionally put after all the estimator I/O
  PropertyCache[cpuInd] = block_timer_.elapsed();
  writeScalarH5();
  RecordCount++;
}

void EstimatorManagerNew::collectScalarEstimators(const RefVector<ScalarEstimatorBase>& estimators)
{
  AverageCache = 0.0;
  for (ScalarEstimatorBase& est : estimators)
    est.addAccumulated(AverageCache.begin());
}

void EstimatorManagerNew::collectOperatorEstimators(const std::vector<RefVector<OperatorEstBase>>& crowd_op_ests)
{
  for (int iop = 0; iop < operator_ests_.size(); ++iop)
  {
    RefVector<OperatorEstBase> this_op_est_for_all_crowds;
    for (int icrowd = 0; icrowd < crowd_op_ests.size(); ++icrowd)
      this_op_est_for_all_crowds.emplace_back(crowd_op_ests[icrowd][iop]);
    operator_ests_[iop]->collect(this_op_est_for_all_crowds);
  }
}

void EstimatorManagerNew::makeBlockAverages(unsigned long accepts, unsigned long rejects)
{
  // accumulate unsigned long counters over ranks.
  // Blocks ends are infrequent, two MPI transfers to preserve type
  // these could be replaced with a singple call MPI_struct_type some packing scheme or even
  // a pack into and out of an fp type that can be assured to hold the integral type exactly
  // IMHO they should not be primarily stored in a vector with magic indexes
  std::vector<unsigned long> accepts_and_rejects(my_comm_->size() * 2, 0);
  accepts_and_rejects[my_comm_->rank()]                    = accepts;
  accepts_and_rejects[my_comm_->size() + my_comm_->rank()] = rejects;
  my_comm_->allreduce(accepts_and_rejects);
  unsigned long total_block_accept =
      std::accumulate(accepts_and_rejects.begin(), accepts_and_rejects.begin() + my_comm_->size(), 0);
  unsigned long total_block_reject = std::accumulate(accepts_and_rejects.begin() + my_comm_->size(),
                                                     accepts_and_rejects.begin() + my_comm_->size() * 2, 0);

  //Transfer FullPrecisionRead data
  const size_t n1 = AverageCache.size();
  const size_t n2 = n1 + PropertyCache.size();

  // This is a hack but it needs to be the correct size

  std::vector<double> reduce_buffer(n2, 0.0);
  {
    auto cur = reduce_buffer.begin();
    copy(AverageCache.begin(), AverageCache.end(), cur);
    copy(PropertyCache.begin(), PropertyCache.end(), cur + n1);
  }

  // This is necessary to use mpi3's C++ style reduce
#ifdef HAVE_MPI
  my_comm_->comm.reduce_in_place_n(reduce_buffer.begin(), reduce_buffer.size(), std::plus<>{});
#endif
  if (my_comm_->rank() == 0)
  {
    auto cur = reduce_buffer.begin();
    copy(cur, cur + n1, AverageCache.begin());
    copy(cur + n1, cur + n2, PropertyCache.begin());
    const RealType invTotWgt = 1.0 / PropertyCache[weightInd];
    AverageCache *= invTotWgt;
    //do not weight weightInd i.e. its index 0!
    for (int i = 1; i < PropertyCache.size(); i++)
      PropertyCache[i] *= invTotWgt;
  }

  // now we put the correct accept ratio in
  PropertyCache[acceptRatioInd] = static_cast<FullPrecRealType>(total_block_accept) /
      static_cast<FullPrecRealType>(total_block_accept + total_block_reject);

  //add the block average to summarize
  energyAccumulator(AverageCache[0]);
  varAccumulator(AverageCache[1]);
}

void EstimatorManagerNew::writeScalarH5()
{
  //Do not assume h_file is valid
  if (h_file)
  {
    for (int o = 0; o < h5desc.size(); ++o)
      // cheating here, remove SquaredAverageCache from API
      h5desc[o].write(AverageCache.data(), AverageCache.data());
    H5Fflush(h_file->getFileID(), H5F_SCOPE_LOCAL);
  }

  if (Archive)
  {
    *Archive << std::setw(10) << RecordCount;
    int maxobjs = std::min(BlockAverages.size(), max4ascii);
    for (int j = 0; j < maxobjs; j++)
      *Archive << std::setw(FieldWidth) << AverageCache[j];
    for (int j = 0; j < PropertyCache.size(); j++)
      *Archive << std::setw(FieldWidth) << PropertyCache[j];
    *Archive << std::endl;
  }
}

void EstimatorManagerNew::reduceOperatorEstimators()
{
  if (operator_ests_.size() > 0)
  {
    std::vector<size_t> operator_data_sizes(operator_ests_.size());
    RefVector<OperatorEstBase> ref_op_ests = convertUPtrToRefVector(operator_ests_);
    for (int iop = 0; iop < operator_data_sizes.size(); ++iop)
    {
      operator_data_sizes[iop] = operator_ests_[iop]->get_data().size();
    }
    // 1 larger because we put the weight in to avoid dependence of the Scalar estimators being reduced firt.
    size_t nops = *(std::max_element(operator_data_sizes.begin(), operator_data_sizes.end())) + 1;
    std::vector<RealType> operator_send_buffer;
    std::vector<RealType> operator_recv_buffer;
    operator_send_buffer.reserve(nops);
    operator_recv_buffer.reserve(nops);
    for (int iop = 0; iop < operator_ests_.size(); ++iop)
    {
      auto& estimator      = *operator_ests_[iop];
      auto& data           = estimator.get_data();
      size_t adjusted_size = data.size() + 1;
      operator_send_buffer.resize(adjusted_size, 0.0);
      operator_recv_buffer.resize(adjusted_size, 0.0);
      std::copy_n(data.begin(), data.size(), operator_send_buffer.begin());
      operator_send_buffer[data.size()] = estimator.get_walkers_weight();
      // This is necessary to use mpi3's C++ style reduce
#ifdef HAVE_MPI
      my_comm_->comm.reduce_n(operator_send_buffer.begin(), adjusted_size, operator_recv_buffer.begin(), std::plus<>{},
                              0);
#else
      operator_recv_buffer = operator_send_buffer;
#endif
      if (my_comm_->rank() == 0)
      {
        std::copy_n(operator_recv_buffer.begin(), data.size(), data.begin());
        size_t reduced_walker_weights = operator_recv_buffer[data.size()];
        RealType invTotWgt            = 1.0 / static_cast<QMCT::RealType>(reduced_walker_weights);
        operator_ests_[iop]->normalize(invTotWgt);
      }
    }
  }
}

void EstimatorManagerNew::writeOperatorEstimators()
{
  if (my_comm_->rank() == 0)
  {
    if (h_file)
    {
      for (auto& op_est : operator_ests_)
        op_est->write();
      H5Fflush(h_file->getFileID(), H5F_SCOPE_LOCAL);
    }
  }
}

void EstimatorManagerNew::zeroOperatorEstimators()
{
  for (auto& op_est : operator_ests_)
    op_est->zero();
}

void EstimatorManagerNew::getApproximateEnergyVariance(RealType& e, RealType& var)
{
  RealType tmp[3];
  tmp[0] = energyAccumulator.count();
  tmp[1] = energyAccumulator.result();
  tmp[2] = varAccumulator.result();
  my_comm_->bcast(tmp, 3);
  e   = tmp[1] / tmp[0];
  var = tmp[2] / tmp[0] - e * e;
}

EstimatorManagerNew::EstimatorType* EstimatorManagerNew::getEstimator(const std::string& a)
{
  std::map<std::string, int>::iterator it = EstimatorMap.find(a);
  if (it == EstimatorMap.end())
    return nullptr;
  else
    return Estimators[(*it).second].get();
}

bool EstimatorManagerNew::put(QMCHamiltonian& H, const ParticleSet& pset, const TrialWaveFunction& twf, const WaveFunctionFactory& wf_factory, xmlNodePtr cur)
{
  std::vector<std::string> extra_types;
  std::vector<std::string> extra_names;
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "estimator")
    {
      std::string est_type("none");
      std::string est_name(MainEstimatorName);
      std::string use_hdf5("yes");
      OhmmsAttributeSet hAttrib;
      hAttrib.add(est_type, "type");
      hAttrib.add(est_name, "name");
      hAttrib.add(use_hdf5, "hdf5");
      hAttrib.put(cur);
      if ((est_name == MainEstimatorName) || (est_name == "elocal"))
      {
        max4ascii = H.sizeOfObservables() + 3;
        add(std::make_unique<LocalEnergyEstimator>(H, use_hdf5 == "yes"), MainEstimatorName);
      }
      else if (est_name == "RMC")
      {
        int nobs(20);
        OhmmsAttributeSet hAttrib;
        hAttrib.add(nobs, "nobs");
        hAttrib.put(cur);
        max4ascii = nobs * H.sizeOfObservables() + 3;
        add(std::make_unique<RMCLocalEnergyEstimator>(H, nobs), MainEstimatorName);
      }
      else if (est_name == "CSLocalEnergy")
      {
        OhmmsAttributeSet hAttrib;
        int nPsi = 1;
        hAttrib.add(nPsi, "nPsi");
        hAttrib.put(cur);
        add(std::make_unique<CSEnergyEstimator>(H, nPsi), MainEstimatorName);
        app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
      }
      else if (est_name == "SpinDensityNew")
      {
        SpinDensityInput spdi;
        spdi.readXML(cur);
        DataLocality dl = DataLocality::crowd;
        if (spdi.get_save_memory())
          dl = DataLocality::rank;
        if (spdi.get_cell().explicitly_defined)
          operator_ests_.emplace_back(std::make_unique<SpinDensityNew>(std::move(spdi), pset.mySpecies, dl));
        else
          operator_ests_.emplace_back(
              std::make_unique<SpinDensityNew>(std::move(spdi), pset.Lattice, pset.mySpecies, dl));
      }
      else if (est_type == "MomentumDistribution")
      {
        MomentumDistributionInput mdi;
        mdi.readXML(cur);
        DataLocality dl = DataLocality::crowd;
        operator_ests_.emplace_back(
          std::make_unique<MomentumDistribution>(std::move(mdi), 
            pset.getTotalNum(), pset.getTwist(), pset.Lattice, dl));
      }
      else if (est_type == "OneBodyDensityMatrices")
      {
        OneBodyDensityMatricesInput obdmi(cur);
        // happens once insures golden particle set is not abused.
        ParticleSet pset_target(pset);
        operator_ests_.emplace_back(
          std::make_unique<OneBodyDensityMatrices>(std::move(obdmi), 
                                                   pset.Lattice, pset.getSpeciesSet(), wf_factory, pset_target));
      }
      else
      {
        extra_types.push_back(est_type);
        extra_names.push_back(est_name);
      }
    }
    cur = cur->next;
  }
  if (Estimators.empty())
  {
    app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
    max4ascii = H.sizeOfObservables() + 3;
    add(std::make_unique<LocalEnergyEstimator>(H, true), MainEstimatorName);
  }
  //Collectables is special and should not be added to Estimators
  if (Collectables == 0 && H.sizeOfCollectables())
  {
    app_log() << "  Using CollectablesEstimator for collectables, e.g. sk, gofr, density " << std::endl;
    Collectables = new CollectablesEstimator(H);
  }
  // Unrecognized estimators are not allowed
  if (!extra_types.empty())
  {
    app_log() << "\nUnrecognized estimators in input:" << std::endl;
    for (int i=0; i<extra_types.size(); i++)
    {
      app_log() << "  type: "<<extra_types[i]<<"     name: "<<extra_names[i]<<std::endl;
    }
    app_log() << std::endl;
    throw UniformCommunicateError("Unrecognized estimators encountered in input.  See log message for more details.");
  }
  return true;
}

int EstimatorManagerNew::add(std::unique_ptr<EstimatorType> newestimator, const std::string& aname)
{
  std::map<std::string, int>::iterator it = EstimatorMap.find(aname);
  int n                                   = Estimators.size();
  if (it == EstimatorMap.end())
  {
    Estimators.push_back(std::move(newestimator));
    EstimatorMap[aname] = n;
  }
  else
  {
    n = (*it).second;
    app_log() << "  EstimatorManagerNew::add replace " << aname << " estimator." << std::endl;
    Estimators[n] = std::move(newestimator);
  }
  return n;
}

} // namespace qmcplusplus
