
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>

#include "Configuration.h"
#include "Utilities/FairDivide.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/FairDivide.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "PropagatorFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"

namespace qmcplusplus
{
namespace afqmc
{
Propagator PropagatorFactory::buildAFQMCPropagator(TaskGroup_& TG,
                                                   xmlNodePtr cur,
                                                   Wavefunction& wfn,
                                                   RandomBase<RealType>& rng)
{
  // allocator for local memory
  using allocator = device_allocator<ComplexType>;
  using pointer   = device_ptr<ComplexType>;
  using CVector   = boost::multi::array<ComplexType, 1, allocator>;
  using CMatrix   = boost::multi::array<ComplexType, 2, allocator>;

  // read options from xml
  if (cur == NULL)
    APP_ABORT(" Error: Null xml ptr in PropagatorFactory::buildAFQMCPropagator.\n");

  std::string info("info0");
  std::string name("");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(info, "info");
  oAttrib.add(name, "name");
  oAttrib.put(cur);

  if (InfoMap.find(info) == InfoMap.end())
  {
    app_error() << "ERROR: Undefined info in WavefunctionFactory. \n";
    APP_ABORT("ERROR: Undefined info in WavefunctionFactory. \n");
  }
  AFQMCInfo& AFinfo = InfoMap[info];

  RealType vbias_bound = 50.0;
  std::string sub("yes");
  std::string printP1eV("no");
  ParameterSet m_param;
  m_param.add(sub, "substractMF");
  m_param.add(printP1eV, "printP1eigval");
  m_param.add(vbias_bound, "vbias_bound");
  m_param.put(cur);

  bool substractMF = true;
  std::transform(sub.begin(), sub.end(), sub.begin(), (int (*)(int))tolower);
  if (sub == "no" || sub == "false")
    substractMF = false;

  if (substractMF)
    app_log() << " Using mean-field subtraction in propagator: " << name << "\n";

  // buld mean field expectation value of the Cholesky matrix
  CVector vMF(iextensions<1u>{wfn.local_number_of_cholesky_vectors()}, allocator{});
  using std::fill_n;
  fill_n(vMF.origin(), vMF.num_elements(), ComplexType(0));
  if (substractMF)
  {
    wfn.vMF(vMF);
    // if (not distribution_over_cholesky_vectors()), vMF needs to be reduced over TG
    if (not wfn.distribution_over_cholesky_vectors())
    {
      if (not TG.TG_local().root())
        fill_n(vMF.origin(), vMF.num_elements(), ComplexType(0));
      TG.TG().all_reduce_in_place_n(to_address(vMF.origin()), vMF.num_elements(), std::plus<>());
    }
  }

  boost::multi::array<ComplexType, 1> vMF_(vMF.extensions());
  using std::copy_n;
  copy_n(vMF.origin(), vMF.num_elements(), vMF_.origin());
  RealType vmax = 0, v_ = 0;
  for (int i = 0; i < vMF_.size(); i++)
    v_ = std::max(v_, std::abs(vMF_[i]));
  TG.Global().all_reduce_n(&v_, 1, &vmax, boost::mpi3::max<>());
  app_log() << " Largest component of Mean-field subtraction potential: " << vmax << std::endl;
  if (vmax > vbias_bound)
    app_log() << " WARNING: Mean-field subtraction potential has components outside vbias_bound.\n"
              << "          Consider increasing vbias_bound. max(vMF[n]), vbias_bound: " << vmax << " " << vbias_bound
              << std::endl;

  // assemble H1(i,j) = h(i,j) + vn0(i,j) + sum_n vMF[n]*Spvn(i,j,n)
  auto H1(wfn.getHamiltonianOperations()->getOneBodyPropagatorMatrix(TG, vMF_));
  // assemble H1(i,j) = h(i,j) + vn0(i,j) + sum_n vMF[n]*Spvn(i,j,n)

  if (TG.getNGroupsPerTG() == 1)
    return Propagator(AFQMCBasePropagator(AFinfo, cur, TG, wfn, std::move(H1), std::move(vMF), rng));
  else
  {
    if (wfn.distribution_over_cholesky_vectors())
      // use specialized distributed algorithm for case
      // when vbias doesn't need reduction over TG
      return Propagator(AFQMCDistributedPropagatorDistCV(AFinfo, cur, TG, wfn, std::move(H1), std::move(vMF), rng));
    else
      return Propagator(AFQMCDistributedPropagator(AFinfo, cur, TG, wfn, std::move(H1), std::move(vMF), rng));
  }
}


} // namespace afqmc

} // namespace qmcplusplus
