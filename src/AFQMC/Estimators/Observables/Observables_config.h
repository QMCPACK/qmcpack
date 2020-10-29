#ifndef QMCPLUSPLUS_AFQMC_OBSERVABLES_CONFIG_HPP
#define QMCPLUSPLUS_AFQMC_OBSERVABLES_CONFIG_HPP

namespace qmcplusplus
{
namespace afqmc
{
enum observables
{
  GFockOpa, // numerator of Extended Koopman's Theorem approach
  GFockOpb, // numerator of Extended Koopman's Theorem approach
  OneRDMFull,
  TwoRDMFull,
  OneRDMc,
  TwoRDMc,
  OneRSDMFull,
  TwoRSDMFull,
  OneRSDMc,
  TwoRSDMc,
  energy,
  force
}; // 12 observables
std::array<std::string, 12> hdf_ids = {"GFockOpa_",
                                       "GFockOpb_",
                                       "one_rdm_",
                                       "two_rdm_",
                                       "contracted_one_rdm_",
                                       "contracted_two_rdm_",
                                       "one_rsdm_",
                                       "two_rsdm_",
                                       "contracted_one_rsdm_",
                                       "contracted_two_rsdm_",
                                       "bp_energy_",
                                       "bp_force_"};


} // namespace afqmc
} // namespace qmcplusplus


#endif
