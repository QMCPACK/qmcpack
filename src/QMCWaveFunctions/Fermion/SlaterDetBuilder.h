//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LCORBITALSETBUILDER_H
#define QMCPLUSPLUS_LCORBITALSETBUILDER_H

#include <vector>
#include "Configuration.h"
#include "WaveFunctionComponentBuilder.h"
#include <hdf/hdf_archive.h>

namespace qmcplusplus
{
class TrialWaveFunction;
class BackflowTransformation;
class DiracDeterminantBase;
class MultiDiracDeterminant;
class SPOSetBuilder;
class SPOSetBuilderFactory;
struct ci_configuration;

/** derived class from WaveFunctionComponentBuilder
 *
 * Builder SlaterDeterminant with LCOrbitalSet
 */
class SlaterDetBuilder : public WaveFunctionComponentBuilder
{
public:
  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  SlaterDetBuilder(Communicate* comm,
                   SPOSetBuilderFactory& factory,
                   ParticleSet& els,
                   TrialWaveFunction& psi,
                   const PSetMap& psets);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   *
   */
  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  /// reference to the sposet_builder_factory, should be const once the legacy input style is removed
  SPOSetBuilderFactory& sposet_builder_factory_;
  ///reference to TrialWaveFunction, should go away as the CUDA code.
  TrialWaveFunction& targetPsi;
  ///reference to a PSetMap
  const PSetMap& ptclPool;

  /** process a determinant element
   * @param cur xml node
   * @param spin_group the spin group of the created determinant
   * @return legacy_input_sposet_builder an sposet builder to handle legacy input
   * @return BFTrans backflow transformations
   */
  std::unique_ptr<DiracDeterminantBase> putDeterminant(
      xmlNodePtr cur,
      int spin_group,
      const std::unique_ptr<SPOSetBuilder>& legacy_input_sposet_builder,
      const std::unique_ptr<BackflowTransformation>& BFTrans);

  bool createMSDFast(std::vector<std::unique_ptr<MultiDiracDeterminant>>& Dets,
                     std::vector<std::vector<size_t>>& C2node,
                     std::vector<ValueType>& C,
                     std::vector<ValueType>& CSFcoeff,
                     std::vector<size_t>& DetsPerCSF,
                     std::vector<RealType>& CSFexpansion,
                     bool& usingCSF,
                     opt_variables_type& myVars,
                     bool& Optimizable,
                     bool& CI_Optimizable,
                     xmlNodePtr cur) const;


  bool readDetList(xmlNodePtr cur,
                   std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                   std::vector<std::vector<size_t>>& C2nodes,
                   std::vector<std::string>& CItags,
                   std::vector<ValueType>& coeff,
                   bool& optimizeCI,
                   std::vector<int>& nptcls,
                   std::vector<ValueType>& CSFcoeff,
                   std::vector<size_t>& DetsPerCSF,
                   std::vector<RealType>& CSFexpansion,
                   bool& usingCSF) const;

  bool readDetListH5(xmlNodePtr cur,
                     std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                     std::vector<std::vector<size_t>>& C2nodes,
                     std::vector<std::string>& CItags,
                     std::vector<ValueType>& coeff,
                     bool& optimizeCI,
                     std::vector<int>& nptcls) const;

  template<typename VT,
           std::enable_if_t<(std::is_same<VT, ValueType>::value) && (std::is_floating_point<VT>::value), int> = 0>
  void readCoeffs(hdf_archive& hin, std::vector<VT>& ci_coeff, size_t n_dets, int ext_level) const
  {
    ///Determinant coeffs are stored in Coeff for the ground state and Coeff_N
    ///for the Nth excited state.
    ///The Ground State is always stored in Coeff
    ///Backward compatibility is insured
    std::string extVar;
    if (ext_level == 0)
      extVar = "Coeff";
    else
      extVar = "Coeff_" + std::to_string(ext_level);

    if (!hin.readEntry(ci_coeff, extVar))
      throw std::runtime_error("Could not read CI coefficients from HDF5");
  }

  template<typename VT,
           std::enable_if_t<(std::is_same<VT, ValueType>::value) &&
                                (std::is_same<VT, std::complex<typename VT::value_type>>::value),
                            int> = 0>
  void readCoeffs(hdf_archive& hin, std::vector<VT>& ci_coeff, size_t n_dets, int ext_level) const
  {
    std::string extVar;
    std::vector<double> CIcoeff_real;
    std::vector<double> CIcoeff_imag;
    CIcoeff_imag.resize(n_dets);
    CIcoeff_real.resize(n_dets);
    fill(CIcoeff_imag.begin(), CIcoeff_imag.end(), 0.0);
    ///Determinant coeffs are stored in Coeff_N where N is Nth excited state.
    ///The Ground State is always stored in Coeff.
    ///Backward compatibility is insured

    std::string ext_var;
    if (ext_level == 0)
      extVar = "Coeff";
    else
      extVar = "Coeff_" + std::to_string(ext_level);


    if (!hin.readEntry(CIcoeff_real, extVar))
      throw std::runtime_error("Could not read CI coefficients from HDF5");

    extVar = extVar + "_imag";
    if (!hin.readEntry(CIcoeff_imag, extVar))
      app_log() << "Coeff_imag not found in h5. Set to zero." << std::endl;

    for (size_t i = 0; i < n_dets; i++)
      ci_coeff[i] = VT(CIcoeff_real[i], CIcoeff_imag[i]);
  }
};

} // namespace qmcplusplus
#endif
