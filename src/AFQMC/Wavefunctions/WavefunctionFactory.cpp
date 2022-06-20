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

#include <random>

#include "hdf/hdf_archive.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Wavefunctions/NOMSD.hpp"
#include "AFQMC/Wavefunctions/PHMSD.hpp"
#include "AFQMC/HamiltonianOperations/HamOpsIO.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
Wavefunction WavefunctionFactory::fromASCII(TaskGroup_& TGprop,
                                            TaskGroup_& TGwfn,
                                            xmlNodePtr cur,
                                            WALKER_TYPES walker_type,
                                            Hamiltonian& h,
                                            RealType cutvn,
                                            int targetNW)
{
  if (cur == NULL)
    APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

  std::string type("msd");
  std::string info("info0");
  std::string init_type("");
  std::string name("");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(type, "type");
  oAttrib.add(info, "info");
  oAttrib.add(init_type, "init");
  oAttrib.add(name, "name");
  oAttrib.put(cur);

  std::transform(type.begin(), type.end(), type.begin(), (int (*)(int))tolower);
  std::transform(init_type.begin(), init_type.end(), init_type.begin(), (int (*)(int))tolower);

  if (InfoMap.find(info) == InfoMap.end())
  {
    app_error() << "ERROR: Undefined info in WavefunctionFactory. \n";
    APP_ABORT("ERROR: Undefined info in WavefunctionFactory. \n");
  }

  RealType cutv2(0.);
  int ndets_to_read(-1); // if not set, read the entire file
  int initial_configuration = 0;
  double randomize_guess(0.0);
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  std::string write_trial_density_matrix("");
  ParameterSet m_param;
  m_param.add(filename, "filename");
  m_param.add(restart_file, "restart_file");
  m_param.add(write_trial_density_matrix, "trial_density_matrix");
  m_param.add(cutv2, "cutoff");
  m_param.add(ndets_to_read, "ndet");
  m_param.add(initial_configuration, "initial_configuration");
  m_param.add(randomize_guess, "randomize_guess");
  m_param.put(cur);

  AFQMCInfo& AFinfo = InfoMap[info];
  auto NCE          = h.getNuclearCoulombEnergy();

  int NMO  = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;
  int NPOL = (walker_type == NONCOLLINEAR) ? 2 : 1;
  if ((walker_type == NONCOLLINEAR) && (NAEB != 0))
    APP_ABORT(" Error in Wavefunctions/WavefunctionFactory::fromASCII: noncollinear && NAEB!=0. \n\n\n ");

  std::ifstream in;
  in.open(filename.c_str());
  if (in.fail())
  {
    app_error() << "Problems opening file:  " << filename << std::endl;
    APP_ABORT("Problems opening file. \n");
  }
  std::string wfn_type = getWfnType(in);

  using Alloc = shared_allocator<ComplexType>;
  if (type == "msd" || type == "nomsd")
  {
    app_log() << " Wavefunction type: NOMSD\n";

    std::vector<ComplexType> ci;
    std::vector<PsiT_Matrix> PsiT;
    if (wfn_type == "occ" || wfn_type == "mixed")
    {
      std::vector<PsiT_Matrix> PsiT_MO; // read_ph_wavefunction returns the full MO matrix
      // careful here!!!!
      // number of terms in PsiT_MO depend on RHF/UHF type, not on walker_type!!!
      ph_excitations<int, ComplexType> abij =
          read_ph_wavefunction(in, ndets_to_read, walker_type, TGwfn.Node(), NMO, NAEA, NAEB, PsiT_MO);
      assert(abij.number_of_configurations() == ndets_to_read);
      int NEL = (walker_type == NONCOLLINEAR) ? (NAEA + NAEB) : NAEA;
      int N_  = NPOL * NMO;
      ComplexType one(1.0, 0.0);
      if (walker_type == COLLINEAR)
        PsiT.reserve(2 * ndets_to_read);
      else
        PsiT.reserve(ndets_to_read);
      ci.reserve(ndets_to_read);
      auto refc = abij.reference_configuration();
      // add reference
      ci.emplace_back(std::get<2>(*abij.configurations_begin()));
      if (wfn_type == "occ")
        PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL, N_}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
      else
        PsiT.emplace_back(
            PsiT_Matrix(tp_ul_ul{NEL, N_}, tp_ul_ul{0, 0}, get_nnz(PsiT_MO[0], refc, NEL, 0), Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        if (wfn_type == "occ")
        {
          for (int k = 0; k < NEL; k++)
            PsiT.back().emplace_back({k, *(refc + k)}, one);
        }
        else
        {
          for (int k = 0; k < NEL; k++)
          {
            size_t ki = *(refc + k); // occupied state #k
            auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
            auto val  = PsiT_MO[0].non_zero_values_data(ki);
            for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
              PsiT.back().emplace_back({k, *col}, *val);
          }
        }
      }
      if (walker_type == COLLINEAR)
      {
        if (wfn_type == "occ")
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB, NMO}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
        else
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB, NMO}, tp_ul_ul{0, 0},
                                        get_nnz(PsiT_MO.back(), refc + NAEA, NAEB, NMO), Alloc(TGwfn.Node())));
        if (TGwfn.Node().root())
        {
          if (wfn_type == "occ")
          {
            for (int k = 0; k < NAEB; k++)
              PsiT.back().emplace_back({k, *(refc + NAEA + k) - NMO}, one);
          }
          else
          {
            for (int k = 0; k < NAEB; k++)
            {
              size_t ki = *(refc + NAEA + k) - NMO; // occupied state #k
              // this should be correct
              auto col = PsiT_MO.back().non_zero_indices2_data(ki);
              auto val = PsiT_MO.back().non_zero_values_data(ki);
              for (size_t ic = 0, icend = PsiT_MO.back().num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
                PsiT.back().emplace_back({k, *col}, *val);
            }
          }
        }
      }
      // work array
      std::vector<int> iwork(NAEA);
      auto configurations = abij.configurations_begin() + 1;
      for (; configurations < abij.configurations_end(); ++configurations)
      {
        ci.emplace_back(std::get<2>(*configurations));
        abij.get_alpha_configuration(std::get<0>(*configurations), iwork);
        if (wfn_type == "occ")
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL, N_}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
        else
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL, N_}, tp_ul_ul{0, 0}, get_nnz(PsiT_MO[0], iwork.data(), NEL, 0),
                                        Alloc(TGwfn.Node())));
        if (TGwfn.Node().root())
        {
          // add excited configuration
          if (wfn_type == "occ")
          {
            for (int k = 0; k < NEL; k++)
              PsiT.back().emplace_back({k, iwork[k]}, one);
          }
          else
          {
            for (int k = 0; k < NEL; k++)
            {
              size_t ki = iwork[k]; // occupied state #k
              auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
              auto val  = PsiT_MO[0].non_zero_values_data(ki);
              for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
                PsiT.back().emplace_back({k, *col}, *val);
            }
          }
        }
        if (walker_type == COLLINEAR)
        {
          abij.get_beta_configuration(std::get<1>(*configurations), iwork);
          if (wfn_type == "occ")
            PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB, NMO}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
          else
            PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB, NMO}, tp_ul_ul{0, 0},
                                          get_nnz(PsiT_MO.back(), iwork.data(), NAEB, NMO), Alloc(TGwfn.Node())));
          if (TGwfn.Node().root())
          {
            if (wfn_type == "occ")
            {
              for (int k = 0; k < NAEB; k++)
                PsiT.back().emplace_back({k, iwork[k] - NMO}, one);
            }
            else
            {
              for (int k = 0; k < NAEB; k++)
              {
                size_t ki = iwork[k] - NMO; // occupied state #k
                auto col  = PsiT_MO.back().non_zero_indices2_data(ki);
                auto val  = PsiT_MO.back().non_zero_values_data(ki);
                for (size_t ic = 0, icend = PsiT_MO.back().num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
                  PsiT.back().emplace_back({k, *col}, *val);
              }
            }
          }
        }
      }
    }
    else if (wfn_type == "matrix")
    {
      read_general_wavefunction(in, ndets_to_read, walker_type, TGwfn.Node(), NMO, NAEA, NAEB, PsiT, ci);
    }
    else
    {
      APP_ABORT("Error: Unknown wfn_type in WavefunctionFactory with MSD wavefunction.\n");
    }
    TGwfn.node_barrier();
    // multideterminant NOMSD needs
    auto HOps(getHamOps(restart_file, walker_type, NMO, NAEA, NAEB, PsiT,
      TGprop, TGwfn, cutvn, cutv2, ndets_to_read, h));
    TGwfn.node_barrier();
    // add initial_guess
    auto guess = initial_guess.find(name);
    if (guess == initial_guess.end())
    {
      auto newg =
          initial_guess.insert(std::make_pair(name, boost::multi::array<ComplexType, 3>({2, NPOL * NMO, NAEA})));
      int iC = (walker_type != COLLINEAR ? initial_configuration : 2 * initial_configuration);
      if (iC >= PsiT.size())
        APP_ABORT(" Error: initial_configuration > ndets_to_read \n");
      if (!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n");
      using ma::conj;
      std::fill_n((newg.first)->second.origin(), 2 * NPOL * NMO * NAEA, ComplexType(0.0, 0.0));
      {
        auto pbegin = PsiT[iC].pointers_begin();
        auto pend   = PsiT[iC].pointers_end();
        auto p0     = pbegin[0];
        auto v0     = PsiT[iC].non_zero_values_data();
        auto c0     = PsiT[iC].non_zero_indices2_data();
        for (int i = 0; i < PsiT[iC].size(0); i++)
          for (int ip = pbegin[i]; ip < pend[i]; ip++)
          {
            ((newg.first)->second)[0][c0[ip - p0]][i] = ma::conj(v0[ip - p0]);
          }
      }
      if (walker_type == COLLINEAR)
      {
        auto pbegin = PsiT[iC + 1].pointers_begin();
        auto pend   = PsiT[iC + 1].pointers_end();
        auto p0     = pbegin[0];
        auto v0     = PsiT[iC + 1].non_zero_values_data();
        auto c0     = PsiT[iC + 1].non_zero_indices2_data();
        for (int i = 0; i < PsiT[iC + 1].size(0); i++)
          for (int ip = pbegin[i]; ip < pend[i]; ip++)
          {
            ((newg.first)->second)[1][c0[ip - p0]][i] = ma::conj(v0[ip - p0]);
          }
      }
    }
    else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");

    if (TGwfn.TG_local().size() > 1)
    {
      SlaterDetOperations SDetOp(SlaterDetOperations_shared<ComplexType>(NPOL * NMO, NAEA));
      return Wavefunction(NOMSD<devcsr_Matrix>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                               std::move(PsiT), walker_type, NCE, targetNW));
    }
    else
    {
      SlaterDetOperations SDetOp(
          SlaterDetOperations_serial<ComplexType, DeviceBufferManager>(NPOL * NMO, NAEA, DeviceBufferManager{}));
      return Wavefunction(NOMSD<devcsr_Matrix>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                               std::move(PsiT), walker_type, NCE, targetNW));
    }
  }
  else if (type == "phmsd")
  {
    app_log() << " Wavefunction type: PHMSD\n";

    /* Implementation notes:
     *  - PsiT: [Nact, NMO] where Nact is the number of active space orbitals,
     *                     those that participate in the ci expansion
     *  - The half rotation is done with respect to the supermatrix PsiT
     *  - Need to calculate Nact and create a mapping from orbital index to actice space index.
     *    Those orbitals in the corresponding virtual space (not in active) map to -1 as a precaution.
     */

    // assuming walker_type==COLLINEAR for now, specialize a type for perfect pairing PHMSD
    // MAM: generatlize for NONCOLLINEAR later!!!
    if (walker_type != COLLINEAR)
      APP_ABORT("Error: PHMSD requires a COLLINEAR calculation.\n");
    std::vector<PsiT_Matrix> PsiT_MO;
    ph_excitations<int, ComplexType> abij =
        read_ph_wavefunction(in, ndets_to_read, walker_type, TGwfn.Node(), NMO, NAEA, NAEB, PsiT_MO);
    int N_ = (walker_type == NONCOLLINEAR) ? 2 * NMO : NMO;
    if (wfn_type == "occ")
    {
      //build PsiT_MO
      ComplexType one(1.0, 0.0);
      // wfn_type == "occ" implies a single reference now, since integrals can't be UHF
      PsiT_MO.reserve(1);
      PsiT_MO.emplace_back(PsiT_Matrix(tp_ul_ul{N_, N_}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
        for (int k = 0; k < N_; k++)
          PsiT_MO.back().emplace_back({k, k}, one);
    }
    else if (wfn_type == "mixed")
    {
      // nothing to do
    }
    else if (wfn_type == "matrix")
    {
      APP_ABORT("Error: wfn_type=matrix not allowed in WavefunctionFactory with PHMSD wavefunction.\n");
    }
    else
    {
      APP_ABORT("Error: Unknown wfn_type in WavefunctionFactory with MSD wavefunction.\n");
    }
    TGwfn.node_barrier();

    // find active space orbitals and create super trial matrix PsiT
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(PsiT_MO.size());
    // expect mapped over range [0-2*NMO], but alpha and beta sectors with 0-based active indexes
    std::map<int, int> mo2active(find_active_space(PsiT_MO.size() == 1, abij, NMO, NAEA, NAEB));
    std::map<int, int> acta2mo;
    std::map<int, int> actb2mo;
    std::vector<int> active_alpha;
    std::vector<int> active_beta;
    std::vector<int> active_combined;
    for (int i = 0; i < NMO; i++)
    {
      if (mo2active[i] >= 0)
      {
        active_alpha.push_back(i);
        acta2mo[mo2active[i]] = i;
      }
      if (mo2active[i + NMO] >= 0)
      {
        active_beta.push_back(i);
        actb2mo[mo2active[i + NMO]] = i + NMO;
      }
      if (mo2active[i] >= 0 || mo2active[i + NMO] >= 0)
        active_combined.push_back(i);
    }
    if (PsiT_MO.size() == 1)
    {
      // RHF reference
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_combined.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[0], active_combined.data(), active_combined.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_combined.size(); k++)
        {
          size_t ki = active_combined[k]; // occupied state #k
          auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[0].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT[0].emplace_back({k, *col}, *val);
        }
      }
    }
    else
    {
      // UHF reference
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_alpha.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[0], active_alpha.data(), active_alpha.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_alpha.size(); k++)
        {
          size_t ki = active_alpha[k]; // occupied state #k
          auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[0].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT.back().emplace_back({k, *col}, *val);
        }
      }
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_beta.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[1], active_beta.data(), active_beta.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_beta.size(); k++)
        {
          size_t ki = active_beta[k]; // occupied state #k
          auto col  = PsiT_MO[1].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[1].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[1].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT[1].emplace_back({k, *col}, *val);
        }
      }
    }
    // now that mappings have been constructed, map indexes of excited state orbitals
    // to the corresponding active space indexes
    if (TGwfn.Node().root())
    {
      // map reference
      auto refc = abij.reference_configuration();
      for (int i = 0; i < NAEA + NAEB; i++, ++refc)
        *refc = mo2active[*refc];
      for (int n = 1; n < abij.maximum_excitation_number()[0]; n++)
      {
        auto it  = abij.alpha_begin(n);
        auto ite = abij.alpha_end(n);
        for (; it < ite; ++it)
        {
          auto exct = (*it) + n; // only need to map excited state indexes
          for (int np = 0; np < n; ++np, ++exct)
            *exct = mo2active[*exct];
        }
      }
      for (int n = 1; n < abij.maximum_excitation_number()[1]; n++)
      {
        auto it  = abij.beta_begin(n);
        auto ite = abij.beta_end(n);
        for (; it < ite; ++it)
        {
          auto exct = (*it) + n; // only need to map excited state indexes
          for (int np = 0; np < n; ++np, ++exct)
            *exct = mo2active[*exct];
        }
      }
    }
    // is PureSD actually faster??? CHECK!!!
    // NOTE: For UHF reference, treat HOps as a 2 determinant wavefunction of a
    //        CLOSED walker type. This way you can treat alpha/beta sectors independently in
    //        HOps through the index corresponding 0/1.
    // never add coulomb to half rotated v2 tensor in PHMSD
    TGwfn.node_barrier();
    auto HOps(getHamOps(restart_file, CLOSED, NMO, NAEA, NAEB, PsiT,
      TGprop, TGwfn, cutvn, cutv2, false, h));
    TGwfn.node_barrier();
    // add initial_guess
    // when propagating Nact states, change this here
    auto guess = initial_guess.find(name);
    if (guess == initial_guess.end())
    {
      auto newg = initial_guess.insert(std::make_pair(name, boost::multi::array<ComplexType, 3>({2, NMO, NAEA})));
      if (!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n");
      auto& Psi0((newg.first)->second);
      randomize_guess = std::abs(randomize_guess);
      if (randomize_guess > 1e-12)
        app_log() << " Randomizing initial guess with uniform distribution: " << randomize_guess << std::endl;
      std::default_random_engine generator(777);
      std::uniform_real_distribution<double> distribution(-randomize_guess, randomize_guess);
      int iC = initial_configuration;
      if (iC >= abij.number_of_configurations())
        APP_ABORT(" Error: initial_configuration > ndets \n");
      using ma::conj;
      std::fill_n((newg.first)->second.origin(), 2 * NMO * NAEA, ComplexType(0.0, 0.0));
      //auto refc = abij.reference_configuration();
      {
        std::vector<int> alphaC(NAEA);
        abij.get_alpha_configuration(std::get<0>(*(abij.configurations_begin() + iC)), alphaC);
        auto pbegin = PsiT[0].pointers_begin();
        auto pend   = PsiT[0].pointers_end();
        auto p0     = pbegin[0];
        auto v0     = PsiT[0].non_zero_values_data();
        auto c0     = PsiT[0].non_zero_indices2_data();
        // only takinf NAEA states, increase later if super SM is needed
        for (int i = 0; i < NAEA; i++)
        {
          //int ik = *(refc+i);
          int ik = alphaC[i];
          for (int ip = pbegin[ik]; ip < pend[ik]; ip++)
            Psi0[0][c0[ip - p0]][i] = ma::conj(v0[ip - p0]);
        }
        if (randomize_guess > 1e-12)
          for (int i = 0; i < NMO; i++)
            for (int a = 0; a < NAEA; a++)
              Psi0[0][i][a] += distribution(generator);
      }
      if (walker_type == COLLINEAR)
      {
        std::vector<int> betaC(NAEB);
        abij.get_beta_configuration(std::get<1>(*(abij.configurations_begin() + iC)), betaC);
        auto pbegin = PsiT.back().pointers_begin();
        auto pend   = PsiT.back().pointers_end();
        auto p0     = pbegin[0];
        auto v0     = PsiT.back().non_zero_values_data();
        auto c0     = PsiT.back().non_zero_indices2_data();
        // only takinf NAEB states, increase later if super SM is needed
        for (int i = 0; i < NAEB; i++)
        {
          //int ik = *(refc+NAEA+i);
          int ik = betaC[i];
          for (int ip = pbegin[ik]; ip < pend[ik]; ip++)
            Psi0[1][c0[ip - p0]][i] = ma::conj(v0[ip - p0]);
        }
        if (randomize_guess > 1e-12)
          for (int i = 0; i < NMO; i++)
            for (int a = 0; a < NAEB; a++)
              Psi0[1][i][a] += distribution(generator);
      }
    }
    else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");

    // setup configuration coupligs
    using index_aos = ma::sparse::array_of_sequences<int, int, shared_allocator<int>, ma::sparse::is_root>;
    //    std::allocator<ComplexType> alloc_{}; //shared_allocator<ComplexType>;
    shared_allocator<int> alloc_{TGwfn.Node()};

    // alpha
    std::vector<int> counts_alpha(abij.number_of_unique_excitations()[0]);
    std::vector<int> counts_beta(abij.number_of_unique_excitations()[1]);
    if (TGwfn.Node().root())
    {
      for (auto it = abij.configurations_begin(); it < abij.configurations_end(); ++it)
      {
        ++counts_alpha[std::get<0>(*it)];
        ++counts_beta[std::get<1>(*it)];
      }
    }
    TGwfn.Node().broadcast_n(counts_alpha.begin(), counts_alpha.size());
    TGwfn.Node().broadcast_n(counts_beta.begin(), counts_beta.size());
    index_aos beta_coupled_to_unique_alpha(counts_alpha.size(), counts_alpha, alloc_);
    index_aos alpha_coupled_to_unique_beta(counts_beta.size(), counts_beta, alloc_);
    if (TGwfn.Node().root())
    {
      int ni = 0;
      for (auto it = abij.configurations_begin(); it < abij.configurations_end(); ++it, ++ni)
      {
        beta_coupled_to_unique_alpha.emplace_back(std::get<0>(*it), ni);
        alpha_coupled_to_unique_beta.emplace_back(std::get<1>(*it), ni);
      }
    }
    TGwfn.Node().barrier();

    //return Wavefunction{};
    return Wavefunction(PHMSD(AFinfo, cur, TGwfn, std::move(HOps), std::move(acta2mo), std::move(actb2mo),
                              std::move(abij), std::move(beta_coupled_to_unique_alpha),
                              std::move(alpha_coupled_to_unique_beta), std::move(PsiT), walker_type, NCE, targetNW));
  }
  else if (type == "generalmsd")
  {
    app_error() << " Error: Wavefunction type GeneralMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type GeneralMSD not yet implemented. \n");
    return Wavefunction{};
  }
  else
  {
    app_error() << " Error: Unknown wave-function type: " << type << std::endl;
    APP_ABORT(" Error: Unknown wave-function type. \n");
    return Wavefunction{};
  }
}

Wavefunction WavefunctionFactory::fromHDF5(TaskGroup_& TGprop,
                                           TaskGroup_& TGwfn,
                                           xmlNodePtr cur,
                                           WALKER_TYPES walker_type,
                                           Hamiltonian& h,
                                           RealType cutvn,
                                           int targetNW)
{
  if (cur == NULL)
    APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

  std::string type("MSD");
  std::string info("info0");
  std::string init_type("");
  std::string name("");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(type, "type");
  oAttrib.add(info, "info");
  oAttrib.add(init_type, "init");
  oAttrib.add(name, "name");
  oAttrib.put(cur);

  std::transform(type.begin(), type.end(), type.begin(), (int (*)(int))tolower);
  std::transform(init_type.begin(), init_type.end(), init_type.begin(), (int (*)(int))tolower);

  if (InfoMap.find(info) == InfoMap.end())
  {
    app_error() << "ERROR: Undefined info in WavefunctionFactory. \n";
    APP_ABORT("ERROR: Undefined info in WavefunctionFactory. \n");
  }

  RealType cutv2(0.);
  int ndets_to_read(-1); // if not set, read the entire file
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  std::string rediag("");
  std::string dense_trial("");
  ParameterSet m_param;
  m_param.add(filename, "filename");
  m_param.add(restart_file, "restart_file");
  m_param.add(cutv2, "cutoff");
  m_param.add(rediag, "rediag");
  m_param.add(ndets_to_read, "ndet");
  m_param.add(dense_trial, "dense_trial");
  m_param.put(cur);
  bool recompute_ci = false;
  std::transform(rediag.begin(), rediag.end(), rediag.begin(), (int (*)(int))tolower);
  if (rediag == "yes" || rediag == "true")
    recompute_ci = true;

  AFQMCInfo& AFinfo = InfoMap[info];
  ValueType NCE     = 0.0;

  int NMO  = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;
  int NPOL = (walker_type == NONCOLLINEAR) ? 2 : 1;
  if ((walker_type == NONCOLLINEAR) && (NAEB != 0))
    APP_ABORT(" Error in Wavefunctions/WavefunctionFactory::fromASCII: noncollinear && NAEB!=0. \n\n\n ");

  std::vector<int> excitations;

  using Alloc = shared_allocator<ComplexType>;
  // HOps, ci, PsiT, NCE
  hdf_archive dump(TGwfn.Global());
  if (!dump.open(filename, H5F_ACC_RDONLY))
  {
    app_error() << " Error hdf5 file in WavefunctionFactory. \n";
    APP_ABORT("");
  }
  if (!dump.push("Wavefunction", false))
  {
    app_error() << " Error in WavefunctionFactory: Group Wavefunction not found. \n";
    APP_ABORT("");
  }

  if (type == "msd" || type == "nomsd")
  {
    app_log() << " Wavefunction type: NOMSD" << std::endl;
    if (!dump.push("NOMSD", false))
    {
      app_error() << " Error in WavefunctionFactory: Group NOMSD not found.\n";
      APP_ABORT("");
    }
    std::vector<ComplexType> ci;

    // Read common trial wavefunction input options.
    WALKER_TYPES input_wtype;
    getCommonInput(dump, NMO, NAEA, NAEB, ndets_to_read, ci, input_wtype, TGwfn.Global().root());
    NCE = h.getNuclearCoulombEnergy();

    TGwfn.Global().broadcast_n(ci.data(), ci.size());
    TGwfn.Global().broadcast_value(NCE);

    // Create Trial wavefunction.
    int nd     = (walker_type == COLLINEAR ? 2 * ndets_to_read : ndets_to_read);
    int ndread = nd;
    if (walker_type == COLLINEAR and input_wtype == CLOSED)
      ndread = ndets_to_read;
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(nd);
    using Alloc = shared_allocator<ComplexType>;
    for (int i = 0; i < ndread; ++i)
    {
      if (!dump.push(std::string("PsiT_") + std::to_string(i), false))
      {
        app_error() << " Error in WavefunctionFactory: Group PsiT not found. \n";
        APP_ABORT("");
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, Alloc>(dump, TGwfn.Node())); //,Alloc(TGwfn.Node())));
      dump.pop();
      if (walker_type == COLLINEAR and input_wtype == CLOSED)
      {
        if (NAEA != NAEB)
          APP_ABORT(" Error: NAEA!=NAEB when initializing collinear wfn from closed shell file.\n");
        // read them again
        if (!dump.push(std::string("PsiT_") + std::to_string(i), false))
        {
          app_error() << " Error in WavefunctionFactory: Group PsiT not found. \n";
          APP_ABORT("");
        }
        PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, Alloc>(dump, TGwfn.Node())); //,Alloc(TGwfn.Node())));
        dump.pop();
      }
    }

    // Set initial walker's Slater matrix.
    getInitialGuess(dump, name, NMO, NAEA, NAEB, walker_type);

    // Restore Hamiltonian Operations Object from file if it exists.
    TGwfn.node_barrier();
    auto HOps(getHamOps(restart_file, walker_type, NMO, NAEA, NAEB, PsiT,
      TGprop, TGwfn, cutvn, cutv2, ndets_to_read, h));
    TGwfn.node_barrier();

    // if not set, get default based on HamTYpe
    // use sparse trial only on KP runs
    if (dense_trial == "")
    {
      dense_trial = std::string("yes");
      if (HOps.getHamType() == KPFactorized || HOps.getHamType() == KPTHC)
        dense_trial = std::string("no");
    }

    if (TGwfn.TG_local().size() > 1)
    {
      SlaterDetOperations SDetOp(SlaterDetOperations_shared<ComplexType>(NPOL * NMO, NAEA));
      if (dense_trial == "yes")
      {
        using MType = ComplexMatrix<node_allocator<ComplexType>>;
        std::vector<MType> PsiT_;
        PsiT_.reserve(PsiT.size());
        auto alloc_shared_(make_node_allocator<ComplexType>(TGwfn));
        for (auto& v : PsiT)
        {
          PsiT_.emplace_back(MType({static_cast<MType::size_type>(v.size(0)), static_cast<MType::size_type>(v.size(1))}, alloc_shared_));
          ma::Matrix2MAREF('N', v, PsiT_.back());
        }
        return Wavefunction(NOMSD<MType>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                         std::move(PsiT_), walker_type, NCE, targetNW));
      }
      else
      {
        return Wavefunction(NOMSD<devcsr_Matrix>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                                 std::move(PsiT), walker_type, NCE, targetNW));
      }
    }
    else
    {
      SlaterDetOperations SDetOp(
          SlaterDetOperations_serial<ComplexType, DeviceBufferManager>(NPOL * NMO, NAEA, DeviceBufferManager{}));
      if (dense_trial == "yes")
      {
        using MType = ComplexMatrix<node_allocator<ComplexType>>;
        std::vector<MType> PsiT_;
        PsiT_.reserve(PsiT.size());
        auto alloc_shared_(make_node_allocator<ComplexType>(TGwfn));
        for (auto& v : PsiT)
        {
          PsiT_.emplace_back(MType({static_cast<MType::size_type>(v.size(0)), static_cast<MType::size_type>(v.size(1))}, alloc_shared_));
          ma::Matrix2MAREF('N', v, PsiT_.back());
        }
        return Wavefunction(NOMSD<MType>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                         std::move(PsiT_), walker_type, NCE, targetNW));
      }
      else
      {
        return Wavefunction(NOMSD<devcsr_Matrix>(AFinfo, cur, TGwfn, std::move(SDetOp), std::move(HOps), std::move(ci),
                                                 std::move(PsiT), walker_type, NCE, targetNW));
      }
    }
  }
  else if (type == "phmsd")
  {
    app_log() << " Wavefunction type: PHMSD" << std::endl;

    /* Implementation notes:
     *  - PsiT: [Nact, NMO] where Nact is the number of active space orbitals,
     *                     those that participate in the ci expansion
     *  - The half rotation is done with respect to the supermatrix PsiT
     *  - Need to calculate Nact and create a mapping from orbital index to actice space index.
     *    Those orbitals in the corresponding virtual space (not in active) map to -1 as a precaution.
     */

    // assuming walker_type==COLLINEAR for now, specialize a type for perfect pairing PHMSD
    if (walker_type != COLLINEAR)
      APP_ABORT("Error: PHMSD requires a COLLINEAR calculation.\n");
    std::vector<PsiT_Matrix> PsiT_MO;
    std::string wfn_type;
    if (!dump.push("PHMSD", false))
    {
      app_error() << " Error in WavefunctionFactory: Group PHMSD not found. \n";
      APP_ABORT("");
    }
    std::vector<int> occbuff;
    std::vector<ComplexType> coeffs;
    // 1. Read occupancies and coefficients.
    app_log() << " Reading PHMSD wavefunction from " << filename << std::endl;
    read_ph_wavefunction_hdf(dump, coeffs, occbuff, ndets_to_read, walker_type, TGwfn.Node(), NMO, NAEA, NAEB, PsiT_MO,
                             wfn_type);
    app_log() << " Finished reading PHMSD wavefunction " << std::endl;
    boost::multi::array_ref<int, 2> occs(to_address(occbuff.data()), {ndets_to_read, NAEA + NAEB});
    // 2. Compute Variational Energy / update coefficients
    app_log() << " Computing variational energy of trial wavefunction." << std::endl;
    computeVariationalEnergyPHMSD(TGwfn, h, occs, coeffs, ndets_to_read, NAEA, NAEB, NMO, recompute_ci);
    app_log() << " Finished computing variational energy of trial wavefunction." << std::endl;
    // 3. Construct Structures.
    ph_excitations<int, ComplexType> abij = build_ph_struct(coeffs, occs, ndets_to_read, TGwfn.Node(), NMO, NAEA, NAEB);
    int N_                                = (walker_type == NONCOLLINEAR) ? 2 * NMO : NMO;
    if (wfn_type == "occ")
    {
      //build PsiT_MO
      ComplexType one(1.0, 0.0);
      // wfn_type == "occ" implies a single reference now, since integrals can't be UHF
      PsiT_MO.reserve(1);
      PsiT_MO.emplace_back(PsiT_Matrix(tp_ul_ul{N_, N_}, tp_ul_ul{0, 0}, 1, Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
        for (int k = 0; k < N_; k++)
          PsiT_MO.back().emplace_back({k, k}, one);
    }
    else if (wfn_type == "mixed")
    {
      // nothing to do
    }
    else if (wfn_type == "matrix")
    {
      APP_ABORT("Error: wfn_type=matrix not allowed in WavefunctionFactory with PHMSD wavefunction.\n");
    }
    else
    {
      APP_ABORT("Error: Unknown wfn_type in WavefunctionFactory with MSD wavefunction.\n");
    }
    TGwfn.node_barrier();

    // find active space orbitals and create super trial matrix PsiT
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(PsiT_MO.size());
    // expect mapped over range [0-2*NMO], but alpha and beta sectors with 0-based active indexes
    std::map<int, int> mo2active(find_active_space(PsiT_MO.size() == 1, abij, NMO, NAEA, NAEB));
    std::map<int, int> acta2mo;
    std::map<int, int> actb2mo;
    std::vector<int> active_alpha;
    std::vector<int> active_beta;
    std::vector<int> active_combined;
    for (int i = 0; i < NMO; i++)
    {
      if (mo2active[i] >= 0)
      {
        active_alpha.push_back(i);
        acta2mo[mo2active[i]] = i;
      }
      if (mo2active[i + NMO] >= 0)
      {
        active_beta.push_back(i);
        actb2mo[mo2active[i + NMO]] = i + NMO;
      }
      if (mo2active[i] >= 0 || mo2active[i + NMO] >= 0)
        active_combined.push_back(i);
    }
    if (PsiT_MO.size() == 1)
    {
      // RHF reference
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_combined.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[0], active_combined.data(), active_combined.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_combined.size(); k++)
        {
          size_t ki = active_combined[k]; // occupied state #k
          auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[0].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT[0].emplace_back({k, *col}, *val);
        }
      }
    }
    else
    {
      // UHF reference
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_alpha.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[0], active_alpha.data(), active_alpha.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_alpha.size(); k++)
        {
          size_t ki = active_alpha[k]; // occupied state #k
          auto col  = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[0].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[0].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT.back().emplace_back({k, *col}, *val);
        }
      }
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_beta.size(), NMO}, tp_ul_ul{0, 0},
                                    get_nnz(PsiT_MO[1], active_beta.data(), active_beta.size(), 0),
                                    Alloc(TGwfn.Node())));
      if (TGwfn.Node().root())
      {
        for (int k = 0; k < active_beta.size(); k++)
        {
          size_t ki = active_beta[k]; // occupied state #k
          auto col  = PsiT_MO[1].non_zero_indices2_data(ki);
          auto val  = PsiT_MO[1].non_zero_values_data(ki);
          for (size_t ic = 0, icend = PsiT_MO[1].num_non_zero_elements(ki); ic < icend; ic++, ++col, ++val)
            PsiT[1].emplace_back({k, *col}, *val);
        }
      }
    }
    // now that mappings have been constructed, map indexes of excited state orbitals
    // to the corresponding active space indexes
    if (TGwfn.Node().root())
    {
      // map reference
      auto refc = abij.reference_configuration();
      for (int i = 0; i < NAEA + NAEB; i++, ++refc)
        *refc = mo2active[*refc];
      for (int n = 1; n < abij.maximum_excitation_number()[0]; n++)
      {
        auto it  = abij.alpha_begin(n);
        auto ite = abij.alpha_end(n);
        for (; it < ite; ++it)
        {
          auto exct = (*it) + n; // only need to map excited state indexes
          for (int np = 0; np < n; ++np, ++exct)
            *exct = mo2active[*exct];
        }
      }
      for (int n = 1; n < abij.maximum_excitation_number()[1]; n++)
      {
        auto it  = abij.beta_begin(n);
        auto ite = abij.beta_end(n);
        for (; it < ite; ++it)
        {
          auto exct = (*it) + n; // only need to map excited state indexes
          for (int np = 0; np < n; ++np, ++exct)
            *exct = mo2active[*exct];
        }
      }
    }

    // is PureSD actually faster??? CHECK!!!
    // NOTE: For UHF reference, treat HOps as a 2 determinant wavefunction of a
    //        CLOSED walker type. This way you can treat alpha/beta sectors independently in
    //        HOps through the index corresponding 0/1.
    // never add coulomb to half rotated v2 tensor in PHMSD
    getInitialGuess(dump, name, NMO, NAEA, NAEB, walker_type);

    TGwfn.node_barrier();
    auto HOps(getHamOps(restart_file, CLOSED, NMO, NAEA, NAEB, PsiT,
      TGprop, TGwfn, cutvn, cutv2, false, h));
    TGwfn.node_barrier();
    // setup configuration couplings
    using index_aos = ma::sparse::array_of_sequences<int, int, shared_allocator<int>, ma::sparse::is_root>;
    shared_allocator<int> alloc_{TGwfn.Node()};

    // alpha
    std::vector<int> counts_alpha(abij.number_of_unique_excitations()[0]);
    std::vector<int> counts_beta(abij.number_of_unique_excitations()[1]);
    if (TGwfn.Node().root())
    {
      for (auto it = abij.configurations_begin(); it < abij.configurations_end(); ++it)
      {
        ++counts_alpha[std::get<0>(*it)];
        ++counts_beta[std::get<1>(*it)];
      }
    }
    TGwfn.Node().broadcast_n(counts_alpha.begin(), counts_alpha.size());
    TGwfn.Node().broadcast_n(counts_beta.begin(), counts_beta.size());
    index_aos beta_coupled_to_unique_alpha(counts_alpha.size(), counts_alpha, alloc_);
    index_aos alpha_coupled_to_unique_beta(counts_beta.size(), counts_beta, alloc_);
    if (TGwfn.Node().root())
    {
      int ni = 0;
      for (auto it = abij.configurations_begin(); it < abij.configurations_end(); ++it, ++ni)
      {
        beta_coupled_to_unique_alpha.emplace_back(std::get<0>(*it), ni);
        alpha_coupled_to_unique_beta.emplace_back(std::get<1>(*it), ni);
      }
    }
    TGwfn.Node().barrier();

    return Wavefunction(PHMSD(AFinfo, cur, TGwfn, std::move(HOps), std::move(acta2mo), std::move(actb2mo),
                              std::move(abij), std::move(beta_coupled_to_unique_alpha),
                              std::move(alpha_coupled_to_unique_beta), std::move(PsiT), walker_type, NCE, targetNW));
    app_error() << " Error: Wavefunction type PHMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type PHMSD not yet implemented. \n");
    return Wavefunction();
  }
  else if (type == "generalmsd")
  {
    app_error() << " Error: Wavefunction type GeneralMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type GeneralMSD not yet implemented. \n");
    return Wavefunction{};
  }
  else
  {
    app_error() << " Error: Unknown wave-function type: " << type << std::endl;
    APP_ABORT(" Error: Unknown wave-function type. \n");
    return Wavefunction{};
  }
}

/*
 * Read Initial walker from file.
*/
void WavefunctionFactory::getInitialGuess(hdf_archive& dump,
                                          std::string& name,
                                          int NMO,
                                          int NAEA,
                                          int NAEB,
                                          WALKER_TYPES walker_type)
{
  int NPOL = (walker_type == NONCOLLINEAR) ? 2 : 1;
  std::vector<int> dims(5);
  if (!dump.readEntry(dims, "dims"))
  {
    app_error() << " Error in getCommonInput(): Problems reading dims. \n";
    APP_ABORT("");
  }
  WALKER_TYPES wtype(initWALKER_TYPES(dims[3]));
  auto guess = initial_guess.find(name);
  if (guess == initial_guess.end())
  {
    auto newg = initial_guess.insert(std::make_pair(name, boost::multi::array<ComplexType, 3>({2, NPOL * NMO, NAEA})));
    if (!newg.second)
      APP_ABORT(" Error: Problems adding new initial guess. \n");
    using ma::conj;
    std::fill_n((newg.first)->second.origin(), 2 * NPOL * NMO * NAEA, ComplexType(0.0, 0.0));
    {
      boost::multi::array<ComplexType, 2> Psi0Alpha({NPOL * NMO, NAEA});
      if (!dump.readEntry(Psi0Alpha, "Psi0_alpha"))
      {
        app_error() << " Error in WavefunctionFactory: Initial wavefunction Psi0_alpha not found. \n";
        APP_ABORT("");
      }
      for (int i = 0; i < NPOL * NMO; i++)
        for (int j = 0; j < NAEA; j++)
          ((newg.first)->second)[0][i][j] = Psi0Alpha[i][j];
    }
    if (walker_type == COLLINEAR)
    {
      if (wtype == COLLINEAR)
      {
        boost::multi::array<ComplexType, 2> Psi0Beta({NMO, NAEB});
        if (!dump.readEntry(Psi0Beta, "Psi0_beta"))
        {
          app_error() << " Error in WavefunctionFactory: Initial wavefunction Psi0_beta not found. \n";
          APP_ABORT("");
        }
        for (int i = 0; i < NMO; i++)
          for (int j = 0; j < NAEB; j++)
            ((newg.first)->second)[1][i][j] = Psi0Beta[i][j];
      }
      else if (wtype == CLOSED)
      {
        boost::multi::array<ComplexType, 2> Psi0Beta({NMO, NAEA});
        assert(NAEA == NAEB);
        if (!dump.readEntry(Psi0Beta, "Psi0_alpha"))
        {
          app_error() << " Error in WavefunctionFactory: Initial wavefunction Psi0_beta not found. \n";
          APP_ABORT("");
        }
        for (int i = 0; i < NMO; i++)
          for (int j = 0; j < NAEB; j++)
            ((newg.first)->second)[1][i][j] = Psi0Beta[i][j];
      }
      else
        APP_ABORT(" Error: Unknown wtype. \n");
    }
  }
  else
    APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");
}


void WavefunctionFactory::computeVariationalEnergyPHMSD(TaskGroup_& TG,
                                                        Hamiltonian& ham,
                                                        boost::multi::array_ref<int, 2>& occs,
                                                        std::vector<ComplexType>& coeff,
                                                        int ndets,
                                                        int NAEA,
                                                        int NAEB,
                                                        int NMO,
                                                        bool recompute_ci)
{
  // CI coefficients can in general be complex and want to avoid two mpi communications so
  // keep everything complex even if Hamiltonian matrix elements are real.
  // Allocate H in Node's shared memory, but use as a raw array with proper synchronization
  int dim((recompute_ci ? ndets : 0));
  boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>> H({dim, dim}, TG.Node());
  boost::multi::array<ComplexType, 1> energy(iextensions<1u>{2});
  using std::fill_n;
  fill_n(H.origin(), H.num_elements(), ComplexType(0.0));           // this call synchronizes
  fill_n(energy.origin(), energy.num_elements(), ComplexType(0.0)); // this call synchronizes
  ValueType enuc = ham.getNuclearCoulombEnergy();
  for (int idet = 0; idet < ndets; idet++)
  {
    // These should already be sorted.
    boost::multi::array_ref<int, 1> deti(occs[idet].origin(), {NAEA + NAEB});
    ComplexType cidet = coeff[idet];
    for (int jdet = idet; jdet < ndets; jdet++)
    {
      // Compute <Di|H|Dj>
      if ((idet * ndets + jdet) % TG.Global().size() == TG.Global().rank())
      {
        if (idet == jdet)
        {
          ComplexType Hii(0.0);
          Hii = slaterCondon0(ham, deti, NMO) + enuc;
          energy[0] += ma::conj(cidet) * cidet * Hii;
          energy[1] += ma::conj(cidet) * cidet;
          if (recompute_ci)
            H[idet][idet] = Hii;
        }
        else
        {
          ComplexType Hij(0.0);
          boost::multi::array_ref<int, 1> detj(occs[jdet].origin(), {NAEA + NAEB});
          ComplexType cjdet = coeff[jdet];
          int perm          = 1;
          std::vector<int> excit;
          int nexcit = getExcitation(deti, detj, excit, perm);
          if (nexcit == 1)
          {
            Hij = ComplexType(perm) * slaterCondon1(ham, excit, detj, NMO);
          }
          else if (nexcit == 2)
          {
            Hij = ComplexType(perm) * slaterCondon2(ham, excit, NMO);
          }
          energy[0] += ma::conj(cidet) * cjdet * Hij + ma::conj(cjdet) * cidet * ma::conj(Hij);
          if (recompute_ci)
          {
            H[idet][jdet] = Hij;
            H[jdet][idet] = ma::conj(Hij);
          }
        }
      }
    }
  }
  TG.Node().barrier();
  if (TG.Node().root() && recompute_ci)
    TG.Cores().all_reduce_in_place_n(to_address(H.origin()), H.num_elements(), std::plus<>());
  TG.Global().all_reduce_in_place_n(energy.origin(), 2, std::plus<>());
  app_log() << " - Variational energy of trial wavefunction: " << std::setprecision(16) << energy[0] / energy[1]
            << "\n";
  if (recompute_ci)
  {
    app_log() << " - Diagonalizing CI matrix.\n";
    using RVector = boost::multi::array<RealType, 1>;
    using CMatrix = boost::multi::array<ComplexType, 2>;
    // Want a "unique" solution for all cores/nodes.
    if (TG.Global().rank() == 0)
    {
      std::pair<RVector, CMatrix> Sol = ma::symEigSelect<RVector, CMatrix>(H, 1);
      app_log() << " - Updating CI coefficients. \n";
      app_log() << " - Recomputed coe  // only makes Sparse ham. 2e rotation more effiecientfficient of first determinant: " << Sol.second[0][0] << "\n";
      for (int idet = 0; idet < ndets; idet++)
      {
        ComplexType ci = Sol.second[0][idet];
        // Do we want this much output?
        //app_log() << idet << " old: " << coeff[idet] << " new: " << ci << "\n";
        coeff[idet] = ci;
      }
      app_log() << " - Recomputed variational energy of trial wavefunction: " << Sol.first[0] << "\n";
    }
    TG.Global().broadcast_n(to_address(coeff.data()), coeff.size(), 0);
  }
}
/*
 * Helper function to get HamOps object from file or from scratch.
*/
HamiltonianOperations WavefunctionFactory::getHamOps(const std::string& restart_file,
                                                     WALKER_TYPES type,
                                                     int NMO,
                                                     int NAEA,
                                                     int NAEB,
                                                     std::vector<PsiT_Matrix>& PsiT,
                                                     TaskGroup_& TGprop,
                                                     TaskGroup_& TGwfn,
                                                     RealType cutvn,
                                                     RealType cutv2,
                                                     int ndets_to_read,
                                                     Hamiltonian& h)
{
  // if requested, create restart file
  // Will use phdf5 in the future, for now only head node writes
  hdf_archive restart(TGwfn.Global());
  if (restart_file != "")
  {
    if (TGwfn.Global().root())
    {
      app_log() << " Open " << restart_file << " for HamOps\n";
      // first, try open existing restart file
      if (!restart.open(restart_file, H5F_ACC_RDONLY))
      { // make new restart file if cannot open existing
        app_log() << " No restart_file create anew.\n";
        if (!restart.create(restart_file, H5F_ACC_EXCL))
        {
          app_error() << " Error in WavefunctionFactory: Failed to create restart_file" << restart_file << "\n";
          APP_ABORT("");
        }
      }
    }
  }

  bool read = restart.is_group("HamiltonianOperations");
  if (read)
  {
    app_log() << " getHamOps using restart file\n";
    return loadHamOps(restart, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2);
  }
  else
  {
    app_log() << " getHamOps from scratch\n";
    bool pureSD = false;  // make Sparse ham. 2e rotation more effiecient
    // pureSD only matters during setup of a large run. Avoid Sparse in that case.
    return h.getHamiltonianOperations(pureSD, ndets_to_read > 1, type, PsiT, cutvn, cutv2, TGprop, TGwfn, restart);
  }

//  if(restart_file != "")
//  {
//    if(TGwfn.Global().root()) {
//      restart.close();
//    }
//  }
}
/**
 * Compute the excitation level between two determinants.
 */
int WavefunctionFactory::getExcitation(boost::multi::array_ref<int, 1>& deti,
                                       boost::multi::array_ref<int, 1>& detj,
                                       std::vector<int>& excit,
                                       int& perm)
{
  std::vector<int> from_orb, to_orb;
  // Work out which orbitals are excited from / to.
  std::set_difference(detj.begin(), detj.end(), deti.begin(), deti.end(), std::inserter(from_orb, from_orb.begin()));
  std::set_difference(deti.begin(), deti.end(), detj.begin(), detj.end(), std::inserter(to_orb, to_orb.begin()));
  int nexcit = from_orb.size();
  if (nexcit <= 2)
  {
    for (int i = 0; i < from_orb.size(); i++)
      excit.push_back(from_orb[i]);
    for (int i = 0; i < to_orb.size(); i++)
      excit.push_back(to_orb[i]);
    int nperm = 0;
    int nmove = 0;
    for (auto o : from_orb)
    {
      auto it = std::find(detj.begin(), detj.end(), o);
      int loc = std::distance(detj.begin(), it);
      nperm += loc - nmove;
      nmove += 1;
    }
    nmove = 0;
    for (auto o : to_orb)
    {
      auto it = std::find(deti.begin(), deti.end(), o);
      int loc = std::distance(deti.begin(), it);
      nperm += loc - nmove;
      nmove += 1;
    }
    perm = nperm % 2 == 1 ? -1 : 1;
  }
  return nexcit;
}

ComplexType WavefunctionFactory::slaterCondon0(Hamiltonian& ham, boost::multi::array_ref<int, 1>& det, int NMO)
{
  ValueType one_body = ValueType(0.0);
  ValueType two_body = ValueType(0.0);
  for (int i = 0; i < det.size(); i++)
  {
    int oi = det[i];
    one_body += ham.H(oi, oi);
    for (int j = i + 1; j < det.size(); j++)
    {
      int oj = det[j];
      two_body += ham.H(oi, oj, oi, oj) - ham.H(oi, oj, oj, oi);
    }
  }
  return ComplexType(one_body + two_body);
}

ComplexType WavefunctionFactory::slaterCondon1(Hamiltonian& ham,
                                               std::vector<int>& excit,
                                               boost::multi::array_ref<int, 1>& det,
                                               int NMO)
{
  int i              = excit[0];
  int a              = excit[1];
  ValueType one_body = ham.H(i, a);
  ValueType two_body = ValueType(0.0);
  for (auto j : det)
  {
    two_body += ham.H(i, j, a, j) - ham.H(i, j, j, a);
  }
  return ComplexType(one_body + two_body);
}

ComplexType WavefunctionFactory::slaterCondon2(Hamiltonian& ham, std::vector<int>& excit, int NMO)
{
  int i = excit[0];
  int j = excit[1];
  int a = excit[2];
  int b = excit[3];
  return ComplexType(ham.H(i, j, a, b) - ham.H(i, j, b, a));
}

//ComplexType WavefunctionFactory::contractOneBody(std::vector<int>& det, std::vector<int>& excit, boost::multi::array_ref<ComplexType,2>& HSPot, int NMO)
//{
//ComplexType oneBody = ComplexType(0.0);
//int spini, spina;
//if(excit.size()==0) {
//for(auto i : det) {
//int oi = decodeSpinOrbital(i, spini, NMO);
//oneBody += HSPot[oi][oi];
//}
//} else {
//int oi = decodeSpinOrbital(excit[0], spini, NMO);
//int oa = decodeSpinOrbital(excit[1], spina, NMO);
//oneBody = HSPot[oi][oa];
//}
//return oneBody;
//}


} // namespace afqmc
} // namespace qmcplusplus
