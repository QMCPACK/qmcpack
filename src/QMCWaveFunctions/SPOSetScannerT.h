//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_SCANNERT_H
#define QMCPLUSPLUS_SPOSET_SCANNERT_H

#include "OhmmsData/AttributeSet.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/SPOSetT.h"

namespace qmcplusplus
{
template<typename T>
struct OutputReportMakerBase
{
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using GradVector  = typename SPOSetT<T>::GradVector;

  const ValueVector& SPO_v_avg;
  const ValueVector& SPO_l_avg;
  const GradVector& SPO_g_avg;
  int nknots;
};

template<typename T>
struct OutputReportMaker : OutputReportMakerBase<T>
{
  using RealType = typename SPOSetT<T>::RealType;

  void operator()(std::ofstream& output_report) const
  {
    output_report << "#   Report: Orb   Value_avg   Gradients_avg   Laplacian_avg" << std::endl;
    for (int iorb = 0; iorb < this->SPO_v_avg.size(); iorb++)
    {
      auto one_over_nknots = static_cast<RealType>(1.0 / this->nknots);
      output_report << "\t" << iorb << "    " << std::scientific << this->SPO_v_avg[iorb] * one_over_nknots << "   "
                    << this->SPO_g_avg[iorb][0] * one_over_nknots << "   " << this->SPO_g_avg[iorb][1] * one_over_nknots
                    << "   " << this->SPO_g_avg[iorb][2] * one_over_nknots << "   "
                    << this->SPO_l_avg[iorb] * one_over_nknots << std::fixed << std::endl;
    }
  }
};

template<typename T>
struct OutputReportMaker<std::complex<T>> : OutputReportMakerBase<std::complex<T>>
{
  using RealType = typename SPOSetT<T>::RealType;

  void operator()(std::ofstream& output_report) const
  {
    output_report << "#   Report: Orb   Value_avg I/R  Gradients_avg Laplacian_avg" << std::endl;
    for (int iorb = 0; iorb < this->SPO_v_avg.size(); iorb++)
    {
      auto one_over_nknots = static_cast<RealType>(1.0 / this->nknots);
      output_report << "\t" << iorb << "    " << std::scientific << this->SPO_v_avg[iorb] * one_over_nknots << "   "
                    << this->SPO_v_avg[iorb].imag() / this->SPO_v_avg[iorb].real() << "   "
                    << this->SPO_g_avg[iorb][0] * one_over_nknots << "   " << this->SPO_g_avg[iorb][1] * one_over_nknots
                    << "   " << this->SPO_g_avg[iorb][2] * one_over_nknots << "   "
                    << this->SPO_l_avg[iorb] * one_over_nknots << std::fixed << std::endl;
    }
  }
};

/** a scanner for all the SPO sets.
 */
template<typename T>
class SPOSetScannerT
{
public:
  using PtclPool    = std::map<std::string, const std::unique_ptr<ParticleSetT<T>>>;
  using SPOSetMap   = typename SPOSetT<T>::SPOMap;
  using RealType    = typename SPOSetT<T>::RealType;
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using HessVector  = typename SPOSetT<T>::HessVector;

  RealType myfabs(RealType s) { return std::fabs(s); }
  template<typename U>
  std::complex<U> myfabs(std::complex<U>& s)
  {
    return std::complex<U>(myfabs(s.real()), myfabs(s.imag()));
  }
  template<typename U>
  TinyVector<U, OHMMS_DIM> myfabs(TinyVector<U, OHMMS_DIM>& s)
  {
    return TinyVector<U, OHMMS_DIM>(myfabs(s[0]), myfabs(s[1]), myfabs(s[2]));
  }

  const SPOSetMap& sposets;
  ParticleSetT<T>& target;
  const PtclPool& ptcl_pool_;
  ParticleSetT<T>* ions;

  // construction/destruction
  SPOSetScannerT(const SPOSetMap& sposets_in, ParticleSetT<T>& targetPtcl, const PtclPool& psets)
      : sposets(sposets_in), target(targetPtcl), ptcl_pool_(psets), ions(0){};
  //~SPOSetScannerT(){};

  // processing scanning
  void put(xmlNodePtr cur)
  {
    app_log() << "Entering the SPO set scanner!" << std::endl;
    // check in the source particle set and search for it in the pool.
    std::string sourcePtcl("ion0");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourcePtcl, "source");
    aAttrib.put(cur);
    auto pit(ptcl_pool_.find(sourcePtcl));
    if (pit == ptcl_pool_.end())
      app_log() << "Source particle set not found. Can not be used as "
                   "reference point."
                << std::endl;
    else
      ions = pit->second.get();

    // scanning the SPO sets
    xmlNodePtr cur_save = cur;
    for (const auto& [name, sposet] : sposets)
    {
      app_log() << "  Processing SPO " << sposet->getName() << std::endl;
      // scanning the paths
      cur = cur_save->children;
      while (cur != NULL)
      {
        std::string trace_name("no name");
        OhmmsAttributeSet aAttrib;
        aAttrib.add(trace_name, "name");
        aAttrib.put(cur);
        std::string cname(getNodeName(cur));
        std::string prefix(sposet->getName() + "_" + cname + "_" + trace_name);
        if (cname == "path")
        {
          app_log() << "    Scanning a " << cname << " called " << trace_name << " and writing to "
                    << prefix + "_v/g/l/report.dat" << std::endl;
          auto spo = sposet->makeClone();
          scan_path(cur, *spo, prefix);
        }
        else
        {
          if (cname != "text" && cname != "comment")
            app_log() << "    Unknown type of scanning " << cname << std::endl;
        }
        cur = cur->next;
      }
    }
    app_log() << "Exiting the SPO set scanner!" << std::endl << std::endl;
  }

  // scanning a path
  void scan_path(xmlNodePtr cur, SPOSetT<T>& sposet, std::string prefix)
  {
    std::string file_name;
    file_name = prefix + "_v.dat";
    std::ofstream output_v(file_name.c_str());
    file_name = prefix + "_g.dat";
    std::ofstream output_g(file_name.c_str());
    file_name = prefix + "_l.dat";
    std::ofstream output_l(file_name.c_str());
    file_name = prefix + "_report.dat";
    std::ofstream output_report(file_name.c_str());

    int nknots(2);
    int from_atom(-1);
    int to_atom(-1);
    TinyVector<double, OHMMS_DIM> from_pos(0.0, 0.0, 0.0);
    TinyVector<double, OHMMS_DIM> to_pos(0.0, 0.0, 0.0);

    OhmmsAttributeSet aAttrib;
    aAttrib.add(nknots, "nknots");
    aAttrib.add(from_atom, "from_atom");
    aAttrib.add(to_atom, "to_atom");
    aAttrib.add(from_pos, "from_pos");
    aAttrib.add(to_pos, "to_pos");
    aAttrib.put(cur);

    // sanity check
    if (nknots < 2)
      nknots = 2;
    // check out the reference atom coordinates
    if (ions)
    {
      if (from_atom >= 0 && from_atom < ions->R.size())
        from_pos = ions->R[from_atom];
      if (to_atom >= 0 && to_atom < ions->R.size())
        to_pos = ions->R[to_atom];
    }

    // prepare a fake particle set
    ValueVector SPO_v, SPO_l, SPO_v_avg, SPO_l_avg;
    GradVector SPO_g, SPO_g_avg;
    int OrbitalSize(sposet.size());
    SPO_v.resize(OrbitalSize);
    SPO_g.resize(OrbitalSize);
    SPO_l.resize(OrbitalSize);
    SPO_v_avg.resize(OrbitalSize);
    SPO_g_avg.resize(OrbitalSize);
    SPO_l_avg.resize(OrbitalSize);
    SPO_v_avg      = 0.0;
    SPO_g_avg      = 0.0;
    SPO_l_avg      = 0.0;
    double Delta   = 1.0 / (nknots - 1);
    int elec_count = target.R.size();
    auto R_saved   = target.R;
    typename ParticleSetT<T>::SingleParticlePos zero_pos(0.0, 0.0, 0.0);
    for (int icount = 0, ind = 0; icount < nknots; icount++, ind++)
    {
      if (ind == elec_count)
        ind = 0;
      target.R[ind][0] = (to_pos[0] - from_pos[0]) * Delta * icount + from_pos[0];
      target.R[ind][1] = (to_pos[1] - from_pos[1]) * Delta * icount + from_pos[1];
      target.R[ind][2] = (to_pos[2] - from_pos[2]) * Delta * icount + from_pos[2];
      target.makeMove(ind, zero_pos);
      sposet.evaluateVGL(target, ind, SPO_v, SPO_g, SPO_l);
      std::ostringstream o;
      o << "x_y_z  " << std::fixed << std::setprecision(7) << target.R[ind][0] << " " << target.R[ind][1] << " "
        << target.R[ind][2];
      output_v << o.str() << " : " << std::scientific << std::setprecision(12);
      output_g << o.str() << " : " << std::scientific << std::setprecision(12);
      output_l << o.str() << " : " << std::scientific << std::setprecision(12);
      for (int iorb = 0; iorb < OrbitalSize; iorb++)
      {
        SPO_v_avg[iorb] += myfabs(SPO_v[iorb]);
        SPO_g_avg[iorb] += myfabs(SPO_g[iorb]);
        SPO_l_avg[iorb] += myfabs(SPO_l[iorb]);
        output_v << SPO_v[iorb] << "  ";
        output_g << SPO_g[iorb][0] << "  " << SPO_g[iorb][1] << "  " << SPO_g[iorb][2] << "  ";
        output_l << SPO_l[iorb] << "  ";
      }
      output_v << std::endl;
      output_g << std::endl;
      output_l << std::endl;
    }
    // restore the whole target.
    target.R = R_saved;
    target.update();
    OutputReportMaker<T>{SPO_v_avg, SPO_l_avg, SPO_g_avg, nknots}(output_report);
    output_v.close();
    output_g.close();
    output_l.close();
    output_report.close();
  }
};
} // namespace qmcplusplus

#endif
