//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "EnergyDensityEstimator.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Particle/DistanceTable.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Utilities/string_utils.h"
#include <string>
#include <vector>

namespace qmcplusplus
{
EnergyDensityEstimator::EnergyDensityEstimator(const PSPool& PSP, const std::string& defaultKE)
    : psetpool(PSP), Pdynamic(0), Pstatic(0), w_trace(0), Td_trace(0), Vd_trace(0), Vs_trace(0)
{
  update_mode_.set(COLLECTABLE, 1);
  defKE      = defaultKE;
  nsamples   = 0;
  ion_points = false;
  nions      = -1;
  request_.request_scalar("weight");
  request_.request_array("Kinetic");
  request_.request_array("LocalPotential");
}


EnergyDensityEstimator::~EnergyDensityEstimator() { delete_iter(spacegrids.begin(), spacegrids.end()); }


bool EnergyDensityEstimator::put(xmlNodePtr cur, ParticleSet& Pdyn)
{
  Pdynamic = &Pdyn;
  return put(cur);
}


/** check xml elements
 */
bool EnergyDensityEstimator::put(xmlNodePtr cur)
{
  input_xml = cur;
  //initialize simple xml attributes
  name_ = "EnergyDensity";
  std::string dyn, stat = "";
  ion_points = false;
  OhmmsAttributeSet attrib;
  attrib.add(name_, "name");
  attrib.add(dyn, "dynamic");
  attrib.add(stat, "static");
  attrib.add(ion_points, "ion_points");
  attrib.put(cur);
  //collect particle sets
  if (!Pdynamic)
    Pdynamic = get_particleset(dyn);
  if (Pdynamic->hasSK())
    Pdynamic->turnOnPerParticleSK();
  nparticles = Pdynamic->getTotalNum();
  std::vector<ParticleSet*> Pref;
  if (stat == "")
  {
    Pstatic      = 0;
    dtable_index = -1;
  }
  else
  {
    Pstatic = get_particleset(stat);
    if (Pstatic->hasSK())
      Pstatic->turnOnPerParticleSK();
    dtable_index = Pdynamic->addTable(*Pstatic);
    Pref.resize(1);
    Pref[0] = Pstatic;
    if (!ion_points)
      nparticles += Pstatic->getTotalNum();
    else
      nions = Pstatic->getTotalNum();
  }
  //size arrays
  R.resize(nparticles);
  EDValues.resize(nparticles, nEDValues);
  if (ion_points)
  {
    EDIonValues.resize(nions, nEDValues);
    Rion.resize(nions, DIM);
    for (int i = 0; i < nions; i++)
      for (int d = 0; d < DIM; d++)
        Rion(i, d) = Pstatic->R[i][d];
  }
  particles_outside.resize(nparticles);
  fill(particles_outside.begin(), particles_outside.end(), true);
  //read xml element contents
  xmlNodePtr element;
  bool stop = false;
  //initialize reference points
  app_log() << "Initializing reference points" << std::endl;
  bool has_ref = false;
  element      = cur->children;
  while (element != NULL)
  {
    std::string name((const char*)element->name);
    if (name == "reference_points")
    {
      if (has_ref)
      {
        APP_ABORT("EnergyDensityEstimator::put: EDE can only have one instance of reference_points.");
      }
      else
      {
        bool ref_succeeded = ref.put(element, *Pdynamic, Pref);
        stop               = stop || !ref_succeeded;
        has_ref            = true;
      }
    }
    element = element->next;
  }
  if (!has_ref)
  {
    bool ref_succeeded = ref.put(*Pdynamic, Pref);
    stop               = stop || !ref_succeeded;
  }
  //initialize grids or other cell partitions
  bool periodic = Pdynamic->getLattice().SuperCellEnum != SUPERCELL_OPEN;
  bool grid_succeeded;
  element     = cur->children;
  int nvalues = (int)nEDValues;
  while (element != NULL)
  {
    std::string name = (const char*)element->name;
    if (name == "spacegrid")
    {
      SpaceGrid* sg = new SpaceGrid(nvalues);
      spacegrids.push_back(sg);
      if (Pstatic)
      {
        set_ptcl();
        grid_succeeded = sg->put(element, ref.points, Rptcl, Zptcl, Pdynamic->getTotalNum(), periodic, false);
        unset_ptcl();
      }
      else
        grid_succeeded = sg->put(element, ref.points, periodic, false);
      stop = stop || !grid_succeeded;
    }
    element = element->next;
  }
  if (stop == true)
  {
    APP_ABORT("EnergyDensityEstimator::put");
  }
  return true;
}


void EnergyDensityEstimator::set_ptcl()
{
  ParticleSet& P = *Pstatic;
  SpeciesSet& species(P.getSpeciesSet());
  int ChargeAttribIndx = species.addAttribute("charge");
  int nspecies         = species.TotalNum;
  int nps              = P.getTotalNum();
  std::vector<RealType> Zspec;
  Zspec.resize(nspecies);
  Zptcl.resize(nps);
  for (int spec = 0; spec < nspecies; spec++)
    Zspec[spec] = species(ChargeAttribIndx, spec);
  for (int i = 0; i < nps; i++)
    Zptcl[i] = Zspec[P.GroupID[i]];
  Rptcl.resize(P.R.size());
  for (int i = 0; i < P.R.size(); i++)
    Rptcl[i] = P.R[i];
  if (P.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    P.applyMinimumImage(Rptcl);
}

void EnergyDensityEstimator::unset_ptcl()
{
  Zptcl.clear();
  Rptcl.clear();
}


ParticleSet* EnergyDensityEstimator::get_particleset(std::string& psname)
{
  auto pit(psetpool.find(psname));
  if (pit == psetpool.end())
  {
    app_log() << "  ParticleSet " << psname << " does not exist" << std::endl;
    APP_ABORT("EnergyDensityEstimator::put");
  }
  return pit->second.get();
}


void EnergyDensityEstimator::getRequiredTraces(TraceManager& tm)
{
  bool write = omp_get_thread_num() == 0;
  w_trace    = tm.get_real_trace("weight");
  Td_trace   = tm.get_real_trace(*Pdynamic, "Kinetic");
  Vd_trace   = tm.get_real_combined_trace(*Pdynamic, "LocalPotential");
  if (Pstatic)
    Vs_trace = tm.get_real_combined_trace(*Pstatic, "LocalPotential");
  have_required_traces_ = true;
}


void EnergyDensityEstimator::write_description(std::ostream& os)
{
  os << "EnergyDensityEstimator::write_description" << std::endl;
  os << std::endl;
  os << "  EnergyDensityEstimator details" << std::endl;
  os << std::endl;
  std::string indent = "    ";
  os << indent + "nparticles  = " << nparticles << std::endl;
  os << indent + "nspacegrids = " << spacegrids.size() << std::endl;
  os << std::endl;
  ref.write_description(os, indent);
  os << std::endl;
  for (int i = 0; i < spacegrids.size(); i++)
  {
    spacegrids[i]->write_description(os, indent);
  }
  os << std::endl;
  os << "  end EnergyDensityEstimator details" << std::endl;
  os << std::endl;
  os << "end EnergyDensityEstimator::write_description" << std::endl;
  return;
}


bool EnergyDensityEstimator::get(std::ostream& os) const
{
  os << "EDM replace this " << std::endl;
  APP_ABORT("EnergyDensityEstimator::get");
  return true;
}


void EnergyDensityEstimator::resetTargetParticleSet(ParticleSet& P)
{
  //remains empty
}


//#define ENERGYDENSITY_CHECK

EnergyDensityEstimator::Return_t EnergyDensityEstimator::evaluate(ParticleSet& P)
{
  if (have_required_traces_)
  {
    Pdynamic = &P;
    //Collect positions from ParticleSets
    int p = 0;
    {
      const ParticlePos& Rs = Pdynamic->R;
      for (int i = 0; i < Rs.size(); i++)
      {
        R[p] = Rs[i];
        p++;
      }
    }
    if (Pstatic && !ion_points)
    {
      const ParticlePos& Rs = Pstatic->R;
      for (int i = 0; i < Rs.size(); i++)
      {
        R[p] = Rs[i];
        p++;
      }
    }
    if (P.getLattice().SuperCellEnum != SUPERCELL_OPEN)
      P.applyMinimumImage(R);
    //Convert information accumulated in ParticleSets into EnergyDensity quantities
    RealType w = w_trace->sample[0];
    p          = 0;
    {
      Vd_trace->combine();
      const ParticleSet& Ps = *Pdynamic;
      const auto& Ts        = Td_trace->sample;
      const auto& Vs        = Vd_trace->sample;
      for (int i = 0; i < Ps.getTotalNum(); i++)
      {
        EDValues(p, W) = w;
        EDValues(p, T) = w * Ts[i];
        EDValues(p, V) = w * Vs[i];
        p++;
      }
    }
    if (Pstatic)
    {
      Vs_trace->combine();
      const ParticleSet& Ps = *Pstatic;
      const auto& Vs        = Vs_trace->sample;
      if (!ion_points)
        for (int i = 0; i < Ps.getTotalNum(); i++)
        {
          EDValues(p, W) = w;
          EDValues(p, T) = 0.0;
          EDValues(p, V) = w * Vs[i];
          p++;
        }
      else
        for (int i = 0; i < Ps.getTotalNum(); i++)
        {
          EDIonValues(i, W) = w;
          EDIonValues(i, T) = 0.0;
          EDIonValues(i, V) = w * Vs[i];
        }
    }
    //Accumulate energy density in spacegrids
    const auto& dtab(P.getDistTableAB(dtable_index));
    fill(particles_outside.begin(), particles_outside.end(), true);
    for (int i = 0; i < spacegrids.size(); i++)
    {
      SpaceGrid& sg = *spacegrids[i];
      sg.evaluate(R, EDValues, P.Collectables, particles_outside, dtab);
    }
    //Accumulate energy density of particles outside any spacegrid
    int bi, v;
    const int bimax = outside_buffer_offset + (int)nEDValues;
    for (int p = 0; p < particles_outside.size(); p++)
    {
      if (particles_outside[p])
      {
        for (bi = outside_buffer_offset, v = 0; bi < bimax; bi++, v++)
        {
          P.Collectables[bi] += EDValues(p, v);
        }
      }
    }
    if (ion_points)
    {
      // Accumulate energy density for ions at a point field
      bi = ion_buffer_offset;
      for (int i = 0; i < nions; i++)
        for (v = 0; v < (int)nEDValues; v++, bi++)
        {
          P.Collectables[bi] += EDIonValues(i, v);
        }
    }
    nsamples++;
#if defined(ENERGYDENSITY_CHECK)
    int thread    = omp_get_thread_num();
    using WP      = WalkerProperties::Indexes;
    RealType Eref = P.PropertyList[WP::LOCALENERGY];
    RealType Vref = P.PropertyList[WP::LOCALPOTENTIAL];
    RealType Tref = Eref - Vref;
#pragma omp critical(edcheck)
    {
      RealType Dsum = 0.0;
      RealType Tsum = 0.0;
      RealType Vsum = 0.0;
      RealType Esum = 0.0;
      for (int p = 0; p < nparticles; p++)
      {
        Dsum += EDValues(p, W);
        Tsum += EDValues(p, T);
        Vsum += EDValues(p, V);
      }
      if (ion_points)
        for (int i = 0; i < nions; i++)
        {
          Dsum += EDIonValues(i, W);
          Tsum += EDIonValues(i, T);
          Vsum += EDIonValues(i, V);
        }
      Esum           = Tsum + Vsum;
      static int cnt = 0;
      //app_log()<<"eval ED Dsum"<<cnt<<" "<<Dsum<< std::endl;
      app_log() << thread << " eval ED " << cnt << " " << Tsum << " " << Vsum << " " << Esum << std::endl;
      int nvals = (int)nEDValues;
      RealType edvals[nvals];
      RealType edtmp[nvals];
      for (int v = 0; v < nvals; v++)
        edvals[v] = 0.0;
      for (int i = 0; i < spacegrids.size(); i++)
      {
        SpaceGrid& sg = *spacegrids[i];
        sg.sum(P.Collectables, edtmp);
        for (int v = 0; v < nvals; v++)
          edvals[v] += edtmp[v];
      }
      for (int v = 0; v < nvals; v++)
        edvals[v] += P.Collectables[outside_buffer_offset + v];
      if (ion_points)
      {
        for (int i = 0; i < nions; i++)
          for (int v = 0; v < nvals; v++)
            edvals[v] += P.Collectables[ion_buffer_offset + i * nvals + v];
      }
      //app_log()<<"eval ES Dsum"<<cnt<<" "<<edvals[W]<< std::endl;
      app_log() << thread << " eval ES " << cnt << " " << edvals[T] << " " << edvals[V] << " " << edvals[T] + edvals[V]
                << std::endl;
      app_log() << thread << " ref  E  " << cnt << " " << Tref << " " << Vref << " " << Eref << std::endl;
      cnt++;
    }
#endif
    //APP_ABORT("EnergyDensityEstimator::evaluate");
  }

  return 0.0;
}


void EnergyDensityEstimator::write_Collectables(std::string& label, int& cnt, ParticleSet& P)
{
  //for(int v=0;v<nEDValues;v++){
  int ii    = spacegrids[0]->buffer_offset;
  int io    = outside_buffer_offset;
  double Ti = P.Collectables[ii + 1] / P.Collectables[ii] * 12.0;
  double To = P.Collectables[io + 1] / P.Collectables[io] * 12.0;
  app_log() << "EDcoll " << label << cnt << " " << Ti << " " << To << std::endl;
  //}
}


void EnergyDensityEstimator::write_EDValues(void)
{
  app_log() << "EDValues" << std::endl;
  for (int p = 0; p < nparticles; p++)
    fprintf(stdout, "  %d %e %e %e\n", p, EDValues(p, 0), EDValues(p, 1), EDValues(p, 2));
}

void EnergyDensityEstimator::write_nonzero_domains(const ParticleSet& P)
{
  app_log() << "Nonzero domains" << std::endl;
  int nd = 1;
  for (int i = 0; i < spacegrids.size(); i++)
    nd += spacegrids[i]->nDomains();
  for (int i = 0; i < nd; i++)
  {
    bool nonzero = false;
    int n        = outside_buffer_offset + i * nEDValues;
    for (int v = 0; v < nEDValues; v++)
    {
      nonzero = nonzero || std::abs(P.Collectables[n + v]) > 1e-8;
    }
    if (nonzero)
    {
      //      fprintf(stdout,"  %d %e %e %e %e %e %e\n",i,P.Collectables[n],
      fprintf(stdout, "  %d %e %e %e \n", i, P.Collectables[n], P.Collectables[n + 1], P.Collectables[n + 2]);
    }
  }
}


void EnergyDensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.size();
  //allocate space for energy density outside of any spacegrid
  outside_buffer_offset = collectables.size();
  int nvalues           = (int)nEDValues;
  std::vector<RealType> tmp(nvalues);
  collectables.add(tmp.begin(), tmp.end());
  //allocate space for spacegrids
  for (int i = 0; i < spacegrids.size(); i++)
  {
    spacegrids[i]->allocate_buffer_space(collectables);
  }
  if (ion_points)
  {
    ion_buffer_offset = collectables.size();
    nvalues           = nions * ((int)nEDValues);
    std::vector<RealType> tmp2(nvalues);
    collectables.add(tmp2.begin(), tmp2.end());
  }
}


void EnergyDensityEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  hdf_path hdf_name{name_};
  h5desc.emplace_back(hdf_name / "variables");
  auto& oh = h5desc.back();
  oh.addProperty(const_cast<int&>(nparticles), "nparticles", file);
  int nspacegrids = spacegrids.size();
  oh.addProperty(const_cast<int&>(nspacegrids), "nspacegrids", file);
  oh.addProperty(const_cast<int&>(nsamples), "nsamples", file);
  if (ion_points)
  {
    oh.addProperty(const_cast<int&>(nions), "nions", file);
    oh.addProperty(const_cast<Matrix<RealType>&>(Rion), "ion_positions", file);
  }

  ref.save(h5desc, file);
  h5desc.emplace_back(hdf_name / "outside");
  auto& ohOutside = h5desc.back();
  std::vector<int> ng(1);
  ng[0] = (int)nEDValues;
  ohOutside.set_dimensions(ng, outside_buffer_offset);
  for (int i = 0; i < spacegrids.size(); i++)
  {
    SpaceGrid& sg = *spacegrids[i];
    sg.registerCollectables(h5desc, file, i);
  }
  if (ion_points)
  {
    std::vector<int> ng2(2);
    ng2[0] = nions;
    ng2[1] = (int)nEDValues;

    h5desc.emplace_back(hdf_name / "ions");
    auto& ohIons = h5desc.back();
    ohIons.set_dimensions(ng2, ion_buffer_offset);
  }
}

void EnergyDensityEstimator::setObservables(PropertySetType& plist)
{
  //remains empty
}

void EnergyDensityEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  //remains empty
}


std::unique_ptr<OperatorBase> EnergyDensityEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  bool write = omp_get_thread_num() == 0;
  if (write)
    app_log() << "EnergyDensityEstimator::makeClone" << std::endl;

  std::unique_ptr<EnergyDensityEstimator> edclone = std::make_unique<EnergyDensityEstimator>(psetpool, defKE);
  edclone->put(input_xml, qp);
  //int thread = omp_get_thread_num();
  //app_log()<<thread<<"make edclone"<< std::endl;
  //edclone->Pdynamic = Pdynamic->get_clone(thread);
  //edclone->Pstatic  = Pstatic->get_clone(thread);
  return edclone;
}


} // namespace qmcplusplus
