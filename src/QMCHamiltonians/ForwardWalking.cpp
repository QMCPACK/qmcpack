//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ForwardWalking.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCDrivers/WalkerProperties.h"
#include <cstdlib>
#include <set>
#include <string>

namespace qmcplusplus
{
ForwardWalking::ForwardWalking() { update_mode_.set(OPTIMIZABLE, 1); }

ForwardWalking::~ForwardWalking() = default;

ForwardWalking::Return_t ForwardWalking::rejectedMove(ParticleSet& P)
{
  for (int i = 0; i < nobservables_; i++)
  {
    int lastindex = t_walker_->PHindex[p_ids_[i]] - 1;
    if (lastindex < 0)
      lastindex += walker_lengths_[i][2];
    t_walker_->addPropertyHistoryPoint(p_ids_[i], t_walker_->PropertyHistory[p_ids_[i]][lastindex]);
  }
  calculate(P);
  return 0.0;
}

ForwardWalking::Return_t ForwardWalking::calculate(ParticleSet& P)
{
  std::vector<RealType>::iterator Vit = values_.begin();
  for (int i = 0; i < nobservables_; i++)
  {
    int j       = 0;
    int FWindex = t_walker_->PHindex[p_ids_[i]] - 1;
    while (j < walker_lengths_[i][1])
    {
      FWindex -= walker_lengths_[i][0];
      if (FWindex < 0)
        FWindex += walker_lengths_[i][2];
      (*Vit) = t_walker_->PropertyHistory[p_ids_[i]][FWindex];
      j++;
      Vit++;
    }
  }
  copy(values_.begin(), values_.end(), t_walker_->getPropertyBase() + first_hamiltonian_ + my_index_);
  return 0.0;
}

ForwardWalking::Return_t ForwardWalking::evaluate(ParticleSet& P)
{
  for (int i = 0; i < nobservables_; i++)
    t_walker_->addPropertyHistoryPoint(p_ids_[i], P.PropertyList[h_ids_[i]]);
  calculate(P);
  return 0.0;
}

bool ForwardWalking::put(xmlNodePtr cur) { return true; }

bool ForwardWalking::putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P)
{
  using WP = WalkerProperties::Indexes;

  first_hamiltonian_ = h.startIndex();
  nobservables_      = 0;
  nvalues_           = 0;
  xmlNodePtr tcur    = cur->children;
  while (tcur != NULL)
  {
    std::string cname((const char*)tcur->name);
    if (cname == "Observable")
    {
      std::string tagName("none");
      int Hindex(-100);
      int blockSeries(0);
      int blockFreq(0);
      OhmmsAttributeSet Tattrib;
      Tattrib.add(tagName, "name");
      Tattrib.add(blockSeries, "max");
      Tattrib.add(blockFreq, "frequency");
      Tattrib.put(tcur);
      if (tagName.find("*") == std::string::npos)
      {
        //Single Observable case
        int numProps = P.PropertyList.Names.size();
        Hindex       = h.getObservable(tagName) + WP::NUMPROPERTIES;
        if (tagName == "LocalPotential")
        {
          Hindex  = WP::LOCALPOTENTIAL;
          tagName = "LocPot";
        }
        else if (tagName == "LocalEnergy")
        {
          Hindex  = WP::LOCALENERGY;
          tagName = "LocEn";
        }
        else if (Hindex == (WP::NUMPROPERTIES - 1))
        {
          app_log() << "Not a valid H element(" << Hindex << ") Valid names are:";
          for (int jk = 0; jk < h.sizeOfObservables(); jk++)
            app_log() << "  " << h.getObservableName(jk);
          app_log() << std::endl;
          exit(-1);
        }
        names_.push_back(tagName);
        h_ids_.push_back(Hindex);
        app_log() << " Hamiltonian Element " << tagName << " was found at " << Hindex << std::endl;
        int numT = blockSeries / blockFreq;
        nobservables_ += 1;
        nvalues_ += numT;
        app_log() << "   " << numT << " values will be calculated every " << blockFreq << "*tau H^-1" << std::endl;
        std::vector<int> pms(3);
        pms[0] = blockFreq;
        pms[1] = numT;
        pms[2] = blockSeries + 2;
        walker_lengths_.push_back(pms);
        int maxWsize = blockSeries + 2;
        int pindx    = P.addPropertyHistory(maxWsize);
        p_ids_.push_back(pindx);
      }
      else
      {
        bool FOUNDH(false);
        // 	    Multiple observables for this tag
        int found = tagName.rfind("*");
        tagName.replace(found, 1, "");
        int numProps = P.PropertyList.Names.size();
        for (int j = 0; j < h.sizeOfObservables(); j++)
        {
          std::string Hname = h.getObservableName(j);
          if (Hname.find(tagName) != std::string::npos)
          {
            FOUNDH = true;
            app_log() << " Hamiltonian Element " << Hname << " was found at " << j << std::endl;
            names_.push_back(Hname);
            Hindex = j + WP::NUMPROPERTIES;
            h_ids_.push_back(Hindex);
            int numT = blockSeries / blockFreq;
            nobservables_ += 1;
            nvalues_ += numT;
            app_log() << "   " << numT << " values will be calculated every " << blockFreq << "*tau H^-1" << std::endl;
            std::vector<int> pms(3);
            pms[0] = blockFreq;
            pms[1] = numT;
            pms[2] = blockSeries + 2;
            walker_lengths_.push_back(pms);
            int maxWsize = blockSeries + 2;
            int pindx    = P.addPropertyHistory(maxWsize);
            p_ids_.push_back(pindx);
          }
        }
        //handle FOUNDH
        if (FOUNDH)
        {
          nobservables_ += 1;
        }
        else
        {
          app_log() << "Not a valid H element(" << Hindex << ") Valid names are:";
          for (int jk = 0; jk < h.sizeOfObservables(); jk++)
            app_log() << "  " << h.getObservableName(jk);
          app_log() << std::endl;
          APP_ABORT("ForwardWalking::put");
        }
      }
    }
    tcur = tcur->next;
  }
  app_log() << "Total number of observables calculated:" << nobservables_ << std::endl;
  app_log() << "Total number of values calculated:" << nvalues_ << std::endl;
  values_.resize(nvalues_);
  return true;
}

bool ForwardWalking::get(std::ostream& os) const
{
  os << "ForwardWalking";
  return true;
}

std::unique_ptr<OperatorBase> ForwardWalking::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //nothing to worry, default copy constructor will do
  return std::make_unique<ForwardWalking>(*this);
}

void ForwardWalking::addObservables(PropertySetType& plist)
{
  //not used
}

void ForwardWalking::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = plist.size();
  int nc    = 0;
  for (int i = 0; i < nobservables_; ++i)
    for (int j = 0; j < walker_lengths_[i][1]; ++j, ++nc)
    {
      std::stringstream sstr;
      sstr << "FWE_" << names_[i] << "_" << j * walker_lengths_[i][0];
      int id = plist.add(sstr.str());
    }
  app_log() << "ForwardWalking::Observables [" << my_index_ << ", " << my_index_ + nc << ")" << std::endl;
}

void ForwardWalking::setObservables(PropertySetType& plist)
{
  copy(values_.begin(), values_.end(), plist.begin() + my_index_);
}

void ForwardWalking::setParticlePropertyList(PropertySetType& plist, int offset)
{
  copy(values_.begin(), values_.end(), plist.begin() + my_index_ + offset);
}
} // namespace qmcplusplus
