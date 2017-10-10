//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_LOCAL_MOMENT_HAMILTONIAN_H
#define QMCPLUSPLUS_LOCAL_MOMENT_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{

/** local moment estimator
 *
 * Compute local moment for up - down spins and all ion types in some particleset
 */
class LocalMomentEstimator: public QMCHamiltonianBase
{
public:

  /** constructor
   * @param elns target particle set
   * @param sources source particle set
   *
   * Use the elns.DistTables to evaluate the pair correlation functions.
   */
  LocalMomentEstimator(ParticleSet& elns, ParticleSet& srcs);

  void resetTargetParticleSet(ParticleSet& P);

  /* evaluate the pair correlation functions */
  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.size();
    int dm=std::floor(Dmax*100);
    for (int i=0; i<lm.size1()*lm.size2(); i++)
    {
      std::stringstream sstr;
      sstr << "lm" <<dm<<names[i];
      plist.add(sstr.str());
    }
  }
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(std::vector<observable_helper*>& h5list, hid_t gid) const;
  void setObservables(PropertySetType& plist)
  {
    int k(0);
    for(int i(0); i<lm.size1(); i++)
      for(int j(0); j<lm.size2(); j++,k++)
        plist[myIndex+k]=lm(i,j);
  }
  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    int k(0);
    for(int i(0); i<lm.size1(); i++)
      for(int j(0); j<lm.size2(); j++,k++)
        plist[myIndex+k+offset]=lm(i,j);
  }
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

private:
  /// maximum distance
  RealType Dmax;
  Matrix<RealType> lm;
  std::vector<int> ion_id, el_id;
  std::vector<RealType> el_nrm;
  int nag;
  DistanceTableData* d_table;
  std::vector<std::string> names;
  ParticleSet& ions;
};

}
#endif

