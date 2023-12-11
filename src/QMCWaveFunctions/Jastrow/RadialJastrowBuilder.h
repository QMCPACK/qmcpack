//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H
#define QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include <string>


namespace qmcplusplus
{
/** JastrowBuilder using an analytic 1d functor
 * Should be able to eventually handle all one and two body jastrows
 * although spline based ones will come later
 */


class RadialJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  enum detail
  {
    CPU,
    CUDA,
    OMPTARGET
  };

  // one body constructor
  RadialJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source);
  // two body constructor
  RadialJastrowBuilder(Communicate* comm, ParticleSet& target);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  /// \xmla{jastrow,name}
  std::string NameOpt;
  /// \xmla{jastrow,type}
  std::string TypeOpt;
  /// \xmla{jastrow,function}
  std::string Jastfunction;
  /// \xmla{jastrow,spin}
  std::string SpinOpt;
  ///particle set for source particle
  ParticleSet* SourcePtcl;

  // has a specialization for RPAFunctor in cpp file
  template<class RadFuncType, bool SPIN = false, unsigned Implementation = detail::CPU>
  std::unique_ptr<WaveFunctionComponent> createJ1(xmlNodePtr cur);

  template<class RadFuncType, unsigned Implementation = detail::CPU>
  std::unique_ptr<WaveFunctionComponent> createJ2(xmlNodePtr cur);

  template<class RadFuncType>
  void initTwoBodyFunctor(RadFuncType& functor, double fac);

  template<class RadFuncType>
  void computeJ2uk(const std::vector<RadFuncType*>& functors);

  void guardAgainstOBC();
  void guardAgainstPBC();
};

} // namespace qmcplusplus

#endif
