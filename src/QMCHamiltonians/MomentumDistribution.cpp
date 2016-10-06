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
    
    
#include "QMCHamiltonians/MomentumDistribution.h"

namespace qmcplusplus
{
typedef QMCTraits::IndexType    IndexType;
typedef QMCTraits::PosType      PosType;
typedef QMCTraits::RealType     RealType;
typedef QMCTraits::ValueType    ValueType;
typedef QMCTraits::ComplexType  ComplexType;

template<> IndexType FFTMomentumDist<1>::getIndex(const TinyVector<IndexType, 1>& LocInDist)
{
  return LocInDist[0];
}

template<> IndexType FFTMomentumDist<2>::getIndex(const TinyVector<IndexType, 2>& LocInDist)
{
  return LocInDist[0] * NumPts[1] + LocInDist[1];
}

template<> IndexType FFTMomentumDist<3>::getIndex(const TinyVector<IndexType, 3>& LocInDist)
{
  return (LocInDist[0] * NumPts[1] + LocInDist[1]) * NumPts[2] + LocInDist[2];
}

template<> void FFTMomentumDist<1>::PlaceCloudInBin(const TinyVector<RealType, 1>& FracDisp,
    ValueType RatioOfWvFcn)
{
  TinyVector<IndexType, 1> IndBelow;
  TinyVector<IndexType, 1> IndAbove;
  RealType ScaledLoc = FracDisp[1] * NumPts[0];
  IndBelow[0] = IndexType(std::floor(ScaledLoc));
  RealType FracBelowR = ScaledLoc - RealType(IndBelow[0]);
  RealType FracAboveR = 1.0 - FracBelowR;
  if (IndBelow[0] < (NumPts[0] - 1))
  {
    IndAbove[0] = IndBelow[0] + 1;
  }
  else
  {
    IndAbove[0] = 0;
  }
  placeIntsInBin(IndAbove, RatioOfWvFcn * FracBelowR);
  placeIntsInBin(IndBelow, RatioOfWvFcn * FracAboveR);
  totalNumSamples++;
}

template<> void FFTMomentumDist<2>::PlaceCloudInBin(const TinyVector<RealType, 2>& FracDisp,
    ValueType RatioOfWvFcn)
{
  TinyVector<RealType, 2> FracBelowR;
  TinyVector<RealType, 2> FracAboveR;
  TinyVector<IndexType, 2> IndBelow;
  TinyVector<IndexType, 2> IndAbove;
  for (IndexType i = 0; i < 2; i++)
  {
    RealType ScaledLoc = FracDisp[i] * NumPts[i];
    IndBelow[i] = IndexType(std::floor(ScaledLoc));
    FracBelowR[i] = ScaledLoc - RealType(IndBelow[i]);
    FracAboveR[i] = 1.0 - FracBelowR[i];
    if (IndBelow[i] < (NumPts[i] - 1))
    {
      IndAbove[i] = IndBelow[i] + 1;
    }
    else
    {
      IndAbove[i] = 0;
    }
  }
  placeIntsInBin(IndBelow, RatioOfWvFcn * FracAboveR[0] * FracAboveR[1]);
  placeIntsInBin(TinyVector<IndexType, 2>(IndBelow[0], IndAbove[1]), RatioOfWvFcn * FracAboveR[0] * FracBelowR[1]);
  placeIntsInBin(TinyVector<IndexType, 2>(IndAbove[0], IndBelow[1]), RatioOfWvFcn * FracBelowR[0] * FracAboveR[1]);
  placeIntsInBin(TinyVector<IndexType, 2>(IndAbove[0], IndAbove[1]), RatioOfWvFcn * FracBelowR[0] * FracBelowR[1]);
}

template<> void FFTMomentumDist<3>::PlaceCloudInBin(const TinyVector<RealType, 3>& FracDisp,
    ValueType RatioOfWvFcn)
{
  TinyVector<RealType, 3> FracBelowR;
  TinyVector<RealType, 3> FracAboveR;
  TinyVector<IndexType, 3> IndBelow;
  TinyVector<IndexType, 3> IndAbove;
  for (IndexType i = 0; i < 3; i++)
  {
    RealType ScaledLoc = FracDisp[i] * NumPts[i];
    IndBelow[i] = IndexType(std::floor(ScaledLoc));
    FracBelowR[i] = ScaledLoc - RealType(IndBelow[i]);
    FracAboveR[i] = 1.0 - FracBelowR[i];
    if (IndBelow[i] < (NumPts[i] - 1))
    {
      IndAbove[i] = IndBelow[i] + 1;
    }
    else
    {
      IndAbove[i] = 0;
    }
  }
  placeIntsInBin(IndBelow,
                 RatioOfWvFcn * FracAboveR[0] * FracAboveR[1] * FracAboveR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndBelow[0], IndBelow[1], IndAbove[2]),
                 RatioOfWvFcn * FracAboveR[0] * FracAboveR[1] * FracBelowR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndBelow[0], IndAbove[1], IndBelow[2]),
                 RatioOfWvFcn * FracAboveR[0] * FracBelowR[1] * FracAboveR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndBelow[0], IndAbove[1], IndAbove[2]),
                 RatioOfWvFcn * FracAboveR[0] * FracBelowR[1] * FracBelowR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndAbove[0], IndBelow[1], IndBelow[2]),
                 RatioOfWvFcn * FracBelowR[0] * FracAboveR[1] * FracAboveR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndAbove[0], IndBelow[1], IndAbove[2]),
                 RatioOfWvFcn * FracBelowR[0] * FracAboveR[1] * FracBelowR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndAbove[0], IndAbove[1], IndBelow[2]),
                 RatioOfWvFcn * FracBelowR[0] * FracBelowR[1] * FracAboveR[2]);
  placeIntsInBin(TinyVector<IndexType, 3>(IndAbove[0], IndAbove[1], IndAbove[2]),
                 RatioOfWvFcn * FracBelowR[0] * FracBelowR[1] * FracBelowR[2]);
  totalNumSamples++;
}


/*
 * Three Dimensional Momentum Distribution.  This will take displacements in three dimensions
 * and produce a three dimensional n(k)
 */
ThreeDimMomDist::ThreeDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes,
                                 PtclChoiceBase* pcb) : FFTMomentumDist<3>(pcb, PtclSet.getTotalNum())
{
  TinyVector<IndexType, 3> inNumPts;
  TinyVector<RealType, 3> InvSpacing;
  for (IndexType i = 0; i < 3; i++)
  {
    inNumPts[i] = DimSizes[i];
    InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
  }
  initialize(inNumPts, InvSpacing);
}

void ThreeDimMomDist::updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi,
    IndexType NumCycles)
{
  TinyVector<IndexType, 3> Indexes;
  for (IndexType cycle = 0; cycle < NumCycles; cycle++)
  {
    pcp->NewWalker();
    PosType dr;
    for (Indexes[0] = 0; Indexes[0] < NumPts[0]; Indexes[0]++)
    {
      dr[1] = 0.0;
      for (Indexes[1] = 0; Indexes[1] < NumPts[1]; Indexes[1]++)
      {
        dr[2] = 0.0;
        for (Indexes[2] = 0; Indexes[2] < NumPts[2]; Indexes[2]++)
        {
          IndexType partToDisplace = (*pcp)();
          PtclSet.makeMove(partToDisplace, dr);
          if (Indexes[0] == 0 && Indexes[1] == 0 && Indexes[2] == 0)
          {
            placeIntsInBin(Indexes, 1.0);
          }
          else
          {
            placeIntsInBin(Indexes, Psi.ratio(PtclSet, partToDisplace));
          }
          totalNumSamples++;
          PtclSet.rejectMove(partToDisplace);
          dr[2] += Spacing[2];
        }
        dr[1] += Spacing[1];
      }
      dr[0] += Spacing[0];
    }
  }
}

/*
  void ThreeDimMomDist::updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumCycles) {

   for (IndexType i = 0; i < NumCycles; i++) {
      IndexType PtclToDisplace = (*pcp)();
      PosType fracDisp;
      for (IndexType j = 0; j < 3; j++) {
	fracDisp[j] = Random();
      }
      PosType dr = PtclSet.Lattice.toCart(fracDisp);

      PtclSet.makeMove(PtclToDisplace, dr);
      PlaceCloudInBin(fracDisp, Psi.ratio(PtclSet, PtclToDisplace));
      totalNumSamples++;
      PtclSet.rejectMove(PtclToDisplace);
    }
  }
 */


/*
 * One Dimensional Momentum Distribution.  This will take displacements in one dimension
 * and produce a one dimensional n(k)
 */

OneDimMomDist::OneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes,
                             PtclChoiceBase* pcb) : FFTMomentumDist<1>(pcb, PtclSet.getTotalNum())
{
  TinyVector<IndexType, 1> inIndSize;
  TinyVector<RealType, 1> InvSpacing;
  for (IndexType i = 0; i < 1; i++)
  {
    inIndSize[i] = DimSizes[i];
    InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
  }
  initialize(inIndSize, InvSpacing);
}

void OneDimMomDist::updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi,
                                       IndexType NumCycles)
{
  for (IndexType cycle = 0; cycle < NumCycles; cycle++)
  {
    pcp->NewWalker();
    TinyVector<IndexType, 1> Indexes(0.0);
    PosType dr(0,0,0);
    placeIntsInBin(Indexes, 1.0);
    totalNumSamples++;
    for (Indexes[0] = 1; Indexes[0] < NumPts[0]; Indexes[0]++)
    {
      dr[0] = dr[1] = dr[2] += Spacing[0];
//        dr[0] += Spacing[0];
      IndexType partToDisplace = (*pcp)();
      PtclSet.makeMove(partToDisplace, dr);
      placeIntsInBin(Indexes, Psi.ratio(PtclSet, partToDisplace));
      totalNumSamples++;
      PtclSet.rejectMove(partToDisplace);
    }
  }
}

/*
 * One Dimensional Averaged Momentum Distribution.  This will take displacements in three dimensions
 * and produce a one dimensional n(k)
 */
AveragedOneDimMomDist::AveragedOneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes,
    PtclChoiceBase* pcb) : FFTMomentumDist<1>(pcb, PtclSet.getTotalNum())
{
  TinyVector<IndexType, 1> inIndSize;
  TinyVector<RealType, 1> InvSpacing;
  for (IndexType i = 0; i < 1; i++)
  {
    inIndSize[i] = DimSizes[i];
    InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
  }
  initialize(inIndSize, InvSpacing);
}

void AveragedOneDimMomDist::updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi,
    IndexType NumCycles)
{
  IndexType targetNumSamples = totalNumSamples + NumCycles;
  while(totalNumSamples < targetNumSamples)
  {
    PosType dr;
    for (IndexType i = 0; i < 3; i++)
    {
      dr += Random() * PtclSet.Lattice.a(i);
    }
    RealType DispVecLength = std::sqrt(dot(dr, dr));
    RealType FracDispVecLength = InvRsLatSize[0] * DispVecLength / 2.0;
    IndexType PtclToDisplace = (*pcp)();
    PtclSet.makeMove(PtclToDisplace, dr);
    PlaceCloudInBin(FracDispVecLength, Psi.ratio(PtclSet, PtclToDisplace));
    totalNumSamples++;
    PtclSet.rejectMove(PtclToDisplace);
  }
}

/*
 * Base class implementation for all Direct Momentum Distributions
 */

Vector<RealType>& DirectMomDist::getNofK()
{
  RealType factor = static_cast<RealType>(NumParticles) / static_cast<RealType>(totalNumSamples);
  for (IndexType GVecNum = 0; GVecNum < GVectors.size(); GVecNum++)
  {
    MomDist[GVecNum] = NofK[GVecNum].real() * factor;
//      std::cout << MomDist[GVecNum] << std::endl;
  }
//    std::cout << std::endl;
  return MomDist;
}

/*
 * Random Momentum Distribution.  This will take displacements in three dimensions
 * and produce a one dimensional n(k) by directly evaluating e^iGr for a list of G
 * vectors that is given to the constructor
 */

RandomMomDist::RandomMomDist(ParticleSet& PtclSet, const Vector<PosType>& inGVectors,
                             PtclChoiceBase* pcb) : DirectMomDist(pcb, PtclSet.getTotalNum())
{
  GVectors.resize(inGVectors.size());
  for (IndexType i = 0; i < inGVectors.size(); i++)
  {
    GVectors[i] = inGVectors[i];
  }
  NofK.resize(inGVectors.size());
  MomDist.resize(inGVectors.size());
}

void RandomMomDist::updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples)
{
  pcp->NewWalker();
  for (IndexType i = 0; i < NumSamples; i++)
  {
    PosType dr = PtclSet.Lattice.toCart(PosType(Random(), Random(), Random()));
    IndexType PtclToDisplace = (*pcp)();
    PtclSet.makeMove(PtclToDisplace, dr);
    ComplexType ratio = Psi.ratio(PtclSet, PtclToDisplace);
    PtclSet.rejectMove(PtclToDisplace);
    totalNumSamples++;
    for (IndexType GVecNum = 0; GVecNum < GVectors.size(); GVecNum++)
    {
      NofK[GVecNum] += ratio * exp(ComplexType(0,-1.0) * dot(GVectors[GVecNum], dr));
    }
  }
}

}
