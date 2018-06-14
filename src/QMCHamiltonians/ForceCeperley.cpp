//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/ForceCeperley.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

ForceCeperley::ForceCeperley(ParticleSet& ions, ParticleSet& elns):
  ForceBase(ions, elns)
{
  ReportEngine PRE("ForceCeperley","ForceCeperley");
  myName = "Ceperley_Force_Base";
  prefix="HFCep";
  // Defaults
  Rcut = 0.4;
  m_exp = 2;
  N_basis = 4;
  forces = 0.0;
  forces_ShortRange.resize(Nnuc);
  forces_ShortRange = 0.0;
  ///////////////////////////////////////////////////////////////
}

void ForceCeperley::InitMatrix()
{
  Sinv.resize(N_basis, N_basis);
  h.resize(N_basis);
  c.resize(N_basis);
  for(int k=0; k<N_basis; k++)
  {
    h[k] = std::pow(Rcut, (k+2))/static_cast<RealType>(k+2);
    for(int j=0; j<N_basis; j++)
    {
      Sinv(k,j) = std::pow(Rcut, (m_exp+k+j+3))/static_cast<RealType>(m_exp+k+j+3);
    }
  }
  // in Numerics/DeterminantOperators.h
  invert_matrix(Sinv, false);
  // in Numerics/MatrixOperators.h
  MatrixOperators::product(Sinv, h.data(), c.data());
}

ForceCeperley::Return_t
ForceCeperley::evaluate(ParticleSet& P)
{
  forces = forces_IonIon;
  forces_ShortRange = forces_IonIon;
  const DistanceTableData* d_ab=P.DistTables[myTableIndex];
  const ParticleScalar_t* restrict Zat=Ions.Z.first_address();
  const ParticleScalar_t* restrict Qat=P.Z.first_address();

  if(d_ab->DTType == DT_SOA)
  {
    for(int jat=0; jat<Nel; jat++)
    {
      for(int iat=0; iat<Nnuc; iat++)
      {
      // electron contribution (special treatment if distance is inside cutoff!)
        RealType r = d_ab->Distances[jat][iat];
        RealType zoverr3=Qat[jat]*Zat[iat]/(r*r*r);
        if(r>=Rcut)
        {
          forces[iat] += zoverr3*d_ab->Displacements[jat][iat];
        }
        else
        {
          RealType g_q=0.0;
          for(int q=0; q<N_basis; q++)
          {
            g_q += c[q] * std::pow(r, m_exp+q+1);
          }
          g_q *= zoverr3;
          // negative sign accounts for definition of target as electrons
          forces[iat] += g_q*d_ab->Displacements[jat][iat];
          forces_ShortRange[iat] += g_q*d_ab->Displacements[jat][iat];
        }
      }
    }
  }
  else
  {
    for(int iat=0; iat<Nnuc; iat++)
    {
      RealType esum = 0.0;
      // electron contribution (special treatment if distance is inside cutoff!)
      for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++)
      {
        RealType zoverr3=Qat[jat]*Zat[iat]*std::pow(d_ab->rinv(nn),3);
        if(d_ab->r(nn)>=Rcut)
        {
          forces[iat] -= zoverr3*d_ab->dr(nn);
        }
        else
        {
          RealType r = d_ab->r(nn);
          RealType g_q=0.0;
          for(int q=0; q<N_basis; q++)
          {
            g_q += c[q] * std::pow(r, m_exp+q+1);
          }
          g_q *= zoverr3;
          // negative sign accounts for definition of target as electrons
          forces[iat] -= g_q*d_ab->dr(nn);
          forces_ShortRange[iat] -= g_q*d_ab->dr(nn);
        }
      }
    }
  }
  return 0.0;
}

bool ForceCeperley::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix, "name");
  attr.add(ionionforce, "addionion");
  attr.put(cur);
  addionion = (ionionforce=="yes") || (ionionforce == "true");
  app_log() << "ionionforce = "<<ionionforce<< std::endl;
  app_log() << "addionion="<<addionion<< std::endl;
  app_log() << "FirstTime= "<<FirstTime<< std::endl;
  ParameterSet fcep_param_set;
  fcep_param_set.add(Rcut, "rcut","real");
  fcep_param_set.add(N_basis, "nbasis", "int");
  fcep_param_set.add(m_exp, "weight_exp", "int");
  fcep_param_set.put(cur);
  app_log() <<"    ForceCeperley Parameters"<< std::endl;
  app_log() <<"        ForceCeperley::Rcut="<<Rcut<< std::endl;
  app_log() <<"        ForceCeperley::N_basis="<<N_basis<< std::endl;
  app_log() <<"        ForceCeperley::m_exp="<<m_exp<< std::endl;
  InitMatrix();
  return true;
}

QMCHamiltonianBase* ForceCeperley::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return new ForceCeperley(*this);
}
}

//  void ForceCeperley::addObservables(PropertySetType& plist) {
//    //cerr << "ForceBase::addObs sound off" << std::endl;
//    //obsName << myName << "0_x";
//    //myIndex = plist.add(obsName.str());
//    //obsName.clear();
//    mySize = Nnuc*OHMMS_DIM;
//    //cerr << "ForceBase mySize is " << Nnuc << " * " << OHMMS_DIM << " = " << mySize << std::endl;
//    checkInit = true;
//    if(myIndex<0) myIndex=plist.size();
//    int tmpIndex;
//    bool firstTime = true;
//    for(int iat=0; iat<Nnuc; iat++) {
//      for(int x=0; x<OHMMS_DIM; x++) {
//        std::ostringstream obsName;
//        obsName << "HFCep_" << iat << "_" << x;
//        tmpIndex = plist.add(obsName.str());
//        //if(firstTime) {
//        //  firstTime = false;
//        //  myIndex = tmpIndex;// + 3;
//        //}
//        std::cerr << iat << ", " << x << " stored at " << tmpIndex << std::endl;
//      }
//    }
//    std::cerr << "AddObs myIndex is " << myIndex << " last " << tmpIndex << std::endl;
//  }
//  // debugging version only z component
//  //void ForceCeperley::addObservables(PropertySetType& plist) {
//  //  mySize = Nnuc*OHMMS_DIM;
//  //  checkInit = true;
//  //  int tmpIndex;
//  //  bool firstTime = true;
//  //  for(int iat=0; iat<Nnuc; iat++) {
//  //    //for(int x=0; x<OHMMS_DIM; x++) {
//  //      std::ostringstream obsName;
//  //      obsName << "HFCep_" << iat << "_Z_sr";
//  //      tmpIndex = plist.add(obsName.str());
//  //      if(firstTime) {
//  //        firstTime = false;
//  //        myIndex = tmpIndex;
//  //        std::cerr << "ForceCeperley addObs setting myindex " << myIndex << std::endl;
//  //      }
//  //      std::ostringstream obsName2;
//  //      obsName2 << "HFCep_" << iat << "_Z_lr";
//  //      tmpIndex = plist.add(obsName2.str());
//  //      //cerr << iat << ", " << x << " stored at " << tmpIndex << std::endl;
//  //    //}
//  //  }
//  //}
//
//  //// overriding base class definition to print out short and long range contribs
//  //void ForceCeperley::setObservables(PropertySetType& plist) {
//  //  //cerr << "ForceBase::setObs storing forces ";
//  //  int index = myIndex;
//  //  for(int iat=0; iat<Nnuc; iat++) {
//  //    // HACK only z component
//  //    //for(int x=0; x<OHMMS_DIM; x++) {
//  //    plist[index] = forces(iat,2);
//  //    //cerr << index << ": " << plist[index] << "; ";
//  //    index++;
//  //    plist[index] = storeForces(iat,2);
//  //    //cerr << index << ": " << plist[index] << "; ";
//  //    index++;
//  //    //}
//  //  }
//  //  //cerr << std::endl;
//  //}
//  void ForceCeperley::setObservables(PropertySetType& plist) {
//    ///cerr << "ForceBase::setObs storing forces";
//    int index = myIndex;
//    for(int iat=0; iat<Nnuc; iat++) {
//      for(int x=0; x<OHMMS_DIM; x++) {
//        plist[index] = forces(iat,x) + storeForces(iat,x);
//        //cerr << " " << index << ":" << plist[index];
//        index++;
//      }
//    }
//    //cerr << std::endl;
//  }

