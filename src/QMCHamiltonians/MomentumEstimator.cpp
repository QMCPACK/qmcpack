//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCHamiltonians/MomentumEstimator.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <Numerics/e2iphi.h>
#include <Numerics/OhmmsBlas.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <Particle/DistanceTableData.h>
#include <Numerics/DeterminantOperators.h>
#include <set>

namespace qmcplusplus
{

MomentumEstimator::MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi)
  :M(1), refPsi(psi), kgrid(4), Lattice(elns.Lattice), norm_nofK(1), hdf5_out(false)
{
  UpdateMode.set(COLLECTABLE,1);
  psi_ratios.resize(elns.getTotalNum());
  twist=elns.getTwist();
}

void MomentumEstimator::resetTargetParticleSet(ParticleSet& P)
{
}

MomentumEstimator::Return_t MomentumEstimator::evaluate(ParticleSet& P)
{
  const int np=P.getTotalNum();
  const int nk=kPoints.size();
  for (int s=0; s<M; ++s)
  {
    PosType newpos;
    for (int i=0; i<OHMMS_DIM; ++i)
      newpos[i]=myRNG();
    //make it cartesian
    vPos[s]=Lattice.toCart(newpos);
    P.makeVirtualMoves(vPos[s]);
    refPsi.evaluateRatiosAlltoOne(P,psi_ratios);
    for (int i=0; i<np; ++i)
      psi_ratios_all[s][i] = psi_ratios[i];

    for (int ik=0; ik<nk; ++ik)
       kdotp[ik] = -dot(kPoints[ik], vPos[s]);
    eval_e2iphi(nk, kdotp.data(), phases_vPos[s].data(0), phases_vPos[s].data(1));
  }

  std::fill_n(nofK.begin(),nk,RealType(0));
  for (int i=0; i<np; ++i)
  {
    for (int ik=0; ik<nk; ++ik)
      kdotp[ik] = dot(kPoints[ik], P.R[i]);
    eval_e2iphi(nk, kdotp.data(), phases.data(0), phases.data(1));
    for (int s=0; s<M; ++s)
    {
      const ComplexType one_ratio(psi_ratios_all[s][i]);
      const RealType ratio_c = one_ratio.real();
      const RealType ratio_s = one_ratio.imag();
      const RealType *restrict phases_c = phases.data(0);
      const RealType *restrict phases_s = phases.data(1);
      const RealType *restrict phases_vPos_c = phases_vPos[s].data(0);
      const RealType *restrict phases_vPos_s = phases_vPos[s].data(1);
      RealType *restrict nofK_here=nofK.data();
      #pragma omp simd aligned(nofK_here,phases_c,phases_s,phases_vPos_c,phases_vPos_s)
      for (int ik=0; ik<nk; ++ik)
        nofK_here[ik] += ( phases_c[ik]*phases_vPos_c[ik] - phases_s[ik]*phases_vPos_s[ik] ) * ratio_c
                       - ( phases_s[ik]*phases_vPos_c[ik] + phases_c[ik]*phases_vPos_s[ik] ) * ratio_s ;
    }
  }

  compQ=0.0;
  for (int iq=0; iq<compQ.size(); ++iq)
  {
    for (int i=0; i<mappedQtonofK[iq].size(); ++i)
      compQ[iq] += nofK[mappedQtonofK[iq][i]];
    compQ[iq] *= mappedQnorms[iq];
  }
  for (int ik=0; ik<nofK.size(); ++ik)
    nofK[ik] *= norm_nofK;
  if (hdf5_out)
  {
    RealType w=tWalker->Weight;
    int j=myIndex;
    for (int ik=0; ik<nofK.size(); ++ik,++j)
      P.Collectables[j]+= w*nofK[ik];
    for (int iq=0; iq<compQ.size(); ++iq,++j)
      P.Collectables[j]+= w*compQ[iq];
  }
  return 0.0;
}

void MomentumEstimator::registerCollectables(std::vector<observable_helper*>& h5desc
    , hid_t gid) const
{
  if (hdf5_out)
  {
    //descriptor for the data, 1-D data
    std::vector<int> ng(1);
    //add nofk
    ng[0]=nofK.size();
    observable_helper* h5o=new observable_helper("nofk");
    h5o->set_dimensions(ng,myIndex);
    h5o->open(gid);
    h5o->addProperty(const_cast<std::vector<PosType>&>(kPoints),"kpoints");
    h5o->addProperty(const_cast<std::vector<int>&>(kWeights),"kweights");
    h5desc.push_back(h5o);
    //add compQ
    ng[0]=Q.size();
    h5o=new observable_helper("compQ");
    h5o->set_dimensions(ng,myIndex+nofK.size());
    h5o->open(gid);
    h5o->addProperty(const_cast<std::vector<RealType>&>(Q),"q");
    h5desc.push_back(h5o);
  }
}


void MomentumEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (hdf5_out)
  {
    myIndex=collectables.size();
    collectables.add(nofK.begin(),nofK.end());
    collectables.add(compQ.begin(),compQ.end());
  }
  else
  {
    myIndex=plist.size();
    for (int i=0; i<nofK.size(); i++)
    {
      std::stringstream sstr;
      sstr << "nofk_" <<i;
      int id=plist.add(sstr.str());
    }
    for (int i=0; i<Q.size(); i++)
    {
      std::stringstream sstr;
      sstr << "Q_" <<i;
      int id=plist.add(sstr.str());
    }
  }
}



void MomentumEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex);
    copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size());
  }
}

void MomentumEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex+offset);
    copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size()+offset);
  }
}

bool MomentumEstimator::putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode)
{
  OhmmsAttributeSet pAttrib;
  std::string hdf5_flag="yes";
  pAttrib.add(hdf5_flag,"hdf5");
  pAttrib.add(kgrid,"grid");
  pAttrib.add(M,"samples");
  pAttrib.put(cur);
  hdf5_out = (hdf5_flag=="yes");
#if OHMMS_DIM==3
  int numqtwists(6*kgrid+3);
  std::vector<int> qk(0);
  mappedQtonofK.resize(numqtwists,qk);
  compQ.resize(numqtwists);
  RealType qn(4.0*M_PI*M_PI*std::pow(Lattice.Volume,-2.0/3.0));
  mappedQnorms.resize(numqtwists,qn*0.5/RealType(M));
  if (twist[0]==0)
    mappedQnorms[kgrid]=qn/RealType(M);
  if (twist[1]==0)
    mappedQnorms[3*kgrid+1]=qn/RealType(M);
  if (twist[2]==0)
    mappedQnorms[5*kgrid+2]=qn/RealType(M);
  Q.resize(numqtwists);
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    PosType kpt;
    kpt[0]=i-twist[0];
    kpt[1]=i-twist[1];
    kpt[2]=i-twist[2];
    kpt=Lattice.k_cart(kpt);
    Q[i+kgrid]=std::abs(kpt[0]);
    Q[i+kgrid+(2*kgrid+1)]=std::abs(kpt[1]);
    Q[i+kgrid+(4*kgrid+2)]=std::abs(kpt[2]);
  }
  app_log()<<" Using all k-space points with (nx^2+ny^2+nz^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<< std::endl;
  app_log()<<"  My twist is:"<<twist[0]<<"  "<<twist[1]<<"  "<<twist[2]<< std::endl;
  int indx(0);
  int kgrid_squared=kgrid*kgrid;
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    for (int j=-kgrid; j<(kgrid+1); j++)
    {
      for (int k=-kgrid; k<(kgrid+1); k++)
      {
        if (i*i+j*j+k*k<=kgrid_squared) //if (std::sqrt(i*i+j*j+k*k)<=kgrid)
        {
          PosType kpt;
          kpt[0]=i-twist[0];
          kpt[1]=j-twist[1];
          kpt[2]=k-twist[2];
          //convert to Cartesian: note that 2Pi is multiplied
          kpt=Lattice.k_cart(kpt);
          kPoints.push_back(kpt);
          mappedQtonofK[i+kgrid].push_back(indx);
          mappedQtonofK[j+kgrid+(2*kgrid+1)].push_back(indx);
          mappedQtonofK[k+kgrid+(4*kgrid+2)].push_back(indx);
          indx++;
        }
      }
    }
  }
#endif
#if OHMMS_DIM==2
  int numqtwists(4*kgrid+2);
  std::vector<int> qk(0);
  mappedQtonofK.resize(numqtwists,qk);
  compQ.resize(numqtwists);
  RealType qn(2.0*M_PI/std::sqrt(Lattice.Volume));
  mappedQnorms.resize(numqtwists,qn*0.5/RealType(M));
  if (twist[0]==0)
    mappedQnorms[kgrid]=qn/RealType(M);
  if (twist[1]==0)
    mappedQnorms[3*kgrid+1]=qn/RealType(M);
  Q.resize(numqtwists);
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    PosType kpt;
    kpt[0]=i-twist[0];
    kpt[1]=i-twist[1];
    kpt=Lattice.k_cart(kpt);
    Q[i+kgrid]=std::abs(kpt[0]);
    Q[i+kgrid+(2*kgrid+1)]=std::abs(kpt[1]);
  }
  app_log()<<" Using all k-space points with (nx^2+ny^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<< std::endl;
  app_log()<<"  My twist is:"<<twist[0]<<"  "<<twist[1]<< std::endl;
  int indx(0);
  int kgrid_squared=kgrid*kgrid;
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    for (int j=-kgrid; j<(kgrid+1); j++)
    {
      if (i*i+j*j<=kgrid_squared) //if (std::sqrt(i*i+j*j+k*k)<=kgrid)
      {
        PosType kpt;
        kpt[0]=i-twist[0];
        kpt[1]=j-twist[1];
        //convert to Cartesian: note that 2Pi is multiplied
        kpt=Lattice.k_cart(kpt);
        kPoints.push_back(kpt);
        mappedQtonofK[i+kgrid].push_back(indx);
        mappedQtonofK[j+kgrid+(2*kgrid+1)].push_back(indx);
        indx++;
      }
    }
  }
#endif
  if (rootNode)
  {
    std::stringstream sstr;
    sstr<<"Kpoints";
    for(int i(0); i<OHMMS_DIM; i++)
      sstr<<"_"<<round(100.0*twist[i]);
    sstr<<".dat";
    std::ofstream fout(sstr.str().c_str());
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout << "# mag_k        ";
    for(int i(0); i<OHMMS_DIM; i++)
      fout << "k_"<<i<<"           ";
    fout << std::endl;
    for (int i=0; i<kPoints.size(); i++)
    {
      float khere(std::sqrt(dot(kPoints[i],kPoints[i])));
      fout<<khere;
      for(int j(0); j<OHMMS_DIM; j++)
        fout<<"   "<<kPoints[i][j];
      fout<< std::endl;
    }
    fout.close();
    sstr.str("");
    sstr<<"Qpoints";
    for(int i(0); i<OHMMS_DIM; i++)
      sstr<<"_"<<round(100.0*twist[i]);
    sstr<<".dat";
    std::ofstream qout(sstr.str().c_str());
    qout.setf(std::ios::scientific, std::ios::floatfield);
    qout << "# mag_q" << std::endl;
    for (int i=0; i<Q.size(); i++)
    {
      qout<<Q[i]<< std::endl;
    }
    qout.close();
  }
  nofK.resize(kPoints.size());
  kdotp.resize(kPoints.size());
  vPos.resize(M);
  phases.resize(kPoints.size());
  phases_vPos.resize(M);
  for(int im=0; im<M; im++)
    phases_vPos[im].resize(kPoints.size());
  psi_ratios_all.resize(M,psi_ratios.size());
  norm_nofK=1.0/RealType(M);
  return true;
}

bool MomentumEstimator::get(std::ostream& os) const
{
  return true;
}

QMCHamiltonianBase* MomentumEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  MomentumEstimator* myclone=new MomentumEstimator(qp,psi);
  myclone->resize(kPoints,Q,M);
  myclone->myIndex=myIndex;
  myclone->kgrid=kgrid;
  myclone->norm_nofK=norm_nofK;
  myclone->mappedQtonofK.resize(mappedQtonofK.size());
  for(int i=0; i<mappedQtonofK.size(); i++)
    myclone->mappedQtonofK[i]=mappedQtonofK[i];
  myclone->hdf5_out=hdf5_out;
  myclone->mappedQnorms=mappedQnorms;
  return myclone;
}

void MomentumEstimator::resize(const std::vector<PosType>& kin, const std::vector<RealType>& qin, const int Min)
{
  //copy kpoints
  kPoints=kin;
  nofK.resize(kin.size());
  kdotp.resize(kPoints.size());
  phases.resize(kPoints.size());
  //copy q
  Q=qin;
  compQ.resize(qin.size());
  //M
  M=Min;
  vPos.resize(M);
  psi_ratios_all.resize(M,psi_ratios.size());
  phases_vPos.resize(M);
  for(int im=0; im<M; im++)
    phases_vPos[im].resize(kPoints.size());
}

void MomentumEstimator::setRandomGenerator(RandomGenerator_t* rng)
{
  //simply copy it
  myRNG=*rng;
}
}

