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
  :M(40), refPsi(psi), Lattice(elns.Lattice), norm_nofK(1), hdf5_out(false)
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
  if (hdf5_out)
  {
    RealType w=tWalker->Weight*norm_nofK;
    int j=myIndex;
    for (int ik=0; ik<nofK.size(); ++ik,++j)
      P.Collectables[j]+= w*nofK[ik];
  }
  else
  {
    for (int ik=0; ik<nofK.size(); ++ik)
      nofK[ik] *= norm_nofK; 
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
  }
}


void MomentumEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (hdf5_out)
  {
    myIndex=collectables.size();
    collectables.add(nofK.begin(),nofK.end());
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
  }
}



void MomentumEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex);
  }
}

void MomentumEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex+offset);
  }
}

bool MomentumEstimator::putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode)
{
  OhmmsAttributeSet pAttrib;
  std::string hdf5_flag="yes";
  ///dims of a grid for generating k points (obtained below)
  int kgrid=0;
  //maximum k-value in the k-grid in cartesian coordinates
  RealType kmax=0.0;
  //maximum k-values in the k-grid along the reciprocal cell axis
  RealType kmax0=0.0;
  RealType kmax1=0.0;
  RealType kmax2=0.0;
  pAttrib.add(hdf5_flag,"hdf5");
  //pAttrib.add(kgrid,"grid");
  pAttrib.add(kmax,"kmax");
  pAttrib.add(kmax0,"kmax0");
  pAttrib.add(kmax1,"kmax1");
  pAttrib.add(kmax2,"kmax2");
  pAttrib.add(M,"samples"); // default value is 40 (in the constructor)
  pAttrib.put(cur);
  hdf5_out = (hdf5_flag=="yes");
  // minimal length as 2 x WS radius.
  RealType min_Length = elns.Lattice.WignerSeitzRadius_G*4.0*M_PI;
  PosType vec_length;
  //length of reciprocal lattice vector
  for (int i=0; i<OHMMS_DIM; i++)
    vec_length[i]=2.0*M_PI*std::sqrt(dot(elns.Lattice.Gv[i],elns.Lattice.Gv[i]));
#if OHMMS_DIM==3
  PosType kmaxs(0);
  kmaxs[0]=kmax0;
  kmaxs[1]=kmax1;
  kmaxs[2]=kmax2;
  RealType sum_kmaxs=kmaxs[0]+kmaxs[1]+kmaxs[2];
  RealType sphere_kmax;
  bool sphere = kmax>0.0 ? true:false;
  bool directional= sum_kmaxs>0.0 ? true:false;
  if (!sphere && !directional)
  {
    // default: kmax = 2 x k_F of polarized non-interacting electron system
    kmax = 2.0*std::pow(6.0*M_PI*M_PI*elns.getTotalNum()/elns.Lattice.Volume,1.0/3);
    sphere = true;
  }
  sphere_kmax=kmax;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    if (kmaxs[i]>kmax)
      kmax=kmaxs[i];
  }
  kgrid = int(kmax/min_Length)+1;
  RealType kgrid_squared[OHMMS_DIM];
  for (int i=0; i<OHMMS_DIM; i++)
    kgrid_squared[i]=kmaxs[i]*kmaxs[i]/vec_length[i]/vec_length[i];
  RealType kmax_squared=sphere_kmax*sphere_kmax;
  std::vector<int> kcount0;
  std::vector<int> kcount1;
  std::vector<int> kcount2;
  kcount0.resize(2*kgrid+1,0);
  kcount1.resize(2*kgrid+1,0);
  kcount2.resize(2*kgrid+1,0);
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    for (int j=-kgrid; j<(kgrid+1); j++)
    {
      for (int k=-kgrid; k<(kgrid+1); k++)
      {
        PosType ikpt,kpt;
        ikpt[0]=i-twist[0];
        ikpt[1]=j-twist[1];
        ikpt[2]=k-twist[2];
        //convert to Cartesian: note that 2Pi is multiplied
        kpt=Lattice.k_cart(ikpt);
        bool not_recorded=true;
        // This collects the k-points within the parallelepiped (if enabled)
        if (directional && ikpt[0]*ikpt[0]<=kgrid_squared[0] && ikpt[1]*ikpt[1]<=kgrid_squared[1] && ikpt[2]*ikpt[2]<=kgrid_squared[2])
        {
          kPoints.push_back(kpt);
          kcount0[kgrid+i]=1;
          kcount1[kgrid+j]=1;
          kcount2[kgrid+k]=1;
          not_recorded=false;
        }
        // This collects the k-points within a sphere (if enabled and the k-point has not been recorded yet)
        if (sphere && not_recorded && kpt[0]*kpt[0]+kpt[1]*kpt[1]+kpt[2]*kpt[2]<=kmax_squared) //if (std::sqrt(kx*kx+ky*ky+kz*kz)<=sphere_kmax)
        {
          kPoints.push_back(kpt);
        }
      }
    }
  }
  if (sphere && !directional)
  {
    app_log()<<"    Using all k-space points with (kx^2+ky^2+kz^2)^0.5 < "<< sphere_kmax <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution is "<< kPoints.size() << std::endl;
  }
  else if (directional && !sphere)
  {
    int sums[3];
    sums[0]=0;
    sums[1]=0;
    sums[2]=0;
    for (int i=0; i<2*kgrid+1; i++)
    {
      sums[0]+=kcount0[i];
      sums[1]+=kcount1[i];
      sums[2]+=kcount2[i];
    }
    app_log()<<"    Using all k-space points within cut-offs "<< kmax0 << ", " << kmax1 << ", " << kmax2 <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution: "<< kPoints.size() << std::endl;
    app_log()<<"      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log()<<"      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
    app_log()<<"      Number of grid points in kmax2 direction: " << sums[2] << std::endl;
  }
  else
  {
    int sums[3];
    sums[0]=0;
    sums[1]=0;
    sums[2]=0;
    for (int i=0; i<2*kgrid+1; i++)
    {
      sums[0]+=kcount0[i];
      sums[1]+=kcount1[i];
      sums[2]+=kcount2[i];
    }
    app_log()<<"    Using all k-space points with (kx^2+ky^2+kz^2)^0.5 < "<< sphere_kmax <<", and"<< std::endl;
    app_log()<<"    within the cut-offs "<< kmax0 << ", " << kmax1 << ", " << kmax2 <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution is "<< kPoints.size() << std::endl;
    app_log()<<"    The number of k-points within the cut-off region: "<< sums[0]*sums[1]*sums[2] << std::endl;
    app_log()<<"      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log()<<"      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
    app_log()<<"      Number of grid points in kmax2 direction: " << sums[2] << std::endl;
  }
  app_log()<<"    Number of samples: "<< M << std::endl;
  app_log()<<"    My twist is: "<<twist[0]<<"  "<<twist[1]<<"  "<<twist[2]<< std::endl;
#endif
#if OHMMS_DIM==2
  PosType kmaxs(0);
  kmaxs[0]=kmax0;
  kmaxs[1]=kmax1;
  RealType sum_kmaxs=kmaxs[0]+kmaxs[1];
  RealType disk_kmax;
  bool disk = kmax>0.0 ? true:false;
  bool directional = sum_kmaxs>0.0 ? true:false;
  if (!disk && !directional)
  {
    // default: kmax = 2 x k_F of polarized non-interacting electron system
    kmax = 2.0*std::pow(4.0*pi*elns.getTotalNum()/elns.Lattice.Volume,0.5);
    disk=true;
  }
  disk_kmax=kmax;
  for (int i=0; i<OHMMS_DIM; i++)
  {
    if (kmaxs[i]>kmax)
      kmax=kmaxs[i];
  }
  kgrid = int(kmax/min_Length)+1;
  RealType kgrid_squared[OHMMS_DIM];
  for (int i=0; i<OHMMS_DIM; i++)
    kgrid_squared[i]=kmaxs[i]*kmaxs[i]/vec_length[i]/vec_length[i];
  RealType kmax_squared=disk_kmax*disk_kmax;
  std::vector<int> kcount0;
  std::vector<int> kcount1;
  kcount0.resize(2*kgrid+1,0);
  kcount1.resize(2*kgrid+1,0);
  for (int i=-kgrid; i<(kgrid+1); i++)
  {
    for (int j=-kgrid; j<(kgrid+1); j++)
    {
      PosType ikpt,kpt;
      ikpt[0]=i-twist[0];
      ikpt[1]=j-twist[1];
      //convert to Cartesian: note that 2Pi is multiplied
      kpt=Lattice.k_cart(ikpt);
      bool not_recorded=true;
      if (directional && ikpt[0]*ikpt[0]<=kgrid_squared[0] && ikpt[1]*ikpt[1]<=kgrid_squared[1])
      {
        kPoints.push_back(kpt);
        kcount0[kgrid+i]=1;
        kcount1[kgrid+j]=1;
        not_recorded=false;
      }
      if (disk && not_recorded && kpt[0]*kpt[0]+kpt[1]*kpt[1]<=kmax_squared) //if (std::sqrt(kx*kx+ky*ky)<=disk_kmax)
      {
        kPoints.push_back(kpt);
      }
    }
  }
  if (disk && !directional)
  {
    app_log()<<"    Using all k-space points with (kx^2+ky^2)^0.5 < "<< disk_kmax <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution is "<< kPoints.size() << std::endl;
  }
  else if (directional && !disk)
  {
    int sums[2];
    sums[0]=0;
    sums[1]=0;
    for (int i=0; i<2*kgrid+1; i++)
    {
      sums[0]+=kcount0[i];
      sums[1]+=kcount1[i];
    }
    app_log()<<"    Using all k-space points within cut-offs "<< kmax0 << ", " << kmax1 <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution: "<< kPoints.size() << std::endl;
    app_log()<<"      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log()<<"      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
  }
  else
  {
    int sums[2];
    sums[0]=0;
    sums[1]=0;
    for (int i=0; i<2*kgrid+1; i++)
    {
      sums[0]+=kcount0[i];
      sums[1]+=kcount1[i];
    }
    app_log()<<"    Using all k-space points with (kx^2+ky^2)^0.5 < "<< disk_kmax <<", and"<< std::endl;
    app_log()<<"    within the cut-offs "<< kmax0 << ", " << kmax1 <<" for Momentum Distribution."<< std::endl;
    app_log()<<"    Total number of k-points for Momentum Distribution is "<< kPoints.size() << std::endl;
    app_log()<<"    The number of k-points within the cut-off region: "<< sums[0]*sums[1] << std::endl;
    app_log()<<"      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log()<<"      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
  }
  app_log()<<"    Number of samples: "<< M << std::endl;
  app_log()<<"    My twist is: "<<twist[0]<<"  "<<twist[1]<<"  "<<twist[2]<< std::endl;
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
  myclone->resize(kPoints,M);
  myclone->myIndex=myIndex;
  myclone->norm_nofK=norm_nofK;
  myclone->hdf5_out=hdf5_out;
  return myclone;
}

void MomentumEstimator::resize(const std::vector<PosType>& kin, const int Min)
{
  //copy kpoints
  kPoints=kin;
  nofK.resize(kin.size());
  kdotp.resize(kPoints.size());
  phases.resize(kPoints.size());
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

