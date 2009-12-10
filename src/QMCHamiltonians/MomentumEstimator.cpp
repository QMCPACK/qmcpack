//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
      :M(4), refPsi(psi), kgrid(4), Lattice(elns.Lattice), norm_nofK(1), hdf5_out(false)
  {
    UpdateMode.set(COLLECTABLE,1);
    psi_ratios.resize(elns.getTotalNum());
    kdotp.resize(elns.getTotalNum());
    phases.resize(elns.getTotalNum());
    twist=elns.getTwist();
  }

  void MomentumEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  MomentumEstimator::Return_t MomentumEstimator::evaluate(ParticleSet& P)
  {
    const int np=P.getTotalNum();
    nofK=0.0;
    compQ=0.0;

    //will use temp[i].r1 for the Compton profile
    const vector<DistanceTableData::TempDistType>& temp(P.DistTables[0]->Temp);
    for (int s=0; s<M; ++s)
      {
        PosType newpos;

        for (int i=0; i<OHMMS_DIM;++i) newpos[i]=myRNG();
        //make it cartesian
        newpos=Lattice.toCart(newpos);
        P.makeVirtualMoves(newpos); //updated: temp[i].r1=|newpos-P.R[i]|, temp[i].dr1=newpos-P.R[i]
        refPsi.get_ratios(P,psi_ratios);
//         for (int i=0; i<np; ++i) app_log()<<i<<" "<<psi_ratios[i].real()<<" "<<psi_ratios[i].imag()<<endl;
        P.rejectMove(0); //restore P.R[0] to the orginal position

        for (int ik=0; ik < kPoints.size(); ++ik)
          {
            for (int i=0; i<np; ++i) kdotp[i]=dot(kPoints[ik],temp[i].dr1_nobox);
            eval_e2iphi(np,kdotp.data(),phases.data());
            RealType nofk_here(std::real(BLAS::dot(np,phases.data(),psi_ratios.data())));
            nofK[ik]+= nofk_here;
          }
          
        

        for (int iq=0; iq < Q.size(); ++iq)
          for (int i=0; i<mappedQtonofK[iq].size(); ++i)
              compQ[iq] += nofK[mappedQtonofK[iq][i]];
            
       


      }
      for (int ik=0; ik<nofK.size(); ++ik) nofK[ik] *= norm_nofK;
      for (int iq=0; iq<compQ.size(); ++iq) compQ[iq] *= norm_nofK*mappedQnorms[iq];
      
    if (hdf5_out)
      {
        int j=myIndex;
        for (int ik=0; ik<nofK.size(); ++ik,++j) P.Collectables[j]+= nofK[ik];
        for (int iq=0; iq<compQ.size(); ++iq,++j) P.Collectables[j]+= compQ[iq];
      }
    return 0.0;
  }

  void MomentumEstimator::registerCollectables(vector<observable_helper*>& h5desc
      , hid_t gid) const
    {
      if (hdf5_out)
        {
          //descriptor for the data, 1-D data
          vector<int> ng(1);

          //add nofk
          ng[0]=nofK.size();
          observable_helper* h5o=new observable_helper("nofk");
          h5o->set_dimensions(ng,myIndex);
          h5o->open(gid);
          h5o->addProperty(const_cast<vector<PosType>&>(kPoints),"kpoints");
          h5o->addProperty(const_cast<vector<int>&>(kWeights),"kweights");
          h5desc.push_back(h5o);

          //add compQ
          ng[0]=Q.size();
          h5o=new observable_helper("compQ");
          h5o->set_dimensions(ng,myIndex+nofK.size());
          h5o->open(gid);
          h5o->addProperty(const_cast<vector<RealType>&>(Q),"q");
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
        for (int i=0;i<nofK.size();i++)
          {
            std::stringstream sstr;
            sstr << "nofk_" <<i;
            int id=plist.add(sstr.str());
          }
        for (int i=0;i<Q.size();i++)
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
        std::copy(nofK.begin(),nofK.end(),plist.begin()+myIndex);
        std::copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size());
      }

  }

  void MomentumEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    if (!hdf5_out)
      {
        std::copy(nofK.begin(),nofK.end(),plist.begin()+myIndex+offset);
        std::copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size()+offset);
      }
  }

bool MomentumEstimator::putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode)
  {
    string ctype("scalar");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(ctype,"mode");
    pAttrib.add(M,"samples");
    pAttrib.put(cur);
    if (ctype=="hdf5") hdf5_out=true;
    else hdf5_out=false;
    
    mappedQtonofK.resize(kgrid+1);
    compQ.resize(kgrid+1);
    mappedQnorms.resize(kgrid+1,1.0/6.0);
    mappedQnorms[0]=1.0/3.0;
    Q.resize(kgrid+1);
    
    xmlNodePtr kids=cur->children;
    while (kids!=NULL)
      {
        string cname((const char*)(kids->name));
        if (cname=="kpoints")
          {
            string ctype("manual");
            OhmmsAttributeSet pAttrib;
            pAttrib.add(ctype,"mode");
            pAttrib.add(kgrid,"grid");
            pAttrib.put(kids);

            app_log()<<" Using all k-space points with (nx^2+ny^2+nz^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<<endl;
            app_log()<<"  My twist is:"<<twist[0]<<"  "<<twist[1]<<"  "<<twist[2]<<endl;

            RealType cof = -1.0;
            int indx(0);
            for (int i=-kgrid;i<(kgrid+1);i++)
            {
              PosType kpt;
              for (int j=-kgrid;j<(kgrid+1);j++) for (int k=-kgrid;k<(kgrid+1);k++)
                  {
                    if (std::sqrt(i*i+j*j+k*k)<=kgrid)
                      {
                        kpt[0]=i+cof*twist[0];
                        kpt[1]=j+cof*twist[1];
                        kpt[2]=k+cof*twist[2];

                        //convert to Cartesian: note that 2Pi is multiplied
                        kpt=Lattice.k_cart(kpt);
                        kPoints.push_back(kpt);
                        
                        mappedQtonofK[abs(i)].push_back(indx);
                        mappedQtonofK[abs(j)].push_back(indx);
                        mappedQtonofK[abs(k)].push_back(indx);
                        indx++;
                        }
                   }
                Q[abs(i)]=kpt[0];
              }
           }
          
        kids=kids->next;
      }
if (rootNode)
{
    RealType KF(std::pow(3*M_PI*M_PI*elns.R.size()/elns.Lattice.Volume,1.0/3.0));
    string fname="Kpoints.dat";
    ofstream fout(fname.c_str());
    fout.setf(ios::scientific, ios::floatfield);
    fout << "# mag_k        kx           ky            kz            koverk_f" << endl;
    for (int i=0;i<kPoints.size();i++)
      {
        float khere(std::sqrt(dot(kPoints[i],kPoints[i])));
        fout<<khere<<"   "<<kPoints[i][0]<<"    "<<kPoints[i][1]<<" "<<kPoints[i][2] 
        << " "<< khere/KF<<endl;
      }
    fout.close();

    fname="Qpoints.dat";
    ofstream qout(fname.c_str());
    qout.setf(ios::scientific, ios::floatfield);
    qout << "# mag_q" << endl;
    for (int i=0;i<Q.size();i++)
      {
        qout<<Q[i]<<endl;
      }
    qout.close();
}
    
    nofK.resize(kPoints.size());
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
    myclone->resize(kPoints,Q);
    myclone->myIndex=myIndex;
    myclone->kgrid=kgrid;
    myclone->norm_nofK=norm_nofK;
    
    myclone->mappedQtonofK.resize(mappedQtonofK.size());
    for(int i=0;i<mappedQtonofK.size();i++) myclone->mappedQtonofK[i]=mappedQtonofK[i];
    myclone->hdf5_out=hdf5_out;
    myclone->mappedQnorms=mappedQnorms;
    return myclone;
  }

  void MomentumEstimator::resize(const vector<PosType>& kin, const vector<RealType>& qin)
  {
    //copy kpoints
    kPoints=kin;
    nofK.resize(kin.size());

    //copy q
    Q=qin;
    compQ.resize(qin.size());

  }

  void MomentumEstimator::setRandomGenerator(RandomGenerator_t* rng)
  {
    //simply copy it
    myRNG=*rng;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
