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
        P.rejectMove(0); //restore P.R[0] to the orginal position

        for (int ik=0; ik < kPoints.size(); ++ik)
          {
            for (int i=0; i<np; ++i) kdotp[i]=dot(kPoints[ik],temp[i].dr1);
            eval_e2iphi(np,kdotp.data(),phases.data());
            RealType nofk_here(real(BLAS::dot(np,phases.data(),psi_ratios.data())));
            nofK[ik]+= nofk_here;
          }


        for (int iq=0; iq < Q.size(); ++iq)
          for (int i=0; i<mappedQtonofK[iq].size(); ++i)
            for (int j=0; j<mappednofKtoK[mappedQtonofK[iq][i]].size(); ++j)
            {
              compQ[iq] += nofK[ mappednofKtoK[ mappedQtonofK[iq][i] ] [j] ] * mappedKnorms[mappedQtonofK[iq][i]];
            }


      }
    if (hdf5_out)
      {
//need normalization factor
        int j=myIndex;
        for (int ik=0; ik<nofK.size(); ++ik,++j) P.Collectables[j]+= norm_nofK*nofK[ik];
        for (int iq=0; iq<compQ.size(); ++iq,++j) P.Collectables[j]+= compQ[iq]*mappedQnorms[iq];
      }
    else
      {
        for (int ik=0; ik<nofK.size(); ++ik) nofK[ik] *= norm_nofK;
        for (int iq=0; iq<compQ.size(); ++iq) compQ[iq] *= mappedQnorms[iq];
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

            Matrix<RealType> BasisMatrix(3);

            ParticleSet::ParticlePos_t R_cart(1);
            R_cart.setUnit(PosUnit::CartesianUnit);
            ParticleSet::ParticlePos_t R_unit(1);
            R_unit.setUnit(PosUnit::LatticeUnit);

            for (int i=0;i<3;i++)
              {
                R_unit[0][0]=0;
                R_unit[0][1]=0;
                R_unit[0][2]=0;
                R_unit[0][i]=1;
                elns.convert2Cart(R_unit,R_cart);
                for (int j=0;j<3;j++) BasisMatrix(i,j)= R_cart[0][j];
              }
            invert_matrix(BasisMatrix, false);
//             app_log()<<" Reciprocal Matrix "<<endl;
//             for (int i=0;i<3;i++)
//               {
//                 for (int j=0;j<3;j++)
//                   {
//                     app_log()<<BasisMatrix(i,j)<<" ";
//                   }
//                 app_log()<<endl;
//               }

            app_log()<<" Using all k-space points with (nx^2+ny^2+nz^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<<endl;
            vector<int> primes(8);
            primes[0]=2; primes[1]=3; primes[2]=5; primes[3]=7; primes[4]=11; primes[5]=13; primes[6]=15;
            primes[7]=17;
            assert(kgrid < 19);//if not then add more primes
            for (int i=-kgrid;i<(kgrid+1);i++) for (int j=-kgrid;j<(kgrid+1);j++) for (int k=-kgrid;k<(kgrid+1);k++)
                  {
                    if (std::sqrt(i*i+j*j+k*k)<=kgrid)
                      {
                        PosType kpt;
                        kpt[0]=RealType(i)*BasisMatrix(0,0)+RealType(j)*BasisMatrix(1,0)+RealType(k)*BasisMatrix(2,0);
                        kpt[1]=RealType(i)*BasisMatrix(0,1)+RealType(j)*BasisMatrix(1,1)+RealType(k)*BasisMatrix(2,1);
                        kpt[2]=RealType(i)*BasisMatrix(0,2)+RealType(j)*BasisMatrix(1,2)+RealType(k)*BasisMatrix(2,2);
                        for (int l=0;l<3;l++) kpt[l] *= 2.0*M_PI;
                        kPoints.push_back(kpt);
                        
                        int nzero(0);
                        for(int l=0;l<3;l++) if(std::abs(kpt[l])<1e-3) nzero++;
                        if (nzero==3)
                        {
                          mappedKnorms.push_back(8.0*M_PI*M_PI*M_PI/elns.Lattice.Volume);
                        }
                        else
                        {
                          vector<vector<int> > primefactorization(3);
                          primefactorization[0].push_back(std::abs(i));
                          primefactorization[1].push_back(std::abs(j));
                          primefactorization[2].push_back(std::abs(k));
                          //break into prime numbers, find common ones
                          for(int l=0;l<3;l++)
                          {
                            int pindx(0);
                            while (primefactorization[l][0] > 1)
                            {
                              if(primefactorization[l][0]%primes[pindx]==0)
                              {
                                primefactorization[l][0] /= primes[pindx];
                                primefactorization[l].push_back(primes[pindx]);
                              }
                              else
                                pindx++;
                            }
                          }
                          //now factorized, remove common elements
                          if (nzero==2)
                          {
                            vector<vector<int> > poppedprime(3);
                            for(int l=0;l<3;l++) poppedprime[l].push_back(0);
                            for(int l=0;l<3;l++) if(kpt[l]!=0) poppedprime[l][0]=1;
                            for(int l=0;l<3;l++) primefactorization[l]=poppedprime[l];
                          }
                          else if (nzero==1)
                          {
                            vector<vector<int> > poppedprime;
                            for(int l=0;l<3;l++) if(kpt[l]!=0) poppedprime.push_back(primefactorization[l]);
                            for(int l=0;l<poppedprime[0].size();l++)
                              for(int m=0;m<poppedprime[1].size();m++)
                                if (poppedprime[0][l]==poppedprime[1][m])
                                {
                                  poppedprime[0][l]=1;
                                  poppedprime[1][m]=1;
                                }
                             int ppi(0);
                             for(int l=0;l<3;l++) if(kpt[l]!=0) primefactorization[l]=poppedprime[ppi++];
                          }
                          else
                          {
                            for(int l=0;l<primefactorization[0].size();l++)
                              for(int m=0;m<primefactorization[1].size();m++)
                                for(int n=0;n<primefactorization[2].size();n++)
                                if ((primefactorization[0][l]==primefactorization[1][m])&&(primefactorization[0][l]==primefactorization[2][n]))
                                {
                                  primefactorization[0][l]=1;
                                  primefactorization[1][m]=1;
                                  primefactorization[2][n]=1;
                                }
                          }
                          
                          vector<int> newkpt(3,1);
                          for(int l=0;l<3;l++) for(int m=0;m<primefactorization[l].size();m++) newkpt[l] *= primefactorization[l][m];
                          
                          kpt[0]=RealType(newkpt[0])*BasisMatrix(0,0)+RealType(newkpt[1])*BasisMatrix(1,0)+RealType(newkpt[2])*BasisMatrix(2,0);
                          kpt[1]=RealType(newkpt[0])*BasisMatrix(0,1)+RealType(newkpt[1])*BasisMatrix(1,1)+RealType(newkpt[2])*BasisMatrix(2,1);
                          kpt[2]=RealType(newkpt[0])*BasisMatrix(0,2)+RealType(newkpt[1])*BasisMatrix(1,2)+RealType(newkpt[2])*BasisMatrix(2,2);
                          for (int l=0;l<3;l++) kpt[l] *= 2.0*M_PI;
                          mappedKnorms.push_back(dot(kpt,kpt));
//                           app_log()<<i<<" "<<j<<" "<<k<<" "<<newkpt[0]<<" "<<newkpt[1]<<" "<<newkpt[2]<<"   "<<dot(kpt,kpt)<<endl;
                        }
                      }
                  }
          }
        kids=kids->next;
      }

    vector<std::pair<float,int> > mappedKpoints;
    for (int i=0;i<kPoints.size();i++)
      {
        float khere(std::sqrt(dot(kPoints[i],kPoints[i])));
        std::pair<float,int> kh(khere,i);
        mappedKpoints.push_back(kh);
      }
    std::sort(mappedKpoints.begin(),mappedKpoints.end());

    set<float> qvalues;
    for (int i=0;i<kPoints.size();i++) qvalues.insert(mappedKpoints[i].first);
    mappedQnorms.resize(qvalues.size());

    vector<int> ktoqmap(kPoints.size());
    int qindx(0);
    set<float>::iterator q_it(qvalues.begin());
    for (int i=0;i<kPoints.size();i++)
      {
        if (mappedKpoints[i].first == *q_it)
          {
            ktoqmap[mappedKpoints[i].second]=qindx;
            mappedQnorms[qindx]++;
          }
        else
          {
            q_it++;
            qindx++;
            ktoqmap[mappedKpoints[i].second]=qindx;
            mappedQnorms[qindx]++;
          }
      }
    
    mappedQtonofK.resize(qvalues.size());
    for (int i=0;i<kPoints.size();i++) mappedQtonofK[ktoqmap[i]].push_back(i);
    for (int i=0;i<mappedQnorms.size();i++) mappedQnorms[i]=1.0/(RealType(M)*mappedQnorms[i]);
    



    kids=cur->children;
    while (kids!=NULL)
      {
        string cname((const char*)(kids->name));
        if (cname=="q")
          {
            string ctype("auto");
            OhmmsAttributeSet pAttrib;
            pAttrib.add(ctype,"mode");
            pAttrib.put(kids);
            if (ctype=="auto")
              {
                Q.resize(qvalues.size());
                qindx=0;
                for (q_it=qvalues.begin(); q_it!=qvalues.end(); q_it++) Q[qindx++] = *q_it;
                
                for (int i=0;i<kPoints.size();i++)
                  {
                    vector<int> includedk;
                    RealType knrm=dot(kPoints[i],kPoints[i]);
                    for (int j=0;j<kPoints.size();j++) if (dot(kPoints[i],kPoints[j])==knrm) includedk.push_back(j);
                    mappednofKtoK.push_back(includedk);
                  }
                
//                 mappedQtoK.resize(qvalues.size());
//                 for (int i=0;i<allincludedk.size();i++) mappedQtoK[ktoqmap[i]].insert(mappedQtoK[ktoqmap[i]].end(),allincludedk[i].begin(),allincludedk[i].end());

                if (rootNode)
                  {
                    string QKname="QPoints.dat";
                    ofstream dout(QKname.c_str());
                    dout.setf(ios::scientific, ios::floatfield);
                    dout << "#q       qnorm       qmag      k" << endl;
                    for (int i=0;i<qvalues.size();i++)
                      {
                        dout<<i<<", "<<mappedQnorms[i]<<", "<<Q[i]<<" :  ";
                        for (int j=0;j<mappedQtonofK[i].size();j++) dout<<mappedQtonofK[i][j]<<" ";
                        dout<<endl;
                      }
                    dout.close();
                    
                    QKname="nofKPoints.dat";
                    ofstream nout(QKname.c_str());
                    nout.setf(ios::scientific, ios::floatfield);
                    nout << "#nofk       q       knorm      ks" << endl;
                    for (int i=0;i<kPoints.size();i++)
                      {
                        nout<<i<<", "<<ktoqmap[i]<<" "<<mappedKnorms[i]<<", ";
                        for (int j=0;j<mappednofKtoK[i].size();j++) nout<<mappednofKtoK[i][j]<<" ";
                        nout<<endl;
                      }
                    nout.close();
       
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
                  }

              }
          }
        kids=kids->next;
      }
    
    nofK.resize(kPoints.size());
    compQ.resize(Q.size());
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
    myclone->mappednofKtoK.resize(mappednofKtoK.size());
    for(int i=0;i<mappednofKtoK.size();i++) myclone->mappednofKtoK[i]=mappednofKtoK[i];
    myclone->hdf5_out=hdf5_out;
    myclone->mappedQnorms=mappedQnorms;
    myclone->mappedKnorms=mappedKnorms;

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
