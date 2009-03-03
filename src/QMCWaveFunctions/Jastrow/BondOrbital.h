#ifndef QMCPLUSPLUS_BONDJASTROW_H
#define QMCPLUSPLUS_BONDJASTROW_H
#include "Configuration.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
  {

  /** @ingroup Bond Jastrow
   * @brief Asymmetric Gaussian Bond Jastrow. centered at a bond location with an optimizable offset.
   *  J_b(r = \vec{r-B}) = 1 + c_0*Exp(-r\bar{C}r)
   *  \bar{C} is a matrix with components c_xx, c_xy, ... c_zz. All components must be greater than zero.
   *  cutoff is determined by Max[2*c_ii,c_ij]*r_c - ln(|c_0|) > 14
   */
  class BondOrbital: public OrbitalBase
    {

      const ParticleSet& CenterRef;
      const DistanceTableData* d_table;

      RealType curVal, curLap;
      PosType curGrad;
      ParticleAttrib<RealType> U,d2U;
      ParticleAttrib<PosType> dU;
      RealType *FirstAddressOfdU, *LastAddressOfdU; 
      ///B is the bond center offset from the Bond particleset.
      ///C is the coefficients for the gaussian c_ij and c_0
      vector<vector<RealType> > B,C;
//       Optimize location as well as components?
      bool LocationOptimize;
      int NumBonds,NumVars,perBond,startIndex;
      
///       for derivatives
    vector<GradVectorType*> gradLogPsi;
    vector<ValueVectorType*> lapLogPsi;
    Vector<RealType> dLogPsi,dHLogPsi;

    public: 

      ///constructor
      BondOrbital(ParticleSet& els, const ParticleSet& centers)
          : CenterRef(centers), d_table(0), FirstAddressOfdU(0), LastAddressOfdU(0)
      {
        U.resize(els.getTotalNum());
        d_table = DistanceTable::add(CenterRef,els);
        NumBonds=d_table->size(SourceIndex);
        B.resize(d_table->size(SourceIndex),vector<RealType>(3,0.0));
        C.resize(d_table->size(SourceIndex),vector<RealType>(7,0.0));
        //allocate vector of proper size  and set them to 0 
      }

      ~BondOrbital() { }

      //evaluate the distance table with P
      void resetTargetParticleSet(ParticleSet& P)
      {
        d_table = DistanceTable::add(CenterRef,P);
      }

      /** check in an optimizable parameter
       */
      void checkInVariables(opt_variables_type& active)
      {
      active.insertFrom(myVars);
      }

      /** check out optimizable variables
       */
      void checkOutVariables(const opt_variables_type& active)
      {
      myVars.getIndex(active);  
//       myVars.print(std::cout);
//         reset( );
      }

      ///reset the value of all the unique Two-Body Jastrow functions
      void resetParameters(const opt_variables_type& active)
      {
      int ia=myVars.where(0); 
      if(ia>-1) 
      {
        for (int i=0;i<NumVars;i++) myVars[i]=active[i+ia];
      } 
      reset( );
      }
      
      void reset()
      {
        for (int i=0;i<NumBonds;i++)
        {
          int indx = i*perBond;
          if (myVars[indx]>-1) C[i][0] = myVars[indx];
          else  C[i][0] = myVars[indx]=-1;
          for(int j=1;j<perBond;j++)
          { 
            if (j<7) 
            {
              if (myVars[indx+j]>1e-1) C[i][j] = myVars[indx+j];
              else myVars[indx+j] = C[i][j] = 1e-1;
            }
            else if (j>6) B[i][j-7] = myVars[indx+j];
          }
        }
      }
      
    bool put(xmlNodePtr cur){
      ReportEngine PRE("MolecularBond","put(cur)");
      
      OrbitalName="MolBonds"; 
      string movable="yes";
      OhmmsAttributeSet bb;
      bb.add(OrbitalName,"name");
      bb.add(movable,"movable");
      bb.put(cur);
      LocationOptimize = (movable=="yes");
      xmlNodePtr tcur=cur->children;
      int nBondsReadIn(0);
      while(tcur != NULL)
      {
        string cname((const char*)tcur->name);
        if (cname == "coefficients") 
        {
//           string type("0"), id("0");
//           OhmmsAttributeSet cAttrib;
//           cAttrib.add(id, "id");
//           cAttrib.add(type, "type");
//           cAttrib.put(xmlCoefs);
          int number=-1;
          string type("0");
          string Cname("0");
          OhmmsAttributeSet aa;
          aa.add(Cname,"name");
          aa.add(number,"bond");
          aa.add(type, "type");
          aa.put(tcur); 
          
          if (type != "Array") 
          {
            PRE.error( "Unknown correlation type " + type + " in BsplineFunctor." + "Resetting to \"Array\"");
            xmlNewProp (tcur, (const xmlChar*) "type", (const xmlChar*) "Array");
          }
          
//           if (Cname=="0")
//           {
//             cout<<Cname<<endl;
//             std::stringstream sstr;
//             sstr<<OrbitalName<<"_"<<number;
//             char* Cn(sstr.str());
//             xmlSetProp (tcur, (const xmlChar*) "name", (const xmlChar*)(&Cn));
//             
//           }
          
          
          vector<RealType> inputs;
          if ((number>-1) & (number<NumBonds))
          {
            putContent(inputs,tcur); 
            if (inputs.size()==7) C[number]=inputs;
            else if ((LocationOptimize)&(inputs.size()==10))
            {
              for(int i=0;i<7;i++) C[number][i]=inputs[i];
              for(int j=0;j<3;j++) B[number][j]=inputs[7+j];
            } 
            else if (inputs.size()>6)
            {
              for(int i=0;i<7;i++) C[number][i]=inputs[i];
            }
          }
        }
        tcur = tcur->next;
      }
      myVars.clear();
      std::stringstream sstr;
      for (int i=0;i<NumBonds;i++)
      {
        for(int j=0;j<7;j++)
        {
          sstr.str("");
          sstr << OrbitalName<<"_"<<i<<"_"<<j; 
          if ((C[i][j]<0)&(j!=0)) C[i][j]=0;
          else if ((C[i][j]<-1)&(j==0)) C[i][j]=-1;
          myVars.insert(sstr.str(),C[i][j],true);
        }
        if (LocationOptimize)
        {
          for(int j=0;j<3;j++)
          {
            sstr.str("");
            sstr << OrbitalName<<"_"<<i<<"_r_"<<j; 
            myVars.insert(sstr.str(),B[i][j],true);
          }
        }
      }
      NumVars=myVars.size();    
      if (LocationOptimize) perBond=10;
      else perBond=7;
      reset();
      
      startIndex = myVars.where(0);
      if(NumVars && dLogPsi.size()==0)
      {
        dLogPsi.resize(NumVars);
        gradLogPsi.resize(NumVars,0);
        lapLogPsi.resize(NumVars,0);
        for(int i=0; i<NumVars; ++i)
        {
          gradLogPsi[i]=new GradVectorType(U.size());
          lapLogPsi[i]=new ValueVectorType(U.size());
        }
        dHLogPsi.resize(NumVars); 
      } 
      return true;
    }

    void reportStatus(ostream& os)
    {
          myVars.print(os);
    }

      /**
       *@param P input configuration containing N particles
       *@param G a vector containing N gradients
       *@param L a vector containing N laplacians
       *@return The wavefunction value  \f$exp(-J({\bf R}))\f$
       *
       *Upon exit, the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
       *and the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$ are accumulated.
       *While evaluating the value of the Jastrow for a set of
       *particles add the gradient and laplacian contribution of the
       *Jastrow to G(radient) and L(aplacian) for local energy calculations
       *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f]
       *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
       */
      RealType evaluateLog(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L)
      {
        LogValue=0.0;
        U=0.0; 

        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij,du;
        for (int i=0; i<d_table->size(SourceIndex); i++)
          {
            for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
              {
                j = d_table->J[nn];
                rij = d_table->dr(nn);
                ///Correct for bond offsets from particleset
                rij[0] -= B[i][0];
                rij[1] -= B[i][1];
                rij[2] -= B[i][2];
                ///calculate F(r) value
                f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                             /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                             /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2];
                expfr = C[i][0]*std::exp(-1.0*f);
                /// and gradient terms
                du[0] = -2.0*(C[i][1]*rij[0] + C[i][2]*rij[1] + C[i][3]*rij[2]) ;
                du[1] = -2.0*(C[i][2]*rij[0] + C[i][4]*rij[1] + C[i][5]*rij[2]) ;
                du[2] = -2.0*(C[i][3]*rij[0] + C[i][5]*rij[1] + C[i][6]*rij[2]) ;
                
                orbVal = 1.0 + expfr;
                uij= std::log(orbVal);
                orbVal = 1.0/orbVal;
                expfr *= orbVal;
 
//                 d2f = (-2.0*(C[i][1] + C[i][4] + C[i][6] ) )*expfr ;
                d2f = -2.0*(C[i][1] + C[i][4] + C[i][6] )*expfr + (du[0]*du[0] +du[1]*du[1] +du[2]*du[2] )*expfr*orbVal ;
 
                du[0] *= expfr;
                du[1] *= expfr;
                du[2] *= expfr;
//                 d2f -= (du[0]*du[0] +du[1]*du[1] +du[2]*du[2]);
                
                LogValue += uij;
                U[j] += uij;
                G[j] += du;
                L[j] += d2f;
              }
          }

        return LogValue;
      }

      ValueType evaluate(ParticleSet& P,
                         ParticleSet::ParticleGradient_t& G,
                         ParticleSet::ParticleLaplacian_t& L)
      {
        return std::exp(evaluateLog(P,G,L));
      }

      /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
       * @param P active particle set
       * @param iat particle that has been moved.
       */
      inline ValueType ratio(ParticleSet& P, int iat)
      {
        curVal=0.0;
        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij;
        int n=d_table->size(VisitorIndex);
        for (int i=0, nn=iat; i<d_table->size(SourceIndex); ++i,nn+= n)
        { 
          rij = d_table->Temp[i].dr1;
          ///Correct for bond offsets from particleset
          rij[0] -= B[i][0];
          rij[1] -= B[i][1];
          rij[2] -= B[i][2];
          ///calculate F(r) value
          f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                      /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                      /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2];
          expfr = C[i][0]*std::exp(-1.0*f);
          orbVal = 1.0 + expfr;
          uij= std::log(orbVal);
          curVal += uij; 
        } 
        return std::exp(curVal-U[iat]);
      }

      /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$ and fill-in the differential gradients/laplacians
       * @param P active particle set
       * @param iat particle that has been moved.
       * @param dG partial gradients
       * @param dL partial laplacians
       */
      inline ValueType ratio(ParticleSet& P, int iat,
                             ParticleSet::ParticleGradient_t& dG,
                             ParticleSet::ParticleLaplacian_t& dL)
      {
        return std::exp(logRatio(P,iat,dG,dL));
      }

      inline GradType evalGrad(ParticleSet& P, int iat)
      {
        int n=d_table->size(VisitorIndex);
        curGrad = 0.0;
        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij,du;
        for (int i=0, nn=iat; i<d_table->size(SourceIndex); ++i,nn+= n)
          { 
              rij = d_table->dr(nn);
              ///Correct for bond offsets from particleset
              rij[0] -= B[i][0];
              rij[1] -= B[i][1];
              rij[2] -= B[i][2];
              ///calculate F(r) value
              f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                          /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                          /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2];
              expfr = C[i][0]*std::exp(-1.0*f);
              /// and gradient terms
              du[0] = -2.0*(C[i][1]*rij[0] + C[i][2]*rij[1] + C[i][3]*rij[2]) ;
              du[1] = -2.0*(C[i][2]*rij[0] + C[i][4]*rij[1] + C[i][5]*rij[2]) ;
              du[2] = -2.0*(C[i][3]*rij[0] + C[i][5]*rij[1] + C[i][6]*rij[2]) ;
              
              orbVal = 1.0 + expfr;
//                   uij= std::log(orbVal);
              orbVal = 1.0/orbVal;
              expfr *= orbVal;
              ///grad log Orbital
              du[0] *= expfr;
              du[1] *= expfr;
              du[2] *= expfr;
              curGrad += du; 
            }
        return curGrad;
      }

      inline GradType evalGradSource(ParticleSet& P,
                                     ParticleSet& source, int isrc)
      {
//         if (&source != &CenterRef)
          return GradType();
//         FT* func=Fs[isrc];
//         if (func == 0) return GradType();
//         GradType G(0.0);
//         RealType dudr, d2udr2;
//         for (int nn=d_table->M[isrc]; nn<d_table->M[isrc+1]; nn++)  {
// 	  RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2);
// 	  dudr *= d_table->rinv(nn);
// 	  G += dudr*d_table->dr(nn);
// 	}
//         return G;
      }


      inline GradType 
      evalGradSource(ParticleSet& P, ParticleSet& source, int isrc,
		     TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
		     TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
      {
//         if (&source != &CenterRef)
          return GradType();
//         FT* func=Fs[isrc];
//         if (func == 0) return GradType();
//         GradType G(0.0);
//         RealType dudr, d2udr2, d3udr3;
//         for (int nn=d_table->M[isrc],iel=0; nn<d_table->M[isrc+1]; nn++,iel++) {
// 	  RealType rinv = d_table->rinv(nn); 
// 	  RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2, d3udr3);
// 	  dudr *= rinv;
// 	  d2udr2 *= rinv * rinv;
// 	  G += dudr*d_table->dr(nn);
// 	  for (int dim_ion=0; dim_ion < OHMMS_DIM; dim_ion++) {
// 	    for (int dim_el=0; dim_el < OHMMS_DIM; dim_el++) 
// 	      grad_grad[dim_ion][iel][dim_el] += d2udr2 * d_table->dr(nn)[dim_ion] * d_table->dr(nn)[dim_el]
// 		- dudr * rinv * rinv * d_table->dr(nn)[dim_ion] * d_table->dr(nn)[dim_el];
// 	    grad_grad[dim_ion][iel][dim_ion] += dudr;
// 	    lapl_grad[dim_ion][iel] += (d3udr3*rinv + 2.0*d2udr2 - 2.0*rinv*rinv*dudr)*d_table->dr(nn)[dim_ion];
// 	  }
// 	}
//         return G;
      }

      inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
      {
        int n=d_table->size(VisitorIndex);
        curVal=0.0;
        curGrad = 0.0;
        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij,du;
        for (int i=0, nn=iat; i<d_table->size(SourceIndex); i++,nn+= n)
          { 
//             rij = d_table->dr(nn);
            rij = d_table->Temp[i].dr1;
            ///Correct for bond offsets from particleset
            rij[0] -= B[i][0];
            rij[1] -= B[i][1];
            rij[2] -= B[i][2];
            ///calculate F(r) value
            f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                        /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                        /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2];
            
            /// and gradient terms
            du[0] = -2.0*(C[i][1]*rij[0] + C[i][2]*rij[1] + C[i][3]*rij[2]) ;
            du[1] = -2.0*(C[i][2]*rij[0] + C[i][4]*rij[1] + C[i][5]*rij[2]) ;
            du[2] = -2.0*(C[i][3]*rij[0] + C[i][5]*rij[1] + C[i][6]*rij[2]) ;
            
            expfr = std::exp(-1.0*f)*C[i][0];
            orbVal = 1.0 + expfr;
            uij= std::log(orbVal);
            orbVal = 1.0/orbVal;
            expfr *= orbVal;
            
            du[0] *= expfr;
            du[1] *= expfr;
            du[2] *= expfr;
            curGrad += du;
            curVal += uij; 
          }
        grad_iat += curGrad;
        return std::exp(curVal-U[iat]);
      }

      inline ValueType logRatio(ParticleSet& P, int iat,
                                ParticleSet::ParticleGradient_t& dG,
                                ParticleSet::ParticleLaplacian_t& dL)
      {
        int n=d_table->size(VisitorIndex);
        curVal=0.0;
        curLap=0.0;
        curGrad = 0.0;

        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij,du;
        for (int i=0, nn=iat; i<d_table->size(SourceIndex); i++,nn+= n)
          { 
//             rij = d_table->dr(nn);
            rij = d_table->Temp[i].dr1;
            ///Correct for bond offsets from particleset
            rij[0] -= B[i][0];
            rij[1] -= B[i][1];
            rij[2] -= B[i][2];
            ///calculate F(r) value
            f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                        /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                        /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2]; 
            /// and gradient terms
            du[0] = -2.0*(C[i][1]*rij[0] + C[i][2]*rij[1] + C[i][3]*rij[2]) ;
            du[1] = -2.0*(C[i][2]*rij[0] + C[i][4]*rij[1] + C[i][5]*rij[2]) ;
            du[2] = -2.0*(C[i][3]*rij[0] + C[i][5]*rij[1] + C[i][6]*rij[2]) ;

            expfr = std::exp(-1.0*f)*C[i][0];
            orbVal = 1.0 + expfr;
            uij= std::log(orbVal);
            orbVal = 1.0/orbVal;
            expfr *= orbVal;
            d2f = -2.0*(C[i][1] + C[i][4] + C[i][6] )*expfr + (du[0]*du[0] +du[1]*du[1] +du[2]*du[2] )*expfr*orbVal ;
            du[0] *= expfr;
            du[1] *= expfr;
            du[2] *= expfr; 
            curGrad += du;
            curVal += uij;
            curLap += d2f;
          }
        dG[iat] += curGrad-dU[iat];
        dL[iat] += curLap-d2U[iat];
        return curVal-U[iat];
      }

      inline void restore(int iat) {}

      void acceptMove(ParticleSet& P, int iat)
      {
        U[iat] = curVal;
        dU[iat]=curGrad;
        d2U[iat]=curLap;
      }


      void update(ParticleSet& P,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL,
                  int iat)
      {
        dG[iat] += curGrad-dU[iat];
        dU[iat]  = curGrad;
        dL[iat] += curLap-d2U[iat];
        d2U[iat] = curLap;
        U[iat]   = curVal;
      }

      void evaluateLogAndStore(ParticleSet& P,
                               ParticleSet::ParticleGradient_t& dG,
                               ParticleSet::ParticleLaplacian_t& dL)
      {
        LogValue = 0.0;
        U=0.0;
        dU=0.0;
        d2U=0.0;        
        RealType d2f,f,expfr,orbVal,uij;
        int j;
        PosType rij ,du;
        for (int i=0; i<d_table->size(SourceIndex); i++)
          { 
            for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
              {
                
                j = d_table->J[nn];
                rij = d_table->dr(nn);
//                 rij = d_table->Temp[i].dr1;
                ///Correct for bond offsets from particleset
                rij[0] -= B[i][0];
                rij[1] -= B[i][1];
                rij[2] -= B[i][2];
                ///calculate F(r) value
                f = C[i][1]*rij[0]*rij[0] + 2.0*C[i][2]*rij[0]*rij[1] + 2.0*C[i][3]*rij[0]*rij[2]
                             /*+ C[i][2]*rij[1]*rij[0]*/ + C[i][4]*rij[1]*rij[1] + 2.0*C[i][5]*rij[1]*rij[2]
                             /*+ C[i][3]*rij[2]*rij[0] + C[i][5]*rij[2]*rij[1]*/ + C[i][6]*rij[2]*rij[2];
                expfr = C[i][0]*std::exp(-1.0*f);
                /// and gradient terms
                du[0] = -2.0*(C[i][1]*rij[0] + C[i][2]*rij[1] + C[i][3]*rij[2]) ;
                du[1] = -2.0*(C[i][2]*rij[0] + C[i][4]*rij[1] + C[i][5]*rij[2]) ;
                du[2] = -2.0*(C[i][3]*rij[0] + C[i][5]*rij[1] + C[i][6]*rij[2]) ;
                
                orbVal = 1.0 + expfr;
                uij= std::log(orbVal);
                orbVal = 1.0/orbVal;
                expfr *= orbVal;
 
                d2f = -2.0*(C[i][1] + C[i][4] + C[i][6] )*expfr + (du[0]*du[0] +du[1]*du[1] +du[2]*du[2] )*expfr*orbVal ;
                du[0] *= expfr;
                du[1] *= expfr;
                du[2] *= expfr;
//                 d2f += (du[0]*du[0] +du[1]*du[1] +du[2]*du[2]);
                LogValue += uij;
                U[j] += uij;
                dG[j] += du;
                dL[j] += d2f;
                
                dU[j] += du;
                d2U[j] += d2f;
              }
          }
      }

      /** equivalent to evalaute with additional data management */
      RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
      {

        //U.resize(d_table->size(VisitorIndex));
        d2U.resize(d_table->size(VisitorIndex));
        dU.resize(d_table->size(VisitorIndex));
        FirstAddressOfdU = &(dU[0][0]);
        LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;

        evaluateLogAndStore(P,P.G,P.L);

        //add U, d2U and dU. Keep the order!!!
        buf.add(U.begin(), U.end());
        buf.add(d2U.begin(), d2U.end());
        buf.add(FirstAddressOfdU,LastAddressOfdU);

        return LogValue;
      }

      RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false)
      {
        evaluateLogAndStore(P,P.G,P.L);
        //LogValue = 0.0;
        //U=0.0; dU=0.0; d2U=0.0;
        //RealType uij, dudr, d2udr2;
        //for(int i=0; i<d_table->size(SourceIndex); i++) {
        //  FT* func=Fs[i];
        //  if(func == 0) continue;
        //  for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
        //    int j = d_table->J[nn];
        //    //uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        //    uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
        //    LogValue-=uij;
        //    U[j]+=uij;
        //    dudr *= d_table->rinv(nn);
        //    dU[j] -= dudr*d_table->dr(nn);
        //    d2U[j] -= d2udr2+2.0*dudr;

        //    //add gradient and laplacian contribution
        //    P.G[j] -= dudr*d_table->dr(nn);
        //    P.L[j] -= d2udr2+2.0*dudr;
        //  }
        //}

        //FirstAddressOfdU = &(dU[0][0]);
        //LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;

        buf.put(U.first_address(), U.last_address());
        buf.put(d2U.first_address(), d2U.last_address());
        buf.put(FirstAddressOfdU,LastAddressOfdU);

        return LogValue;
      }

      /** copy the current data from a buffer
       *@param P the ParticleSet to operate on
       *@param buf PooledData which stores the data for each walker
       *
       *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
       */
      void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
      {
        buf.get(U.first_address(), U.last_address());
        buf.get(d2U.first_address(), d2U.last_address());
        buf.get(FirstAddressOfdU,LastAddressOfdU);
      }

      /** return the current value and copy the current data to a buffer
       *@param P the ParticleSet to operate on
       *@param buf PooledData which stores the data for each walker
       */
      inline RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
      {
        RealType sumu = 0.0;
        for (int i=0; i<U.size(); i++) sumu+=U[i];

        buf.put(U.first_address(), U.last_address());
        buf.put(d2U.first_address(), d2U.last_address());
        buf.put(FirstAddressOfdU,LastAddressOfdU);
        return sumu;
        //return std::exp(-sumu);
      }

      OrbitalBasePtr makeClone(ParticleSet& tqp) const
        {
          BondOrbital* j1copy=new BondOrbital(tqp, CenterRef);
          j1copy->myVars.clear();
          j1copy->myVars.insertFrom(myVars);
          j1copy->LocationOptimize=LocationOptimize;
          j1copy->NumVars=NumVars;
          j1copy->startIndex=startIndex;
          j1copy->perBond=perBond;
          j1copy->NumBonds=NumBonds;
          j1copy->dLogPsi.resize(NumVars);
          j1copy->dHLogPsi.resize(NumVars); 
          j1copy->reset();
          return j1copy;
        }

      void copyFrom(const OrbitalBase& old)
      {
        //nothing to do
      }
      
      void evaluateDerivatives(ParticleSet& P, RealType ke0, const opt_variables_type& active, vector<RealType>& dlogpsi, vector<RealType>& dhpsioverpsi)
  { 
    int curIndex = 0;
    dLogPsi=0.0;
    dHLogPsi=0.0;
    PosType TG;
    TG[0]=0; TG[1]=0; TG[2]=0;
    vector<PosType> GLogPsi(NumVars,TG); 
    vector<TinyVector<RealType,3> > derivs(10);
    vector<PosType> gradderivs(10,TG);
//     RealType d2f,f,expfr,cexpfr,orbVal,uij;
//     int j;
    /*PosType rij,du;*/ 
    for(int i=0; i<d_table->size(SourceIndex); i++) 
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) 
      {
        std::fill(derivs.begin(),derivs.end(),0.0); 
        int j = d_table->J[nn];

        PosType rij = d_table->dr(nn); 
//         rij =d_table->Temp[i].dr1;
//         rij = -1.0*CenterRef.R[i]+P.R[j]; 
//         cout<<rij[0]<<"  "<<rij[1]<<"  "<<rij[2]<<endl;
//         cout<<rij[0]-(P.R[j][0]-CenterRef.R[i][0])<<"  "<<rij[1]-(P.R[j][1]-CenterRef.R[i][1])<<"  "<<rij[2]-(P.R[j][2]-CenterRef.R[i][2])<<endl;
        
        ///Correct for bond offsets from particleset
        rij[0] -= B[i][0];
        rij[1] -= B[i][1];
        rij[2] -= B[i][2];
        
        RealType x = rij[0];
        RealType y = rij[1];
        RealType z = rij[2];
        RealType x2 = rij[0]*rij[0];
        RealType y2 = rij[1]*rij[1];
        RealType z2 = rij[2]*rij[2];
        RealType xy = rij[0]*rij[1];
        RealType xz = rij[0]*rij[2];
        RealType yz = rij[1]*rij[2];
        
        RealType FF = C[i][1]*x2 + 2.0*C[i][2]*xy + 2.0*C[i][3]*xz + C[i][4]*y2 + 2.0*C[i][5]*yz + C[i][6]*z2;
//         cout<<C[i][1]*x2<<"  "<<2.0*C[i][2]*xy<<"  "<< 2.0*C[i][3]*xz<<"  "<<C[i][4]*y2 <<"  "<<2.0*C[i][5]*yz<<"  "<< C[i][5]<<"  "<<C[i][6]*z2<<endl;
//         cout<<x<<"  "<<x2<<"  "<<y<<"  "<<y2<<"  "<<z<<"  "<<z2<<"  "<<FF<<endl;

        
        RealType c12  = -2.0*(C[i][1] + C[i][4] + C[i][6]);
        RealType c12x = -2*( C[i][1]*x +C[i][2]*y +C[i][3]*z);
        RealType c12y = -2*( C[i][2]*x +C[i][4]*y +C[i][5]*z);
        RealType c12z = -2*( C[i][3]*x +C[i][5]*y +C[i][6]*z);
        RealType c12x2 = c12x*c12x;
        RealType c12y2 = c12y*c12y;
        RealType c12z2 = c12z*c12z;
        
//         du[0] = c12x ;
//         du[1] = c12y ;
//         du[2] = c12z ;
        
        RealType expfr = std::exp(-1.0*FF);
        RealType cexpfr = expfr*C[i][0]; 
        RealType orbVal = 1.0/(1.0 + cexpfr);
//         uij= std::log(orbVal);
//         orbVal = 1.0/orbVal; 
        RealType cexpOrbVal = cexpfr*orbVal;
        RealType expOrbVal = expfr*orbVal;

        RealType d2f =  c12*cexpOrbVal + (c12x2 + c12y2 + c12z2)*cexpOrbVal*orbVal;

        derivs[0][0]= expfr*orbVal;
        derivs[1][0]= -1.0*x2*cexpOrbVal;
        derivs[2][0]= -2.0*xy*cexpOrbVal;
        derivs[3][0]= -2.0*xz*cexpOrbVal;
        derivs[4][0]= -1.0*y2*cexpOrbVal;
        derivs[5][0]= -2.0*yz*cexpOrbVal;
        derivs[6][0]= -1.0*z2*cexpOrbVal;
        
        derivs[7][0]= -1.0*c12x*cexpOrbVal;
        derivs[8][0]= -1.0*c12y*cexpOrbVal;
        derivs[9][0]= -1.0*c12z*cexpOrbVal;
//         RealType dxoff =  x*d_table->rinv(nn)*derivs[7][0];
//         RealType dyoff =  y*d_table->rinv(nn)*derivs[8][0];
//         RealType dzoff =  z*d_table->rinv(nn)*derivs[9][0];
//         derivs[7][0] -= dyoff+dzoff;
//         derivs[8][0] -= dxoff+dzoff;
//         derivs[9][0] -= dyoff+dxoff;
//         PosType gradexpfrOrbVal = du*expOrbVal;
//         du*=cexpOrbVal;
//         gradderivs[0] = gradexpfrOrbVal - du*derivs[0][0];
//         
//         gradderivs[1] = 1.0*x2*du- du*derivs[1][0];
//         gradderivs[1][0] += -2.0*x*cexpOrbVal;
//         
//         gradderivs[2] = +2.0*xy*du- du*derivs[2][0];
//         gradderivs[2][0] += -2.0*x*cexpOrbVal;
//         gradderivs[2][1] += -2.0*y*cexpOrbVal;
//         
//         gradderivs[3] = +2.0*xz*du- du*derivs[3][0];
//         gradderivs[3][0] += -2.0*x*cexpOrbVal;
//         gradderivs[3][2] += -2.0*z*cexpOrbVal;
//         
//         gradderivs[4] = 1.0*y2*du- du*derivs[4][0];
//         gradderivs[4][1] += -2.0*y*cexpOrbVal;
//         
//         gradderivs[5] = 2.0*yz*du- du*derivs[5][0];
//         gradderivs[5][1] += -2.0*y*cexpOrbVal;
//         gradderivs[5][2] += -2.0*z*cexpOrbVal;
//         
//         gradderivs[6] = 1.0*z2*du- du*derivs[6][0];
//         gradderivs[6][2] += -2.0*z*cexpOrbVal;
//         
//         gradderivs[7] = 1.0*c12x*du - du*derivs[7][0];
//         gradderivs[7][0] += -2*C[i][1]*cexpOrbVal;
//         gradderivs[7][1] += -2*C[i][2]*cexpOrbVal;
//         gradderivs[7][2] += -2*C[i][3]*cexpOrbVal;
//         
//         gradderivs[8] = 1.0*c12y*du - du*derivs[8][0];
//         gradderivs[8][0] += -2*C[i][2]*cexpOrbVal;
//         gradderivs[8][1] += -2*C[i][4]*cexpOrbVal;
//         gradderivs[8][2] += -2*C[i][5]*cexpOrbVal;
//         
//         gradderivs[9] = 1.0*c12z*du - du*derivs[9][0];
//         gradderivs[9][0] += -2*C[i][3]*cexpOrbVal;
//         gradderivs[9][1] += -2*C[i][5]*cexpOrbVal;
//         gradderivs[9][2] += -2*C[i][6]*cexpOrbVal;

        
//         derivs[0][2]= 
//         derivs[1][2]= 
//         derivs[2][2]= 
//         derivs[3][2]= 
//         derivs[4][2]= 
//         derivs[5][2]= 
//         derivs[6][2]= 
//         derivs[7][2]= 
//         derivs[8][2]= 
//         derivs[9][2]=  
        


//         RealType LapA = c12 + (c12x2 + c12y2 + c12z2);
//         derivs[0][2]=  (2*expfr*(-(C[i][0]*(C[i][1] + C[i][4] + C[i][6] + 2*std::pow(C[i][1],2)*std::pow(x,2) + 2*std::pow(C[i][2],2)*std::pow(x,2) + 2*std::pow(C[i][3],2)*std::pow(x,2) + 4*C[i][1]*C[i][2]*x*y + 4*C[i][2]*C[i][4]*x*y + 
//             4*C[i][3]*C[i][5]*x*y + 2*std::pow(C[i][2],2)*std::pow(y,2) + 2*std::pow(C[i][4],2)*std::pow(y,2) + 2*std::pow(C[i][5],2)*std::pow(y,2) + 4*C[i][1]*C[i][3]*x*z + 4*C[i][2]*C[i][5]*x*z + 4*C[i][3]*C[i][6]*x*z + 
//             4*C[i][2]*C[i][3]*y*z + 4*C[i][4]*C[i][5]*y*z + 4*C[i][5]*C[i][6]*y*z + 2*std::pow(C[i][3],2)*std::pow(z,2) + 2*std::pow(C[i][5],2)*std::pow(z,2) + 2*std::pow(C[i][6],2)*std::pow(z,2))) + 
//        expfr*(-C[i][6] + 2*std::pow(C[i][1],2)*std::pow(x,2) + 2*std::pow(C[i][2],2)*std::pow(x,2) + 2*std::pow(C[i][3],2)*std::pow(x,2) + 4*C[i][3]*C[i][5]*x*y + 2*std::pow(C[i][2],2)*std::pow(y,2) + 
//           2*std::pow(C[i][4],2)*std::pow(y,2) + 2*std::pow(C[i][5],2)*std::pow(y,2) + 4*C[i][2]*C[i][5]*x*z + 4*C[i][3]*C[i][6]*x*z + 4*C[i][2]*C[i][3]*y*z + 4*C[i][5]*C[i][6]*y*z + 2*std::pow(C[i][3],2)*std::pow(z,2) + 
//           2*std::pow(C[i][5],2)*std::pow(z,2) + 2*std::pow(C[i][6],2)*std::pow(z,2) + C[i][1]*(-1 + 4*C[i][2]*x*y + 4*C[i][3]*x*z) + C[i][4]*(-1 + 4*C[i][2]*x*y + 4*C[i][5]*y*z))))/std::pow(C[i][0] + expfr,3);
//         RealType Ci0 = C[i][0];
//         RealType Ci1 = C[i][1];
//         RealType Ci2 = C[i][2];
//         RealType Ci3 = C[i][3];
//         RealType Ci4 = C[i][4];
//         RealType Ci5 = C[i][5];
//         RealType Ci6 = C[i][6];
//         RealType Exppart = cexpfr;
//         derivs[1][2]= Exppart*(-2 + 2*(Ci1 + Ci4 + Ci6)*std::pow(x,2)) - (4*Exppart*x*(-2*Ci1*x - 2*Ci2*y - 2*Ci3*z))/(1 + Exppart) + 
//    (C[i][0]*Exppart*std::pow(x,2)*(std::pow(-2*Ci1*x - 2*Ci2*y - 2*Ci3*z,2) + std::pow(-2*Ci2*x - 2*Ci4*y - 2*Ci5*z,2) + std::pow(-2*Ci3*x - 2*Ci5*y - 2*Ci6*z,2)))/
//     std::pow(1 + Exppart,2) - (Exppart*std::pow(x,2)*(std::pow(-2*Ci1*x - 2*Ci2*y - 2*Ci3*z,2) + std::pow(-2*Ci2*x - 2*Ci4*y - 2*Ci5*z,2) + 
//         std::pow(-2*Ci3*x - 2*Ci5*y - 2*Ci6*z,2)))/(1 + Exppart);
//         derivs[3][2]= 
//         derivs[4][2]= 
//         derivs[5][2]= 
//         derivs[6][2]= 
//         derivs[7][2]= 
//         derivs[8][2]= 
//         derivs[9][2]=  
        RealType LapA = c12 + (c12x2 + c12y2 + c12z2);
        derivs[0][2]=  LapA*expOrbVal;
        derivs[1][2]= ((-2.0 - 4.0*x*c12x) - x2*LapA)*cexpOrbVal;
        derivs[2][2]= (-4.0*( x*c12y+y*c12x) - 2.0*xy*LapA)*cexpOrbVal;
        derivs[3][2]= (-4.0*( z*c12x+x*c12z) - 2.0*xz*LapA)*cexpOrbVal;
        derivs[4][2]= ((-2.0 - 4.0*y*c12y) - y2*LapA)*cexpOrbVal;
        derivs[5][2]= (-4.0*( z*c12y+y*c12z) - 2.0*yz*LapA)*cexpOrbVal;
        derivs[6][2]= ((-2.0 - 4.0*z*c12z) - z2*LapA)*cexpOrbVal;
        derivs[7][2]= -cexpOrbVal*(-4.0*(C[i][1]*c12x +C[i][2]*c12y +C[i][3]*c12z) + LapA*c12x);
        derivs[8][2]= -cexpOrbVal*(-4.0*(C[i][2]*c12x +C[i][4]*c12y +C[i][5]*c12z) + LapA*c12y);
        derivs[9][2]= -cexpOrbVal*(-4.0*(C[i][3]*c12x +C[i][5]*c12y +C[i][6]*c12z) + LapA*c12z);
               
        for(int ip=0; ip<perBond; ip++ )
        {
          dLogPsi[curIndex+ip] += derivs[ip][0] ;
//           GLogPsi[curIndex+ip] += gradderivs[ip];
          dHLogPsi[curIndex+ip] += derivs[ip][2] - (P.L[j]+dot(P.G[j],P.G[j]))*derivs[ip][0]; 
        }
//         cout<<"New e"<<endl;
//         cout<<j<<"  "<<x<<"  "<<y<<"  "<<z<<"  "<<FF<<"  "<<orbVal<<endl;
//         for(int ip=0; ip<perBond; ip++ ) cout<<ip<<"  "<<derivs[ip][0]<<"  "<<derivs[ip][2]<<endl;
      }
      curIndex+=perBond;
    }

/*      for(int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if(kk<0) continue;
        dlogpsi[kk]=dLogPsi[k]; 
        dhpsioverpsi[kk]= -0.5* dHLogPsi[k] - ke0* dLogPsi[k];
      }      */  
      int k0=myVars.where(0);
      for(int ip=0; ip<NumVars; ip++ )
      {
        dlogpsi[k0+ip]  = dLogPsi[ip]; 
        dhpsioverpsi[k0+ip]  = -0.5* dHLogPsi[ip] ;
//         dhpsioverpsi[k0+ip]  = -0.5* dHLogPsi[ip] - ke0*dLogPsi[ip];
//         dhpsioverpsi[k0+ip]  = -0.5* (dHLogPsi[ip] - Dot(GLogPsi[ip],P.G));
      }
//       for(int ip=0; ip<NumVars; ip++ ) cout<<ip<<"  "<<dlogpsi[k0+ip]<<"  "<<dhpsioverpsi[k0+ip]<<endl;
  }
      
      
    };
 

}
#endif


