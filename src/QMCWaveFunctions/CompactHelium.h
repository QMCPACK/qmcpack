#ifndef COMPACT_HELIUM_H
#define COMPACT_HELIUM_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {


///Class for the Helium wavefunction using two electrons and one Ion. Bad version
class CompactHeliumTwo: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    CompactHeliumTwo(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      CompactHeliumTwo* cloned = new CompactHeliumTwo(tqp,CenterRef);
//       cloned->nameA = nameA;
//       cloned->nameB = nameB;
//       cloned->nameC = nameC;
//       cloned->nameD = nameD;
      cloned->myVars.insert(nameA,pA,true);
      cloned->myVars.insert(nameB,pB,true);
      cloned->myVars.insert(nameC,pC,true);
      cloned->myVars.insert(nameD,pD,true);
//       cloned->OrbitalName=OrbitalName;
      cloned->pA=pA;
      cloned->pB=pB;
      cloned->pC=pC;
      cloned->pD=pD;
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-9.85875926e-02;
      pB=6.42881288e-02 ;
      pC=7.47029314e-02 ;
      pD=-8.07044762e-02;
      
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_B";
      nameB = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_C";
      nameC = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_D";
      nameD = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	if(pname[0]=='B') putContent(pB,cur);
	if(pname[0]=='C') putContent(pC,cur);
	if(pname[0]=='D') putContent(pD,cur);
	cur=cur->next;
      }
      
      
      
      myVars.insert(nameA,pA,true);
      myVars.insert(nameB,pB,true);
      myVars.insert(nameC,pC,true);
      myVars.insert(nameD,pD,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      int ib=myVars.where(1); if(ib>-1) pB=active[ib];
      int ic=myVars.where(2); if(ic>-1) pC=active[ic];
      int id=myVars.where(3); if(id>-1) pD=active[id];
      reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D){
      pA=A;
      pB=B;
      pC=C;
      pD=D;
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<"WF parameters: A="<<pA<<"  B="<<pB<<"  C="<<pC<<"  D="<<pD<<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
      ValueType r1 = Ie_table->r(1);
      ValueType r01 = ee_table->r(0);
      
      ValueType s = r0+r1;
      ValueType t = r0-r1;
      ValueType u = r01;
      
      ValueType expm2s = std::exp(-2*s);
      ValueType expAu = std::exp(pA*u);
      ValueType part1 = (1.0+0.5*u)*expAu;
      ValueType part2 = 1.0 + pB*s*u + pC*t*t + pD*u*u;
      ValueType mpart1 = 1.0/part1;
      ValueType mpart2 = 1.0/part2;
      
//       Gradients
      ValueType bu2ct  = (pB*u+2.0*pC*t)*mpart2;
      ValueType bs2du  = (pB*s+2.0*pD*u)*mpart2;
      ValueType bu2mct = (pB*u-2.0*pC*t)*mpart2;
      ValueType mupart = 1.0/(1.0+0.5*u);

      ValueType F01 = (-2.0 + bu2ct );
      ValueType F02 = (pA + 0.5*mupart + bs2du );
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);
      
      ValueType F11 = (-2.0 + bu2mct );
      PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] = J0;
      G[1] = J1;  
      
      ValueType L0 = -0.25*mupart*mupart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2ct*bu2ct;
      L0 -= dot(Ie_table->dr(0),ee_table->dr(0))*Ie_table->rinv(0)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2ct);
      L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;

      ValueType L1= -0.25*mupart*mupart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2mct*bu2mct;
      L1 += dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2mct);
      L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02; 
      
      L[0]= L0;
      L[1]= L1; 
      
       return expm2s*part1*part2;
    }


};

///Class for the Helium wavefunction using two electrons and one Ion. Re-corrected form. This one works the best.
class NewCompactHeliumTwo: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    NewCompactHeliumTwo(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      NewCompactHeliumTwo* cloned = new NewCompactHeliumTwo(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->nameB = nameB;
      cloned->nameC = nameC;
      cloned->nameD = nameD;
      cloned->myVars.insert(nameA,pA,true);
      cloned->myVars.insert(nameB,pB,true);
      cloned->myVars.insert(nameC,pC,true);
      cloned->myVars.insert(nameD,pD,true);
//       cloned->OrbitalName=OrbitalName;
      cloned->reset(pA,pB,pC,pD);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-1.013;
      pB=0.2119256858;
      pC=0.1416353426;
      pD=-0.0113074102;
      
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_B";
      nameB = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_C";
      nameC = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_D";
      nameD = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	if(pname[0]=='B') putContent(pB,cur);
	if(pname[0]=='C') putContent(pC,cur);
	if(pname[0]=='D') putContent(pD,cur);
	cur=cur->next;
      }
      reset(pA,pB,pC,pD);
      
      
      myVars.insert(nameA,pA,true);
      myVars.insert(nameB,pB,true);
      myVars.insert(nameC,pC,true);
      myVars.insert(nameD,pD,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      int ib=myVars.where(1); if(ib>-1) pB=active[ib];
      int ic=myVars.where(2); if(ic>-1) pC=active[ic];
      int id=myVars.where(3); if(id>-1) pD=active[id];
      reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D){
      pA=A;
      pB=B;
      pC=C;
      pD=D;
      
//       myVars.clear();
//       myVars.insert(nameA,pA,true);
//       myVars.insert(nameB,pB,true);
//       myVars.insert(nameC,pC,true);
//       myVars.insert(nameD,pD,true);
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<"WF parameters: A="<<pA<<"  B="<<pB<<"  C="<<pC<<"  D="<<pD<<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
      ValueType r1 = Ie_table->r(1);
      ValueType r01 = ee_table->r(0);
      
      ValueType s = r0+r1;
      ValueType t = r0-r1;
      ValueType u = r01;
      
      ValueType expm2s = std::exp(-2*s);
      ValueType expAu = std::exp(pA*u);
      ValueType part1 = 1.0+0.5*u*expAu;
      ValueType part2 = 1.0 + pB*s*u + pC*t*t + pD*u*u;
      ValueType mpart1 = 1.0/part1;
      ValueType mpart2 = 1.0/part2;
      
//       Gradients
      ValueType bu2ct  = (pB*u+2.0*pC*t)*mpart2;
      ValueType bs2du  = (pB*s+2.0*pD*u)*mpart2;
      ValueType bu2mct = (pB*u-2.0*pC*t)*mpart2;
      ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (-2.0 + bu2ct );
      ValueType F02 = (upart + bs2du );
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);
      
      ValueType F11 = (-2.0 + bu2mct );
//       app_log()<<F11<<endl;
      PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] = J0;
      G[1] = J1;  
      
      ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2ct*bu2ct;
      L0 -= dot(Ie_table->dr(0),ee_table->dr(0))*Ie_table->rinv(0)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2ct);
      L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;

      ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2mct*bu2mct;
      L1 += dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2mct);
      L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02; 
      
      L[0]= L0;
      L[1]= L1; 
      
       return expm2s*part1*part2;
    }


};

class NewCuspCompactHeliumTwo: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    NewCuspCompactHeliumTwo(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      NewCuspCompactHeliumTwo* cloned = new NewCuspCompactHeliumTwo(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->nameB = nameB;
      cloned->nameC = nameC;
      cloned->nameD = nameD;
      cloned->myVars.insert(nameA,pA,true);
      cloned->myVars.insert(nameB,pB,true);
//       cloned->myVars.insert(nameC,pC,true);
      cloned->myVars.insert(nameD,pD,true);
//       cloned->OrbitalName=OrbitalName;
      cloned->reset(pA,pB,pC,pD);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-1.013;
      pB=0.2119256858;
      pC=pB*0.5;
      pD=-0.0113074102;
      
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_B";
      nameB = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_C";
      nameC = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_D";
      nameD = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	if(pname[0]=='B') putContent(pB,cur);
	if(pname[0]=='C') putContent(pC,cur);
	if(pname[0]=='D') putContent(pD,cur);
	cur=cur->next;
      }
      reset(pA,pB,pC,pD);
      
      
      myVars.insert(nameA,pA,true);
      myVars.insert(nameB,pB,true);
//       myVars.insert(nameC,pC,true);
      myVars.insert(nameD,pD,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      int ib=myVars.where(1); if(ib>-1) pB=active[ib];
//       int ic=myVars.where(2); if(ic>-1) pC=active[ic];
      int id=myVars.where(2); if(id>-1) pD=active[id];
      reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D){
      pA=A;
      pB=B;
      pC=0.5*B;
      pD=D;
      
//       myVars.clear();
//       myVars.insert(nameA,pA,true);
//       myVars.insert(nameB,pB,true);
//       myVars.insert(nameC,pC,true);
//       myVars.insert(nameD,pD,true);
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA<<endl;
      os<<nameB<<"  "<<pB<<endl;
      os<<nameC<<"  "<<pC<<endl;
      os<<nameD<<"  "<<pD<<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
      ValueType r1 = Ie_table->r(1);
      ValueType r01 = ee_table->r(0);
      
      ValueType s = r0+r1;
      ValueType t = r0-r1;
      ValueType u = r01;
      
      ValueType expm2s = std::exp(-2*s);
      ValueType expAu = std::exp(pA*u);
      ValueType part1 = 1.0+0.5*u*expAu;
      ValueType part2 = 1.0 + pB*s*u + pC*t*t + pD*u*u;
      ValueType mpart1 = 1.0/part1;
      ValueType mpart2 = 1.0/part2;
      
//       Gradients
      ValueType bu2ct  = (pB*u+2.0*pC*t)*mpart2;
      ValueType bs2du  = (pB*s+2.0*pD*u)*mpart2;
      ValueType bu2mct = (pB*u-2.0*pC*t)*mpart2;
      ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (-2.0 + bu2ct );
      ValueType F02 = (upart + bs2du );
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);
      
      ValueType F11 = (-2.0 + bu2mct );
//       app_log()<<F11<<endl;
      PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] = J0;
      G[1] = J1;  
      
      ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2ct*bu2ct;
      L0 -= dot(Ie_table->dr(0),ee_table->dr(0))*Ie_table->rinv(0)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2ct);
      L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;

      ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2mct*bu2mct;
      L1 += dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2mct);
      L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02; 
      
      L[0]= L0;
      L[1]= L1; 
      
       return expm2s*part1*part2;
    }


};

class SimpleCompactHelium: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    SimpleCompactHelium(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      SimpleCompactHelium* cloned = new SimpleCompactHelium(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->myVars.insert(nameA,pA,true);
      cloned->reset(pA);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-0.08;
     
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	cur=cur->next;
      }
      reset(pA);
      
      
      myVars.insert(nameA,pA,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      reset(pA);
    }
    void reset(ValueType A){
      pA=A;
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA <<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
      ValueType r1 = Ie_table->r(1);
      ValueType r01 = ee_table->r(0);
      
      ValueType s = r0+r1;
      ValueType t = r0-r1;
      ValueType u = r01;
      
      ValueType expm2s = std::exp(-2*s);
      ValueType expAu = std::exp(pA*u);
      ValueType part1 = 1.0+0.5*u*expAu;
      ValueType mpart1 = 1.0/part1;

      
//       Gradients
      ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (-2.0);
      ValueType F02 = (upart);
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);
      
      ValueType F11 = (-2.0 );
//       app_log()<<F11<<endl;
      PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] += J0;
      G[1] += J1;  
      
      ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart ;
      L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;

      ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart;
      L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02; 
      
      L[0] += L0;
      L[1] += L1; 
      
       return expm2s*part1;
    }


};

class SimpleCompactHeliumElectronCorrelation: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    SimpleCompactHeliumElectronCorrelation(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      SimpleCompactHeliumElectronCorrelation* cloned = new SimpleCompactHeliumElectronCorrelation(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->myVars.insert(nameA,pA,true);
      cloned->reset(pA);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-0.08;
     
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	cur=cur->next;
      }
      reset(pA);
      
      
      myVars.insert(nameA,pA,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      reset(pA);
    }
    void reset(ValueType A){
      pA=A;
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA <<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
//       ValueType r0 = Ie_table->r(0);
//       ValueType r1 = Ie_table->r(1);
      ValueType r01 = ee_table->r(0);
      
//       ValueType s = r0+r1;
//       ValueType t = r0-r1;
      ValueType u = r01;
      
//       ValueType expm2s = std::exp(-2*s);
      ValueType expAu = std::exp(pA*u);
      ValueType part1 = 1.0+0.5*u*expAu;
      ValueType mpart1 = 1.0/part1;
// 
//       
// //       Gradients
      ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;
// 
//       ValueType F01 = (-2.0);
      ValueType F02 = (upart);
      PosType J0 = -1.0* ee_table->dr(0)*F02*ee_table->rinv(0);
//       
//       ValueType F11 = (-2.0 );
// //       app_log()<<F11<<endl;
      PosType J1 = ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] += J0;
      G[1] += J1;  
//       
      ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart ;
      L0 += 2.0*ee_table->rinv(0)*F02;
// 
      ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart;
      L1 += 2.0*ee_table->rinv(0)*F02; 
//       
      L[0] += L0;
      L[1] += L1; 
      
       return part1;
    }


};


class SimpleCompactHeliumOrbitalPart: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    SimpleCompactHeliumOrbitalPart(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
      ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      SimpleCompactHeliumOrbitalPart* cloned = new SimpleCompactHeliumOrbitalPart(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->myVars.insert(nameA,pA,true);
      cloned->nameB = nameB;
      cloned->myVars.insert(nameB,pB,true);
      cloned->reset(pA,pB);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-0.08; pB=-2.0;
     
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_B";
      nameB = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	if(pname[0]=='B') putContent(pB,cur);
	cur=cur->next;
      }
      reset(pA,pB);
      
      
      myVars.insert(nameA,pA,true);
      myVars.insert(nameB,pB,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
      ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      int ib=myVars.where(1); if(ib>-1) pB=active[ia];
      reset(pA,pB);
    }
    void reset(ValueType A,ValueType B){
      pA=A;
      pB=B;
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA<<endl;
      os<<nameB<<"  "<<pB<<endl;

    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
      ValueType r1 = Ie_table->r(1);
//       ValueType r01 = ee_table->r(0);
      
      ValueType s = r0+r1;
//       ValueType t = r0-r1;
//       ValueType u = r01;
      
      ValueType expm2s = std::exp(pB*s);
//       ValueType expAu = std::exp(pA*u);
//       ValueType part1 = 1.0+0.5*u*expAu;
//       ValueType mpart1 = 1.0/part1;

      
//       Gradients
//       ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (pB);
//       ValueType F02 = (upart);
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0);
//       - ee_table->dr(0)*F02*ee_table->rinv(0);
      
      ValueType F11 = (pB );
//       app_log()<<F11<<endl;
      PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1);
//       + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] += J0;
      G[1] += J1;  
      
//       ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart ;
      ValueType L0  = 2.0*Ie_table->rinv(0)*F01;
//       + 2.0*ee_table->rinv(0)*F02;

//       ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart;
      ValueType L1  = 2.0*Ie_table->rinv(1)*F11;
//       + 2.0*ee_table->rinv(0)*F02; 
      
      L[0] += L0;
      L[1] += L1; 
      
       return expm2s;
//        *part1;
    }


};

class SingleSlaterOrbital: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;

    SingleSlaterOrbital(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
//       ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      SingleSlaterOrbital* cloned = new SingleSlaterOrbital(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->myVars.insert(nameA,pA,true);
      cloned->nameB = nameB;
      cloned->myVars.insert(nameB,pB,true);
      cloned->reset(pA,pB);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-0.08; pB=-1.0;
     
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
      sstr.str("");
      sstr << HEprefix << "_B";
      nameB = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
	if(pname[0]=='B') putContent(pB,cur);
	cur=cur->next;
      }
      reset(pA,pB);
      
      
      myVars.insert(nameA,pA,true);
      myVars.insert(nameB,pB,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
//       ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
      int ib=myVars.where(1); if(ib>-1) pB=active[ia];
      reset(pA,pB);
    }
    void reset(ValueType A,ValueType B){
      pA=A;
      pB=B;
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA<<endl;
      os<<nameB<<"  "<<pB<<endl;

    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
//       ValueType r1 = Ie_table->r(1);
//       ValueType r01 = ee_table->r(0);
      
      ValueType s = r0 ;
//       ValueType t = r0-r1;
//       ValueType u = r01;
      
      ValueType expm2s = std::exp(pB*s);
//       ValueType expAu = std::exp(pA*u);
//       ValueType part1 = 1.0+0.5*u*expAu;
//       ValueType mpart1 = 1.0/part1;

      
//       Gradients
//       ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (pB);
//       ValueType F02 = (upart);
      PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0);
//       - ee_table->dr(0)*F02*ee_table->rinv(0);
      
//       ValueType F11 = (pB );
//       app_log()<<F11<<endl;
//       PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1);
//       + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] += J0;
//       G[1] += J1;  
      
//       ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart ;
      ValueType L0  = 2.0*Ie_table->rinv(0)*F01;
//       + 2.0*ee_table->rinv(0)*F02;

//       ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart;
//       ValueType L1  = 2.0*Ie_table->rinv(1)*F11;
//       + 2.0*ee_table->rinv(0)*F02; 
      
      L[0] += L0;
//       L[1] += L1; 
      
       return expm2s;
//        *part1;
    }


};
// class JeremyCompactHelium: public OrbitalBase
// {
//   public:
//     ParticleSet& CenterRef;
//     DistanceTableData* Ie_table;
//     DistanceTableData* ee_table;
//     ValueType pA, pB, pC, pD, pE;
//     string nameA,nameB,nameC,nameD,nameE;
// 
//     JeremyCompactHelium(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
//     {
//       Ie_table = DistanceTable::add(CenterRef,electrons);
//       ee_table = DistanceTable::add(electrons);
//       
//       ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
// //       pA=-1.013;
// //       pB=0.2119;
// //       pC=0.1406;
// //       pD=-0.003;
//     }
//     
//     OrbitalBase* makeClone(ParticleSet& tqp) const
//     {
//       JeremyCompactHelium* cloned = new JeremyCompactHelium(tqp,CenterRef);
//       cloned->nameA = nameA;
//       cloned->nameB = nameB;
//       cloned->nameC = nameC;
//       cloned->nameD = nameD;
//       cloned->nameE = nameE;
//       cloned->myVars.insert(nameA,pA,true);
// //       cloned->myVars.insert(nameC,pC,true);
//       cloned->myVars.insert(nameD,pD,true);
//       cloned->myVars.insert(nameE,pE,true);
// //       cloned->OrbitalName=OrbitalName;
//       cloned->reset(pA,pB,pC,pD,pE);
//       
//       return cloned;
//     }
//     
//     
//     bool put(xmlNodePtr cur){
//       pA=-0.08;
//       pB=0.0;
//       pC=0.0;
//       pD=0.0;
//       pE=0.0;
//       
//       
//       OrbitalName="CHe4";
//       string HEprefix("HE");
//       OhmmsAttributeSet bb;
//       bb.add(HEprefix,"id");
//       bb.add(OrbitalName,"name");
//       bb.put(cur);
//       
//       std::stringstream sstr;
//       sstr << HEprefix << "_A";
//       nameA = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_B";
//       nameB = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_C";
//       nameC = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_D";
//       nameD = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_E";
//       nameE = sstr.str();
//       
//       cur=cur->children;
//       while(cur != NULL)
//       {
// 	string pname="0";
// 	OhmmsAttributeSet aa;
// 	aa.add(pname,"name");
// 	aa.put(cur);
// 	if(pname[0]=='A') putContent(pA,cur);
// // 	if(pname[0]=='C') putContent(pC,cur);
// 	if(pname[0]=='D') putContent(pD,cur);
// 	if(pname[0]=='E') putContent(pE,cur);
// 	cur=cur->next;
//       }
//       reset(pA,pB,pC,pD,pE);
//       
//       
//       myVars.insert(nameA,pA,true);
// //       myVars.insert(nameC,pC,true);
//       myVars.insert(nameD,pD,true);
//       myVars.insert(nameE,pE,true);
//       reportStatus(app_log());
//       return true;
//     }
//     
//     void resetTargetParticleSet(ParticleSet& P) 
//     {
//       Ie_table = DistanceTable::add(CenterRef,P);
//       ee_table = DistanceTable::add(P);
//     }
//     
//     void checkInVariables(opt_variables_type& active)
//     {
//       active.insertFrom(myVars);
//     }
//     
//     void checkOutVariables(const opt_variables_type& active)
//     {
//       myVars.getIndex(active);
//       myVars.print(std::cout);
//     }
//     
//     void resetParameters(const opt_variables_type& active)
//     {
//       int ia=myVars.where(0); if(ia>-1) pA=active[ia];
// //       int ic=myVars.where(1); if(ic>-1) pC=active[ic];
//       int id=myVars.where(1); if(id>-1) pD=active[id];
//       int ie=myVars.where(2); if(ie>-1) pE=active[ie];
//       reset(pA,pB,pC,pD,pE);
//     }
//     void reset(ValueType A, ValueType B, ValueType C, ValueType D, ValueType E){
//       pA=A;
//       pB=B;
//       pC=C;
//       pD=D;
//       pE=E;
//       
// //       myVars.clear();
// //       myVars.insert(nameA,pA,true);
// //       myVars.insert(nameB,pB,true);
// //       myVars.insert(nameC,pC,true);
// //       myVars.insert(nameD,pD,true);
//     }
//     
//     RealType
//     evaluateLog(ParticleSet& P, 
//         ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
// 	{
// 	 return std::log(evaluate(P,G,L)); 
// 	}
//     
//     ValueType ratio(ParticleSet& P, int iat,
// 			    ParticleSet::ParticleGradient_t& dG,
// 			    ParticleSet::ParticleLaplacian_t& dL)
// 			    {return 0; }
//     
//     void acceptMove(ParticleSet& P, int iat){}
//     void restore(int iat){}
//     ValueType ratio(ParticleSet& P, int iat){return 0;}
//     void update(ParticleSet& P, 
// 			ParticleSet::ParticleGradient_t& dG, 
// 			ParticleSet::ParticleLaplacian_t& dL,
// 			int iat) {}
//     RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
//     
//     RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
//     RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
//     void copyFromBuffer(ParticleSet& P, BufferType& buf){}
//     
//      
//     void reportStatus(ostream& os)
//     {
//       os<<"WF parameters: A="<<pA<<"  B="<<pB<<"  C="<<pC<<"  D="<<pD<<"  E="<<pE<<endl;
//     }
//     
//     ValueType evaluate(ParticleSet& P, 
// 			ParticleSet::ParticleGradient_t& G, 
// 			ParticleSet::ParticleLaplacian_t& L)
//     {
//       ValueType r0 = Ie_table->r(0);
//       ValueType r1 = Ie_table->r(1);
//       ValueType r01 = ee_table->r(0);
//       
//       ValueType s = r0+r1;
//       ValueType t = r0-r1;
//       ValueType u = r01;
//       ValueType v = r1*r1+r0*r0;
//       ValueType w = r0*r0-r1*r1;
//       ValueType w2 = w*w;
//       
//       ValueType vEr1r0m2 = std::pow(v,pE*r0*r1-2);
//       ValueType vEr1r0m1 = std::pow(v,pE*r0*r1-1);
//       ValueType vEr1r0 = v*vEr1r0m1;
//       ValueType logv = std::log(v);
//       
//       ValueType expm2s = std::exp(-2*s);
//       ValueType expAu = std::exp(pA*u);
//       ValueType part1 = 1.0+0.5*u*expAu ;
//       ValueType part2 = 1.0 + pC*w2*u*u + pD*vEr1r0;
//       ValueType mpart1 = 1.0/part1;
//       ValueType mpart2 = 1.0/part2;
//       
// //       Gradients
//       ValueType dr1P2A = 4.0*pC*w*u*u*r0;
//       ValueType dr1P2B = pD*pE*vEr1r0m1*r1*(2.0*r0*r0 + v*logv);
//       ValueType dr1P2 =  dr1P2A + dr1P2B;
//       ValueType dr12P1 = 0.5*(1.0+pA*u)*expAu;
//       ValueType dr12P2 = 2.0*pC*w2*u;
//       
//       ValueType F01 = (-2.0 + dr1P2*mpart2 );
//       ValueType F02 = (dr12P1*mpart1 + dr12P2*mpart2 );
//       PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);
// 
//       ValueType dr2P2 = -4.0*pC*w*u*u*r1 + pD*pE*vEr1r0m1*r0*(2.0*r1*r1 + v*logv);
//       ValueType F11 = (-2.0 + dr2P2*mpart2);
//       ValueType F12 = F02;
//       PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F12*ee_table->rinv(0);      
//       G[0] = J0;
//       G[1] = J1;  
//       
//       ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - dr12P1*dr12P1*mpart1*mpart1;
//       L0 += 2.0*pC*w2*mpart2 - dr12P2*dr12P2*mpart2*mpart2;
//       L0 += pD*pE*r1*vEr1r0m2*(2.0*r0*(r0*r0 + 2.0*pE*r0*r0*r0*r1 + 3.0*r1*r1)+ pE*r1*v*logv*(4.0*r0*r0 + v*logv))*mpart2 - dr1P2*dr1P2*mpart2*mpart2;
//       L0 += 4.0*pC*u*u*(3.0*r0*r0-r1*r1)*mpart2;
// //       L0 += pD*pE*vEr1r0m1*(4*r0*r1*(2.0 + pE*r0*r1) + pE*logv*(8.0*r0*r0*r1*r1 + v*v*logv))*mpart2 - dr1P2*dr1P2*mpart2*mpart2;
//       L0 -= 2.0* dot(Ie_table->dr(0),ee_table->dr(0))*Ie_table->rinv(0)*ee_table->rinv(0)* (8.0*pC*u*w*r0*mpart2 - dr1P2*dr12P2*mpart2*mpart2);
//       L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;
// // 
//       ValueType L1 = (0.5*pA*pA*u+pA)*expAu*mpart1 - dr12P1*dr12P1*mpart1*mpart1;
//       L1 += 2.0*pC*w2*mpart2 - dr12P2*dr12P2*mpart2*mpart2;
//       
//       L1 += pD*pE*r0*vEr1r0m2*(2.0*r1*(r1*r1 + 2.0*pE*r1*r1*r1*r0 + 3.0*r0*r0)+ pE*r0*v*logv*(4.0*r1*r1 + v*logv))*mpart2 - dr2P2*dr2P2*mpart2*mpart2;
//       L1 += -4.0*pC*u*u*(r0*r0-3.0*r1*r1)*mpart2;
// //       L0 += pD*pE*vEr1r0m1*(4*r0*r1*(2.0 + pE*r0*r1) + pE*logv*(8.0*r0*r0*r1*r1 + v*v*logv))*mpart2 - dr1P2*dr1P2*mpart2*mpart2;
//       L1 += 2.0* dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* (-8.0*pC*u*w*r1*mpart2 - dr2P2*dr12P2*mpart2*mpart2);
//       L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F12;
//       
//       L[0]= L0;
//       L[1]= L1; 
//       
// //       app_log()<<L0<<"  "<<L1<<"  "<<dot(J0,J0)<<"  "<<dot(J1,J1)<<endl;
//       
//        return expm2s*part1*part2;
//     }
// 
// 
// };


















/// One electron Helium WF. Other electron is frozen at center.
class NewCompactHeliumOne: public OrbitalBase
{
  public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
//     DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD, pK;
    string nameA,nameB,nameC,nameD;

    NewCompactHeliumOne(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0)
//     , ee_table(0) 
    {
      Ie_table = DistanceTable::add(CenterRef,electrons);
//       ee_table = DistanceTable::add(electrons);
      
      ///Default parameterization is from  CW David, PRA 74 059904(E) (2006)
//       pA=-1.013;
//       pB=0.2119;
//       pC=0.1406;
//       pD=-0.003;
    }
    
    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
      NewCompactHeliumOne* cloned = new NewCompactHeliumOne(tqp,CenterRef);
      cloned->nameA = nameA;
      cloned->nameB = nameB;
      cloned->nameC = nameC;
      cloned->nameD = nameD;
      cloned->myVars.insert(nameA,pA,true);
//       cloned->myVars.insert(nameB,pB,true);
//       cloned->myVars.insert(nameC,pC,true);
//       cloned->myVars.insert(nameD,pD,true);
//       cloned->OrbitalName=OrbitalName;
      cloned->reset(pA);
      
      return cloned;
    }
    
    
    bool put(xmlNodePtr cur){
      pA=-0.08;
//       pB= 0.2119256858;
//       pC=0.1416353426;
//       pD=-0.0113074102;
      
      
      OrbitalName="CHe4";
      string HEprefix("HE");
      OhmmsAttributeSet bb;
      bb.add(HEprefix,"id");
      bb.add(OrbitalName,"name");
      bb.put(cur);
      
      std::stringstream sstr;
      sstr << HEprefix << "_A";
      nameA = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_B";
//       nameB = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_C";
//       nameC = sstr.str();
//       sstr.str("");
//       sstr << HEprefix << "_D";
//       nameD = sstr.str();
      
      cur=cur->children;
      while(cur != NULL)
      {
	string pname="0";
	OhmmsAttributeSet aa;
	aa.add(pname,"name");
	aa.put(cur);
	if(pname[0]=='A') putContent(pA,cur);
// 	if(pname[0]=='B') putContent(pB,cur);
// 	if(pname[0]=='C') putContent(pC,cur);
// 	if(pname[0]=='D') putContent(pD,cur);
	cur=cur->next;
      }
      reset(pA );   
      
      
      myVars.insert(nameA,pA,true);
//       myVars.insert(nameB,pB,true);
//       myVars.insert(nameC,pC,true);
//       myVars.insert(nameD,pD,true);
      reportStatus(app_log());
      return true;
    }
    
    void resetTargetParticleSet(ParticleSet& P) 
    {
      Ie_table = DistanceTable::add(CenterRef,P);
//       ee_table = DistanceTable::add(P);
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }
    
    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      myVars.print(std::cout);
    }
    
    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) pA=active[ia];
//       int ib=myVars.where(1); if(ib>-1) pB=active[ib];
//       int ic=myVars.where(2); if(ic>-1) pC=active[ic];
//       int id=myVars.where(3); if(id>-1) pD=active[id];
      reset(pA );
    }
    void reset(ValueType A ){
      pA=A;
//       pB=B;
//       pC=C;
//       pD=D;
//       pK=B+C+D;
      
//       myVars.clear();
//       myVars.insert(nameA,pA,true);
//       myVars.insert(nameB,pB,true);
//       myVars.insert(nameC,pC,true);
//       myVars.insert(nameD,pD,true);
    }
    
    RealType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
	{
	 return std::log(evaluate(P,G,L)); 
	}
    
    ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL)
			    {return 0; }
    
    void acceptMove(ParticleSet& P, int iat){}
    void restore(int iat){}
    ValueType ratio(ParticleSet& P, int iat){return 0;}
    void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}
    
    RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
    void copyFromBuffer(ParticleSet& P, BufferType& buf){}
    
     
    void reportStatus(ostream& os)
    {
      os<<nameA<<"  "<<pA <<endl;
    }
    
    ValueType evaluate(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L)
    {
      ValueType r0 = Ie_table->r(0);
//       ValueType r1 = 0.0;
//       ValueType r01 = r0;
//       
//       ValueType s = r0+r1;
//       ValueType t = r0-r1;
      ValueType u = r0;
      
      ValueType expm2s = std::exp(-2*r0);
      ValueType expAu = std::exp(pA*r0);
      ValueType part1 = 1.0+0.5*r0*expAu;
//       ValueType part2 = 1.0 + pB*s*u + pC*t*t + pD*u*u;
      ValueType mpart1 = 1.0/part1;
//       ValueType mpart2 = 1.0/part2;
      
//       Gradients
//       ValueType bu2ct  = (pB*u+2.0*pC*t)*mpart2;
//       ValueType bs2du  = (pB*s+2.0*pD*u)*mpart2;
//       ValueType bu2mct = (pB*u-2.0*pC*t)*mpart2;
      ValueType upart = 0.5*(1+pA*u)*expAu*mpart1;

      ValueType F01 = (-2.0 );
      ValueType F02 = (upart );
      PosType J0 = Ie_table->dr(0)*(F01+F02)*Ie_table->rinv(0);
      
//       ValueType F11 = (-2.0 + bu2mct );
//       PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);      
      G[0] += J0;
//       G[1] = J1;  
      
     ValueType L0 = (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart ; 
      L0 += 2.0*Ie_table->rinv(0)*(F01+F02) ;

//       ValueType L1= (0.5*pA*pA*u+pA)*expAu*mpart1 - upart*upart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2mct*bu2mct;
//       L1 += dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2mct);
//       L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02; 
      
      L[0] += L0;
//       L[1]= L1; 
      
       return expm2s*part1 ;
    }


};



}
#endif