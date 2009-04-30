#ifndef COMPACT_HELIUM_H
#define COMPACT_HELIUM_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"

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


    bool put(xmlNodePtr cur) {
        pA=-9.85875926e-02;
        pB=6.42881288e-02 ;
        pC=7.47029314e-02 ;
        pD=-8.07044762e-02;


        OrbitalName="CHe4";
        OhmmsAttributeSet bb;
        bb.add(OrbitalName,"name");
        bb.put(cur);

        cur=cur->children;
        while (cur != NULL)
        {
            string pname="0";
            string pid="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.add(pid,"id");
            aa.put(cur);
            if (pname[0]=='A')
            {
                putContent(pA,cur);
                nameA=pid;
            }
            if (pname[0]=='B')
            {
                putContent(pB,cur);
                nameB=pid;
            }
            if (pname[0]=='C')
            {
                putContent(pC,cur);
                nameC=pid;
            }
            if (pname[0]=='D')
            {
                putContent(pD,cur);
                nameD=pid;
            }
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
//       myVars.print(std::cout);
    }

    void resetParameters(const opt_variables_type& active)
    {
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        int ib=myVars.where(1);
        if (ib>-1) myVars[1]=pB=active[ib];
        int ic=myVars.where(2);
        if (ic>-1) myVars[2]=pC=active[ic];
        int id=myVars.where(3);
        if (id>-1) myVars[3]=pD=active[id];
        reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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


    bool put(xmlNodePtr cur) {
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
        while (cur != NULL)
        {
            string pname="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.put(cur);
            if (pname[0]=='A') putContent(pA,cur);
            if (pname[0]=='B') putContent(pB,cur);
            if (pname[0]=='C') putContent(pC,cur);
            if (pname[0]=='D') putContent(pD,cur);
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        int ib=myVars.where(1);
        if (ib>-1) myVars[1]=pB=active[ib];
        int ic=myVars.where(2);
        if (ic>-1) myVars[2]=pC=active[ic];
        int id=myVars.where(3);
        if (id>-1) myVars[3]=pD=active[id];
        reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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


    bool put(xmlNodePtr cur) {
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
        while (cur != NULL)
        {
            string pname="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.put(cur);
            if (pname[0]=='A') putContent(pA,cur);
            if (pname[0]=='B') putContent(pB,cur);
            if (pname[0]=='C') putContent(pC,cur);
            if (pname[0]=='D') putContent(pD,cur);
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        int ib=myVars.where(1);
        if (ib>-1) myVars[1]=pB=active[ib];
//       int ic=myVars.where(2); if(ic>-1) pC=active[ic];
        int id=myVars.where(2);
        if (id>-1) myVars[2]=pD=active[id];
        reset(pA,pB,pC,pD);
    }
    void reset(ValueType A, ValueType B, ValueType C, ValueType D) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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

class DSimpleCompactHelium: public DiffOrbitalBase
{
public:
    ParticleSet& CenterRef ;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA, pB, pC, pD;
    string nameA,nameB,nameC,nameD;
    ///variables handled by this orbital
    opt_variables_type myVars;

    DSimpleCompactHelium(ParticleSet& electrons, ParticleSet& Ion):CenterRef(Ion), Ie_table(0), ee_table(0)
    {
        Ie_table = DistanceTable::add(CenterRef,electrons);
        ee_table = DistanceTable::add(electrons);
    }

    DiffOrbitalBase* makeClone(ParticleSet& tqp) const
    {
        DSimpleCompactHelium* cloned = new DSimpleCompactHelium(tqp,CenterRef);
        cloned->nameA = nameA;
        cloned->myVars.insert(nameA,pA,true);
        cloned->reset(pA);

        return cloned;
    }

    void resetTargetParticleSet(ParticleSet& P)
    {
        Ie_table = DistanceTable::add(CenterRef,P);
        ee_table = DistanceTable::add(P);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
        myVars.getIndex(active);
        myVars.print(std::cout);
    }

    void resetParameters(const opt_variables_type& active)
    {
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA);
    }
    void reset(ValueType A) {
        pA=A;
    }
    void reportStatus(ostream& os)
    {
        myVars.print(os);
    }

    void evaluateDerivatives(ParticleSet& P, RealType ke0,
                             const opt_variables_type& active,
                             vector<RealType>& dlogpsi,
                             vector<RealType>& dhpsioverpsi)
    {
//     assert( &els == &P);
        ValueType dLogPsi=0.0;
        ValueType d2rdLogPsi=0.0;
        PosType gradLogPsi;
        ValueType u = ee_table->r(0);
        ValueType expAu = std::exp(pA*u);
        ValueType part1 = 1.0+0.5*u*expAu;
        ValueType mpart1 = 1.0/part1;
        dLogPsi = 0.5*u*u*mpart1*expAu;
//       PosType R01 = P.R[0] - P.R[1];
//       ValueType r01 = std::sqrt(R01[0]*R01[0]+ R01[1]*R01[1]+ R01[2]*R01[2]);
//       cout<<pA<<"  "<<u<<"  "<<r01<<endl;
        ValueType Gval = (expAu*u*(1.0 + 0.5*pA*u + 0.25*u*expAu)*mpart1*mpart1);
        gradLogPsi = ee_table->dr(0)*ee_table->rinv(0)*Gval;
        qmcplusplus::PtclOnLatticeTraits::ParticleGradient_t elgrads(2);
        elgrads[0]=-1.0*gradLogPsi;
        elgrads[1]=gradLogPsi;

        d2rdLogPsi = expAu*(1.0+pA*u*(2.0+pA*u*(0.5-0.25*expAu*u)))*mpart1*mpart1*mpart1;
        d2rdLogPsi += 2.0*ee_table->rinv(0)*Gval;
        int kk=myVars.where(0);
        if (kk>-1)
        {
//       cout<<pA<<"  "<<dLogPsi<<endl;
            dlogpsi[kk]= dLogPsi;
            dhpsioverpsi[kk]= -(d2rdLogPsi + Dot(P.G,elgrads));
        }
    }


};

class SimpleCompactHelium: public OrbitalBase
{
public:
    ParticleSet& CenterRef,els;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA;
    string nameA;
//     DiffOrbitalBase* dPsi;
//     DSimpleCompactHelium Dp;

    SimpleCompactHelium(ParticleSet& electrons, ParticleSet& Ion): els(electrons), CenterRef(Ion), Ie_table(0), ee_table(0)
//     , Dp(electrons,Ion)
    {
        Ie_table = DistanceTable::add(CenterRef,electrons);
        ee_table = DistanceTable::add(electrons);

//         dPsi = &Dp;
//         setDiffOrbital(dPsi);
    }

    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
        SimpleCompactHelium* cloned = new SimpleCompactHelium(tqp,CenterRef);

        cloned->nameA = nameA;
        cloned->myVars=(myVars);
        cloned->reset(pA);

//         (cloned->Dp).myVars=cloned->myVars);
//         (cloned->Dp).nameA=cloned->nameA;
//         (cloned->Dp).reset(pA);
        return cloned;
    }


    bool put(xmlNodePtr cur) {
        pA=-0.08;


        OrbitalName="He2";
        OhmmsAttributeSet bb;
        bb.add(OrbitalName,"name");
        bb.put(cur);

        nameA = "He2_A_dflt";
        cur=cur->children;
        while (cur != NULL)
        {
            string pid="0";
            OhmmsAttributeSet aa;
            aa.add(nameA,"name");
            aa.add(pid,"id");
            aa.put(cur);
            if (pid[0]=='A') putContent(pA,cur);
            cur=cur->next;
        }
        reset(pA);
        myVars.insert(nameA,pA,nameA != "He2_A_dflt");
        reportStatus(app_log());
//         Dp.myVars.insert(nameA,pA,nameA != "He2_A_dflt");
//         Dp.nameA=nameA;
        return true;
    }

    void resetTargetParticleSet(ParticleSet& P)
    {
        Ie_table = DistanceTable::add(CenterRef,P);
        ee_table = DistanceTable::add(P);
//         dPsi->resetTargetParticleSet(P);
    }

    void checkInVariables(opt_variables_type& active)
    {
        active.insertFrom(myVars);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
        myVars.getIndex(active);
        myVars.print(std::cout);
//         dPsi->checkOutVariables(active);
    }

    void resetParameters(const opt_variables_type& active)
    {
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA);
//         dPsi->resetParameters(active);
    }
    void reset(ValueType A) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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
    
        void evaluateDerivatives(ParticleSet& P, RealType ke0,
                             const opt_variables_type& active,
                             vector<RealType>& dlogpsi,
                             vector<RealType>& dhpsioverpsi)
    {
//     assert( &els == &P);
        ValueType dLogPsi=0.0;
        ValueType d2rdLogPsi=0.0;
        PosType gradLogPsi;
        ValueType u = ee_table->r(0);
        ValueType expAu = std::exp(pA*u);
        ValueType part1 = 1.0+0.5*u*expAu;
        ValueType mpart1 = 1.0/part1;
        dLogPsi = 0.5*u*u*mpart1*expAu;
//       PosType R01 = P.R[0] - P.R[1];
//       ValueType r01 = std::sqrt(R01[0]*R01[0]+ R01[1]*R01[1]+ R01[2]*R01[2]);
//       cout<<pA<<"  "<<u<<"  "<<r01<<endl;
        ValueType Gval = (expAu*u*(1.0 + 0.5*pA*u + 0.25*u*expAu)*mpart1*mpart1);
        gradLogPsi = ee_table->dr(0)*ee_table->rinv(0)*Gval;
        qmcplusplus::PtclOnLatticeTraits::ParticleGradient_t elgrads(2);
        elgrads[0]=-1.0*gradLogPsi;
        elgrads[1]=gradLogPsi;

        d2rdLogPsi = expAu*(1.0+pA*u*(2.0+pA*u*(0.5-0.25*expAu*u)))*mpart1*mpart1*mpart1;
        d2rdLogPsi += 2.0*ee_table->rinv(0)*Gval;
        int kk=myVars.where(0);
        if (kk>-1)
        {
//       cout<<pA<<"  "<<dLogPsi<<endl;
            dlogpsi[kk]= dLogPsi;
            dhpsioverpsi[kk]= -(d2rdLogPsi + Dot(P.G,elgrads));
        }
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


    bool put(xmlNodePtr cur) {
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
        while (cur != NULL)
        {
            string pname="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.put(cur);
            if (pname[0]=='A') putContent(pA,cur);
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA);
    }
    void reset(ValueType A) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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


    bool put(xmlNodePtr cur) {
        pA=-0.08;
        pB=-2.0;


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
        while (cur != NULL)
        {
            string pname="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.put(cur);
            if (pname[0]=='A') putContent(pA,cur);
            if (pname[0]=='B') putContent(pB,cur);
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        int ib=myVars.where(1);
        if (ib>-1) myVars[1]=pB=active[ia];
        reset(pA,pB);
    }
    void reset(ValueType A,ValueType B) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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


class DSingleSlaterOrbital: public DiffOrbitalBase
{
public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA;
    string nameA;
    opt_variables_type myVars;

    DSingleSlaterOrbital(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0)
    {
        Ie_table = DistanceTable::add(CenterRef,electrons);
    }

    DiffOrbitalBase* makeClone(ParticleSet& tqp) const
    {
        DSingleSlaterOrbital* cloned = new DSingleSlaterOrbital(tqp,CenterRef);
        cloned->nameA = nameA;
        cloned->myVars.insertFrom(myVars);
        return cloned;
    }


    bool put(xmlNodePtr cur) {
        return true;
    }

    void resetTargetParticleSet(ParticleSet& P)
    {
        Ie_table = DistanceTable::add(CenterRef,P);
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA);
    }
    void reset(ValueType A ) {
        pA=A;
    }

    void reportStatus(ostream& os)
    {
        myVars.print(os);
    }

    void evaluateDerivatives(ParticleSet& P, RealType ke0,
                             const opt_variables_type& active,
                             vector<RealType>& dlogpsi,
                             vector<RealType>& dhpsioverpsi)
    {
        ValueType r = Ie_table->r(0);
        PosType J0 = -1.0*Ie_table->dr(0)* Ie_table->rinv(0);

        int kk=myVars.where(0);
        if (kk>-1)
        {
            dlogpsi[kk]= r;
            dhpsioverpsi[kk]=  -1.0*(Ie_table->rinv(0) + pA);
        }
    }


};

class SingleSlaterOrbital: public OrbitalBase
{
public:
    ParticleSet& CenterRef;
    DistanceTableData* Ie_table;
    DistanceTableData* ee_table;
    ValueType pA;
    string nameA;
//     DSingleSlaterOrbital DP;
//     DiffOrbitalBase* dPsi;

    SingleSlaterOrbital(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0), ee_table(0) 
//     , DP(electrons,Ion)
    {
        Ie_table = DistanceTable::add(CenterRef,electrons);
//         dPsi = &DP;
//         setDiffOrbital(dPsi);
    }

    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
        SingleSlaterOrbital* cloned = new SingleSlaterOrbital(tqp,CenterRef);
        cloned->nameA = nameA;
        cloned->myVars=(myVars);
        cloned->pA=pA;
//         (cloned->DP).myVars.insertFrom(myVars);
//         (cloned->DP).nameA=nameA;
//         (cloned->DP).pA=pA;
        return cloned;
    }


    bool put(xmlNodePtr cur) {
        pA=-1.0;


        OrbitalName="H";
        OhmmsAttributeSet bb;
        bb.add(OrbitalName,"name");
        bb.put(cur);

//       std::stringstream sstr;
//       sstr << OrbitalName << "_A";
        nameA = "H_A_dflt";

        cur=cur->children;
        while (cur != NULL)
        {
            string pname="H_A_dflt";
            string pid="0";
            OhmmsAttributeSet aa;
            aa.add(pname,"name");
            aa.add(pid,"id");
            aa.put(cur);
            if (pid=="A")
            {
                putContent(pA,cur);
                nameA=pname;
            }
            cur=cur->next;
        }
        reset(pA);

        myVars.insert(nameA,pA, nameA !="H_A_dflt");
        reportStatus(app_log());
//         DP.myVars.insertFrom(myVars);
//         DP.nameA=nameA;
//         DP.reset(pA);
        return true;
    }

    void resetTargetParticleSet(ParticleSet& P)
    {
        Ie_table = DistanceTable::add(CenterRef,P);
//         dPsi->resetTargetParticleSet(P);
    }

    void checkInVariables(opt_variables_type& active)
    {
        active.insertFrom(myVars);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
        myVars.getIndex(active);
        myVars.print(std::cout);
//         dPsi->checkOutVariables(active);
    }

    void resetParameters(const opt_variables_type& active)
    {
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA);
//         dPsi->resetParameters(active);
    }
    void reset(ValueType A ) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
    }

    ValueType evaluate(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
    {
        ValueType r = Ie_table->r(0);

        ValueType expm2s = std::exp(pA*r);
        ValueType F01 = (pA);
        PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0);
        G[0] += J0;

        ValueType L0  = 2.0*pA*Ie_table->rinv(0);

        L[0] += L0;
        return expm2s;
    }
        void evaluateDerivatives(ParticleSet& P, RealType ke0,
                             const opt_variables_type& active,
                             vector<RealType>& dlogpsi,
                             vector<RealType>& dhpsioverpsi)
    {
        ValueType r = Ie_table->r(0);
        PosType J0 = -1.0*Ie_table->dr(0)* Ie_table->rinv(0);

        int kk=myVars.where(0);
        if (kk>-1)
        {
            dlogpsi[kk]= r;
            dhpsioverpsi[kk]=  -1.0*(Ie_table->rinv(0) + pA);
        }
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
    ValueType pA;
    string nameA;

    NewCompactHeliumOne(ParticleSet& electrons, ParticleSet& Ion): CenterRef(Ion), Ie_table(0)
//     , ee_table(0)
    {
        Ie_table = DistanceTable::add(CenterRef,electrons);
    }

    OrbitalBase* makeClone(ParticleSet& tqp) const
    {
        NewCompactHeliumOne* cloned = new NewCompactHeliumOne(tqp,CenterRef);
        cloned->nameA = nameA;
        cloned->myVars=myVars;
        cloned->reset(pA);

        return cloned;
    }


    bool put(xmlNodePtr cur) {
        pA=-0.08;

        OrbitalName="CHe4";
        OhmmsAttributeSet bb;
        bb.add(OrbitalName,"name");
        bb.put(cur);
        
        nameA = "He1_A_dflt";

        cur=cur->children;
        while (cur != NULL)
        {
            string pid="0";
            OhmmsAttributeSet aa;
            aa.add(nameA,"name");
            aa.add(pid,"id");
            aa.put(cur);
            if (pid[0]=='A') putContent(pA,cur);
            cur=cur->next;
        }
        reset(pA );
        myVars.insert(nameA,pA,nameA!="He1_A_dflt");
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
        int ia=myVars.where(0);
        if (ia>-1) myVars[0]=pA=active[ia];
        reset(pA );
    }
    void reset(ValueType A ) {
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
    {
        return 0;
    }

    void acceptMove(ParticleSet& P, int iat) {}
    void restore(int iat) {}
    ValueType ratio(ParticleSet& P, int iat) {
        return 0;
    }
    void update(ParticleSet& P,
                ParticleSet::ParticleGradient_t& dG,
                ParticleSet::ParticleLaplacian_t& dL,
                int iat) {}
    RealType evaluateLog(ParticleSet& P,BufferType& buf) {
        return 0;
    }

    RealType registerData(ParticleSet& P, BufferType& buf) {
        return 0;
    }
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) {
        return 0;
    }
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}


    void reportStatus(ostream& os)
    {
        myVars.print(os);
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



//class VariableGaussianBond: public OrbitalBase
//{
//public:
//ValueType pA, pB, pC , RC, RCmax;
//string nameA,nameB,nameC,nameD ;
//vector<ValueType> center1, center2, line, BondCenter;


//VariableGaussianBond(ParticleSet& electrons, vector<ValueType> CenterPosition , vector<ValueType> IonPosition ) : center1(CenterPosition), center2(IonPosition), BondCenter(IonPosition)
//{
//RC=0.0;
//for( int i=0;i<3;i++){
//line.push_back(IonPosition[i]-CenterPosition[i]);
//RC += line[i]*line[i];
//}
//RCmax = RC = std::sqrt(RC) - 1e-2;
//}

//OrbitalBase* makeClone(ParticleSet& tqp) const
//{
//VariableGaussianBond* cloned = new VariableGaussianBond(tqp,center1,center2);
//cloned->nameA = nameA;
//cloned->nameB = nameB;
//cloned->nameC = nameC;
//cloned->nameD = nameD;
//cloned->myVars.insert(nameA,pA,true);
//cloned->myVars.insert(nameB,pB,true);
//cloned->myVars.insert(nameC,pC,true);
//cloned->myVars.insert(nameD,RC,true);
//cloned->reset(pA,pB,pC,RC);

//return cloned;
//}


//bool put(xmlNodePtr cur){
//return true;
//}

//void resetTargetParticleSet(ParticleSet& P)
//{
//}

//void checkInVariables(opt_variables_type& active)
//{
//active.insertFrom(myVars);
//}

//void checkOutVariables(const opt_variables_type& active)
//{
//myVars.getIndex(active);
//}

//void resetParameters(const opt_variables_type& active)
//{
//int ia=myVars.where(0); if(ia>-1) myVars[0]=pA=active[ia];
//int ib=myVars.where(1); if(ib>-1) myVars[1]=pB=active[ib];
//int ic=myVars.where(2); if(ic>-1) myVars[2]=pC=active[ic];
//int id=myVars.where(3); if(id>-1) myVars[3]=RC=active[id];
//reset(pA,pB,pC,RC);
//}

//void reset(ValueType A, ValueType B, ValueType C, ValueType D){
//pA=A; pB=B; pC=C; RC=D;
////safety checks
//pC = std::max( std::min(1.0,pC) , 0.0);
//RC = std::min(RCmax,RC);

////move bond center along line
//for( int i=0;i<3;i++) BondCenter[i] = center2[i]-c*line[i];


//}

//inline RealType expPart(RealType rx, RealType ry, RealType rz)
//{
//RealType dx = rx-BondCenter[0];
//RealType dy = ry-BondCenter[1];
//RealType dz = rz-BondCenter[2];
//return std::exp(-(pA*(dx*dx+dy*dy) + pB*dz*dz) );
//}

//inline RealType logexpPart(RealType rx, RealType ry, RealType rz)
//{
//RealType dx = rx-BondCenter[0];
//RealType dy = ry-BondCenter[1];
//RealType dz = rz-BondCenter[2];
//return -(pA*(dx*dx+dy*dy) + pB*dz*dz) ;
//}



//RealType
//evaluateLog(ParticleSet& P,
//ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
//{
//RealType Value(0.0);



//return Value;
//}

//ValueType ratio(ParticleSet& P, int iat,
//ParticleSet::ParticleGradient_t& dG,
//ParticleSet::ParticleLaplacian_t& dL)
//{return 0; }

//void acceptMove(ParticleSet& P, int iat){}
//void restore(int iat){}
//ValueType ratio(ParticleSet& P, int iat){return 0;}
//void update(ParticleSet& P,
//ParticleSet::ParticleGradient_t& dG,
//ParticleSet::ParticleLaplacian_t& dL,
//int iat) {}
//RealType evaluateLog(ParticleSet& P,BufferType& buf){return 0;}

//RealType registerData(ParticleSet& P, BufferType& buf){return 0;}
//RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false){return 0;}
//void copyFromBuffer(ParticleSet& P, BufferType& buf){}


//void reportStatus(ostream& os)
//{
//myVars.print(os);
//}

//ValueType evaluate(ParticleSet& P,
//ParticleSet::ParticleGradient_t& G,
//ParticleSet::ParticleLaplacian_t& L)
//{
//ValueType r0 = Ie_table->r(0);
//ValueType r1 = Ie_table->r(1);
//ValueType r01 = ee_table->r(0);

//ValueType s = r0+r1;
//ValueType t = r0-r1;
//ValueType u = r01;

//ValueType expm2s = std::exp(-2*s);
//ValueType expAu = std::exp(pA*u);
//ValueType part1 = (1.0+0.5*u)*expAu;
//ValueType part2 = 1.0 + pB*s*u + pC*t*t + pD*u*u;
//ValueType mpart1 = 1.0/part1;
//ValueType mpart2 = 1.0/part2;

////       Gradients
//ValueType bu2ct  = (pB*u+2.0*pC*t)*mpart2;
//ValueType bs2du  = (pB*s+2.0*pD*u)*mpart2;
//ValueType bu2mct = (pB*u-2.0*pC*t)*mpart2;
//ValueType mupart = 1.0/(1.0+0.5*u);

//ValueType F01 = (-2.0 + bu2ct );
//ValueType F02 = (pA + 0.5*mupart + bs2du );
//PosType J0 = Ie_table->dr(0)*F01*Ie_table->rinv(0) - ee_table->dr(0)*F02*ee_table->rinv(0);

//ValueType F11 = (-2.0 + bu2mct );
//PosType J1 = Ie_table->dr(1)*F11*Ie_table->rinv(1) + ee_table->dr(0)*F02*ee_table->rinv(0);
//G[0] = J0;
//G[1] = J1;

//ValueType L0 = -0.25*mupart*mupart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2ct*bu2ct;
//L0 -= dot(Ie_table->dr(0),ee_table->dr(0))*Ie_table->rinv(0)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2ct);
//L0 += 2.0*Ie_table->rinv(0)*F01 + 2.0*ee_table->rinv(0)*F02;

//ValueType L1= -0.25*mupart*mupart + 2.0*(pD+pC)*mpart2 - bs2du*bs2du - bu2mct*bu2mct;
//L1 += dot(Ie_table->dr(1),ee_table->dr(0))*Ie_table->rinv(1)*ee_table->rinv(0)* ( 2.0*pB*mpart2 - 2.0*bs2du*bu2mct);
//L1 += 2.0*Ie_table->rinv(1)*F11 + 2.0*ee_table->rinv(0)*F02;

//L[0]= L0;
//L[1]= L1;

//return expm2s*part1*part2;
//}


//};




}
#endif
