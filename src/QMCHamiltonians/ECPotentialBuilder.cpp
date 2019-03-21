//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/L2Potential.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OneDimNumGridFunctor.h"
#ifdef QMC_CUDA
#include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
#include "QMCHamiltonians/LocalECPotential_CUDA.h"
#include "QMCHamiltonians/NonLocalECPotential_CUDA.h"
#endif

namespace qmcplusplus
{
/** constructor
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\param psi trial wavefunction
 */
ECPotentialBuilder::ECPotentialBuilder(QMCHamiltonian& h,
                                       ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi,
                                       Communicate* c):
  MPIObjectBase(c),
  hasLocalPot(false),hasNonLocalPot(false),hasL2Pot(false),
  targetH(h), IonConfig(ions), targetPtcl(els), targetPsi(psi)
{ }

bool ECPotentialBuilder::put(xmlNodePtr cur)
{
  if(localPot.empty())
  {
    int ng = IonConfig.getSpeciesSet().getTotalNum();
    localZeff.resize(ng,1);
    localPot.resize(ng,0);
    nonLocalPot.resize(ng,0);
    L2Pot.resize(ng,0);
  }
  std::string ecpFormat("table");
  std::string NLPP_algo(omp_get_nested()?"batched":"default");
#ifdef ENABLE_OFFLOAD
  NLPP_algo = "batched"; // set "batched" as the default
#endif
  std::string pbc("yes");
  std::string forces("no");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(ecpFormat,"format");
  pAttrib.add(NLPP_algo,"algorithm");
  pAttrib.add(pbc,"pbc");
  pAttrib.add(forces,"forces");
  pAttrib.put(cur);
  bool doForces = (forces == "yes") || (forces == "true");
  //const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"format");
  //if(t != NULL) {
  //  ecpFormat= (const char*)t;
  //}
  if(ecpFormat == "xml")
  {
    useXmlFormat(cur);
  }
  else
  {
    useSimpleTableFormat();
  }

  ///create LocalECPotential
  bool usePBC =
    !(IonConfig.Lattice.SuperCellEnum == SUPERCELL_OPEN || pbc =="no");


  if(hasLocalPot)
  {
    if(IonConfig.Lattice.SuperCellEnum == SUPERCELL_OPEN || pbc =="no")
    {
#ifdef QMC_CUDA
      LocalECPotential_CUDA* apot =
        new LocalECPotential_CUDA(IonConfig,targetPtcl);
#else
      LocalECPotential* apot = new LocalECPotential(IonConfig,targetPtcl);
#endif
      for(int i=0; i<localPot.size(); i++)
        if(localPot[i]) apot->add(i,localPot[i],localZeff[i]);
      targetH.addOperator(apot,"LocalECP");
    }
    else
    {
      if (doForces)
        app_log() << "  Will compute forces in CoulombPBCAB.\n" << std::endl;
#ifdef QMC_CUDA
      CoulombPBCAB_CUDA* apot=
        new CoulombPBCAB_CUDA(IonConfig,targetPtcl, doForces);
#else
      CoulombPBCAB* apot =
        new CoulombPBCAB(IonConfig,targetPtcl, doForces);
#endif
      for(int i=0; i<localPot.size(); i++)
      {
        if(localPot[i])
          apot->add(i,localPot[i]);
      }
      targetH.addOperator(apot,"LocalECP");
    }
  }
  if(hasNonLocalPot)
  {
    //resize the sphere
    RealType rc2=0.0;
#ifdef QMC_CUDA
    NonLocalECPotential_CUDA* apot =
      new NonLocalECPotential_CUDA(IonConfig,targetPtcl,targetPsi,usePBC,doForces);
#else
    NonLocalECPotential* apot =
      new NonLocalECPotential(IonConfig,targetPtcl,targetPsi, doForces, NLPP_algo=="batched");
#endif
    int nknot_max=0;
    for(int i=0; i<nonLocalPot.size(); i++)
    {
      if(nonLocalPot[i])
      {
        rc2=std::max(rc2,nonLocalPot[i]->Rmax);
        nknot_max=std::max(nknot_max,nonLocalPot[i]->nknot);
        apot->add(i,nonLocalPot[i]);
      }
    }
    app_log() << "\n  Using NonLocalECP potential \n"
              << "    Maximum grid on a sphere for NonLocalECPotential: "
              << nknot_max << std::endl;
    if(NLPP_algo=="batched") app_log() << "    Using batched ratio computing in NonLocalECP" << std::endl;

    targetPtcl.checkBoundBox(2*rc2);

    targetH.addOperator(apot,"NonLocalECP");
  }
  if(hasL2Pot)
  {
    L2Potential* apot = new L2Potential(IonConfig,targetPtcl,targetPsi);
    for(int i=0;i<L2Pot.size();i++)
      if(L2Pot[i])
        apot->add(i,L2Pot[i]);
    app_log()<< "\n  Using L2 potential"<<std::endl;
    targetH.addOperator(apot,"L2");
  }

  app_log().flush();
  return true;
}

void ECPotentialBuilder::useXmlFormat(xmlNodePtr cur)
{
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "pseudo")
    {
      std::string href("none");
      std::string ionName("none");
      std::string format("xml");
      //RealType rc(2.0);//use 2 Bohr
      OhmmsAttributeSet hAttrib;
      hAttrib.add(href,"href");
      hAttrib.add(ionName,"elementType");
      hAttrib.add(ionName,"symbol");
      hAttrib.add(format,"format");
      //hAttrib.add(rc,"cutoff");
      hAttrib.put(cur);
      SpeciesSet& ion_species(IonConfig.getSpeciesSet());
      int speciesIndex=ion_species.findSpecies(ionName);
      int chargeIndex=ion_species.findAttribute("charge");
      int AtomicNumberIndex=ion_species.findAttribute("atomicnumber");
      if(AtomicNumberIndex==-1) AtomicNumberIndex=ion_species.findAttribute("atomic_number");
      bool success=false;
      if(speciesIndex < ion_species.getTotalNum())
      {
        app_log() << std::endl << "  Adding pseudopotential for " << ionName << std::endl;
        RealType rmax=0.0;

        ECPComponentBuilder ecp(ionName,myComm);
        if(format == "xml")
        {
          if(href == "none")
          {
            success=ecp.put(cur);
          }
          else
          {
            success=ecp.parse(href,cur);
          }
        }
        else if(format == "casino")
        {
          //success=ecp.parseCasino(href,rc);
          success=ecp.parseCasino(href,cur);
        }
        if(success)
        {
          if(OHMMS::Controller->rank()==0) ecp.printECPTable();
          if(ecp.pp_loc)
          {
            localPot[speciesIndex]=ecp.pp_loc;
            localZeff[speciesIndex]=ecp.Zeff;
            hasLocalPot=true;
          }
          if(ecp.pp_nonloc)
          {
            hasNonLocalPot=true;
            nonLocalPot[speciesIndex]=ecp.pp_nonloc;
            rmax=std::max(rmax,ecp.pp_nonloc->Rmax);
          }
          if(ecp.pp_L2)
          {
            hasL2Pot=true;
            L2Pot[speciesIndex]=ecp.pp_L2;
            // should this be added or not?
            //rmax=std::max(rmax,ecp.pp_L2->rcut);
          }
          int rcutIndex=ion_species.addAttribute("rmax_core");
          ion_species(rcutIndex,speciesIndex)=rmax;
          if(chargeIndex == -1)
          {
            app_error() << "  Ion species " << ionName << " needs parameter \'charge\'" << std::endl;
            success=false;
          }
          else
          {
            RealType ion_charge = ion_species(chargeIndex,speciesIndex);
            if( std::fabs(ion_charge - ecp.Zeff) > 1e-4 )
            {
              app_error() << "  Ion species " << ionName << " charge " << ion_charge
                          << " pseudopotential charge " << ecp.Zeff << " mismatch!" << std::endl;
              success=false;
            }
          }
          if(AtomicNumberIndex == -1)
          {
            app_error() << "  Ion species " << ionName << " needs parameter \'atomicnumber\'" << std::endl;
            success=false;
          }
          else
          {
            int atomic_number = ion_species(AtomicNumberIndex,speciesIndex);
            if(atomic_number != ecp.AtomicNumber)
            {
              app_error() << "  Ion species " << ionName << " atomicnumber " << atomic_number
                          << " pseudopotential atomic-number " << ecp.AtomicNumber << " mismatch!" << std::endl;
              success=false;
            }
          }
        }
      }
      else
      {
        app_error() << "  Ion species " << ionName << " is not found." << std::endl;
      }
      if(!success)
      {
        app_error() << " Failed to add pseudopotential for element " << ionName << std::endl;
        APP_ABORT("ECPotentialBuilder::useXmlFormat failed!");
      }
    }
    cur=cur->next;
  }
}

/** reimplement simple table format used by NonLocalPPotential
 */
void ECPotentialBuilder::useSimpleTableFormat()
{
  SpeciesSet& Species(IonConfig.getSpeciesSet());
  int ng(Species.getTotalNum());
  int icharge(Species.addAttribute("charge"));
  for(int ig=0; ig<ng; ig++)
  {
    std::vector<RealType> grid_temp, pp_temp;
    std::string species(Species.speciesName[ig]);
    std::string fname = species+".psf";
    std::ifstream fin(fname.c_str(),std::ios_base::in);
    if(!fin)
    {
      ERRORMSG("Could not open file " << fname)
      exit(-1);
    }
    // Read Number of potentials (local and non) for this atom
    int npotentials;
    fin >> npotentials;
    int lmax=-1;
    int numnonloc=0;
    RealType rmax(0.0);
    app_log() << "  ECPotential for " << species << std::endl;
    NonLocalECPComponent* mynnloc=0;
    typedef OneDimCubicSpline<RealType> CubicSplineFuncType;
    for (int ij=0; ij<npotentials; ij++)
    {
      int angmom,npoints;
      fin >> angmom >> npoints;
      OneDimNumGridFunctor<RealType> inFunc;
      inFunc.put(npoints,fin);
      if(angmom < 0)
        //local potential, input is rescale by -r/z
      {
        RealType zinv=-1.0/Species(icharge,ig);
        int ng=npoints-1;
        RealType rf=5.0;
        ng=static_cast<int>(rf*100)+1;//use 1e-2 resolution
        GridType * agrid= new LinearGrid<RealType>;
        agrid->set(0,rf,ng);
        std::vector<RealType> pp_temp(ng);
        pp_temp[0]=0.0;
        for (int j=1; j<ng; j++)
        {
          RealType r((*agrid)[j]);
          pp_temp[j]=r*zinv*inFunc.splint(r);
        }
        pp_temp[ng-1]=1.0;
        RadialPotentialType *app = new RadialPotentialType(agrid,pp_temp);
        app->spline();
        localPot[ig]=app;
        app_log() << "    LocalECP l=" << angmom << std::endl;
        app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << std::endl;
        hasLocalPot=true; //will create LocalECPotential
      }
      else
      {
        hasNonLocalPot=true; //will create NonLocalECPotential
        if(mynnloc == 0)
          mynnloc = new NonLocalECPComponent;
        RealType rf=inFunc.rmax();
        GridType * agrid= new LinearGrid<RealType>;
        int ng=static_cast<int>(rf*100)+1;
        agrid->set(0.0,rf,ng);
        app_log() << "    NonLocalECP l=" << angmom << " rmax = " << rf << std::endl;
        app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << std::endl;
        std::vector<RealType> pp_temp(ng);
        //get the value
        pp_temp[0]=inFunc(0);
        for (int j=1; j<ng; j++)
        {
          pp_temp[j]=inFunc.splint((*agrid)[j]);
        }
        RadialPotentialType *app = new RadialPotentialType(agrid,pp_temp);
        app->spline();
        mynnloc->add(angmom,app);
        lmax=std::max(lmax,angmom);
        rmax=std::max(rmax,rf);
        numnonloc++;
      }
      if(mynnloc)
      {
        mynnloc->lmax=lmax;
        mynnloc->Rmax=rmax;
        app_log() << "    Maximum cutoff of NonLocalECP " << rmax << std::endl;
      }
    }
    fin.close();
    if(mynnloc)
    {
      nonLocalPot[ig]=mynnloc;
      int numsgridpts=0;
      std::string fname = species+".sgr";
      std::ifstream fin(fname.c_str(),std::ios_base::in);
      if(!fin)
      {
        app_error() << "Could not open file " << fname << std::endl;
        exit(-1);
      }
      PosType xyz;
      RealType weight;
      while(fin >> xyz >> weight)
      {
        mynnloc->addknot(xyz,weight);
        numsgridpts++;
      }
      //cout << "Spherical grid : " << numsgridpts << " points" << std::endl;
      mynnloc->resize_warrays(numsgridpts,numnonloc,lmax);
    }
  }//species
}
}
