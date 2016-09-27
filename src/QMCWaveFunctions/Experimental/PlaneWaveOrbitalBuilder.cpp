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
    
    
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/PlaneWaveOrbitalBuilder.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

PlaneWaveOrbitalBuilder::PlaneWaveOrbitalBuilder(ParticleSet& els, TrialWaveFunction& psi):
  OrbitalBuilderBase(els,psi)
{
  PWOSet = NULL;
}

PlaneWaveOrbitalBuilder::~PlaneWaveOrbitalBuilder()
{
  if(PWOSet)
    delete PWOSet;
}

//All data parsing is handled here, outside storage classes.
bool PlaneWaveOrbitalBuilder::put(xmlNodePtr cur)
{
  //Find how many up and down electrons in ParticleSet
  nup=-1;
  SpeciesSet tspecies=targetPtcl.getSpeciesSet();
  int MemberAttrib = tspecies.addAttribute("membersize");
  for(int i=0; i<tspecies.TotalNum; i++)
  {
    if(strcmp(tspecies.speciesName[i].c_str(),"u")==0)
    {
      nup = static_cast<int>(tspecies(MemberAttrib,i));
      upindx = i;
    }
  }
  if(nup<0)
  {
    LOGMSG("Up channel not found in electron particleset.\n");
    OHMMS::Controller->abort();
  }
  ndown = targetPtcl.getTotalNum() - nup;
  //
  //Get wavefunction data and parameters from XML and HDF5
  //
  std::string hdf5file;
  RealType ecut=-1.0;
  //Move through the XML tree and read basis information
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "basisset")
    {
      const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"ecut");
      if(aptr)
        ecut=atof((const char*)aptr);
    }
    else
      if(cname == "coefficients")
      {
        const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"hdata");
        if(!aptr)
        {
          LOGMSG("Error finding HDF5 filename in PlaneWaveOrbitalBuilder");
          OHMMS::Controller->abort();
        }
        hdf5file=(const char*)aptr;
        LOGMSG("\tReading wavefunction data from " << hdf5file);
      }
      else
        if(cname == "slaterdeterminant")
        {
          //OCCUPATION MODE NOT YET READ; ASSUMED GROUND STATE
          //Find which orbitals should be occupied in each determinant
          //Which spin channel? Ground state or excited?
          //Default: use spin channel 0 and ground-state determinants
          updetspinindex = downdetspinindex = 0;
          xmlNodePtr slaterdetnode=cur->children;
          while(slaterdetnode != NULL)
          {
            std::string cname2((const char*)(slaterdetnode->name));
            //Which determinant?
            if(cname2 == "determinant")
            {
              const xmlChar* aptr = xmlGetProp(slaterdetnode,(const xmlChar*)"id");
              if(strcmp((const char*)aptr,"updet")==0)
                //Up determinant
              {
                //Get occupation property
                xmlNodePtr detnode = slaterdetnode->children;
                while(detnode!=NULL)
                {
                  std::string cname3((const char*)(detnode->name));
                  if(cname3 == "occupation")
                  {
                    const xmlChar* aptr2 = xmlGetProp(detnode,(const xmlChar*)"spindataset");
                    updetspinindex = std::abs(atoi((const char*)aptr2));
                    //const xmlChar* aptr2 = xmlGetProp(detnode,(const xmlChar*)"mode");
                  }
                  detnode = detnode->next;
                }
              }
              else
                if(strcmp((const char*)aptr,"downdet")==0)
                  //Down determinant
                {
                  //Get occupation property
                  xmlNodePtr detnode = slaterdetnode->children;
                  while(detnode!=NULL)
                  {
                    std::string cname3((const char*)(detnode->name));
                    if(cname3 == "occupation")
                    {
                      const xmlChar* aptr2 = xmlGetProp(detnode,(const xmlChar*)"spindataset");
                      downdetspinindex = std::abs(atoi((const char*)aptr2));
                      //const xmlChar* aptr2 = xmlGetProp(detnode,(const xmlChar*)"mode");
                    }
                    detnode = detnode->next;
                  }
                }
            }
            slaterdetnode = slaterdetnode->next;
          }
        }
    cur=cur->next;
  }
  //Open the HDF5 file
  hid_t hfile = H5Fopen(hdf5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  //Check the wavefunction version number
  std::vector<int> version;
  HDFAttribIO<std::vector<int> > hdfver(version);
  hdfver.read(hfile,"version");
  LOGMSG("\tWavefunction HDF version: " << version[0] << "." << version[1]);
  //Only version 0.10 (test) supported now.
  if(version[0] == 0 && version[1] == 10)
  {
    ReadHDFWavefunction010(hfile,ecut);
  }
  else
  {
    LOGMSG("\n\tUnsupported wavefunction HDF version");
    OHMMS::Controller->abort();
  }
  H5Fclose(hfile);
  //Finished reading orbital data.
  //Build the Slater determinant from the prescription in the XML input
  typedef SlaterDeterminant<PlaneWaveOrbitalSet> SlaterDeterminant_t;
  typedef DiracDeterminant<PlaneWaveOrbitalSet> Det_t;
  SlaterDeterminant_t *sdet  = new SlaterDeterminant_t;
  //create up determinant
  LOGMSG("  \tBuilding updet with particles "<<targetPtcl.first(upindx)<<"-"<<targetPtcl.last(upindx)-1);
  Det_t *updet = new Det_t(*PWOSet,targetPtcl.first(upindx));
  updet->set(targetPtcl.first(upindx),nup); //first member in ptclset, length
  sdet->add(updet);
  //create down determinant
  LOGMSG("  \tBuilding downdet with particles "<<targetPtcl.first(1-upindx)<<"-"<<targetPtcl.last(1-upindx)-1<<"\n");
  Det_t *downdet = new Det_t(*PWOSet,targetPtcl.first(1-upindx));
  downdet->set(targetPtcl.first(1-upindx),ndown);
  sdet->add(downdet);
  //add Slater determinant to targetPsi
  targetPsi.addOrbital(sdet);
  return true;
}



//
// The read routine - get data from XML and H5. Process it and build orbitals.
//
void PlaneWaveOrbitalBuilder::ReadHDFWavefunction010(hid_t hfile, RealType& ecut)
{
  hid_t grp_id;
  RealType h5ecut;
  int nh5gvecs,nkpts,nbands,nspins;
  bool h5coefsreal;
  //Get the system run parameters and check that ecut and nh5gvecs are consistent
  char GrpName[128];
  sprintf(GrpName,"parameters");
  grp_id = H5Gopen(hfile,GrpName);
  //read the dataset - nkpts, h5ecut
  int idata;
  double ddata;
  HDFAttribIO<int> hdfint(idata);
  HDFAttribIO<double> hdfdbl(ddata);
  hdfint.read(grp_id,"num_twists");
  nkpts = idata;
  if(nkpts > 1)
  {
    LOGMSG("  Warning: HDF5 wavefunction data contains multiple twists");
    LOGMSG("  Only the first twist angle will be used in this code");
  }
  hdfint.read(grp_id,"num_bands");
  nbands = idata;
  if(nbands < nup || nbands < ndown)
  {
    LOGMSG("Not enough bands in h5 file");
    OHMMS::Controller->abort();
  }
  //How many spin channels are in wavefunction file?
  hdfint.read(grp_id,"num_spins");
  nspins = idata;
  if((updetspinindex >= nspins) ||( downdetspinindex >= nspins))
  {
    LOGMSG("  Attempt to use spin channel not present on HDF5 file.");
    OHMMS::Controller->abort();
  }
  //Are the PW coefficients in H5 file real or complex?
  hdfint.read(grp_id,"complex_coefficients");
  idata==0 ? h5coefsreal=true : h5coefsreal=false;
  if(h5coefsreal && nkpts > 1)
  {
    LOGMSG("HDF Wavefunction file is inconsistent");
    LOGMSG("Number of twists > 1 (" << nkpts << ")");
    LOGMSG("But coefficients are real.");
    OHMMS::Controller->abort();
  }
  hdfdbl.read(grp_id,"maximum_ecut");
  h5ecut = ddata;
  if(ecut > h5ecut)
  {
    ecut = h5ecut;
    LOGMSG("  Warning: chosen Ecut exceeds available data.");
    LOGMSG("  Ecut reset to: " << ecut);
  }
  else
    if(ecut < 0.0)
    {
      ecut = h5ecut;
      LOGMSG("  Plane-wave energy cutoff not specified.");
      LOGMSG("  Choosing maximum: " << ecut);
    }
  //Close the parameters group - back to root
  H5Gclose(grp_id);
  sprintf(GrpName,"eigenstates");
  hid_t es_grp_id = H5Gopen(hfile,GrpName);
  //Open the twist group...only 1 right now
  int twistindx = 0;
  sprintf(GrpName,"twist%d",twistindx);
  hid_t twist_grp_id = H5Gopen(es_grp_id,GrpName);
  //Read the twist angle
  TinyVector<double,3> twist;
  HDFAttribIO<TinyVector<double,3> > hdfobj_twist(twist);
  hdfobj_twist.read(twist_grp_id,"twist_angle");
  if(h5coefsreal)
  {
    //Test that twist = gamma if real coefficients
    for(int idim=0; idim<3; idim++)
      if(std::abs(twist[idim]) > 1.e-6)
      {
        LOGMSG("Error: real wavefunction with non-Gamma twist angle");
        OHMMS::Controller->abort();
      }
  }
#if !defined (QMC_COMPLEX)
  else
  {
    //Complex wavefunctions in file with real code
    LOGMSG("This binary doesn't support complex wavefunctions.");
    LOGMSG("Recompile with complex support (set in Configuration.h).");
    OHMMS::Controller->abort();
  }
#endif
  //Eventually need one of these sets per twist.
  //Create a container for the orbitals. This also creates the basis.
  PWOSet = new PlaneWaveOrbitalSet(&targetPtcl,nspins,nbands,twist);
  PlaneWaveBasis *CurBasis = PWOSet->BasisSet;
  //Read the planewave basisset.
  //Note that the same data is opened here for each twist angle-avoids duplication in the
  //h5 file (which may become very large).
  sprintf(GrpName,"basis");
  grp_id = H5Gopen(hfile,GrpName);
  CurBasis->readbasis(grp_id,ecut,nh5gvecs,targetPtcl.Lattice);
  H5Gclose(grp_id); //Close PW group
  //Read each eigenvector and add it to the set.
  //Resize the coefficient storage first
  //TODO: THIS ALL ASSUMES GROUNDSTATE FOR NOW.
  PWOSet->resizeCoefs();
  std::vector<ValueType> coefs;
  coefs.resize(nh5gvecs);
  HDFAttribIO<std::vector<ValueType> > hdfobj_coefs(coefs);
  //Go through all bands in inputfile.
  for(int ib=0; ib<nbands; ib++)
  {
    sprintf  (GrpName,"band%d",ib);
    hid_t band_grp_id = H5Gopen(twist_grp_id,GrpName);
    //Is this band occupied by either the up channel or the down?
    //If so, add it to the Coefs in PWOSet. Technicality: nup and ndown
    //may be different. This means that we need to be consistent when
    //filing the orbitals in Coefs so that the order matches ParticleSet
    //(i.e. 0..nup-1, nup..ndown+nup for upindx==0
    // or   0..ndown-1, ndown..nup+ndown for upindx==1) otherwise
    //low-energy bands may be used for high-energy electrons.
    //Add coefficient for up determinant...if occupied
    if(ib<nup)
    {
      //Open the correct spin channel in HDF5 file for up
      sprintf(GrpName,"spin%d",updetspinindex);
      hid_t spin_grp_id = H5Gopen(band_grp_id,GrpName);
      //When to file in coefs:
      int jorb = targetPtcl.first(upindx) + ib;
      //Read PW coefficients.
      hdfobj_coefs.read(spin_grp_id,"eigenvector");
      //File it
      PWOSet->addVector(coefs,jorb);
      //Close HDF group
      H5Gclose(spin_grp_id);
    }
    //Now add for down determinant...if occupied
    if(ib<ndown)
    {
      //Open the correct spin channel in HDF5 file for down
      sprintf(GrpName,"spin%d",downdetspinindex);
      hid_t spin_grp_id = H5Gopen(band_grp_id,GrpName);
      //When to file in coefs:
      int jorb = targetPtcl.first(1-upindx) + ib;
      //Read PW coefficients.
      hdfobj_coefs.read(spin_grp_id,"eigenvector");
      //File
      PWOSet->addVector(coefs,jorb);
      //Close HDF group
      H5Gclose(spin_grp_id);
    }
    H5Gclose(band_grp_id);
  } // ib
  H5Gclose(twist_grp_id);
  H5Gclose(es_grp_id);
}
}
