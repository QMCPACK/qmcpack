//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "NLPPClass.h"
#include "XMLWriterClass2.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <sstream>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>

void
ChannelPotentialClass::WriteChannelLog (XMLWriterClass &writer, 
					bool writeVl)
{
  string channels[] = {"s", "p", "d", "f", "g", "h", "i", "j"};
  // Use a logarithmic grid:
  // r(i) = a*(exp(b*i) - 1)
  const int numPoints = 2001;
  double end = Vl.grid.End();
  double step  = 0.00625;
  double scale = end/(expm1(step*(numPoints-1)));
  if (writeVl)
    writer.StartElement("vps");
  else
    writer.StartElement("pswf");
  writer.WriteAttribute("principal-n", n_principal);
  writer.WriteAttribute("l", channels[l]);
  writer.WriteAttribute("spin", -1);
  if (writeVl) {
    writer.WriteAttribute("cutoff", Cutoff);
    writer.WriteAttribute("occupation", Occupation);
  }
  writer.StartElement("radfunc");
  writer.StartElement("grid");
  writer.WriteAttribute("type", "log");
  writer.WriteAttribute("units", "bohr");
  writer.WriteAttribute("scale", scale, true);
  writer.WriteAttribute("step", step, true);
  writer.WriteAttribute("npts", numPoints);
  writer.EndElement(); // "grid"
  vector<double> data;
  for (int i=0; i<numPoints; i++) {
    double r  = scale * (expm1(step*i));
    if (writeVl) {
      double rmin = min (r, Vl.grid.End());
      data.push_back(Vlr(rmin));
    }
    else {
      r = min (r, ul.grid.End());
      data.push_back(ul(r));
    }
  }

  writer.WriteElement("data", data);
  writer.EndElement(); // "radfunc"
  writer.EndElement(); // "vps" or "pswf"
}



void
ChannelPotentialClass::WriteChannelLinear (XMLWriterClass &writer, 
					   double dr, double rmax, bool writeVl)
{
  string channels[] = {"s", "p", "d", "f", "g", "h", "i", "j"};
  // Use a linear grid
  // r(i) = dr * i
  const int numPoints = (int)round(rmax/dr) + 1;
  if (writeVl)
    writer.StartElement("vps");
  else
    writer.StartElement("pswf");
  writer.WriteAttribute("principal-n", n_principal);
  writer.WriteAttribute("l", channels[l]);
  writer.WriteAttribute("spin", -1);
  if (writeVl) {
    writer.WriteAttribute("cutoff", Cutoff);
    writer.WriteAttribute("occupation", Occupation);
  }
  writer.StartElement("radfunc");
  writer.StartElement("grid");
  writer.WriteAttribute("type", "linear");
  writer.WriteAttribute("units", "bohr");
  writer.WriteAttribute("ri", 0.0);
  writer.WriteAttribute("rf", dr*(double)(numPoints-1));
  writer.WriteAttribute("npts", numPoints);
  writer.EndElement(); // "grid"
  vector<double> data;
  for (int i=0; i<numPoints; i++) {
    double r  = dr*(double)i;
    if (writeVl) {
      r = min (r, ul.grid.End());
      if (r < Vl.grid.Start())
  	data.push_back(0.0);
      else
  	data.push_back(Vlr(r));
    }
    else {
      r = min (r, ul.grid.End());
      data.push_back(ul(r));
    }
  }

  writer.WriteElement("data", data);
  writer.EndElement(); // "radfunc"
  writer.EndElement(); // "vps" or "pswf"
}


void
PseudoClass::SetupMaps()
{
  UnitToHartreeMap[string("hartree")] = 1.0;
  UnitToHartreeMap[string("hartrees")] = 1.0;
  UnitToHartreeMap[string("Hartrees")] = 1.0;
  UnitToHartreeMap[string("rydberg")] = 0.5;
  UnitToHartreeMap[string("ev")]      = 0.03674932595264097934;
  
  UnitToBohrMap[string("bohr")]     = 1.0;
  UnitToBohrMap[string("atomic")]   = 1.0;
  UnitToBohrMap[string("angstrom")] = 1.8897261;

  XCMap[XC_LDA] ="LDA" ; XCRevMap["LDA"] =XC_LDA;  XCRevMap["lda"] =XC_LDA;
  XCMap[XC_GGA] ="GGA" ; XCRevMap["GGA"] =XC_GGA;  XCRevMap["gga"] =XC_LDA;
  XCMap[XC_HF]  ="HF"  ; XCRevMap["HF"]  =XC_HF;   XCRevMap["hf"]  =XC_LDA;
  XCMap[XC_DF]  ="DF"  ; XCRevMap["DF"]  =XC_DF;   XCRevMap["df"]  =XC_DF;
  XCMap[XC_NONE]="NONE"; XCRevMap["NONE"]=XC_NONE; XCRevMap["none"]=XC_NONE;
  
  ChannelMap[0] = "s";  ChannelRevMap["s"] = 0;
  ChannelMap[1] = "p";  ChannelRevMap["p"] = 1;
  ChannelMap[2] = "d";  ChannelRevMap["d"] = 2;
  ChannelMap[3] = "f";  ChannelRevMap["f"] = 3;
  ChannelMap[4] = "g";  ChannelRevMap["g"] = 4;
  ChannelMap[5] = "h";  ChannelRevMap["h"] = 5;
  ChannelMap[6] = "i";  ChannelRevMap["i"] = 6;
  ZToSymbolMap[1]   = "H";  ZToSymbolMap[2]   = "He";
  ZToSymbolMap[3]   = "Li"; ZToSymbolMap[4]   = "Be";
  ZToSymbolMap[5]   = "B";  ZToSymbolMap[6]   = "C";
  ZToSymbolMap[7]   = "N";  ZToSymbolMap[8]   = "O";
  ZToSymbolMap[9]   = "F";  ZToSymbolMap[10]  = "Ne";
  ZToSymbolMap[11]  = "Na"; ZToSymbolMap[12]  = "Mg";
  ZToSymbolMap[13]  = "Al"; ZToSymbolMap[14]  = "Si";
  ZToSymbolMap[15]  = "P";  ZToSymbolMap[16]  = "S";
  ZToSymbolMap[17]  = "Cl"; ZToSymbolMap[18]  = "Ar";
  ZToSymbolMap[19]  = "K";  ZToSymbolMap[20]  = "Ca";
  ZToSymbolMap[21]  = "Sc"; ZToSymbolMap[22]  = "Ti";
  ZToSymbolMap[23]  = "V";  ZToSymbolMap[24]  = "Cr";
  ZToSymbolMap[25]  = "Mn"; ZToSymbolMap[26]  = "Fe";
  ZToSymbolMap[27]  = "Co"; ZToSymbolMap[28]  = "Ni";
  ZToSymbolMap[29]  = "Cu"; ZToSymbolMap[30]  = "Zn";
  ZToSymbolMap[31]  = "Ga"; ZToSymbolMap[32]  = "Ge";
  ZToSymbolMap[33]  = "As"; ZToSymbolMap[34]  = "Se";
  ZToSymbolMap[35]  = "Br"; ZToSymbolMap[36]  = "Kr";
  ZToSymbolMap[37]  = "Rb"; ZToSymbolMap[38]  = "Sr";
  ZToSymbolMap[39]  = "Y";  ZToSymbolMap[40]  = "Zr";
  ZToSymbolMap[41]  = "Nb"; ZToSymbolMap[42]  = "Mo";
  ZToSymbolMap[43]  = "Tc"; ZToSymbolMap[44]  = "Ru";
  ZToSymbolMap[45]  = "Rh"; ZToSymbolMap[46]  = "Pd";
  ZToSymbolMap[47]  = "Ag"; ZToSymbolMap[48]  = "Cd";
  ZToSymbolMap[49]  = "In"; ZToSymbolMap[50]  = "Sn";
  ZToSymbolMap[51]  = "Sb"; ZToSymbolMap[52]  = "Te";
  ZToSymbolMap[53]  = "I";  ZToSymbolMap[54]  = "Xe";
  ZToSymbolMap[55]  = "Cs"; ZToSymbolMap[56]  = "Ba";
  ZToSymbolMap[57]  = "La"; ZToSymbolMap[58]  = "Ce";
  ZToSymbolMap[59]  = "Pr"; ZToSymbolMap[60]  = "Nd";
  ZToSymbolMap[61]  = "Pm"; ZToSymbolMap[62]  = "Sm";
  ZToSymbolMap[63]  = "Eu"; ZToSymbolMap[64]  = "Gd";
  ZToSymbolMap[65]  = "Tb"; ZToSymbolMap[66]  = "Dy";
  ZToSymbolMap[67]  = "Ho"; ZToSymbolMap[68]  = "Er";
  ZToSymbolMap[69]  = "Tm"; ZToSymbolMap[70]  = "Yb";
  ZToSymbolMap[71]  = "Lu"; ZToSymbolMap[72]  = "Hf";
  ZToSymbolMap[73]  = "Ta"; ZToSymbolMap[74]  = "W";
  ZToSymbolMap[75]  = "Re"; ZToSymbolMap[76]  = "Os";
  ZToSymbolMap[77]  = "Ir"; ZToSymbolMap[78]  = "Pt";
  ZToSymbolMap[79]  = "Au"; ZToSymbolMap[80]  = "Hg";
  ZToSymbolMap[81]  = "Tl"; ZToSymbolMap[82]  = "Pb";
  ZToSymbolMap[83]  = "Bi"; ZToSymbolMap[84]  = "Po";
  ZToSymbolMap[85]  = "At"; ZToSymbolMap[86]  = "Rn";
  ZToSymbolMap[87]  = "Fr"; ZToSymbolMap[88]  = "Ra";
  ZToSymbolMap[89]  = "Ac"; ZToSymbolMap[90]  = "Th";
  ZToSymbolMap[91]  = "Pa"; ZToSymbolMap[92]  = "U";
  ZToSymbolMap[93]  = "Np"; ZToSymbolMap[94]  = "Pu";
  ZToSymbolMap[95]  = "Am"; ZToSymbolMap[96]  = "Cm";
  ZToSymbolMap[97]  = "Bk"; ZToSymbolMap[98]  = "Cf";
  ZToSymbolMap[99]  = "Es"; ZToSymbolMap[100] = "Fm";
  ZToSymbolMap[101] = "Mc"; ZToSymbolMap[102] = "No";
  ZToSymbolMap[103] = "Lw"; 

  ZToMassMap[  1] = 1.00794;
  ZToMassMap[  2] = 4.002602;
  ZToMassMap[  3] = 6.941;
  ZToMassMap[  4] = 9.012182;
  ZToMassMap[  5] = 10.811;
  ZToMassMap[  6] = 12.0107;
  ZToMassMap[  7] = 14.0067;
  ZToMassMap[  8] = 15.9994;
  ZToMassMap[  9] = 18.9984032;
  ZToMassMap[ 10] = 20.1797;
  ZToMassMap[ 11] = 22.98976928;
  ZToMassMap[ 12] = 24.3050;
  ZToMassMap[ 13] = 26.9815386;
  ZToMassMap[ 14] = 28.0855;
  ZToMassMap[ 15] = 30.973762;
  ZToMassMap[ 16] = 32.065;
  ZToMassMap[ 17] = 35.453;
  ZToMassMap[ 18] = 39.948;
  ZToMassMap[ 19] = 39.0983;
  ZToMassMap[ 20] = 40.078;
  ZToMassMap[ 21] = 44.955912;
  ZToMassMap[ 22] = 47.867;
  ZToMassMap[ 23] = 50.9415;
  ZToMassMap[ 24] = 51.9961;
  ZToMassMap[ 25] = 54.938049;
  ZToMassMap[ 26] = 55.845;
  ZToMassMap[ 27] = 58.933200;
  ZToMassMap[ 28] = 58.6934;
  ZToMassMap[ 29] = 63.546;
  ZToMassMap[ 30] = 65.39;
  ZToMassMap[ 31] = 69.723;
  ZToMassMap[ 32] = 72.61;
  ZToMassMap[ 33] = 74.92160;
  ZToMassMap[ 34] = 78.96 ;
  ZToMassMap[ 35] = 79.904;
  ZToMassMap[ 36] = 83.80;
  ZToMassMap[ 37] = 85.4678;
  ZToMassMap[ 38] = 87.62;
  ZToMassMap[ 39] = 88.90585;
  ZToMassMap[ 40] = 91.224;
  ZToMassMap[ 41] = 92.90638;
  ZToMassMap[ 42] = 95.94;
  ZToMassMap[ 43] = 98;
  ZToMassMap[ 44] = 101.07;
  ZToMassMap[ 45] = 102.90550;
  ZToMassMap[ 46] = 106.42;
  ZToMassMap[ 47] = 107.8682;
  ZToMassMap[ 48] = 112.411;
  ZToMassMap[ 49] = 114.818;
  ZToMassMap[ 50] = 118.710;
  ZToMassMap[ 51] = 121.760;
  ZToMassMap[ 52] = 127.60;
  ZToMassMap[ 53] = 126.90447;
  ZToMassMap[ 54] = 131.29;
  ZToMassMap[ 55] = 132.90545;
  ZToMassMap[ 56] = 137.327;
  ZToMassMap[ 57] = 138.9055;
  ZToMassMap[ 58] = 140.116;
  ZToMassMap[ 59] = 140.90765;
  ZToMassMap[ 60] = 144.24;
  ZToMassMap[ 61] = 145;
  ZToMassMap[ 62] = 150.36;
  ZToMassMap[ 63] = 151.964;
  ZToMassMap[ 64] = 157.25;
  ZToMassMap[ 65] = 158.92534;
  ZToMassMap[ 66] = 162.50;
  ZToMassMap[ 67] = 164.93032;
  ZToMassMap[ 68] = 167.26;
  ZToMassMap[ 69] = 168.93421;
  ZToMassMap[ 70] = 173.04;
  ZToMassMap[ 71] = 174.967;
  ZToMassMap[ 72] = 178.49;
  ZToMassMap[ 73] = 180.9479;
  ZToMassMap[ 74] = 183.84;
  ZToMassMap[ 75] = 186.207;
  ZToMassMap[ 76] = 190.23;
  ZToMassMap[ 77] = 192.217;
  ZToMassMap[ 78] = 195.078;
  ZToMassMap[ 79] = 196.96655;
  ZToMassMap[ 80] = 200.59;
  ZToMassMap[ 81] = 204.3833;
  ZToMassMap[ 82] = 207.2;
  ZToMassMap[ 83] = 208.98038;
  ZToMassMap[ 84] = 209;
  ZToMassMap[ 85] = 210;
  ZToMassMap[ 86] = 222;
  ZToMassMap[ 87] = 223;
  ZToMassMap[ 88] = 226;
  ZToMassMap[ 89] = 227;
  ZToMassMap[ 90] = 232.0381;
  ZToMassMap[ 91] = 231.03588;
  ZToMassMap[ 92] = 238.0289;
  ZToMassMap[ 93] = 237;
  ZToMassMap[ 94] = 244;
  ZToMassMap[ 95] = 243;
  ZToMassMap[ 96] = 247;
  ZToMassMap[ 97] = 247;
  ZToMassMap[ 98] = 251;
  ZToMassMap[ 99] = 252;
  ZToMassMap[100] = 257;
  ZToMassMap[101] = 258;
  ZToMassMap[102] = 259;
  ZToMassMap[103] = 262;
  ZToMassMap[104] = 261;
  ZToMassMap[105] = 262;
  ZToMassMap[106] = 263;
  ZToMassMap[107] = 262;
  ZToMassMap[108] = 265;
  ZToMassMap[109] = 266;
  ZToMassMap[110] = 269;
  ZToMassMap[111] = 272;
  ZToMassMap[112] = 277;
  ZToMassMap[113] = 284;
  ZToMassMap[114] = 289;
  ZToMassMap[115] = 288;
  ZToMassMap[116] = 293;
  ZToMassMap[117] = 293; // UNKNOWN
  ZToMassMap[118] = 294;

  for (int i=1; i<103; i++)
    SymbolToZMap[ZToSymbolMap[i]] = i;
}

void
PseudoClass::WriteFPMD(string fname)
{
  const int numPoints = 2000;
  const double dr = 0.005;

  XMLWriterClass writer;
  writer.StartDocument(fname, "1.0", "UTF-8");
  writer.StartElement("fpmd:species");
  writer.WriteAttribute("xmlns:fpmd", 
			"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0");
  writer.WriteAttribute("xmlns:xsi", 
			"http://www.w3.org/2001/XMLSchema-instance");
  writer.WriteAttribute
    ("xsi:schemaLocation",
     "http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 species.xsd");

  writer.StartElement("description");
  writer.WriteData("ppconvert");
  writer.EndElement();
  
  writer.StartElement("symbol");
  writer.WriteData(ZToSymbolMap[AtomicNumber]);
  writer.EndElement(); // "symbol"

  writer.StartElement("atomic_number");
  writer.WriteData(AtomicNumber);
  writer.EndElement(); // "atomic_number"

  writer.StartElement("mass");
  writer.WriteData(ZToMassMap[AtomicNumber]);
  writer.EndElement(); // "mass"

  writer.StartElement("norm_conserving_pseudopotential");

  writer.StartElement("valence_charge");
  writer.WriteData((int)round(PseudoCharge));
  writer.EndElement(); // "valence_charge"

  writer.StartElement("lmax");
  writer.WriteData(ChannelPotentials.size()-1);
  writer.EndElement(); // "lmax"

  writer.StartElement("llocal");
  writer.WriteData(LocalChannel);
  writer.EndElement(); // "llocal"

  writer.StartElement("nquad");
  writer.WriteData(0);
  writer.EndElement();

  writer.StartElement("rquad");
  writer.WriteData(3.0);
  writer.EndElement(); // "rquad"



  writer.StartElement("mesh_spacing");
  writer.WriteData(dr);
  writer.EndElement(); // "rquad"

  for (int il=0; il<ChannelPotentials.size(); il++) {
    ChannelPotentialClass &pot = ChannelPotentials[il];
    writer.StartElement("projector");
    writer.WriteAttribute("l", pot.l);
    writer.WriteAttribute("size", numPoints);
    
    writer.StartElement("radial_potential");
    vector<double> vl(numPoints);
    for (int ir=0; ir<numPoints; ir++) {
      double r = (double)ir*dr;
      vl[ir] = (ir==0) ? pot.Vl(0.0) : pot.Vlr(r)/r;
    }
    writer.WriteData(vl);
    writer.EndElement(); // "radial_potential"

    if (pot.HasProjector && pot.l != LocalChannel) {
      writer.StartElement("radial_function");
      vector<double> ul(numPoints);
      for (int ir=0; ir<numPoints; ir++) {
	double r = (ir==0) ? 1.0e-8 : (double)ir*dr;
	ul[ir] = pot.ul(r)/r;
      }
      ul[0] = 2.0*ul[1] - ul[2];
      writer.WriteData(ul);
      writer.EndElement(); // "radial_function"
    }

    writer.EndElement(); // "projector"
  }



  writer.EndElement(); // "norm_conserving_pseudopotential"
  writer.EndElement(); // "fpmd:species"
  writer.EndDocument();
}

void
PseudoClass::WriteCASINO(string fname)
{
  const int numPoints = 10001;
  const double rMax = 10.0;

  FILE *fout = fopen (fname.c_str(), "w");
  if (!fout) {
    cerr << "Error trying to write " << fname << ".\n";
    abort();
  }

  fprintf (fout, "Pseudopotential in real space for %s\n",
	   ZToSymbolMap[AtomicNumber].c_str());
  fprintf (fout, "Atomic number and pseudo-charge\n%d %1.2f\n",
	   AtomicNumber, PseudoCharge);
  fprintf (fout, "Energy units (rydberg/hartree/ev):\nrydberg\n");
  fprintf (fout, "Angular momentum of local component (0=s,1=p,2=d..)\n%d\n",
	   LocalChannel);
  fprintf (fout, "NLRULE override (1) VMC/DMC (2) config gen "
	   "(0 ==> input/default value)\n0 0\n");
  fprintf (fout, "Number of grid points\n%d\n",
	   numPoints);
  fprintf (fout, "R(i) in atomic units\n");
  for (int i=0; i<numPoints; i++) {
    double r = (double) i/(double)(numPoints-1)*rMax;
    fprintf (fout, "  %22.16e\n", r);
  }
  
  for (int l=0; l<ChannelPotentials.size(); l++) {
    ChannelPotentialClass &pot = (l < ChannelPotentials.size()) ?
      ChannelPotentials[l] : ChannelPotentials[LocalChannel];
    if (l < ChannelPotentials.size() && l != pot.l) {
      cerr << "Channel mismatch in WriteCASINO.\n";
      abort();
    }
    fprintf (fout, "r*potential (L=%d) in Ry\n",
	     l);
    for (int i=0; i<numPoints; i++) {
      double r = (double) i/(double)(numPoints-1)*rMax;
      fprintf (fout, "  %22.16e\n", 
	       (i==0) ? 0.0 : 2.0*pot.Vlr(r));
    }
  }
  fclose(fout);
}

void
PseudoClass::WriteXML(string fname)
{
//   int rc;
//   xmlTextWriterPtr writer = xmlNewTextWriterFilename (fname.c_str(), 0);

//   rc = xmlTextWriterStartDocument (writer, NULL, "UTF-8", NULL);
//   // Start psuedo
//   rc = xmlTextWriterStartElement(writer, (xmlChar*)"pseudo");
//   rc = xmlTextWriterWriteAttribute
//     (writer, (xmlChar*) "version", (xmlChar*)"0.5");
//   rc = xmlTextWriterEndElement(writer); // "pseudo"
//   rc = xmlTextWriterEndDocument (writer);
//   xmlFreeTextWriter(writer);
  for (int l=0; l<ChannelPotentials.size(); l++)
    if (!ChannelPotentials[l].HasProjector) {
      cerr << "Missing projector for l=" << l 
	   << ".  Aborting." << endl;
      exit(1);
    }
  XMLWriterClass writer;
  writer.StartDocument(fname, "1.0", "UTF-8");
  writer.StartElement("pseudo");
  writer.WriteAttribute ("version", "0.5");
  writer.StartElement("header");
  writer.WriteAttribute("symbol", ZToSymbolMap[AtomicNumber]);
  writer.WriteAttribute("atomic-number", (double)AtomicNumber);
  writer.WriteAttribute("zval", PseudoCharge);
  writer.WriteAttribute("relativistic", Relativistic ? "yes" : "no");
  writer.WriteAttribute("polarized", "no");
  writer.WriteAttribute("creator", "ppconvert");
  writer.WriteAttribute("flavor", "Troullier-Martins");

  writer.WriteAttribute("core-corrections", "no");
  //  writer.WriteAttribute("xc-functional-type", XCMap[XC]);
  writer.WriteAttribute("xc-functional-type", "GGA");
  writer.WriteAttribute("xc-functional-parametrization", "Perdew-Burke-Ernzerhof");
  writer.EndElement(); // "header"

  double rmax;
  // Write the grid information:

  // NEVER USE LOG GRID
  // Leads to inaccuracy in qmcPACK
  if (WriteLogGrid && false) {
    const int numPoints = 2001;
    rmax = PotentialGrid.End();
    double step  = 0.00625;
    double scale = rmax/(expm1(step*(numPoints-1)));
    writer.StartElement("grid");
    writer.WriteAttribute("type", "log");
    writer.WriteAttribute("units", "bohr");
    writer.WriteAttribute("scale", scale, true);
    writer.WriteAttribute("step", step, true);
    writer.WriteAttribute("npts", numPoints-1);
    writer.EndElement(); // "grid"
  }
  else {
    rmax = min (10.0, PotentialGrid.End());
    int numPoints = (int)round(rmax/grid_delta) + 1;
    writer.StartElement("grid");
    writer.WriteAttribute("type", "linear");
    writer.WriteAttribute("units", "bohr");
    writer.WriteAttribute("ri", 0.0);
    writer.WriteAttribute("rf", grid_delta*(double)(numPoints-1), true);
    writer.WriteAttribute("npts", numPoints);
    writer.EndElement(); // "grid"
  }

  writer.StartElement("semilocal"); 
  writer.WriteAttribute("units", "hartree");
  writer.WriteAttribute("format", "r*V");
  writer.WriteAttribute("npots-down", (int)ChannelPotentials.size());
  writer.WriteAttribute("npots-up"  , 0);
  writer.WriteAttribute("l-local", LocalChannel);
  for (int l=0; l<ChannelPotentials.size(); l++)
    if (WriteLogGrid && false)
      ChannelPotentials[l].WriteChannelLog (writer, true);
    else
      ChannelPotentials[l].WriteChannelLinear (writer, grid_delta, rmax, true);
  writer.EndElement(); // "semilocal"

  writer.StartElement("pseudowave-functions");
  writer.WriteAttribute("units", "electrons/bohr^(-3/2)");
  writer.WriteAttribute("format", "u_n,l (r) = 1/r R_n,l (r)");
  writer.WriteAttribute("n-pseudowave-functions-down", 
			(int)ChannelPotentials.size());
  writer.WriteAttribute("n-pseudowave-functions-up", 0);
  for (int l=0; l<ChannelPotentials.size(); l++) 
    if (WriteLogGrid && false)
      ChannelPotentials[l].WriteChannelLog (writer, false);
    else
      ChannelPotentials[l].WriteChannelLinear (writer, grid_delta, rmax, false);
  writer.EndElement(); // "pseudowave-functions"

  writer.EndElement(); // "pseudo"
  writer.EndDocument();
}


void
PseudoClass::Write4Block (FILE *fout, vector<double>& data, int indent)
{
  int N = data.size();
  for (int i=0; i<N/4; i++) {
    for (int j=0; j<indent; j++)   fprintf (fout, " ");
    fprintf (fout, "%17.11e %1.17e %1.17e %1.17e\n", 
	     data[4*i+0], data[4*i+1], data[4*i+2], data[4*i+3]);
  }
  int last = N/4;
  for (int i=0; i<indent; i++)   fprintf (fout, " ");
  for (int j=0; j<(N % 4); j++)
    fprintf (fout, "%17.11e ", data[4*last+j]);
  fprintf (fout, "\n");
}


void
PseudoClass::WriteUPF (string fileName)
{
  cerr << "LocalChannel = " << LocalChannel << endl;
  FILE *fout = fopen (fileName.c_str(), "w");
  if (fout == NULL) {
    cerr << "Error writing UPF file " << fileName << endl;
    return;
  }

  int N = PotentialGrid.NumPoints();

  // Write PP_INFO section
  fprintf (fout, "<PP_INFO>\n");
  fprintf (fout, "Generated using ppconvert\n");
  // Get date
  time_t now = time(NULL);
  struct tm &ltime = *(localtime (&now));
  stringstream date;
  int year = ltime.tm_year % 100;
  int mon  = ltime.tm_mon + 1;
  int day  = ltime.tm_mday;
  date << ((mon<10)  ? "0" : "") << mon << "/";
  date << ((day<10)  ? "0" : "") << day << "/";
  date << ((year<10) ? "0" : "") << year;
  fprintf (fout, "Author:  Kenneth P. Esler  Generation date: %s\n",
	   date.str().c_str());

	   
  fprintf (fout, "</PP_INFO>\n");

  // Write PP_HEADER section
  fprintf (fout, "<PP_HEADER>\n");
  fprintf (fout, "   0                  Version Number\n");
  fprintf (fout, "  %2s                  Element\n", 
	   ZToSymbolMap[AtomicNumber].c_str());
  fprintf (fout, "  NC                  Norm - Conserving pseudopotential\n");
  fprintf (fout, "   F                  Nonlinear Core Correction\n");
  string F_EX, F_CORR, F_XC_GRAD, F_CORR_GRAD;
  F_EX = "SLA";
  F_CORR = "PW";
  F_XC_GRAD = "PBE";
  F_CORR_GRAD = "PBE"; 
  fprintf (fout, " %4s %4s %4s %4s  Exchange-Corelation functional\n",
	   F_EX.c_str(), F_CORR.c_str(), F_XC_GRAD.c_str(), F_CORR_GRAD.c_str());
  fprintf (fout, "  %1.10f        Z valence\n",
	   PseudoCharge);
  fprintf (fout, "  %1.10f       Total energy\n", TotalEnergy);
  fprintf (fout, "  %3.6f  %3.6f  Suggested cutoff for wfc and rho\n",
	   0.0, 0.0);
  fprintf (fout, "  %d                   Max angular momentum component\n",
	   ChannelPotentials.size()-1);
  fprintf (fout, "  %4d                Number of points in mesh\n",
	   N);
  fprintf (fout, "  %d   %d               Number of Wavefunctions, Number of Projectors\n", ChannelPotentials.size(), ChannelPotentials.size()-1);
  fprintf (fout, " WaveFunctions      nl   l   occ    \n");
  for (int l=0; l<ChannelPotentials.size(); l++) 
    fprintf (fout, "                    %d%s   %d   %1.4f\n",
	     ChannelPotentials[l].n_principal,
	     ChannelMap[l].c_str(),
	     ChannelPotentials[l].l,
	     ChannelPotentials[l].Occupation);
    
  fprintf (fout, "</PP_HEADER>\n");

  // Write PP_MESH section
  fprintf (fout, "<PP_MESH>\n");
  fprintf (fout, "  <PP_R>\n");
  Write4Block (fout, PotentialGrid.Points());
  fprintf (fout, "  </PP_R>\n");
  fprintf (fout, "  <PP_RAB>\n");
  vector<double> dr(N);
  for (int i=1; i<(N-1); i++)
    dr[i] = 0.5*(PotentialGrid[i+1]-PotentialGrid[i-1]);
  dr[0]   = 2.0*dr[ 1 ]-dr[ 2 ];
  dr[N-1] = 2.0*dr[N-2]-dr[N-3];
  
  Write4Block (fout, dr);
  fprintf (fout, "  </PP_RAB>\n");
  fprintf (fout, "</PP_MESH>\n");

  // Write PP_LOCAL section
  fprintf (fout, "<PP_LOCAL>\n");
  vector<double> vloc(N);
  for (int i=0; i<N; i++)
    vloc[i] = 2.0*ChannelPotentials[LocalChannel].Vl(i); // /PotentialGrid[i];
  //vloc[0] = 2.0*ChannelPotentials[LocalChannel].Vl(1.0e-7)/1.0e-7;
  Write4Block (fout, vloc);
  fprintf (fout, "</PP_LOCAL>\n");
  
  // Write the nonlocal section
  fprintf (fout, "<PP_NONLOCAL>\n");
  int betaNum = 1;
  for (int l=0; l<ChannelPotentials.size(); l++)
    if (l != LocalChannel) {
      vector<double> beta(N);
      for (int i=0; i<N; i++) {
	double Vloc = ChannelPotentials[LocalChannel].Vlr(i);
	double Vl   = ChannelPotentials[l].Vlr(i);
	double ul   = ChannelPotentials[l].ul(i);
	double r    = PotentialGrid[i];
	beta[i] = 2.0*ul*(Vl - Vloc)/r;
	//beta[i] = 2.0*ul*(Vl - Vloc);
	//beta[i] = 2.0*r*ul*(Vl - Vloc);
      }
      beta[0] = 2.0*ChannelPotentials[l].ul(0.0)*
 	(ChannelPotentials[l].Vlr(1.0e-7)-ChannelPotentials[LocalChannel].Vlr(1.0e-7))/1.0e-7;
//       beta[0] = 2.0*ChannelPotentials[l].ul(0.0)*
// 	(ChannelPotentials[l].Vl(1.0e-7)-ChannelPotentials[LocalChannel].Vl(1.0e-7));
      fprintf (fout, "  <PP_BETA>\n");
      fprintf (fout, "    %d %d             Beta L\n", betaNum, l);
      fprintf (fout, "    %d  \n", N);
      Write4Block (fout, beta, 4);
      fprintf (fout, "  </PP_BETA>\n");
      betaNum++;
    }
  // Write D_ij matrix
  fprintf (fout, "  <PP_DIJ>\n");
  fprintf (fout, "    %d         Number of nonzero Dij\n",
	   ChannelPotentials.size()-1);
  betaNum = 1;
  // Compute D_ij matrix.  In our case, with single projectors, D_ij
  // is diagonal.  The elements are
  // D_ll = 1.0/<u_l |Vl - Vloc | u_l>
  // We will numerically integrate with Simpson's rule
  for (int l=0; l<ChannelPotentials.size(); l++) 
    if (l != LocalChannel) {
      const int M = 100000;
      double DijInv = 0.0;
      double norm = 0.0;
      double dr = PotentialGrid[N-1]/(double)(M);
      double third = 1.0/3.0;
      double prefactor;
      double r = 1.0e-8;
      for (int i=0; i<=M; i++) {
	double Vl   = ChannelPotentials[l].Vlr(r)/r;
	double Vloc = ChannelPotentials[LocalChannel].Vlr(r)/r;
	double u = ChannelPotentials[l].ul(r);
	// Prefactors for Simpsons's rule
	if ((i==0) || (i==M))
	  prefactor = third;
	else if ((i&1)==1)
	  prefactor = 4.0* third;
	else
	  prefactor = 2.0*third;
	DijInv += prefactor*dr *u * u * 2.0 * (Vl - Vloc);
	norm += prefactor*u*u*dr;
	r += dr;
      }
      double Dij = 1.0/DijInv;
      fprintf (fout, "    %d %d  %1.20e\n", betaNum, betaNum, Dij);
      fprintf (stderr, "l = %d, betaNum = %d, norm = %1.10e  Dij = %1.10e\n", 
	       l, betaNum, norm, Dij);
      betaNum++;
    }
  fprintf (fout, "  </PP_DIJ>\n");
  fprintf (fout, "</PP_NONLOCAL>\n");
  
  // Write PP_PSWFC section
  fprintf (fout, "<PP_PSWFC>\n");
  for (int l=0; l<ChannelPotentials.size(); l++) {
    ChannelPotentialClass &pot = ChannelPotentials[l];
    fprintf (fout, "%d%s    %d  %6.4f          Wavefunction\n",
	     pot.n_principal, ChannelMap[l].c_str(), l, pot.Occupation);
    vector<double> u(N);
    for (int i=0; i<N; i++)
      u[i] = pot.ul(i);
    Write4Block (fout, u);
  }
  fprintf (fout, "</PP_PSWFC>\n");

  // Write PP_RHOATOM section
  fprintf (fout, "<PP_RHOATOM>\n");
  vector<double> rho_at(N);
  for (int i=0; i<N; i++) {
    if (RhoAtom.Initialized)
      rho_at[i] = RhoAtom(i)*PotentialGrid[i]*PotentialGrid[i];
    else {
      rho_at[i] = 0.0;
      for (int l=0; l<ChannelPotentials.size(); l++) {
	ChannelPotentialClass &pot = ChannelPotentials[l];
	double u = pot.ul(i);
	rho_at[i] += pot.Occupation * u*u;
      }
    }
  }
  Write4Block (fout, rho_at);
  
  fprintf (fout, "</PP_RHOATOM>\n");

  // Close the file
  fclose (fout);
}


void
PseudoClass::WriteFHI (string fileName)
{
  const double amesh =  1.013084867359809;
  const int numPoints = 1118;


  FILE *fout = fopen (fileName.c_str(), "w");
  assert (fout != NULL);
  // Write ABINIT header
  fprintf (fout, "Written by ppconvert:  %s psuedopotential\n",
	   ZToSymbolMap[AtomicNumber].c_str());
  time_t now = time(NULL);
  struct tm &ltime = *(localtime (&now));
  stringstream date;
  int year = ltime.tm_year % 100;
  int mon  = ltime.tm_mon + 1;
  int day  = ltime.tm_mday;
  date << ((year<10) ? "0" : "") << year;
  date << ((mon<10)  ? "0" : "") << mon;
  date << ((day<10)  ? "0" : "") << day;

  fprintf (fout, "%9.5f %9.5f %s"
	   "                    zatom, zion, pspdat\n", (double)AtomicNumber,
	   PseudoCharge, date.str().c_str());

  fprintf (fout, "    6   7  %d  %d  %d  0.0000     "
	   "pspcod,pspxc,lmax,lloc,mmax,r2well\n",
	   ChannelPotentials.size()-1, LocalChannel, numPoints);
  fprintf (fout, "   0.00000   0.00000   0.00000                  "
	   "rchrg,fchrg,qchrg\n");
  fprintf (fout, "5 --- reserved for future features\n"
	   "6 --- reserved for future features\n"
	   "7 --- Here follows the cpi file in the fhi98pp format -\n");

  // Write FHI file proper
  fprintf (fout, "%d %d\n", (int)round(PseudoCharge), 
	   ChannelPotentials.size());
  for (int i=0; i<10; i++)
    fprintf (fout, "0.0 0.0 0.0\n");

  for (int l=0; l<ChannelPotentials.size(); l++) {
    fprintf (fout, "%d %1.18f\n", numPoints, amesh);
    double r = 5.848035476425733e-05;
    double nrm = 0.0;
    for (int i=0; i<numPoints; i++) {
      double Vl = ChannelPotentials[l].Vlr(r)/r;
      double ul = ChannelPotentials[l].ul(r);
      nrm += 0.5*r*(amesh-1.0/amesh)*ul*ul;
    
      fprintf (fout, "%5d  %20.16e %20.16e %20.16e\n",
	       i+1, r, ul, Vl);
      r *= amesh;
    }
    fprintf (stderr, "Norm = %1.8f\n", nrm);
  }
  fclose (fout);  
}


void
PseudoClass::WriteABINIT (string fileName)
{
  if (fileName == "") 
    fileName = ZToSymbolMap[AtomicNumber] + ".psp";

  FILE *fout = fopen (fileName.c_str(), "w");
  assert (fout != NULL);

  time_t now = time(NULL);
  struct tm &ltime = *(localtime (&now));
  stringstream date;
  int year = ltime.tm_year % 100;
  int mon  = ltime.tm_mon + 1;
  int day  = ltime.tm_mday;
  date << ((year<10) ? "0" : "") << year;
  date << ((mon<10)  ? "0" : "") << mon;
  date << ((day<10)  ? "0" : "")  << day;

  fprintf (fout, "%s  %s", ZToSymbolMap[AtomicNumber].c_str(),
	   asctime(&ltime));
  fprintf (fout, "%9.5f %9.5f %s"
	   "                    zatom, zion, pspdat\n", (double)AtomicNumber,
	   PseudoCharge, date.str().c_str());
  const int pspcod = 1;
  const int pspxc  = 0;
  const int lmax   = ChannelPotentials.size()-1;
  // HACK HACK HAK
  const int lloc   = LocalChannel;
  const int mmax  = 2001;
  double r2well = 0.0;
  // The following are informational only
  const double e99    = 0.0;
  const double e999   = 0.0;
  const double rms    = 0.0;
  const double ekb1   = 0.0;
  const double ekb2   = 0.0;
  const double epsatm = 0.0;
  const double fchrg  = 0.0;
  const double rchrg  = 1.0;
  const double qchrg  = 0.0;

  fprintf (fout, "%d   %d   %d   %d   %d   %8.5f"
	   "              pspcod,pspxc,lmax,lloc,mmax,r2well\n",
	   pspcod, pspxc, lmax, lloc, mmax, r2well);

  for (int l=0; l<ChannelPotentials.size(); l++) {
    // The local channel needs no projector
    int nproj = (l == lloc) ? 0 : 1;
    fprintf (fout, "%d %5.8f %5.8f %d %5.8f"
	     "          l,e99.0,e99.9,nproj,rcpsp\n", 
	     l, e99, e999, nproj, ChannelPotentials[l].Cutoff);
    fprintf (fout, "%14.10f %14.10f %14.10f %14.10f"
	     "  rms,ekb1,ekb2,epsatm\n", rms, ekb1, ekb2, epsatm);
  }
  fprintf (fout, "%8.5f %8.5f %8.5f"
	   "                    rchrg,fchrg,qchrg\n", rchrg, fchrg, qchrg);

  // Write out potentials
  for (int l=0; l<ChannelPotentials.size(); l++) {
    fprintf (fout, "%d = l for CASINO pseudopotential\n", l);
    for (int i=0; i<mmax; i++) {
      double x = (double)i/(double)(mmax-1);
      x += 0.01;
      double r = 100.0*x*x*x*x*x - 1.0e-8;
      double Vl;
      if (r >= ChannelPotentials[l].Vl.grid.End()) 
	Vl = -1.0*PseudoCharge/r;
      else
	Vl = ChannelPotentials[l].Vl(r); // /r;
      double Vlold = Vl;
      Vl = ChannelPotentials[l].Vlr(r)/r;
      if (i == 0) 
	Vl = ChannelPotentials[l].Vl(0.0);
      // fprintf (stderr, "r = %1.10f  Interp diff = %1.8e\n", r, Vl-Vlold);
      fprintf (fout, "%23.16e ", Vl);
      if ((i%3)==2)
	fprintf (fout, "\n");
    }
  }
  // Now write out radial wave functions
  for (int l=0; l<ChannelPotentials.size(); l++) {
      fprintf (fout, "%d = l for CASINO wave function\n", l);
      for (int i=0; i<mmax; i++) {
	double x = (double)i/(double)(mmax-1);
	x += 0.01;
	double r = 100.0*x*x*x*x*x - 1.0e-8;
	double ul;
	double rend = ChannelPotentials[l].ul.grid.End();
	if (r >= rend) {
	  int i=ChannelPotentials[l].ul.grid.NumPoints();
	  while (fabs(ChannelPotentials[l].ul(i)) == 0.0) i--;
	  double r2 = ChannelPotentials[l].ul.grid[i];
	  double r1 = ChannelPotentials[l].ul.grid[i-10];
	  double u2 = ChannelPotentials[l].ul(i);
	  double u1 = ChannelPotentials[l].ul(i-10);
	  double alpha = log(u1/u2)/(r2-r1);
	  // rend = min (rend, 60.0);
	  // double uend = ChannelPotentials[l].ul(rend);
	  // double alpha = log (ChannelPotentials[l].ul(rend-1.0)/uend);
	  //cerr << "alpha = " << alpha << " rend = " << r2 << " u2 = " << u2 << endl;
	  // ul = uend * exp(-alpha*(r-rend));
	  ul = u2 * exp(-alpha*(r - r2));
	}
	else
	  ul= ChannelPotentials[l].ul(r);
	fprintf (fout, "%23.16e ", ul);
	if ((i%3)==2)
	  fprintf (fout, "\n");
      }
  }

  fclose (fout);
}

void
PseudoClass::WriteASCII()
{
  FILE *fout = fopen ("pp.dat", "w");
  for (int i=1; i<PotentialGrid.NumPoints(); i++) {
    double r= PotentialGrid[i];
    fprintf (fout, "%24.16e ", r);
    for (int l=0; l<ChannelPotentials.size(); l++)
      fprintf (fout, "%24.16e ", ChannelPotentials[l].Vl(r));// /r);
    fprintf (fout, "\n");
  }
  fclose (fout);

  fout = fopen ("ul.dat", "w");
  for (int i=1; i<PotentialGrid.NumPoints(); i++) {
    double r= PotentialGrid[i];
    fprintf (fout, "%24.16e ", r);
    for (int l=0; l<ChannelPotentials.size(); l++)
      fprintf (fout, "%24.16e ", ChannelPotentials[l].ul(r));// /r);
    fprintf (fout, "\n");
  }
  fclose (fout);


}

bool
PseudoClass::ReadBFD_PP (string fileName)
{
  MemParserClass parser;
  assert(parser.OpenFile (fileName));
  
  assert (parser.FindToken ("Element Symbol:"));
  string symbol;
  assert (parser.ReadWord(symbol));
  AtomicNumber = SymbolToZMap[symbol];
  cerr << "Reading element " << symbol 
       << ", which has Z=" << AtomicNumber << endl;
  int numProtons, numProjectors;
  assert (parser.FindToken ("Number of replaced protons:"));
  assert (parser.ReadInt (numProtons));
  assert (parser.FindToken ("Number of projectors:"));
  assert (parser.ReadInt (numProjectors));
  ChannelPotentials.resize(numProjectors+1);
  assert (parser.FindToken ("Local component:"));
  assert (parser.FindToken ("Exp."));
  double c0, c1, c2, exp0, exp1, exp2;
  int power;
  assert (parser.ReadDouble(c0));
  PseudoCharge = c0;
  assert (parser.ReadInt (power));
  assert (power == -1);
  assert (parser.ReadDouble(exp0));

  assert (parser.ReadDouble(c1));
  assert (parser.ReadInt (power));
  assert (power == 1);
  assert (parser.ReadDouble(exp1));

  assert (parser.ReadDouble(c2));
  assert (parser.ReadInt (power));
  assert (power == 0);
  assert (parser.ReadDouble(exp2));

  // Create local channel spline
  vector<double> gridPoints, Vlocal, Vlocalr;
  FILE *fout = fopen ("Vlocal.dat", "w");
  for (double r=1.0e-8; r<150.0; r+=0.01) {
    gridPoints.push_back (r);
    double Vloc = (c0/r)*(expm1(-exp0*r*r)) +
      c1*r*exp(-exp1*r*r) +
      c2*exp(-exp2*r*r);
    double Vlocr = (c0)*(expm1(-exp0*r*r)) +
      c1*r*r*exp(-exp1*r*r) +
      c2*r*exp(-exp2*r*r);
    Vlocal.push_back(Vloc);
    Vlocalr.push_back(Vlocr);
    fprintf (fout, "%1.16e %1.16e %1.16e\n", r, Vloc, Vlocr);
  }
  fclose (fout);
  int local = numProjectors;
  if (LocalChannel<0) LocalChannel = local;
  PotentialGrid.Init (gridPoints);
  ChannelPotentials[local].Vl.Init(PotentialGrid, Vlocal);
  ChannelPotentials[local].Vlr.Init(PotentialGrid, Vlocalr);
  ChannelPotentials[local].l = local;

  // Now read nonlocal component
  assert (parser.FindToken ("Non-local component:"));
  assert (parser.FindToken ("Proj."));
  for (int l=0; l<numProjectors; l++) {
    char fname[100];
    snprintf (fname, 100, "V_%d.dat", l);
    FILE *fout = fopen (fname, "w");
    double c0, exp0;
    int power, channel;
    assert (parser.ReadDouble (c0));
    assert (parser.ReadInt (power));
    assert (power == 0);
    assert (parser.ReadDouble (exp0));
    assert (parser.FindToken ("|"));
    assert (parser.ReadInt(channel));
    assert (channel == l);
    assert (parser.FindToken("|"));
    vector<double> Vl(gridPoints.size()), Vlr(gridPoints.size());
    for (int i=0; i<gridPoints.size(); i++) {
      double r = gridPoints[i];
      Vl[i]  =  Vlocal[i] +   c0*exp(-exp0*r*r);
      Vlr[i] = Vlocalr[i] + r*c0*exp(-exp0*r*r);
      fprintf (fout, "%1.16e %1.16e %1.16e\n", r, Vl[i], Vlr[i]);
    }
    fclose (fout);
    ChannelPotentials[l].Vl.Init  (PotentialGrid, Vl);
    ChannelPotentials[l].Vlr.Init (PotentialGrid, Vlr);
    ChannelPotentials[l].l = l;
  }


  parser.CloseFile();
  return true;
}


bool
PseudoClass::ReadFHI_PP (string fileName)
{
  MemParserClass parser;
  if (!parser.OpenFile (fileName)) {
    cerr << "Could not open file \"" << fileName << "\".  Exitting.\n";
    exit(-1);
  }
  string temp;
  // Read top comment line
  assert (parser.ReadLine(temp));
  double atomCharge;
  assert (parser.ReadDouble(atomCharge));
  AtomicNumber = (int)round(atomCharge);
  assert (parser.ReadDouble(PseudoCharge));
  int date, pspcode, pspxc, lmax, lloc, mmax;
  assert (parser.ReadInt(date));
  assert (parser.ReadLine(temp));
  assert (parser.ReadInt(pspcode));
  assert (parser.ReadInt(pspxc));
  assert (parser.ReadInt(lmax));
  assert (parser.ReadInt(LocalChannel));
  assert (parser.ReadInt(mmax));
  // Read the rest of the line and skip the next 4
  for (int i=0; i<5; i++)
    assert (parser.FindToken("\n"));
  // Now we are in the FHI file proper
  double chrg;
  assert (parser.ReadDouble(chrg));
  if(fabs(chrg-PseudoCharge) > 1.0e-10) {
    cerr << "PseudoCharge = " << PseudoCharge << endl;
    cerr << "chrg = " << chrg << endl;
  }

  // Read number of channels and resize
  int numChannels;
  assert (parser.ReadInt(numChannels));
  ChannelPotentials.resize(numChannels);

  // Skip the next 10 lines
  for (int i=0; i<11; i++)
    assert (parser.FindToken("\n"));

  // Read the number of grid points
  for (int l=0; l<numChannels; l++) {
    int numPoints;
    double a_ratio;
    assert (parser.ReadInt(numPoints));
    assert (parser.ReadDouble(a_ratio));
    vector<double> points(numPoints), ul(numPoints), Vl(numPoints), 
      Vlr(numPoints);
    for (int m=0; m<numPoints; m++) {
      int mtest;
      assert (parser.ReadInt(mtest));
      assert (mtest == (m+1));
      assert (parser.ReadDouble(points[m]));
      assert (parser.ReadDouble(ul[m]));
      assert (parser.ReadDouble(Vl[m]));
      Vlr[m] = Vl[m]*points[m];
    }
    PotentialGrid.Init (points);
    ChannelPotentials[l].Vl.Init  (PotentialGrid, Vl);
    ChannelPotentials[l].Vlr.Init (PotentialGrid, Vlr);
    ChannelPotentials[l].ul.Init  (PotentialGrid, ul);
    ChannelPotentials[l].l = l;
    ChannelPotentials[l].n_principal = -1;
    ChannelPotentials[l].Cutoff = 0.0;
    ChannelPotentials[l].HasProjector = true;
  }
  // Now, figure out what rcut should be
  double rmax = PotentialGrid.End();
  int n = PotentialGrid.NumPoints()-1;
  bool diff = false;
  while ((n>0) && !diff) {
    double r = PotentialGrid[n];
    rmax = r;
    for (int l1=0; l1<numChannels; l1++) {
      double Vl1 = ChannelPotentials[l1].Vl(r);
      for (int l2=l1+1; l2<numChannels; l2++) {
	double Vl2 = ChannelPotentials[l2].Vl(r);
	if (fabs(Vl1-Vl2) > 1.0e-5)
	  diff = true;
      }
    }
    n--;
  }
  for (int l=0; l<ChannelPotentials.size(); l++)
    ChannelPotentials[l].Cutoff = rmax;

  parser.CloseFile();
  return true;
}

class KBprojector
{
public:
  int l;
  std::vector<double> beta, u, V;
  KBprojector() {}
  KBprojector(int l_, std::vector<double> beta_) :
    l(l_), beta(beta_) {}
};

bool
PseudoClass::ReadUPF_PP (string fileName)
{
  MemParserClass parser;
  if (!parser.OpenFile (fileName)) {
    cerr << "Could not open file \"" << fileName << "\".  Exitting.\n";
    exit(-1);
  }
  string temp;

  assert(parser.FindToken("<PP_HEADER>")); assert(parser.NextLine());
  // Skip version number
  assert(parser.NextLine());
  string element, ppType, coreCorrection;  
  assert(parser.ReadWord(element));  assert(parser.NextLine());
  assert(SymbolToZMap.find(element) != SymbolToZMap.end());
  AtomicNumber = SymbolToZMap[element];
  //  cerr << "Element = \"" << element << "\"\n";
  
  assert(parser.ReadWord(ppType));
  if (ppType != "NC") {
    cerr << "Psuedopotential type \"" << ppType 
	 << "\" is not norm-conserving.\n";
    abort();
  }
  assert (parser.NextLine());
  assert (parser.ReadWord(coreCorrection)); assert(parser.NextLine());
  if (coreCorrection != "F" && coreCorrection != "f" &&
      coreCorrection != "false" && coreCorrection != "False" &&
      coreCorrection != ".F." && coreCorrection != ".f" &&
      coreCorrection != "FALSE") {
    cerr << "It appears that a nonlinear core correction was used.  Cannot convert PP.\n";
    abort();
  }

  // Skip XC function info
  assert(parser.NextLine());
  assert(parser.ReadDouble(PseudoCharge));
  // Skip to end of line
  assert(parser.NextLine());
  // Skip total energy
  assert(parser.NextLine());
  // Skip energy cutoff
  assert(parser.NextLine());
  // Skip maximum angular momentum component
  assert(parser.NextLine());

  int numPoints;
  assert(parser.ReadInt(numPoints)); assert(parser.NextLine());
  int numWF, numProj;
  assert (parser.ReadInt(numWF)); assert(parser.ReadInt(numProj));
  int numChannels = numProj+1;
  
  ChannelPotentials.resize(numChannels);
  

  assert(parser.FindToken("Wavefunctions"));
//   for (int i=0; i<numWF; i++) {
//     assert(parser.FindToken("NL"));
//     int l;
//     double occ;
//     assert(parser.ReadInt(l));
//     assert(parser.ReadDouble(occ));
//   }
  assert(parser.FindToken("</PP_HEADER>"));
  
  assert(parser.FindToken("<PP_MESH>"));
  assert(parser.FindToken("<PP_R>"));
  vector<double> r(numPoints), dr(numPoints);
  for (int i=0; i<numPoints; i++) 
    assert(parser.ReadDouble(r[i]));
  assert(parser.FindToken("</PP_R>"));
  PotentialGrid.Init(r);


  assert(parser.FindToken("<PP_RAB>"));
  for (int i=0; i<numPoints; i++) 
    assert(parser.ReadDouble(dr[i]));
  assert(parser.FindToken("</PP_RAB>"));
  assert(parser.FindToken("</PP_MESH>"));
	 
  assert(parser.FindToken("<PP_LOCAL>"));
  vector<double> Vlocal(numPoints);
  for (int i=0; i<numPoints; i++) {
    assert(parser.ReadDouble(Vlocal[i]));
    Vlocal[i] *= 0.5;
  }
  assert(parser.FindToken("</PP_LOCAL>"));
  
  assert(parser.FindToken("<PP_NONLOCAL>"));

  // Store which channels have nonlocal projectors:
  // Input:  l    Output:  index of projector
  std::map<int,int> lmap;
  std::vector<KBprojector> projectors;
  for (int proj=0; proj<numProj; proj++) {
    assert (parser.FindToken("<PP_BETA>"));
    int Beta;
    assert (parser.ReadInt(Beta));
    assert(Beta == (proj+1));
    int l;
    assert (parser.ReadInt(l));    parser.NextLine();
    lmap[l] = proj;
    //    fprintf (stderr, "Found a projector for l=%d\n", l);

    int npt;
    assert (parser.ReadInt(npt));

    vector<double> beta(numPoints);
    for (int ir=0; ir<numPoints; ir++) 
      assert(parser.ReadDouble(beta[ir]));
    assert (parser.FindToken("</PP_BETA>"));

    projectors.push_back(KBprojector(l, beta));
  }
  assert(parser.FindToken("</PP_NONLOCAL>"));

  int llocal=0;
  while (lmap.find(llocal) != lmap.end()) llocal++;
  LocalChannel = llocal;
  //  fprintf (stderr, "First channel without a projector=%d\n", llocal);
  ChannelPotentials[llocal].Vl.Init(PotentialGrid, Vlocal);
  vector<double> Vlocal_r(numPoints);
  for (int ir=0; ir<numPoints; ir++)
    Vlocal_r[ir] = r[ir]*Vlocal[ir];
  ChannelPotentials[llocal].Vlr.Init(PotentialGrid, Vlocal_r);
  ChannelPotentials[llocal].l = llocal;

  if (numWF) {
    assert (parser.FindToken("PP_PSWFC"));
    assert (parser.NextLine());
    for (int iwf=0; iwf<numWF; iwf++) {
      string NL;
      assert (parser.ReadWord(NL));
      int l;
      assert (parser.ReadInt(l));
      double occ;
      assert (parser.ReadDouble(occ));
      parser.NextLine();
      std::vector<double> wf(numPoints);
      for (int ir=0; ir<numPoints; ir++) 
	assert (parser.ReadDouble(wf[ir]));
      // Now, setup the channel potentials
      if (l != llocal) {
	assert (lmap.find(l) != lmap.end());
	KBprojector &proj = projectors[lmap[l]];
	vector<double> Vl(numPoints), Vlr(numPoints);
	for (int ir=0; ir<numPoints; ir++) {
	  Vl[ir]  = Vlocal[ir];
	  if (wf[ir] != 0.0) 
	    Vl[ir] += proj.beta[ir]/(2.0*wf[ir]);
	  Vlr[ir] = r[ir]*Vl[ir];
	}
	Vl[0] = 2.0*Vl[1] - Vl[2];
	Vlr[0] = r[0] * Vl[0];
	
// 	for (int ir=0; ir<numPoints; ir++) 
// 	  fprintf (stderr, "Vl(%8.3f) = %22.16f\n", r[ir], Vl[ir]);

	ChannelPotentials[l].Vl.Init(PotentialGrid, Vl);
	ChannelPotentials[l].Vlr.Init(PotentialGrid, Vlr);
	ChannelPotentials[l].l = l;
      }
      ChannelPotentials[l].ul.Init(PotentialGrid, wf);
      ChannelPotentials[l].Occupation = occ;
      ChannelPotentials[l].HasProjector = true;
    }
    
    assert (parser.FindToken("/PP_PSWFC"));
  }

  // Now, figure out what rcut should be
  double rmax = PotentialGrid.End();
  int n = PotentialGrid.NumPoints()-1;
  bool diff = false;
  while ((n>0) && !diff) {
    double r = PotentialGrid[n];
    rmax = r;
    for (int l1=0; l1<numChannels; l1++) {
      double Vl1 = ChannelPotentials[l1].Vl(r);
      for (int l2=l1+1; l2<numChannels; l2++) {
	double Vl2 = ChannelPotentials[l2].Vl(r);
	if (fabs(Vl1-Vl2) > 1.0e-5)
	  diff = true;
      }
    }
    n--;
  }
  fprintf (stderr, "rmax = %1.5f  numChannels = %d\n", rmax, numChannels);
  for (int l=0; l<ChannelPotentials.size(); l++)
    ChannelPotentials[l].Cutoff = rmax;

  parser.CloseFile();
	 
  // Read top comment line
//   assert (parser.ReadLine(temp));
//   double atomCharge;
//   assert (parser.ReadDouble(atomCharge));
//   AtomicNumber = (int)round(atomCharge);
//   assert (parser.ReadDouble(PseudoCharge));
//   int date, pspcode, pspxc, lmax, lloc, mmax;
//   assert (parser.ReadInt(date));
//   assert (parser.ReadLine(temp));
//   assert (parser.ReadInt(pspcode));
//   assert (parser.ReadInt(pspxc));
//   assert (parser.ReadInt(lmax));
//   assert (parser.ReadInt(LocalChannel));
//   assert (parser.ReadInt(mmax));
//   // Read the rest of the line and skip the next 4
//   for (int i=0; i<5; i++)
//     assert (parser.FindToken("\n"));
//   // Now we are in the FHI file proper
//   double chrg;
//   assert (parser.ReadDouble(chrg));
//   if(fabs(chrg-PseudoCharge) > 1.0e-10) {
//     cerr << "PseudoCharge = " << PseudoCharge << endl;
//     cerr << "chrg = " << chrg << endl;
//   }

//   // Read number of channels and resize
//   int numChannels;
//   assert (parser.ReadInt(numChannels));
//   ChannelPotentials.resize(numChannels);

//   // Skip the next 10 lines
//   for (int i=0; i<11; i++)
//     assert (parser.FindToken("\n"));

//   // Read the number of grid points
//   for (int l=0; l<numChannels; l++) {
//     int numPoints;
//     double a_ratio;
//     assert (parser.ReadInt(numPoints));
//     assert (parser.ReadDouble(a_ratio));
//     vector<double> points(numPoints), ul(numPoints), Vl(numPoints), 
//       Vlr(numPoints);
//     for (int m=0; m<numPoints; m++) {
//       int mtest;
//       assert (parser.ReadInt(mtest));
//       assert (mtest == (m+1));
//       assert (parser.ReadDouble(points[m]));
//       assert (parser.ReadDouble(ul[m]));
//       assert (parser.ReadDouble(Vl[m]));
//       Vlr[m] = Vl[m]*points[m];
//     }
//     PotentialGrid.Init (points);
//     ChannelPotentials[l].Vl.Init  (PotentialGrid, Vl);
//     ChannelPotentials[l].Vlr.Init (PotentialGrid, Vlr);
//     ChannelPotentials[l].ul.Init  (PotentialGrid, ul);
//     ChannelPotentials[l].l = l;
//     ChannelPotentials[l].n_principal = -1;
//     ChannelPotentials[l].Cutoff = 0.0;
//     ChannelPotentials[l].HasProjector = true;
//   }
//   // Now, figure out what rcut should be
//   double rmax = PotentialGrid.End();
//   int n = PotentialGrid.NumPoints()-1;
//   bool diff = false;
//   while ((n>0) && !diff) {
//     double r = PotentialGrid[n];
//     rmax = r;
//     for (int l1=0; l1<numChannels; l1++) {
//       double Vl1 = ChannelPotentials[l1].Vl(r);
//       for (int l2=l1+1; l2<numChannels; l2++) {
// 	double Vl2 = ChannelPotentials[l2].Vl(r);
// 	if (fabs(Vl1-Vl2) > 1.0e-5)
// 	  diff = true;
//       }
//     }
//     n--;
//   }
//   for (int l=0; l<ChannelPotentials.size(); l++)
//     ChannelPotentials[l].Cutoff = rmax;

//   parser.CloseFile();
   return true;
}

bool
PseudoClass::ReadGAMESS_PP (string fileName)
{
  MemParserClass parser;
  assert (parser.OpenFile (fileName));

  string psp_name, psp_type;
  assert (parser.ReadWord(psp_name));
  assert (parser.ReadWord(psp_type));
  // izcore is the number of core electrons removed
  // Find atomic number.
  string symbol;
  int i=0;
  while ((i<psp_name.size()) && (psp_name[i] != '-')) {
    symbol = symbol + psp_name[i];
    i++;
  }
  cerr << "Atomic symbol = " << symbol << endl;
  int Z = SymbolToZMap[symbol];
  AtomicNumber = Z;
  int izcore, lmax;
  assert (parser.ReadInt(izcore));
  PseudoCharge = (double)(Z-izcore);
  assert (parser.ReadInt(lmax));
  assert (parser.FindToken("\n"));
  if (LocalChannel<0) LocalChannel = lmax;

  ChannelPotentials.resize(lmax+1);

  // Setup the potential grid
  vector<double> rPoints;
  if (WriteLogGrid) {
    const int numPoints = 2001;
    double end = 150.0;
    double step  = 0.00625;
    double scale = end/(expm1(step*(numPoints-1)));
    for (int i=1; i<numPoints; i++) 
      rPoints.push_back(scale * (expm1(step*i)));
  }
  else {
    for (double r=1.0e-10; r<=150.00001; r+=0.005)
      rPoints.push_back(r);
  }
  PotentialGrid.Init(rPoints);
  int N = rPoints.size();
  vector<double> Vlocal(N), Vl(N);
  vector<double> Vlocal_r(N), Vl_r(N);
  cerr << "PseudoCharge = " << PseudoCharge << endl;
  cerr << "lmax = " << lmax << endl;

  for (int index=0; index<=lmax; index++) {
    // lmax is first, followed by 0, 1, etc.
    int l = (index+lmax)%(lmax+1);
    int ngpot;
    assert (parser.ReadInt (ngpot));
    assert (parser.FindToken("\n"));
    vector<double> coefs(ngpot), exponents(ngpot);
    vector<int> powers(ngpot);
    
    // The local channel is given first
    for (int i=0; i<ngpot; i++) {
      assert (parser.ReadDouble(    coefs[i]));
      assert (parser.ReadInt(      powers[i]));
      assert (parser.ReadDouble(exponents[i]));
      assert (parser.FindToken ("\n"));
    }
    // Now, setup potentials
    if (index == 0) {
      for (int n=0; n<N; n++) {
	double r = rPoints[n];
	double val = 0.0;
	for (int i=0; i<ngpot; i++)
	  val += coefs[i] * pow(r,(double)(powers[i]-2)) 
	    * exp(-exponents[i]*r*r);
	Vlocal[n]   =   val -(double)PseudoCharge/r;
	Vlocal_r[n] = r*val - (double)PseudoCharge;
      }
      ChannelPotentials[l].l = l;
      ChannelPotentials[l].Vl.Init (PotentialGrid, Vlocal);
      ChannelPotentials[l].Vlr.Init (PotentialGrid, Vlocal_r);
    }
    else {
      for (int n=0; n<N; n++) {
	double r = rPoints[n];
	double val = 0.0;
	for (int i=0; i<ngpot; i++)
	  val += coefs[i] * pow(r,(double)(powers[i]-2)) 
	    * exp(-exponents[i]*r*r);
	Vl[n]   = Vlocal[n] + val;
	Vl_r[n] = Vlocal_r[n] + r*val;
      }
      ChannelPotentials[l].l = l;
      ChannelPotentials[l].Vl.Init (PotentialGrid, Vl);
      ChannelPotentials[l].Vlr.Init (PotentialGrid, Vl_r);
    }
  }
  // Now, compute cutoff radius
  bool done=false;
  int ir = rPoints.size()-1;
  if (ChannelPotentials.size() > 1) {
    while (ir > 0 && !done) {
      double vlocal = ChannelPotentials[LocalChannel].Vl(ir);
      for (int l=0; l<ChannelPotentials.size(); l++) 
	if (ChannelPotentials[l].l != LocalChannel)
	  if (fabs(ChannelPotentials[l].Vl(ir) - vlocal) > 1.0e-5)
	    done = true;
      ir--;
    }
  }
  else {
    while (ir > 0 && !done) {
      double vlocal = ChannelPotentials[LocalChannel].Vl(ir);
      double r = rPoints[ir];
      done = fabs (vlocal + PseudoCharge/r) > 1.0e-5;
      ir--;
    }
  }
  for (int l=0; l<ChannelPotentials.size(); l++)
    ChannelPotentials[l].Cutoff = rPoints[ir];


  cerr << "Finished reading GAMESS pseudopotential.\n";
  return true;
}

bool
ReadInt (string &s, int &n)
{
  int pos=0;
  string sint;
  while (s[pos] >= '0' && s[pos] <= '9') 
    sint.push_back(s[pos++]);
  s.erase(0,pos);
  n = atoi (sint.c_str());
  return true;
}

bool
ReadDouble (string &s, double &x)
{
  int pos=0;
  string sd;
  while ((s[pos] >= '0' && s[pos] <= '9') || s[pos]=='.') 
    sd.push_back(s[pos++]);
  s.erase(0,pos);
  x = strtod(sd.c_str(), (char**)NULL);
  return true;
}


Config::Config(string s)
{
  ChannelRevMap['s'] = 0;
  ChannelRevMap['p'] = 1;
  ChannelRevMap['d'] = 2;
  ChannelRevMap['f'] = 3;
  ChannelRevMap['g'] = 4;
  while (s.size() > 0) {
    int n, l;
    double occ;
    ReadInt (s, n);
    l = ChannelRevMap[s[0]];
    s.erase(0,1);
    assert (s[0] == '(');
    s.erase(0,1);
    ReadDouble (s, occ);
    assert (s[0] == ')');
    s.erase(0,1);
    States.push_back (State(n,l,occ));
  }
}



bool 
PseudoClass::ReadCASINO_PP (string fileName)
{
  MemParserClass parser;
  parser.OpenFile (fileName);

  assert (parser.FindToken("pseudo-charge"));
  assert (parser.ReadInt (AtomicNumber));
  assert (parser.ReadDouble (PseudoCharge));
  assert (parser.FindToken("(rydberg/hartree/ev):"));
  assert (parser.ReadWord(EnergyUnit));
  cerr << "EnergyUnit = " << EnergyUnit << endl;
  assert (parser.FindToken("(0=s,1=p,2=d..)"));
  assert (parser.ReadInt (LocalChannel));
  assert (parser.FindToken("grid points"));
  int numGridPoints;
  assert (parser.ReadInt(numGridPoints));
  vector<double> gridPoints(numGridPoints), Vl(numGridPoints), 
    Vlr(numGridPoints);
  assert (parser.FindToken ("in"));
  assert (parser.ReadWord(LengthUnit));
  cerr << "LengthUnit = " << LengthUnit << endl;
  assert (parser.FindToken("\n"));
  //  assert (parser.FindToken("atomic units"));
  for (int i=0; i<numGridPoints; i++) {
    assert (parser.ReadDouble(gridPoints[i]));
    gridPoints[i] *= UnitToBohrMap [LengthUnit];
  }
  PotentialGrid.Init (gridPoints);
  int l=0;
  bool done(false);
  while (!done) {
    if (!parser.FindToken("(L="))
      done = true;
    else {
      assert (parser.ReadInt(l));
      assert (parser.FindToken("\n"));
      for (int i=0; i<numGridPoints; i++) {
	assert (parser.ReadDouble(Vl[i]));
	Vl[i] *= UnitToHartreeMap[EnergyUnit];
	Vlr[i] = Vl[i];
	Vl[i] /= PotentialGrid[i];
      }
      Vl[0] = 2.0*Vl[1]-Vl[2];
      ChannelPotentials.push_back(ChannelPotentialClass());
      ChannelPotentials[l].Vl.Init(PotentialGrid, Vl);
      ChannelPotentials[l].Vlr.Init(PotentialGrid, Vlr);
      ChannelPotentials[l].l = l;
    }
  }
  cerr << "Found " << (l+1) << " l-channel potentials.\n";

  // Now read the summary file to get core radii
  if (!parser.OpenFile ("summary.txt")) {
    cerr << "Could not find summary.txt file.  Aborting.\n";
    abort();
  }
  string summaryUnits;
  assert (parser.FindToken("Units:"));
  assert (parser.ReadWord (summaryUnits));
  assert (parser.FindToken("core radii"));
  assert (parser.FindToken("\n"));
  assert (parser.FindToken("\n"));
  // Read core radii
  double rc;
  for (int i=0; i<=l; i++) {
    assert (parser.FindToken (ChannelMap[i]));
    assert (parser.ReadDouble (rc));
    ChannelPotentials[i].Cutoff = rc;
  }

  assert (parser.FindToken("r_loc"));
  assert (parser.FindToken("\n"));
  assert (parser.FindToken("\n"));
  assert (parser.FindToken("(Grid)"));
  double rloc;
  assert (parser.ReadDouble(rloc));
  
  // Overwrite radius with localization radius
  for (int i=0; i<=l; i++)
    ChannelPotentials[i].Cutoff = rloc;

  // Find ground state occupations
  assert (parser.FindToken ("Eigenvalues of the"));
  string GSconfig;
  parser.ReadWord (GSconfig);
  Config c(GSconfig);
  for (int is=0; is<c.size(); is++)
    ChannelPotentials[c[is].l].Occupation = c[is].occ;
  // char cStr[2];
  // cStr[1] = '\0';
  // cerr << "GSconfig = " << GSconfig << endl;
  // cStr[0] = GSconfig[1];
  // l = ChannelRevMap[cStr];
  // assert (GSconfig[2] == '(');
  // ChannelPotentials[l].Occupation = (double)(GSconfig[3]-'0');
  // assert (GSconfig[4] == ')');
  // if (GSconfig.size() > 5) {
  //   cStr[0] = GSconfig[6];
  //   l = ChannelRevMap[cStr];
  //   assert (GSconfig[7] == '(');
  //   ChannelPotentials[l].Occupation = (double)(GSconfig[8]-'0');
  //   assert (GSconfig[9] == ')');
  // }
  // if (GSconfig.size() > 10) {
  //   cStr[0] = GSconfig[11];
  //   l = ChannelRevMap[cStr];
  //   assert (GSconfig[12] == '(');
  //   ChannelPotentials[l].Occupation = (double)(GSconfig[13]-'0');
  //   assert (GSconfig[14] == ')');
  // }    

  // Now read eigenvalues of ground state
  assert (parser.FindToken ("Pseudo HF"));
  assert (parser.FindToken ("\n"));
  assert (parser.FindToken ("\n"));
  string channel;
  parser.ReadWord (channel);
  while ((channel=="s") || (channel=="p") || (channel=="d") || 
	 (channel=="f") || (channel=="g")) {
    int l = ChannelRevMap[channel];
    string tmp;
    parser.ReadWord (tmp);
    assert (parser.ReadDouble (ChannelPotentials[l].Eigenvalue));
    ChannelPotentials[l].Eigenvalue *= UnitToHartreeMap[summaryUnits];
    cerr << "Eigenvalue for " << channel << "-channel = " 
	 << ChannelPotentials[l].Eigenvalue << endl;
    parser.ReadLine (tmp);
    parser.ReadWord (channel);
  }
  // Read total energy
  assert (parser.FindToken ("total energy"));
  assert (parser.FindToken ("="));
  assert (parser.ReadDouble (TotalEnergy));
  TotalEnergy *= UnitToHartreeMap[summaryUnits];
  cerr << "Total energy = " << TotalEnergy << endl;

//   for (int i=0; i<=l; i++)
//     ChannelPotentials[i].Cutoff = rloc;

  return true;
}


/// Note:  the pseudopotentials must be read before the wave
/// functions.  
bool
PseudoClass::ReadCASINO_WF (string fileName, int l)
{
  // Make sure that l is not too high
  assert (l < ChannelPotentials.size());
  MemParserClass parser;
  parser.OpenFile (fileName);

  // Make sure this a wave function file
  assert (parser.FindToken ("wave function"));
  // Find atomic number
  assert (parser.FindToken ("Atomic number"));
  int atomicNumber;
  assert (parser.ReadInt(atomicNumber));
  assert (atomicNumber == AtomicNumber);
  int numOrbitals;
  assert (parser.FindToken("number of orbitals"));
  assert (parser.ReadInt(numOrbitals));
  assert (parser.FindToken("Radial grid"));
  assert (parser.FindToken("\n"));
  int numPoints;
  assert (parser.ReadInt(numPoints));
  assert (PotentialGrid.NumPoints() == numPoints);
  // Make sure we have the same grid as the potential file
  for (int i=0; i<numPoints; i++) {
    double r;
    assert(parser.ReadDouble(r));
    assert (fabs(r-PotentialGrid[i]) < 1.0e-10);
  }
  bool orbFound = false;
  for (int i=0; i<numOrbitals; i++) {
    assert (parser.FindToken ("Orbital #"));
    assert (parser.FindToken ("\n"));
    int spin, n, thisl;
    assert (parser.ReadInt(spin));
    assert (parser.ReadInt(n));
    assert (parser.ReadInt(thisl));
    if (thisl == l) {
      vector<double> ul(numPoints);
      for (int i=0; i<numPoints; i++)
	assert (parser.ReadDouble(ul[i]));
      ChannelPotentials[l].ul.Init(PotentialGrid, ul);
      ChannelPotentials[l].HasProjector = true;
      ChannelPotentials[l].n_principal = n;
      orbFound = true;
    }
  }
  if (!orbFound) {
    cerr << "Could not file orbital with l=" << l 
	 << " in file ""filename""" << endl; 
    exit(1);
  }
  return true;
}

int
PseudoClass::GetNumChannels()
{
  return ChannelPotentials.size();
}

bool
PseudoClass::HaveProjectors()
{
  bool has = true;
  for (int i=0; i<ChannelPotentials.size(); i++)
    has = has && ChannelPotentials[i].HasProjector;
  return has;
}



main(int argc, char **argv)
{
  // Create list of acceptable command-line parameters
  list<ParamClass> argList;
  // input formats
  argList.push_back(ParamClass("fhi_pot",   true));
  argList.push_back(ParamClass("casino_pot", true));
  argList.push_back(ParamClass("bfd_pot",    true));
  argList.push_back(ParamClass("gamess_pot", true));
  argList.push_back(ParamClass("upf_pot", true));

  // Projector parameters
  argList.push_back(ParamClass("casino_us",  true));
  argList.push_back(ParamClass("casino_up",  true));
  argList.push_back(ParamClass("casino_ud",  true));
  argList.push_back(ParamClass("casino_uf",  true));
  argList.push_back(ParamClass("s_ref",  true));
  argList.push_back(ParamClass("p_ref",  true));
  argList.push_back(ParamClass("d_ref",  true));
  argList.push_back(ParamClass("f_ref",  true));
  argList.push_back(ParamClass("g_ref",  true));
  argList.push_back(ParamClass("xml",    true));
  argList.push_back(ParamClass("tm",     true));
  argList.push_back(ParamClass("upf",    true));
  argList.push_back(ParamClass("fhi",    true));
  argList.push_back(ParamClass("fpmd",   true));
  argList.push_back(ParamClass("casino", true));
  argList.push_back(ParamClass("log_grid", false));
  argList.push_back(ParamClass("local_channel", true));

  CommandLineParserClass parser(argList);
  bool success = parser.Parse(argc, argv);
  if (!success || parser.NumFiles()!=0) {
    cerr << "Usage:  ppconvert  options\n"
	 << "  Options include:    \n"
	 << "   --casino_pot fname \n"
	 << "   --fhi_pot    fname \n"
         << "   --upf_pot    fname \n"
	 << "   --bfd_pot    fname \n"
         << "   --gamess_pot fname \n"
	 << "   --casino_us fname  \n"
	 << "   --casino_up fname  \n"
	 << "   --casino_ud fname  \n"
	 << "   --xml  fname.xml   \n"
	 << "   --tm   fname.tm    \n"
	 << "   --upf  fname.upf   \n"
	 << "   --fhi  fname.fhi   \n"
	 << "   --fpmd  fname.xml  \n"
	 << "   --casino fname.xml \n"
	 << "   --log_grid         \n"
	 << "   --local_channel l  \n";
    exit(1);
  }


  string xmlFile, tmFile;
  // Read the PP file
  PseudoClass nlpp;

  if (parser.Found("log_grid"))
    nlpp.WriteLogGrid = true;

  if (parser.Found("local_channel"))
    nlpp.SetLocalChannel(atoi(parser.GetArg("local_channel").c_str()));

  if (parser.Found("casino_pot"))
    nlpp.ReadCASINO_PP(parser.GetArg("casino_pot"));
  else if (parser.Found ("bfd_pot"))
    nlpp.ReadBFD_PP (parser.GetArg ("bfd_pot"));
  else if (parser.Found ("fhi_pot")) 
    nlpp.ReadFHI_PP (parser.GetArg ("fhi_pot"));
  else if (parser.Found ("upf_pot"))
    nlpp.ReadUPF_PP (parser.GetArg ("upf_pot"));
  else if (parser.Found ("gamess_pot"))
    nlpp.ReadGAMESS_PP (parser.GetArg("gamess_pot"));
  else {
    cerr << "Need to specify a potential file with --casino_pot "
	 << "or --bfd_pot or --fhi_pot or --upf_pot.\n";
    abort();
  }

  // Now check how the projectors are specified
  if (!nlpp.HaveProjectors()) {
    int numChannels = nlpp.GetNumChannels();
    if (numChannels > 0) {
      if (parser.Found ("s_ref")) 
	nlpp.CalcProjector (parser.GetArg("s_ref"), 0);
      else if (parser.Found("casino_us"))
	nlpp.ReadCASINO_WF(parser.GetArg("casino_us"), 0);
      else {
	cerr << "Please specify the s-channel projector with either "
	     << " --s_ref or --casino_us.\n";
	// exit(-1);
      }
    }
    if (numChannels > 1) {
      if (parser.Found ("p_ref"))
	nlpp.CalcProjector (parser.GetArg("p_ref"), 1);
      else if (parser.Found("casino_up"))
	nlpp.ReadCASINO_WF(parser.GetArg("casino_up"), 1);
      else {
	cerr << "Please specify the p-channel projector with either "
	     << " --p_ref or --casino_up.\n";
	// exit(-1);
      }
    }
    if (numChannels > 2) {
      if (parser.Found ("d_ref"))
	nlpp.CalcProjector (parser.GetArg("d_ref"), 2);
      else if(parser.Found("casino_ud"))
	nlpp.ReadCASINO_WF(parser.GetArg("casino_ud"), 2);
      else {
	cerr << "Please specify the d-channel projector with either "
	     << " --d_ref or --casino_ud.\n";
	// exit(-1);
      }
    }
    if (numChannels > 3) {
      if (parser.Found ("f_ref"))
	nlpp.CalcProjector (parser.GetArg("f_ref"), 3);
      else if(parser.Found("casino_uf"))
	nlpp.ReadCASINO_WF(parser.GetArg("casino_uf"), 3);
      else {
	cerr << "Please specify the f-channel projector with either "
	     << " --f_ref or --casino_uf.\n";
	// exit(-1);
      }
    }
    if (numChannels > 4) {
      if (parser.Found ("g_ref"))
	nlpp.CalcProjector (parser.GetArg("g_ref"), 4);
      else if(parser.Found("casino_ug"))
	nlpp.ReadCASINO_WF(parser.GetArg("casino_ug"), 4);
      else {
	cerr << "Please specify the g-channel projector with either "
	     << " --g_ref or --casino_ug.\n";
	// exit(-1);
      }
    }
  }

  
  if (parser.Found("xml")) 
    nlpp.WriteXML(parser.GetArg("xml"));
  if (parser.Found("tm"))
    nlpp.WriteABINIT(parser.GetArg("tm"));
  if (parser.Found ("upf"))
    nlpp.WriteUPF (parser.GetArg("upf"));
  if (parser.Found ("fhi"))
    nlpp.WriteFHI (parser.GetArg("fhi"));
  if (parser.Found ("fpmd"))
    nlpp.WriteFPMD (parser.GetArg("fpmd"));
  if (parser.Found ("casino"))
    nlpp.WriteCASINO (parser.GetArg("casino"));

//   nlpp.WriteXML(xmlFile);
//   nlpp.WriteABINIT();
  nlpp.WriteASCII();

//   nlpp.ReadCASINO_PP ("b_pp.data.HF");
//   nlpp.ReadCASINO_WF ("awfn.data_s2p1_2P", 0);
//   nlpp.ReadCASINO_WF ("awfn.data_s2p1_2P", 1);
//   nlpp.ReadCASINO_WF ("awfn.data_s2d1_2D", 2);
//   nlpp.WriteXML("test.xml");
}





#include "common/IO.h"

Array<double,1> vector2Array(vector<double> &vec)
{
  blitz::Array<double,1> array(vec.size());
  for (int i=0; i<vec.size(); i++)
    array(i) = vec[i];
  return array;
}



bool
PseudoClass::GetNextState (string &state, int &n, int &l, double &occ)
{
  if (state.length() <= 0)
    return false;

  if ((state[0]<'1') || (state[0] > '9')) {
    cerr << "Error parsing reference state.\n";
    abort();
  }
  n = state[0] - '0';
  string channel = state.substr (1,1);
  if ((channel != "s") && (channel != "p") &&
      (channel != "d") && (channel != "f") && (channel != "g")) {
    cerr << "Unrecognized angular momentum channel in reference state.\n";
    abort();
  }
  l = ChannelRevMap[channel];
  // cerr << "n=" << n << "  l=" << l << endl;
  if (state[2] != '(') {
    cerr << "Error parsing occupancy in reference state.\n";
    abort();
  }
  state.erase(0,3);
  
  string occString;
  while ((state.length() > 0) && (state[0] != ')')) {
    occString.append (state.substr(0,1));
    state.erase(0,1);
  }
  char *endptr;
  const char *occStr = occString.c_str();
  occ = strtod (occStr, &endptr);
  if (state[0] != ')') {
    cerr << "Expected a ')' in parsing reference state.\n";
    abort();
  }
  state.erase(0,1);
  
  return true;
}


#include "common/DFTAtom.h"
void
PseudoClass::CalcProjector(string refstate, int lchannel)
{
  DFTAtom atom;
  string saveState = refstate;
  // Temp storage for refstate;
  vector<double> occList;
  vector<int> nList, lList;
  bool done = false;
  int n, l;
  double occ;
  int lindex = -1;
  int index = 0;
  while (GetNextState(refstate, n, l, occ)) {
    nList.push_back (n);
    lList.push_back (l);
    occList.push_back (occ);
    if (l == lchannel && lindex == -1)
      lindex = index;
    index++;
  }
  if (lindex == -1) {
    nList.push_back(lchannel+1);
    lList.push_back(lchannel);
    occList.push_back(0.0);
    lindex = nList.size()-1;
  }
    
  atom.RadialWFs.resize(nList.size());
  for (int i=0; i<atom.RadialWFs.size(); i++) {
    // atom.RadialWFs(i).n = lList[i]+1;
    atom.RadialWFs(i).n = nList[i];
    atom.RadialWFs(i).l = lList[i];
    atom.RadialWFs(i).Occupancy = occList[i];
    atom.RadialWFs(i).Energy = -0.5;
  }
  // Define a grid for solving
  OptimalGrid *grid = new OptimalGrid(PseudoCharge, PotentialGrid.End());
  atom.SetGrid (grid);
  // Set the potential
  atom.SetBarePot (this);
  // Solve atom
  atom.NewMix = 0.75;
  cerr << "Solving atom for reference state " << saveState << ":\n";
  atom.Solve();
  
  // Now, initialize the channel u functions
  vector<double> ul(PotentialGrid.NumPoints()),
    rhoAtom(PotentialGrid.NumPoints());
  for (int i=0; i<ul.size(); i++) {
    double r = PotentialGrid[i];
    ul[i] = sqrt(4.0*M_PI)*atom.RadialWFs(lindex).u(r);
    rhoAtom[i] = atom.rho(r);
  }
  ChannelPotentials[lchannel].ul.Init (PotentialGrid, ul);
  ChannelPotentials[lchannel].HasProjector = true;
  ChannelPotentials[lchannel].Occupation = occList[lindex];
  RhoAtom.Init(PotentialGrid, rhoAtom);
}

double
PseudoClass::V(double r)
{ return ChannelPotentials[LocalChannel].Vl(r); }

double
PseudoClass::dVdr(double r)
{ return ChannelPotentials[LocalChannel].Vl.Deriv(r); }

double
PseudoClass::d2Vdr2(double r)
{ return ChannelPotentials[LocalChannel].Vl.Deriv2(r); }


double
PseudoClass::V(int l, double r)
{ 
  return ChannelPotentials[l].Vl(r); 
}

double
PseudoClass::dVdr(int l, double r)
{ return ChannelPotentials[l].Vl.Deriv(r); }

double
PseudoClass::d2Vdr2(int l, double r)
{ return ChannelPotentials[l].Vl.Deriv2(r); }

void
PseudoClass::Write(IOSectionClass &out)
{
}

void
PseudoClass::Read(IOSectionClass &in)
{
}

