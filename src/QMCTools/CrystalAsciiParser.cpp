//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Cody Melton, camelto2@ncsu.edu, North Carolina State University
//
// File created by: Cody Melton, camelto2@ncsu.edu, North Carolina State University
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCTools/CrystalAsciiParser.h"
#include <map>

CrystalAsciiParser::CrystalAsciiParser()
{
    PBC=true;
    ECP=false;
    BohrUnit=true;
    FixValence=true;
    SpinRestricted=true;
}

CrystalAsciiParser::CrystalAsciiParser(int argc, char** argv):
    QMCGaussianParserBase(argc,argv)
{
    PBC=true;
    ECP=false;
    BohrUnit=true;
    FixValence=true;
    SpinRestricted=true;
}

void CrystalAsciiParser::parse(const std::string &fname) 
{
    std::ifstream fin(fname.c_str());
    if (fin.fail())
    {
	std::cerr << "Error when opening file: " << fname << std::endl;
	abort();
    }
    pivot_begin=fin.tellg();


    std::string aline;
    search(fin,"CALCULATION",aline);
    parsewords(aline.c_str(),currentWords);
    if (currentWords[0] != "CRYSTAL")
    {
	std::cerr << "Currently only converts CRYSTAL calculations" << std::endl;
	abort();
    }
    PBC=true;

    search(fin,"N. OF ATOMS PER CELL",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfAtoms = atoi(currentWords[5].c_str());
    search(fin,"NUMBER OF AO",aline);
    parsewords(aline.c_str(),currentWords);
    numAO = atoi(currentWords[3].c_str());
    search(fin,"N. OF ELECTRONS PER CELL",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfEls = atoi(currentWords[5].c_str());

    search(fin,"TYPE OF CALCULATION",aline);
    parsewords(aline.c_str(),currentWords);
    if (currentWords[4] == "RESTRICTED")
    {
	SpinRestricted=true;
	NumberOfAlpha = NumberOfEls/2;
	NumberOfBeta = NumberOfAlpha;
	SpinMultiplicity = 0;
    }
    else
    {
	std::cerr<<"Currently CrystalAsciiParser only supports RHF calculations" << std::endl;
	abort();
    }

    std::cout << "NUMBER OF AOs: " << numAO << std::endl;
    SizeOfBasisSet = numAO;
    std::cout << "Size of Basis Set: " << SizeOfBasisSet << std::endl;
    numMO = numAO;
    std::cout << "NUMBER OF MOs: " << numMO << std::endl;


    BohrUnit=true; //will convert everything to Bohr
    std::cout << "Number of alpha electrons: " << NumberOfAlpha << std::endl;
    std::cout << "Number of beta electrons: " << NumberOfBeta << std::endl;
    std::cout << "Number of electrons: " << NumberOfEls << std::endl;
    std::cout << "SPIN MULTIPLICITY: " << SpinMultiplicity << std::endl;
    std::cout << "NUMBER OF ATOMS: " << NumberOfAtoms << std::endl;

    IonSystem.create(NumberOfAtoms);
    GroupName.resize(NumberOfAtoms);

    getCell(fin);

    getKpts(fin);

    getGeometry(fin);

    getGaussianCenters(fin);

    getMO(fin);

}

struct atzeff {
    int atomicNum;
    int zeff;
    atzeff() {
	atomicNum=0;
	zeff=0;
    }
    atzeff(int a, int b) {
	atomicNum = a;
	zeff = b;
    }
};

void CrystalAsciiParser::getGeometry(std::istream & is)
{

    const double ang_to_bohr = 1.0/0.529177e0;
    tags.clear();
    is.seekg(pivot_begin);
    std::string aline;
    search(is,"NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL",aline);
    parsewords(aline.c_str(),currentWords);
    int num_irr_ats = atoi(currentWords[8].c_str());
    search(is,"ATOM AT. N.",aline);
    std::map<int,atzeff> convAtNum;
    for (int i=0; i<num_irr_ats; i++)
    {
	getwords(currentWords,is);
	int at = atoi(currentWords[1].c_str());
	atzeff x(at,at);
	convAtNum.insert(std::pair<int,atzeff>(at,x));
    }

    std::vector<int> atomic_number;
    std::vector<int> idx(NumberOfAtoms);
    std::map<int,int> AtomIndexmap;
    Matrix<double> IonPos(NumberOfAtoms,3);
    search(is,"CARTESIAN COORDINATES - PRIMITIVE CELL",aline);
    getwords(currentWords,is);
    getwords(currentWords,is);
    getwords(currentWords,is);
    for (int i=0; i<NumberOfAtoms; i++)
    {
	getwords(currentWords,is);
	idx[i] = atoi(currentWords[1].c_str());
	tags.push_back(currentWords[2]);
	for (int d=0; d<3; d++)
	{
	    IonPos[i][d] = atof(currentWords[3+d].c_str())*ang_to_bohr;
	}

    }
    for (int i=0; i<NumberOfAtoms; i++)
    {
	AtomIndexmap.insert(std::pair<int,int>(i,idx[i]));
    }

    //Find which of the unique atoms have ECPs to determine the Zeff and atomic number
    if(lookFor(is,"PSEUDOPOTENTIAL INFORMATION",aline))
    {
	ECP=true;
	for (int i=0; i < num_irr_ats; i++)
	{
	    if(lookFor(is,"ATOMIC NUMBER",aline))
	    {
		parsewords(aline.c_str(),currentWords);
		int atnum = atoi(currentWords[2].c_str());
		int zeff = atoi(currentWords[5].c_str());
		convAtNum.at(atnum+200).atomicNum = atnum;
		convAtNum.at(atnum+200).zeff = zeff;
	    }
	}
    }

    SpeciesSet & species(IonSystem.getSpeciesSet());
    for (int i=0; i<NumberOfAtoms; i++)
    {
	for (int d=0; d<3; d++)
	{
	    IonSystem.R[i][d] = IonPos[i][d];
	}
	GroupName[i] = IonName[convAtNum.at(AtomIndexmap.at(i)).atomicNum];
	int speciesID = species.addSpecies(GroupName[i]);
	IonSystem.GroupID[i] = speciesID;
	species(AtomicNumberIndex,speciesID)=convAtNum.at(AtomIndexmap.at(i)).atomicNum;
	species(IonChargeIndex,speciesID)=convAtNum.at(AtomIndexmap.at(i)).zeff;
    }
}

void CrystalAsciiParser::getCell(std::istream & is)
{
    X.resize(3);
    Y.resize(3);
    Z.resize(3);

    is.seekg(pivot_begin);
    std::string aline;
    search(is,"DIRECT LATTICE VECTORS COMPON. (A.U.)",aline);
    getwords(currentWords,is);
    getwords(currentWords,is);
    X[0] = atof(currentWords[0].c_str());
    X[1] = atof(currentWords[1].c_str());
    X[2] = atof(currentWords[2].c_str());
    getwords(currentWords,is);
    Y[0] = atof(currentWords[0].c_str());
    Y[1] = atof(currentWords[1].c_str());
    Y[2] = atof(currentWords[2].c_str());
    getwords(currentWords,is);
    Z[0] = atof(currentWords[0].c_str());
    Z[1] = atof(currentWords[1].c_str());
    Z[2] = atof(currentWords[2].c_str());
    std::cout<<"Lattice parameters in Bohr:"<<std::endl;
    std::cout<<X[0]<<"  "<<X[1]<<"  "<<X[2]<<std::endl;
    std::cout<<Y[0]<<"  "<<Y[1]<<"  "<<Y[2]<<std::endl;
    std::cout<<Z[0]<<"  "<<Z[1]<<"  "<<Z[2]<<std::endl;
}

void CrystalAsciiParser::getKpts(std::istream& is)
{

    std::cout << "getKpts needs to be checked" << std::endl;
    is.seekg(pivot_begin);
    std::string aline;
    search(is,"NUMBER OF K POINTS IN THE IBZ",aline);
    parsewords(aline.c_str(),currentWords);
    NbKpts = atoi(currentWords[12].c_str());
    search(is,"K POINTS COORDINATES (OBLIQUE COORDINATES IN UNITS OF IS",aline);
    std::string ignore = ")";
    parsewords(aline.c_str(),currentWords,ignore);
    int is_unit = atoi(currentWords[10].c_str());

    Kpoints_Coord.resize(NbKpts);
    for (int i = 0; i < NbKpts; i++)
    {
	getwords(currentWords,is,ignore);
	int nkpts;
	//Has at max 4 kpoints per line
	nkpts = currentWords.size()/4;
	for (int j = 0; j < nkpts; j++)
	{
	    Kpoints_Coord[i].resize(3);
	    int a = atoi(currentWords[4*j+1].c_str());
	    int b = atoi(currentWords[4*j+2].c_str());
	    int c = atoi(currentWords[4*j+3].c_str()); 
	    Kpoints_Coord[i][0]=a/(double)is_unit;
	    Kpoints_Coord[i][1]=b/(double)is_unit;
	    Kpoints_Coord[i][2]=c/(double)is_unit;
	    std::cout << "Kpoint " << i << " " << Kpoints_Coord[i][0]
            << " " << Kpoints_Coord[i][1] << " " << Kpoints_Coord[i][2] << std::endl;
	    i++;
	}
	i--;
    }
}

void CrystalAsciiParser::getMO(std::istream& is)
{
}

void CrystalAsciiParser::getGaussianCenters(std::istream& is)
{
    is.clear();
    is.seekg(pivot_begin);
    gBound.resize(NumberOfAtoms+1);
    int ng,nx;
    std::string aline;
    std::map<std::string,int> basisDataMap;
    int nUniqAt=0;
    for(int i=0; i<NumberOfAtoms;i++)
    {
	std::map<std::string,int>::iterator it(basisDataMap.find(tags[i]));
	if(it==basisDataMap.end())
	{
	    basisDataMap[tags[i]]=nUniqAt++;
	}
    }

    std::vector<std::vector<double> > expo(nUniqAt),coef(nUniqAt),coef2(nUniqAt);
    std::vector<int> nshll(nUniqAt,0);
    std::vector<std::vector<int> > ncoeffpershell(nUniqAt);
    std::vector<std::vector<std::string> > shID(nUniqAt);
    std::map<std::string,int> gsMap;
    gsMap[std::string("S")]=1;
    gsMap[std::string("SP")]=2;
    gsMap[std::string("P")]=3;
    gsMap[std::string("D")]=4;
    gsMap[std::string("F")]=5;
    gsMap[std::string("G")]=6;

    search(is,"LOCAL ATOMIC FUNCTIONS BASIS SET",aline);
    getwords(currentWords,is); 
    getwords(currentWords,is); 
    getwords(currentWords,is); 

    int currPos=-1;
    while(true)
    {
	getwords(currentWords,is);
	if ((currentWords[0] == "INFORMATION") 
           || (currentWords[0] == "*******************************************************************************"))
	{
	    break;
	}
	if(currentWords.size()==5)//found species
	{
	    std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[1]));
	    if(it == basisDataMap.end())
	    {
		std::cerr << "Error in parser.\n";
		abort();
	    }
	    currPos=it->second;
	    bool newgroup=(nshll[currPos]==0);
	    if(newgroup)
	    {
	        ncoeffpershell[currPos].clear();
	        ncoeffpershell[currPos].push_back(0);
	        shID[currPos].clear();
	        shID[currPos].push_back("None");
	    }
            std::string type;
	    while(true)
	    {
		std::streampos pivot=is.tellg();
		getwords(currentWords,is);
		if((currentWords.size() == 5) 
		  || (currentWords[0] == "INFORMATION")
                  || (currentWords[0] == "*******************************************************************************"))
		{
		    is.seekg(pivot);
		    break;
		}
		if (gsMap.find(currentWords[1])!=gsMap.end() || gsMap.find(currentWords[2])!=gsMap.end())
		{
		    if (newgroup)
		    {
		        if(gsMap.find(currentWords[1])!=gsMap.end())
		        {
		            type = currentWords[1];
		        }
		        else
		            type = currentWords[2];

		        if (gsMap[type] == 2)
		        {
		            std::cerr << "Can't handle SP basis states yet. Fix latex. \n";
		            abort();
		        }
		        if(gsMap.find(type) == gsMap.end())
		        {
		            std::cerr << "Unhandled primitive type " << type << std::endl;
		            abort();
		        }

		        nshll[currPos]++;
	    	        ncoeffpershell[currPos].push_back(0);
	    	        shID[currPos].push_back("None");
		    }
		}
		else 
		{
		    if(newgroup)
		    {
			double expt,coeff1,coeff2,coeff3;
			if (currentWords.size()==3)
			{//The primitive has a - sign between coefficients
			    int j;
			    for(int i=0; i < 3; i++)
			    {
				if(currentWords[i].length() == 19)
				{
				    j = i;
				}
			    }
			    std::string str1,str2;
			    for (int i=0; i<19; i++)
			    {
				if (i<9)
				{
				    str1+=currentWords[j][i];
				}
				else
				{
				    str2+=currentWords[j][i];
				}
			    }

			    if (j == 0)
			    {
				expt = atof(str1.c_str());
				coeff1 = atof(str2.c_str());
				coeff2 = atof(currentWords[1].c_str());
				coeff3 = atof(currentWords[2].c_str());
			    }
			    else if (j == 1)
			    {
				expt = atof(currentWords[0].c_str());
				coeff1 = atof(str1.c_str());
				coeff2 = atof(str2.c_str());
				coeff3 = atof(currentWords[2].c_str());
			    }
			    else 
			    {
				expt = atof(currentWords[0].c_str());
				coeff1 = atof(currentWords[1].c_str());
				coeff2 = atof(str1.c_str());
				coeff3 = atof(str2.c_str());
			    }
			}
			else{
			    expt = atof(currentWords[0].c_str());
			    coeff1 = atof(currentWords[1].c_str());
			    coeff2 = atof(currentWords[2].c_str());
			    coeff3 = atof(currentWords[3].c_str());
			}
			expo[currPos].push_back(expt);
			if(gsMap[type]==1)
			{
			    coef[currPos].push_back(coeff1);
			}
			else if (gsMap[type]==3)
			{
			    coef[currPos].push_back(coeff2);
			}
			else
			{
			    coef[currPos].push_back(coeff3);
			}
			ncoeffpershell[currPos][nshll[currPos]-1]++;
			shID[currPos][nshll[currPos]-1] = type;
		    }
		    if (debug)
		    {
		        std::cout << currPos<< ":"<< expo[currPos].back() << " " << coef[currPos].back() << " " 
		    	          << ncoeffpershell[currPos][nshll[currPos]]
		    	          << " " << shID[currPos][nshll[currPos]] << std::endl;
		    }
		}
	    }
	}
    }

    gShell.clear();
    gNumber.clear();
    gExp.clear();
    gC0.clear();
    gC1.clear();
    int gtot = 0;
    for (int i=0; i<NumberOfAtoms; i++)
    {
	std::map<std::string,int>::iterator it(basisDataMap.find(tags[i]));
	if (it == basisDataMap.end())
	{
	    std::cerr << "Error in parser" << std::endl;
	    abort();
	}
	gBound[i] = gtot;
	int indx = it->second;
	gtot+=nshll[indx];
	for(int k = 0; k<nshll[indx]; k++) {
	    gShell.push_back(gsMap[shID[indx][k]]);
	    gNumber.push_back(ncoeffpershell[indx][k]);
	}
	for (int k=0; k<expo[indx].size(); k++) {
	    gExp.push_back(expo[indx][k]);
	    gC0.push_back(coef[indx][k]);
	}
    }
    gBound[NumberOfAtoms] = gtot;

}
