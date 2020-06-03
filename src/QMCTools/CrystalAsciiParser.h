//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TOOLS_CRYSTAL_OUT_H
#define QMCPLUSPLUS_TOOLS_CRYSTAL_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <map>

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

class CrystalAsciiParser: public QMCGaussianParserBase,
    public OhmmsAsciiParser
{

public:

    CrystalAsciiParser();

    CrystalAsciiParser(int argc, char** argv);

    std::streampos pivot_begin;
    std::vector<std::string> tags;
    std::string MOtype;
    int readtype, numAO;
    int NFZC, NEXT, NTOT, NAC;

    void parse(const std::string& fname);

    void getGeometry(std::istream& is);
    void getCell(std::istream& is);
    void getKpts(std::istream& is);
    void getMO(std::istream& is);
    void getGaussianCenters(std::istream& is);

private:

    void getKMO(std::istream & is, std::vector< std::vector< double > > & Mat);
    void getKMO(std::istream & is, std::vector< std::vector< std::complex<double> > > & CMat);
    void scaleBasis();
    void dumpHDF5(const std::string &fname);
    std::string version;
    int NumberOfSpecies;
    std::map<int,atzeff> convAtNum;
    std::map<int,int> AtomIndexmap;
    bool IsComplex;
    std::map<std::string,int> basisDataMap;
    std::vector<std::vector<double> > expo,coef;
    std::vector<int> nshll;
    std::vector<std::vector<int> > ncoeffpershell;
    std::vector<std::vector<std::string> > shID;
    std::vector<std::vector<double> > Kpoints_Coord;
    std::map<std::string,int> gsMap;
    std::vector< std::vector< std::vector< double> > > real_kmos;
    std::vector< std::vector< std::vector< std::complex<double> > > > complex_kmos;
    std::vector< std::vector< double > > eigvals;
    std::map<std::string,int> str_to_l;
    int spin;
};

#endif
