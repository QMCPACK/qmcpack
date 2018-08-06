#ifndef QMCPLUSPLUS_TOOLS_CRYSTAL_OUT_H
#define QMCPLUSPLUS_TOOLS_CRYSTAL_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <complex>

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
};

#endif
