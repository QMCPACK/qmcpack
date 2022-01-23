#include "SkParserASCII.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include <fstream>
#include <iostream>


namespace qmcplusplus
{
using RealType = SkParserASCII::RealType;
using PosType  = SkParserASCII::PosType;

vector<vector<RealType>> SkParserASCII::read_sk_file(const string& fname)
{
  vector<vector<RealType>> skdata(0);

  vector<RealType> tmp(5);

  ifstream f;
  f.open(fname.c_str(), ifstream::in);

  string tmpstring;           //just a sink for getline.
  std::getline(f, tmpstring); //skip the first line.

  while (!f.eof())
  {
    RealType x = 0, y = 0;
    f >> tmp[KX] >> tmp[KY] >> tmp[KZ] >> tmp[SK] >> tmp[SKERR];
    //corresponds to kx, ky, kz, S(k), err
    if (!f.eof())
      skdata.push_back(tmp);
  }

  return skdata;
}

vector<PosType> SkParserASCII::get_grid_from_data(vector<vector<RealType>>& filedata)
{
  vector<PosType> kgrid(filedata.size());

  for (int i = 0; i < filedata.size(); i++)
  {
    kgrid[i][0] = filedata[i][KX];
    kgrid[i][1] = filedata[i][KY];
    kgrid[i][2] = filedata[i][KZ];
  }

  return kgrid;
}

vector<RealType> SkParserASCII::get_sk_from_data(vector<vector<RealType>>& filedata)
{
  vector<RealType> sk(filedata.size());
  for (int i = 0; i < filedata.size(); i++)
    sk[i] = filedata[i][SK];

  return sk;
}

vector<RealType> SkParserASCII::get_skerr_from_data(vector<vector<RealType>>& filedata)
{
  vector<RealType> skerr(filedata.size());
  for (int i = 0; i < filedata.size(); i++)
    skerr[i] = filedata[i][SKERR];

  return skerr;
}

void SkParserASCII::parse(const string& fname)
{
  vector<vector<RealType>> rawdata(0);
  rawdata   = read_sk_file(fname);
  kgridraw  = get_grid_from_data(rawdata);
  skraw     = get_sk_from_data(rawdata);
  skerr_raw = get_skerr_from_data(rawdata);


  //  cout<<"Ok.  In SkParserASCII\n";
  //  cout<<" print kgridraw, skraw, skerr\n";
  //  for(int i=0; i<kgridraw.size();i++) cout<<kgridraw[i][0]<<" "<<kgridraw[i][1]<<" "<<kgridraw[i][2]<<" "<<skraw[i]<<" "<<skerr_raw[i]<<endl;
  hasGrid        = false;
  isNormalized   = true;
  isParseSuccess = true;
}
} // namespace qmcplusplus
