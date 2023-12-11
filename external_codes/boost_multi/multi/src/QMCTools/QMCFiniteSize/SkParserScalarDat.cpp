#include "SkParserScalarDat.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include <fstream>
#include <iostream>


namespace qmcplusplus
{
using RealType = SkParserScalarDat::RealType;
using PosType  = SkParserScalarDat::PosType;

void SkParserScalarDat::read_sk_file(const std::string& fname)
{
  IndexType IDENT = 7; //This is the position of character which denotes difference
                       //between rhok_e_e_X, rhok_e_i_X, and rhok_e_r_X.  Either e,i, or r.
  std::vector<std::vector<RealType>> skdata(0);

  std::vector<RealType> rho_ee(0);
  std::vector<RealType> rho_ee_err(0);
  std::vector<RealType> rho_i(0);
  std::vector<RealType> rho_r(0);
  std::vector<RealType> rho_i_err(0);
  std::vector<RealType> rho_r_err(0);

  std::ifstream f;
  f.open(fname.c_str(), std::ifstream::in);

  std::string obsname(""), eq(""), pm(""); //just a sink for getline.
  RealType val_tmp(0.0), err_tmp(0.0);


  while (!f.eof())
  {
    f >> obsname >> eq >> val_tmp >> pm >> err_tmp;

    if (obsname.find("rhok_e_") != std::string::npos)
    {
      //Great.  Found the sk data.  Now to determine if its rhok*rhok, Re(Rhok), or Im(rhok)
      if (obsname[IDENT] == 'e')
      {
        //        app_log()<<" rho_ee_"<<rho_ee.size()<<" = "<<val_tmp<<" +/- "<<err_tmp<<endl;
        rho_ee.push_back(val_tmp);
        rho_ee_err.push_back(err_tmp);
      }
      else if (obsname[IDENT] == 'i')
      {
        //        app_log()<<"  rho_i_"<<rho_i.size()<<" = "<<val_tmp<<" +/- "<<err_tmp<<endl;
        rho_i.push_back(val_tmp);
        rho_i_err.push_back(err_tmp);
      }
      else if (obsname[IDENT] == 'r')
      {
        //      app_log()<<"  rho_r_"<<rho_r.size()<<" = "<<val_tmp<<" +/- "<<err_tmp<<endl;
        rho_r.push_back(val_tmp);
        rho_r_err.push_back(err_tmp);
      }
      else
        APP_ABORT("ERROR SkParserScalarDat.  Unexpected rhok_e_X... found.  ");

      obsname = "";
    }
    else
      continue;
    //corresponds to kx, ky, kz, S(k), err
  }

  //check to make sure all vectors have the same length;
  IndexType Nk = rho_ee.size();
  if (rho_ee_err.size() != Nk)
    APP_ABORT("ERROR SkParserScalarDat: density data not the same size");
  if (rho_i.size() != Nk || rho_i_err.size() != Nk)
    APP_ABORT("ERROR SkParserScalarDat: density data not the same size");
  if (rho_r.size() != Nk || rho_r_err.size() != Nk)
    APP_ABORT("ERROR SkParserScalarDat: density data not the same size");

  skraw.resize(Nk);
  skerr_raw.resize(Nk);
  //Now we have all the data.  Need to build Ne*dS(k)=<rho-k*rhok> - <rho-k><rhok>
  //  app_log()<<"DEBUG:\n";
  //  app_log()<<"i sk skerr rho_ee rho_ee_err rho_r rho_r_err rho_i rho_i_err\n";
  for (IndexType i = 0; i < Nk; i++)
  {
    RealType r_r(0.0);
    RealType r_i(0.0);
    //The 3.0 in the following statements is the standard deviation.
    //  If function value is greater than 3standard deviations above error,
    //  then we set it.  Otherwise default to zero.
    if (rho_r_err[i] == 0 || std::abs(rho_r[i]) / rho_r_err[i] > 3.0)
      r_r = rho_r[i];
    if (rho_i_err[i] == 0 || std::abs(rho_i[i]) / rho_i_err[i] > 3.0)
      r_i = rho_i[i];

    skraw[i]     = rho_ee[i] - (r_r * r_r + r_i * r_i); //<rho_-k><rho_k> = |rho_k|^2
    skerr_raw[i] = rho_ee_err[i]; //Actual error distribution is a gaussian + a product distribution.
                                  //However, error in rho-k*rhok is drastically larger than rhok, so
                                  //we ignore the rhok contribution.

    //   app_log()<<i<<" "<<skraw[i]<<" "<<skerr_raw[i]<<" "<<rho_ee[i]<<" "<<rho_ee_err[i]<<" "<<rho_r[i]<<" "<<rho_r_err[i]<<" "<<rho_i[i]<<" "<<rho_i_err[i]<<endl;
  }

  return;
}


void SkParserScalarDat::parse(const std::string& fname)
{
  // vector<vector<RealType> > rawdata(0);
  read_sk_file(fname);
  // kgridraw=get_grid_from_data(rawdata);
  // skraw=get_sk_from_data(rawdata);
  //  skerr_raw=get_skerr_from_data(rawdata);


  //  cout<<"Ok.  In SkParserScalarDat\n";
  //  cout<<" print kgridraw, skraw, skerr\n";
  //  for(int i=0; i<kgridraw.size();i++) cout<<kgridraw[i][0]<<" "<<kgridraw[i][1]<<" "<<kgridraw[i][2]<<" "<<skraw[i]<<" "<<skerr_raw[i]<<endl;
  hasGrid        = false;
  isNormalized   = false;
  isParseSuccess = true;
}
} // namespace qmcplusplus
