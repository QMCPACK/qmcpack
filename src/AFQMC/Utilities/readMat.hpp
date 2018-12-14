#ifndef AFQMC_READMAT_HPP
#define AFQMC_READMAT_HPP

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<ctype.h>

#include "Utilities/SimpleParser.h"

#include "AFQMC/config.h"

namespace qmcplusplus
{

bool readMat( std::string fileName, ComplexMatrix& Mat)
{

  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
     return false;
  }

  bool fullMOMat = false;
  bool Cstyle = true;
  bool wfn_type;
  std::string type;

  std::vector<std::string> words;
  getwords(words,in);
  do {
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
    getwords(words,in);
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
  } while((words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));

  int i,j;
  ComplexType v;

  while(!in.eof() && !in.fail()) {

    in>>i >>j >>v;
    if(in.fail() || in.eof()) break;
    
    Mat(i,j) += v; 
        
  }

  in.close();

}

}

#endif

