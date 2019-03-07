#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include "XmlRep.h"
#include "WriteEshdf.h"
using namespace std;

void convertToVecStrings(int argc, char* argv[], std::vector<std::string>& vec) {
  for (int i = 1; i < argc; i++) {
    vec.push_back(argv[i]);
  }
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[]) {
  vector<string> vecParams;
  convertToVecStrings(argc, argv, vecParams);
  string fname;
  for (int i = 0; i < vecParams.size(); i++) {
    // in principle check for other flags here
    fname = vecParams[i];
  }
  if (file_exists(fname) == false) {
    cout << "must specify a valid xml file name as an argument" << endl;
    exit(1);
  }

  ifstream is(fname.c_str());
  xmlNode qboxSampleXml(&is, 0, true);


  string ofname = "qboxEshdf.h5";
  eshdfFile outFile(ofname);
  outFile.writeBoilerPlate("qbox", qboxSampleXml);
  outFile.writeSupercell(qboxSampleXml);
  outFile.writeAtoms(qboxSampleXml);
  outFile.writeElectrons(qboxSampleXml);
}
