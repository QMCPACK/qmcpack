#ifndef WRITE_ESHDF_H
#define WRITE_ESHDF_H
#include "hdf5.h"
#include <string>
#include <vector>
#include <cmath>

class xmlNode;
class fftContainer;

class eshdfFile {
private:
  hid_t file;
  herr_t error;

  int wrapped(int i, int size) const;
  void writeApplication(const std::string& appName, int major, int minor, int sub);
  void writeCreator();
  void writeFormat();
  void writeVersion();

  void readInEigFcn(const xmlNode& nd, fftContainer& cont);
  void handleSpinGroup(const xmlNode* nd, hid_t groupLoc, double& nocc, fftContainer& cont);
  double getOccupation(const xmlNode* nd) const;

  eshdfFile(const eshdfFile& f); // no copy constructor
  eshdfFile& operator=(const eshdfFile& f); // operator= not allowed
public:
  eshdfFile(const std::string& hdfFileName);
  ~eshdfFile();

  void writeBoilerPlate(const std::string& appName, const xmlNode& qboxSample);
  void writeSupercell(const xmlNode& qboxSample);
  void writeAtoms(const xmlNode& qboxSample);
  void writeElectrons(const xmlNode& qboxSample);
};

#endif
