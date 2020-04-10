/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////


#ifndef WRITE_ESHDF_H
#define WRITE_ESHDF_H
#include "io/hdf_archive.h"
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <complex>

class XmlNode;
class FftContainer;
class KPoint;

class EshdfFile
{
  typedef std::map<std::vector<int>, std::complex<double>> momap_t;
  typedef std::pair<std::vector<int>, std::complex<double>> mopair_t;

private:
  qmcplusplus::hdf_archive outfile_;

  int wrapped(int i, int size) const;
  int getIntsOnly(const std::string& str) const;
  void writeApplication(const std::string& appName, int major, int minor, int sub);
  void writeCreator();
  void writeFormat();
  void writeVersion();

  // helper functions meant for qbox
  void readInEigFcn(const XmlNode& nd, FftContainer& cont);
  void handleSpinGroup(const XmlNode* nd, double& nocc, FftContainer& cont);
  double getOccupation(const XmlNode* nd) const;

  // helper functions meant for espresso
  void handleDensity(const XmlNode& qeXml, const std::string& dir_name, int spinpol);
  std::vector<double> getPtvs(const XmlNode& qeXml);
  void processKPts(const XmlNode& band_structure_xml,
                   const std::vector<double>& ptvs,
                   std::vector<std::vector<double>>& eigenvals,
                   std::vector<std::vector<double>>& occupations,
                   std::vector<KPoint>& kpts,
                   std::vector<double>& weights,
                   std::vector<int>& ngvecs);
  void getNumElectrons(std::vector<std::vector<double>>& occupations,
                       std::vector<double>& weights,
                       int& nup,
                       int& ndn,
                       int spinpol,
                       int ncol);
  void handleKpt(int kpt_num,
                 const std::string& dir_name,
                 KPoint& kpt,
                 const std::vector<double>& eigenvalues,
                 double weight,
                 int spinpol,
                 int noncol,
                 const momap_t& moref);

  void readKptGvecs(int kpt_num, const std::string& dir_name, int spinpol, momap_t& morefmap);

  qmcplusplus::hdf_archive openHdfFileForRead(const std::string& fname);

  EshdfFile(const EshdfFile& f);            // no copy constructor
  EshdfFile& operator=(const EshdfFile& f); // operator= not allowed
public:
  EshdfFile(const std::string& hdfFileName);
  ~EshdfFile();

  void writeQEBoilerPlate(const XmlNode& qeXml);
  void writeQboxBoilerPlate(const XmlNode& qboxSample);
  void writeQESupercell(const XmlNode& qeXml);
  void writeQboxSupercell(const XmlNode& qboxSample);
  void writeQEAtoms(const XmlNode& qeXml);
  void writeQboxAtoms(const XmlNode& qboxSample);
  void writeQEElectrons(const XmlNode& qeXml, const std::string& dir_name);
  void writeQboxElectrons(const XmlNode& qboxSample);
};

#endif
