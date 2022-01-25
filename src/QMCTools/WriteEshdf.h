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
#include "hdf/hdf_archive.h"
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <complex>

class XmlNode;
class FftContainer;
class KPoint;

/*! helper class meant to control access ot the eshdf file that is being
    created and provide helper functions so that the proper xml data
    from QBOX or quantum espresso can be handed to it and the data can
    be processed and written to the eshdf file 
*/
class EshdfFile
{
  using momap_t  = std::map<std::vector<int>, std::complex<double>>;
  using mopair_t = std::pair<std::vector<int>, std::complex<double>>;

private:
  qmcplusplus::hdf_archive outfile_;

  /*! helper indexing function to go from the interval [0,size) to (-size/2,size/2)     */
  int wrapped(int i, int size) const;
  /*! helper function that takes a string and returns an int consisting of only the numerals in the string */
  int getIntsOnly(const std::string& str) const;
  /*! writes to the application field of outfile_ (code -> appName), (version -> major,minor,sub) */
  void writeApplication(const std::string& appName, int major, int minor, int sub);
  /*! writes to the creator field of outfile_ (program_name -> convertpw4qmc), (version -> 0,1,0) */
  void writeCreator();
  /*! writes to the format field of outfile_ (format -> ES-HDF) */
  void writeFormat();
  /*! writes to the version field of outfile_ (version -> 2,1,0) */
  void writeVersion();

  // helper functions meant for qbox
  /*! to be handed a grid_function tag of the qbox sample file.  Will then read in the
      eigenfunciton (spo) in real-space and use the FftContainer cont to transform
      them to k-space */
  void readInEigFcn(const XmlNode& nd, FftContainer& cont);
  /*! to be handed a slater_determinant tag of the qbox sample file.  Will read in the
      occupations and call readInEigFcn for each state, also needs an FftContainer to
      pass on to readInEigFcn */
  void handleSpinGroup(const XmlNode* nd, double& nocc, FftContainer& cont);
  /*! to be handed a slater_determinant tag from the qbox sample file.  Will read in 
      the occupations from the density_matrix subnode add the total which will be
      returned as output */
  double getOccupation(const XmlNode* nd) const;

  // helper functions meant for espresso
  /*! to be handed the top level xml node from data-file-schema, also a string naming 
      the directory where the charge-density.hdf5 file is located and an integer
      spinpol which is 0 for a non spin polarized calculation and 1 for one that is */
  void handleDensity(const XmlNode& qeXml, const std::string& dir_name, int spinpol);
  /*! to be handed the top level xml node from data-file-schema, reads the 
      reciprocal lattice vectors from the xml, calculates the primitive translation
      vectors and puts them sequentially in a 9 element long std::vector */
  std::vector<double> getPtvs(const XmlNode& qeXml);
  /*! to be handed a band_structure xml node from data-file-schema as well as the 9 
      element vector of ptvs.  Then will populate the eigenvalues per kpoint, 
      occupations per kpoint, kpoints,
      and the weights and number of g-vectors for each kpoint */
  void processKPts(const XmlNode& band_structure_xml,
                   const std::vector<double>& ptvs,
                   std::vector<std::vector<double>>& eigenvals,
                   std::vector<std::vector<double>>& occupations,
                   std::vector<KPoint>& kpts,
                   std::vector<double>& weights,
                   std::vector<int>& ngvecs);
  /*! will use the occupations and weights to figure out the total number
      of up and down electrons (nup,ndn).  Also needs a flag to see if the
      calculation was spin_polarized (0 for no, 1 for yes) and another to 
      see if noncolinear spins were used (0 for no, 1 for yes)
   */
  void getNumElectrons(std::vector<std::vector<double>>& occupations,
                       std::vector<double>& weights,
                       int& nup,
                       int& ndn,
                       int spinpol,
                       int ncol);
  /*! helper function to read in the data from the wfc(kpt_num).hdf5 files located in
      directory dir_name.  Needs to know the KPoint, its eigenvalues, weight, the
      flags for spinolarized and noncolinear as well as the momap_t that handles
      the global mapping from g_vectors at individual k_points to global ones 
  */
  void handleKpt(int kpt_num,
                 const std::string& dir_name,
                 KPoint& kpt,
                 const std::vector<double>& eigenvalues,
                 double weight,
                 int spinpol,
                 int noncol,
                 const momap_t& moref);
  /*! handles reading in the g-vectors for a particular k-point so that they can be inserted 
      into a momap_t.  Needs to know the number of the kpoint (kpt_num), the directory where
      it is located (dir_name) and whether it is spin polarized (spinpol)
  */
  void readKptGvecs(int kpt_num, const std::string& dir_name, int spinpol, momap_t& morefmap);
  /*! will open the hdf file fname and will make an hdf_archive */
  qmcplusplus::hdf_archive openHdfFileForRead(const std::string& fname);

  EshdfFile(const EshdfFile& f);            // no copy constructor
  EshdfFile& operator=(const EshdfFile& f); // operator= not allowed
public:
  /*! constructor that sets up the hdf_archive data structure that will hold the eshdf
      file to be created.  Needs to be give the name for that file (hdfFileName)
  */
  EshdfFile(const std::string& hdfFileName);
  ~EshdfFile();

  /*! functions to get the boilerplate written for qbox or quantum espresso
      calls writeApplication, writeVersion, writeCreator and writeFormat
      nees to be handed the top level xml node (data-file-schema or qbox.sample)
  */
  void writeQEBoilerPlate(const XmlNode& qeXml);
  void writeQboxBoilerPlate(const XmlNode& qboxSample);

  /*! functions to write the supercell and primitive_vectors sections of the 
      eshdf file.  Needs to be handed the top level xml node (data-file-schema or qbox.sample)
  */
  void writeQESupercell(const XmlNode& qeXml);
  void writeQboxSupercell(const XmlNode& qboxSample);

  /*! functions to write the atoms section of the eshdf file
      needs to be handed the top level xml node (data-file-schema or qbox.sample) 
  */
  void writeQEAtoms(const XmlNode& qeXml);
  void writeQboxAtoms(const XmlNode& qboxSample);

  /*! functions to write the electrons section of the eshdf file
      need to be handed the top level xml node (data-file-schema or qbox.sample)
      in the case of espresso, also needs to be handed the directory name
      where the wfc.hdf5 and charge-density.hdf5 files are located
  */
  void writeQEElectrons(const XmlNode& qeXml, const std::string& dir_name);
  void writeQboxElectrons(const XmlNode& qboxSample);
};

#endif
