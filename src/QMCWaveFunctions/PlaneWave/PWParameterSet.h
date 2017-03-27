//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file PWParameterSet.cpp
 * @brief Utility class to handle hdf5
 */
#ifndef QMCPLUSPLUS_PWPARAMETERSET_H
#define QMCPLUSPLUS_PWPARAMETERSET_H
#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/MPIObjectBase.h"

namespace qmcplusplus
{
/** class to handle various name conventions for hdf5 file
 */
struct PWParameterSet: public MPIObjectBase
{
  ///true if spin channel exists
  bool hasSpin;
  ///version
  TinyVector<int,2> version;
  ///index of the twist angle
  int twistIndex;
  ///number of input bands
  int numBands;
  ///energy cutoff for QMC wavefunction
  double Ecut;
  ///cutoff radius for truncated orbitals
  double Rcut;
  ///radius of buffer layer for truncated orbitals
  double BufferRadius;
  ///cell multiplications
  TinyVector<int,OHMMS_DIM> BoxDup;
  ///tag for the parameters
  std::string paramTag;
  ///tag for the basis
  std::string basisTag;
  ///tag for the planewaves
  std::string pwTag;
  ///tag for the multipliers of the planewaves
  std::string pwMultTag;
  ///tag for eigentstates
  std::string eigTag;
  ///tag for twist angles
  std::string twistTag;
  ///tag for the band
  std::string bandTag;
  ///tag for the spin
  std::string spinTag;
  ///tag for eigvector
  std::string eigvecTag;
  ///xml processor
  ParameterSet m_param;

  PWParameterSet();

  bool put(xmlNodePtr cur)
  {
    return m_param.put(cur);
  }

  void checkVersion(hid_t h);

  double getEcut(double ecut);

  /** get the dimensions of the eigenvectors
   * @param h fileid
   * @return true, if the data is complex
   */
  bool getEigVectorType(hid_t h);

  bool hasComplexData(hid_t h);

  std::string getTwistAngleName();

  std::string getTwistName();

  std::string getTwistName(int i);

  std::string getBandName(int ib, int ispin);

  std::string getBandName(int ib);

  std::string getSpinName(int ispin);

  std::string getEigVectorName(const std::string& hg, int ib, int ispin);
  std::string getEigVectorName(int ib, int ispin);
  std::string getCenterName(const std::string& hg,int ib);
  std::string getOriginName(const std::string& hg,int ib);

  void writeParameters(hid_t gid);

};
}
#endif
