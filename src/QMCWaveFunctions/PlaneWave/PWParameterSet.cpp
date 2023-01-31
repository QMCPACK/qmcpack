//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWParameterSet.cpp
 * @brief Utility class to handle hdf5
 */
#include "PWParameterSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{
PWParameterSet::PWParameterSet(Communicate* comm)
    : MPIObjectBase(comm),
      hasSpin(true),
      twistIndex(0),
      numBands(0),
      Ecut(-1),
      Rcut(-1),
      BufferRadius(-1),
      BoxDup(1),
      paramTag("parameters"),
      basisTag("basis"),
      pwTag("planewaves"),
      pwMultTag("multipliers"),
      eigTag("eigenstates"),
      twistTag("twist"),
      bandTag("band"),
      spinTag("spin"),
      eigvecTag("eigenvector")
{
  m_param.setName("h5tag");
  m_param.add(twistIndex, "twistIndex");
  m_param.add(Rcut, "rcut");
  m_param.add(BufferRadius, "bufferLayer");
  m_param.add(BoxDup, "expand");
  m_param.add(paramTag, "parameters");
  m_param.add(basisTag, "basis");
  m_param.add(pwTag, "planewaves");
  m_param.add(pwMultTag, "multiplers");
  m_param.add(eigTag, "eigenstates");
  m_param.add(twistTag, "twist");
  m_param.add(bandTag, "band");
  m_param.add(spinTag, "spin");
  m_param.add(eigvecTag, "eigenvector");
}

double PWParameterSet::getEcut(double ecut)
{
  if (Ecut < 0 || Ecut >= ecut)
    Ecut = ecut;
  return Ecut;
}

bool PWParameterSet::getEigVectorType(hid_t h)
{
  int rank = 0;
  if (is_manager())
  {
    std::ostringstream oss;
    oss << "/" << eigTag << "/" << twistTag << twistIndex << "/" << bandTag << 0;
    //if(version[1]==10)
    if (hasSpin)
      oss << "/" << spinTag << 0;
    oss << "/eigenvector";
    hsize_t dimTot[4];
    hid_t dataset   = H5Dopen(h, oss.str().c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    rank            = H5Sget_simple_extent_ndims(dataspace);
    int status_n    = H5Sget_simple_extent_dims(dataspace, dimTot, NULL);
  }
  myComm->bcast(rank);
  return rank == 4;
}

bool PWParameterSet::hasComplexData(hdf_archive& h_file)
{
  int iscomplex = 0;
  // Should be the tag "/electrons/psi_r_is_complex", but the test HDF files
  //  don't have this set
#if 0
  if(is_manager())
  {
    std::ostringstream oss;
    oss << paramTag << "/complex_coefficients";
    h_file.read(iscomplex, oss.str());
  }
#endif
  myComm->bcast(iscomplex);
  return iscomplex;
}

std::string PWParameterSet::getTwistAngleName()
{
  std::ostringstream oss;
  oss << eigTag << "/" << twistTag << twistIndex << "/twist_angle";
  return oss.str();
}

std::string PWParameterSet::getTwistName() { return getTwistName(twistIndex); }

std::string PWParameterSet::getTwistName(int i)
{
  std::ostringstream oss;
  oss << twistTag << i;
  return oss.str();
}

std::string PWParameterSet::getBandName(int ib, int ispin)
{
  std::ostringstream oss;
  oss << "spin_" << ispin << "/"
      << "state_" << ib;
  return oss.str();
}

std::string PWParameterSet::getEigVectorName(const std::string& hg, int ib, int ispin)
{
  std::ostringstream oss;
  oss << hg << "/" << bandTag << ib;
  //if(version[1]==10)
  if (hasSpin)
  {
    oss << "/" << spinTag << ispin;
  }
  oss << "/eigenvector";
  return oss.str();
}

std::string PWParameterSet::getCenterName(const std::string& hg, int ib)
{
  std::ostringstream oss;
  oss << hg << "/" << bandTag << ib << "/center";
  return oss.str();
}

std::string PWParameterSet::getOriginName(const std::string& hg, int ib)
{
  std::ostringstream oss;
  oss << hg << "/" << bandTag << ib << "/origin";
  return oss.str();
}

std::string PWParameterSet::getEigVectorName(int ib, int ispin)
{
  std::ostringstream oss;
  oss << "/" << eigTag << "/" << twistTag << twistIndex << "/" << bandTag << ib;
  //if(version[1]==10)
  if (hasSpin)
  {
    oss << "/" << spinTag << ispin;
  }
  oss << "/eigenvector";
  return oss.str();
}

std::string PWParameterSet::getBandName(int ib)
{
  std::ostringstream oss;
  oss << bandTag << ib;
  return oss.str();
}

std::string PWParameterSet::getSpinName(int ispin)
{
  std::ostringstream oss;
  oss << spinTag << ispin;
  return oss.str();
}

void PWParameterSet::checkVersion(hdf_archive& h)
{
  if (is_manager())
  {
    hid_t dataset         = H5Dopen(h.getFileID(), "version", H5P_DEFAULT);
    hid_t datatype        = H5Dget_type(dataset);
    H5T_class_t classtype = H5Tget_class(datatype);
    H5Tclose(datatype);
    H5Dclose(dataset);
    if (classtype == H5T_INTEGER)
    {
      h.read(version, "version");
    }
    else if (classtype == H5T_FLOAT)
    {
      TinyVector<double, 2> vt;
      h.read(vt, "version");
      version[0] = static_cast<int>(vt[0]);
      version[1] = static_cast<int>(vt[1]);
    }
    else
    {
      APP_ABORT("PWParameterSet::checkVersion  The type of version is not integer or double.");
    }
  }
  myComm->bcast(version);
  app_log() << "\tWavefunction HDF version: " << version[0] << "." << version[1] << std::endl;
  if (version[0] == 0)
  {
    if (version[1] == 11)
    {
      hasSpin   = false;
      paramTag  = "parameters_0";
      basisTag  = "basis_1";
      pwTag     = "planewaves";
      pwMultTag = "multipliers";
      eigTag    = "eigenstates_3";
      twistTag  = "twist_";
      bandTag   = "band_";
    }
    else if (version[1] == 10)
    {
      pwMultTag = "planewaves";
      pwTag     = "0";
    }
  }
}
} // namespace qmcplusplus
