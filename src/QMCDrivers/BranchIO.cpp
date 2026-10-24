//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "BranchIO.h"
#include "hdf/HDFVersion.h"
#include "Message/CommOperators.h"
#include "hdf/hdf_archive.h"
#if defined(HAVE_LIBBOOST)
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>
#endif

//#include <boost/archive/text_oarchive.hpp>
//
namespace qmcplusplus
{
#if defined(HAVE_LIBBOOST)
template<typename T>
inline void put_histogram(std::string name, accumulator_set<T> a, boost::property_tree::ptree& pt)
{
  pt.put(name + ".value", a.properties[0]);
  pt.put(name + ".value_squared", a.properties[1]);
  pt.put(name + ".weight", a.properties[2]);
}

template<typename T>
inline void get_histogram(std::string name, accumulator_set<T>& a, boost::property_tree::ptree& pt)
{
  a.properties[0] = pt.get<T>(name + ".value");
  a.properties[1] = pt.get<T>(name + ".value_squared");
  a.properties[2] = pt.get<T>(name + ".weight");
}
#endif

template<typename T>
struct h5data_proxy<accumulator_set<T>> : public h5_space_type<T, 1>
{
  enum
  {
    CAPACITY = accumulator_set<T>::CAPACITY
  };
  using h5_space_type<T, 1>::dims;
  using h5_space_type<T, 1>::get_address;
  using data_type = accumulator_set<T>;

  inline h5data_proxy(const data_type& a) { dims[0] = CAPACITY; }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  { return h5d_read(grp, aname, get_address(ref.properties), xfer_plist); }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  { return h5d_write(grp, aname.c_str(), this->size(), dims, get_address(ref.properties), xfer_plist); }
};

template<class SFNB>
std::vector<std::string> BranchIO<SFNB>::vParamName;
template<class SFNB>
std::vector<std::string> BranchIO<SFNB>::iParamName;

template<class SFNB>
void BranchIO<SFNB>::initAttributes()
{
  if (vParamName.size())
    return;
  vParamName.resize(10);
  vParamName[0] = "tau";
  vParamName[1] = "taueff";
  vParamName[2] = "etrial";
  vParamName[3] = "eref";
  vParamName[4] = "branchmax";
  vParamName[5] = "branchcutoff";
  vParamName[6] = "branchfilter";
  vParamName[7] = "sigma";
  vParamName[8] = "acc_energy";
  vParamName[9] = "acc_samples";

  iParamName.resize(7);
  iParamName[0] = "warmupsteps";
  iParamName[1] = "energyupdateinterval";
  iParamName[2] = "counter";
  iParamName[3] = "targetwalkers";
  iParamName[4] = "maxwalkers";
  iParamName[5] = "minwalkers";
  iParamName[6] = "brnachinterval";
}

template<class SFNB>
bool BranchIO<SFNB>::write(const std::string& fname)
{
  if (myComm.rank())
    return true;

#if defined(HAVE_LIBBOOST)
  initAttributes();
  using boost::property_tree::ptree;
  ptree pt;

  // Put log filename in property tree
  pt.put("state.branchmode", ref.branch_mode_);

  auto it_vparam = ref.v_param_.begin();
  for (auto& name : vParamName)
    pt.put("state.vparam." + name, *(it_vparam++));
  for (int i = 0; i < iParamName.size(); ++i)
    pt.put("state.iparam." + iParamName[i], ref.iParam[i]);

  put_histogram("state.energy", ref.energy_hist_, pt);
  put_histogram("state.variance", ref.variance_hist_, pt);
  put_histogram("state.r2accepted", ref.r2_accepted_, pt);
  put_histogram("state.r2proposed", ref.r2_proposed_, pt);
  std::string xname = fname + ".qmc.xml";
  write_xml(xname, pt);
#else
  //append .qmc.h5 if missing
  std::string h5name(fname);
  if (fname.find("qmc.h5") >= fname.size())
    h5name.append(".qmc.h5");
  hdf_archive dump(myComm);
  hid_t fid = dump.create(h5name);
  dump.push(hdf::main_state);
  dump.push(hdf::qmc_status);
  std::string v_header("tau:taueff:etrial:eref:branchmax:branchcutoff:branchfilter:sigma:acc_energy:acc_samples");
  std::string i_header("warmupsteps:energyupdateinterval:counter:targetwalkers:maxwalkers:minwalkers:branchinterval");
  dump.write(v_header, "vparam_def");
  dump.write(i_header, "iparam_def");
  dump.write(ref.v_param_, "vparam");
  dump.write(ref.iParam, "iparam");
  dump.write(ref.branch_mode_, "branchmode");
  dump.push("histogram");
  dump.write(ref.energy_hist_, "energy");
  dump.write(ref.variance_hist_, "variance");
  dump.write(ref.r2_accepted_, "r2accepted");
  dump.write(ref.r2_proposed_, "r2proposed");
  //PopHist is not being used in 2010-10-19
  //if(ref.branch_mode_[SimpleFixedNodeBranch::B_DMC])
  //{
  //  dump.push("population");
  //  dump.write(ref.PopHist.myData,"histogram");
  //}
#endif
  return true;
}

template<class SFNB>
bool BranchIO<SFNB>::read(const std::string& fname)
{
  int found_config = 0;

  if (myComm.rank() == 0)
  {
    initAttributes();
    using boost::property_tree::ptree;
    ptree pt;
    std::string xname = fname + ".qmc.xml";
    read_xml(xname, pt);
    if (!pt.empty())
    {
      ref.branch_mode_ = pt.get<BranchModeType>("state.branchmode");

      get_histogram("state.energy", ref.energy_hist_, pt);
      get_histogram("state.variance", ref.variance_hist_, pt);
      get_histogram("state.r2accepted", ref.r2_accepted_, pt);
      get_histogram("state.r2proposed", ref.r2_proposed_, pt);

      int i = 0;
      BOOST_FOREACH (const ptree::value_type& v, pt.get_child("state.iparam"))
      {
        ref.iParam[i++] = v.second.get_value<int>();
      }
      auto it_vparam = ref.v_param_.begin();
      BOOST_FOREACH (const ptree::value_type& v, pt.get_child("state.vparam"))
      {
        *(it_vparam++) = v.second.get_value<double>();
      }
      found_config = 1;
    }
  }
  myComm.bcast(found_config);

  if (!found_config)
    return false;

  bcast_state();

  return true;
}

template<class SFNB>
void BranchIO<SFNB>::bcast_state()
{
  int n = ref.v_param_.size() + ref.iParam.size();
  std::vector<RealType> pdata(n + 1 + 16, -1);

  if (myComm.rank() == 0)
  {
    copy(ref.v_param_.begin(), ref.v_param_.end(), pdata.begin());
    copy(ref.iParam.begin(), ref.iParam.end(), pdata.begin() + ref.v_param_.size());
    int offset      = n;
    pdata[offset++] = ref.branch_mode_.to_ulong();
    copy(ref.energy_hist_.properties, ref.energy_hist_.properties + 4, pdata.begin() + offset);
    offset += 4;
    copy(ref.variance_hist_.properties, ref.variance_hist_.properties + 4, pdata.begin() + offset);
    offset += 4;
    copy(ref.r2_accepted_.properties, ref.r2_accepted_.properties + 4, pdata.begin() + offset);
    offset += 4;
    copy(ref.r2_proposed_.properties, ref.r2_proposed_.properties + 4, pdata.begin() + offset);
  }

  //broadcast to the nodes : need to add a namespace mpi::
  myComm.bcast(pdata);

  if (myComm.rank())
  {
    int ii = 0;
    for (auto& vpar : ref.v_param_)
      vpar = pdata[ii++];
    for (int i = 0; i < ref.iParam.size(); ++i, ++ii)
      ref.iParam[i] = static_cast<int>(pdata[ii]);
    ref.branch_mode_ = static_cast<unsigned long>(pdata[ii]);
  }
  {
    //update historgram
    int ii = n + 1;
    ref.energy_hist_.reset(pdata[ii], pdata[ii + 1], pdata[ii + 2]);
    ii += 4;
    ref.variance_hist_.reset(pdata[ii], pdata[ii + 1], pdata[ii + 2]);
    ii += 4;
    ref.r2_accepted_.reset(pdata[ii], pdata[ii + 1], pdata[ii + 2]);
    ii += 4;
    ref.r2_proposed_.reset(pdata[ii], pdata[ii + 1], pdata[ii + 2]);
  }
}

template class BranchIO<SimpleFixedNodeBranch>;
} // namespace qmcplusplus
