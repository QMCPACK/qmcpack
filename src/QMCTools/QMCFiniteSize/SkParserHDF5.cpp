#include "QMCTools/QMCFiniteSize/SkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"
#include <cmath>

namespace qmcplusplus
{
SkParserHDF5::SkParserHDF5(){};
SkParserHDF5::~SkParserHDF5(){};

void SkParserHDF5::parse(const string& fname)
{
  bool result = statfile.open(fname);
  if (!result)
  {
    cout << "SkParserHDF5::parse could not open " << fname << std::endl;
    exit(1);
  }

  vector<int> readShape;
  statfile.getShape<int>("skall/kpoints/value", readShape);
  assert(readShape[1] == 3);
  array<int, 2> kdims{readShape[0], readShape[1]};
  vector<RealType> ktmp;
  statfile.readSlabReshaped(ktmp, kdims, "skall/kpoints/value");
  for (int ik = 0; ik < readShape[0]; ik++)
  {
    PosType k;
    k[0] = ktmp[3 * ik];
    k[1] = ktmp[3 * ik + 1];
    k[2] = ktmp[3 * ik + 2];
    kgridraw.push_back(k);
  }
  int nKpts = kgridraw.size();

  statfile.getShape<int>("skall/rhok_e_e/value", readShape);
  assert(readShape[1] == nKpts);
  vector<RealType> rhok_e_tmp;
  array<int, 2> rhok_e_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_e_tmp, rhok_e_dims, "skall/rhok_e_e/value");
  int nBlocks = readShape[0];

  vector<RealType> rhok_ee;
  vector<RealType> rhok_ee_err;
  for (int ik = 0; ik < nKpts; ik++)
  {
    vector<RealType> block_data;
    for (int ib = 0; ib < nBlocks; ib++)
      block_data.push_back(rhok_e_tmp[nKpts * ib + ik]);
    RealType avg, err;
    getStats(block_data, avg, err);
    rhok_ee.push_back(avg);
    rhok_ee_err.push_back(err);
  }

  statfile.getShape<int>("skall/rhok_e_i/value", readShape);
  vector<RealType> rhok_i_tmp;
  array<int, 2> rhok_i_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_i_tmp, rhok_i_dims, "skall/rhok_e_i/value");
  assert(readShape[1] == nKpts);
  assert(readShape[1] == nBlocks);

  vector<RealType> rhok_ei;
  vector<RealType> rhok_ei_err;
  for (int ik = 0; ik < nKpts; ik++)
  {
    vector<RealType> block_data;
    for (int ib = 0; ib < nBlocks; ib++)
      block_data.push_back(rhok_i_tmp[nKpts * ib + ik]);
    RealType avg, err;
    getStats(block_data, avg, err);
    rhok_ei.push_back(avg);
    rhok_ei_err.push_back(err);
  }

  statfile.getShape<int>("skall/rhok_e_r/value", readShape);
  vector<RealType> rhok_r_tmp;
  array<int, 2> rhok_r_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_r_tmp, rhok_r_dims, "skall/rhok_e_r/value");
  assert(readShape[1] == nKpts);
  assert(readShape[1] == nBlocks);

  vector<RealType> rhok_er;
  vector<RealType> rhok_er_err;
  for (int ik = 0; ik < nKpts; ik++)
  {
    vector<RealType> block_data;
    for (int ib = 0; ib < nBlocks; ib++)
      block_data.push_back(rhok_r_tmp[nKpts * ib + ik]);
    RealType avg, err;
    getStats(block_data, avg, err);
    rhok_er.push_back(avg);
    rhok_er_err.push_back(err);
  }

  for (int ik = 0; ik < nKpts; ik++)
  {
    RealType r_r(0.0);
    RealType r_i(0.0);
    //The 3.0 in the following statements is the standard deviation.
    //  If function value is greater than 3standard deviations above error,
    //  then we set it.  Otherwise default to zero.
    if (rhok_er_err[ik] == 0 || std::abs(rhok_er[ik]) / rhok_er_err[ik] > 3.0)
      r_r = rhok_er[ik];
    if (rhok_ei_err[ik] == 0 || std::abs(rhok_ei[ik]) / rhok_ei_err[ik] > 3.0)
      r_i = rhok_ei[ik];
    skraw.push_back(rhok_ee[ik] - (r_r * r_r + r_i * r_i));
    skerr_raw.push_back(rhok_ee_err[ik]);
  }

  hasGrid        = false;
  isNormalized   = false;
  isParseSuccess = true;

  statfile.close();
}

} // namespace qmcplusplus
