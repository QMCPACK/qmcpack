#include "SkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"
#include <cmath>

namespace qmcplusplus
{
void SkParserHDF5::parse(const string& fname)
{
  bool result = statfile.open(fname);
  if (!result)
  {
    cout << "SkParserHDF5::parse could not open " << fname << std::endl;
    exit(1);
  }

  //read the kpoints
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

  //read the <rho_-k*rho_k> term
  statfile.getShape<int>("skall/rhok_e_e/value", readShape);
  assert(readShape[1] == nKpts);
  vector<RealType> rhok_e_tmp;
  array<int, 2> rhok_e_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_e_tmp, rhok_e_dims, "skall/rhok_e_e/value");
  int nBlocks = readShape[0];

  //read the Im(rho_k) term
  statfile.getShape<int>("skall/rhok_e_i/value", readShape);
  vector<RealType> rhok_i_tmp;
  array<int, 2> rhok_i_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_i_tmp, rhok_i_dims, "skall/rhok_e_i/value");
  assert(readShape[1] == nKpts);
  assert(readShape[1] == nBlocks);

  //read the Re(rho_k)
  statfile.getShape<int>("skall/rhok_e_r/value", readShape);
  vector<RealType> rhok_r_tmp;
  array<int, 2> rhok_r_dims{readShape[0], readShape[1]};
  statfile.readSlabReshaped(rhok_r_tmp, rhok_r_dims, "skall/rhok_e_r/value");
  assert(readShape[1] == nKpts);
  assert(readShape[1] == nBlocks);

  //For each k, evaluate the flucuating S(k) for all blocks.
  //Then perform a simple equilibration estimate for this particular S(k) value
  //Store the  average and error after throwing out the equilibration
  for (int ik = 0; ik < nKpts; ik++)
  {
    vector<RealType> block_data;
    for (int ib = 0; ib < nBlocks; ib++)
    {
      RealType re, rr, ri;
      re = rhok_e_tmp[nKpts * ib + ik];
      rr = rhok_r_tmp[nKpts * ib + ik];
      ri = rhok_i_tmp[nKpts * ib + ik];
      block_data.push_back(re - (rr * rr + ri * ri));
    }
    int ieq = estimateEquilibration(block_data);
    RealType avg, err;
    getStats(block_data, avg, err, ieq);
    skraw.push_back(avg);
    skerr_raw.push_back(err);
  }

  hasGrid        = false;
  isNormalized   = false;
  isParseSuccess = true;

  statfile.close();
}

} // namespace qmcplusplus
