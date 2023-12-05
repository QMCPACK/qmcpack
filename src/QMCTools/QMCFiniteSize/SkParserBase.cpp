#include "SkParserBase.h"
#include "Message/Communicate.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"

namespace qmcplusplus
{
SkParserBase::SkParserBase() : isParseSuccess(false), isGridComputed(false), isSkComputed(false), skname("SkAll")
{
  skraw.resize(0);
  sk.resize(0);
  skerr.resize(0);
  skerr_raw.resize(0);
  kgridraw.resize(0);
  kgrid.resize(0);
}

void SkParserBase::compute_grid()
{
  if (!isParseSuccess)
    APP_ABORT("SkParserBase::compute_grid(..) : Initial parse failed");

  //  cout<<" We're about to get grid info...\n";
  RealType lx(0), rx(0);
  RealType ly(0), ry(0);
  RealType lz(0), rz(0);

  RealType dx(0), dy(0), dz(0);
  IndexType Nx(0), Ny(0), Nz(0);
  get_gridinfo_from_posgrid(kgridraw, 0, lx, rx, dx, Nx);
  get_gridinfo_from_posgrid(kgridraw, 1, ly, ry, dy, Ny);
  get_gridinfo_from_posgrid(kgridraw, 2, lz, rz, dz, Nz);

  // cout<<" Done with grid info...\n";
  kgrid.resize(Nx * Ny * Nz);

  xgrid.set(lx, rx, Nx);
  ygrid.set(ly, ry, Ny);
  zgrid.set(lz, rz, Nz);


  isGridComputed = true;
}

void SkParserBase::set_grid(const std::vector<TinyVector<IndexType, OHMMS_DIM>>& kgridraw1)
{
  if (skraw.size() != kgridraw1.size())
    APP_ABORT("SkParserBase::set_grid:  S(k) and k-grid don't match");
  kgridraw.resize(kgridraw1.size());
  for (IndexType i = 0; i < kgridraw.size(); i++)
    for (IndexType j = 0; j < OHMMS_DIM; j++)
      kgridraw[i][j] = RealType(kgridraw1[i][j]);
  compute_grid();
}

void SkParserBase::set_grid(const std::vector<PosType>& kgridraw1)
{
  if (skraw.size() != kgridraw1.size())
    APP_ABORT("SkParserBase::set_grid:  S(k) and k-grid don't match");
  kgridraw = kgridraw1;
  compute_grid();
}

void SkParserBase::get_grid(Grid_t& xgrid_i, Grid_t& ygrid_i, Grid_t& zgrid_i)
{
  //  cout<<"In get_grid(..)\n";
  if (!isGridComputed)
    compute_grid();
  //  cout<<"done with compute_grid()\n";
  xgrid_i.set(xgrid.rmin(), xgrid.rmax(), xgrid.size());
  ygrid_i.set(ygrid.rmin(), ygrid.rmax(), ygrid.size());
  zgrid_i.set(zgrid.rmin(), zgrid.rmax(), zgrid.size());
}

void SkParserBase::compute_sk()
{
  if (!isParseSuccess)
    APP_ABORT("SkParserBase::compute_sk() : Initial parse failed");
  //  cout<<"In compute_sk()\n";

  if (kgridraw.size() != skraw.size())
    APP_ABORT("SkParserBase::compute_sk() : Kgrid and SK not the same size");

  if (!isGridComputed)
    compute_grid();

  IndexType nx(0), ny(0), nz(0);
  IndexType Nx(0), Ny(0), Nz(0);
  IndexType newindex(0);

  Nx = xgrid.size();
  Ny = ygrid.size();
  Nz = zgrid.size();

  sk.resize(Nx * Ny * Nz);
  skerr.resize(Nx * Ny * Nz);

  //set k=(0,0,0), S(0)=0, Serr(0)=0

  nx = xgrid.getIndex(0);
  ny = ygrid.getIndex(0);
  nz = zgrid.getIndex(0);

  newindex = nx * Ny * Nz + ny * Nz + nz;

  kgrid[newindex] = 0;
  sk[newindex]    = 0.0;
  skerr[newindex] = 0.0;

  for (IndexType i = 0; i < kgridraw.size(); i++)
  {
    nx = xgrid.getIndex(kgridraw[i][0]);
    ny = ygrid.getIndex(kgridraw[i][1]);
    nz = zgrid.getIndex(kgridraw[i][2]);

    newindex        = nx * Ny * Nz + ny * Nz + nz;
    kgrid[newindex] = kgridraw[i];
    sk[newindex]    = skraw[i];
    skerr[newindex] = skerr_raw[i];
  }

  isSkComputed = true;
}

void SkParserBase::get_sk(std::vector<RealType>& sk_i, std::vector<RealType>& skerr_i)
{
  //  cout<<"In get_sk(..)\n";
  if (!isSkComputed)
    compute_sk();

  sk_i    = sk;
  skerr_i = skerr;
}

} // namespace qmcplusplus
