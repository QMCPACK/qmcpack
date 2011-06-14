//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////
/** helper functions for EinsplineSetBuilder
 */
#ifndef QMCPLUSPLUS_EINSPLINEBUILDER_HELPER_H
#define QMCPLUSPLUS_EINSPLINEBUILDER_HELPER_H
#include <complex>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsArray.h>
namespace qmcplusplus
{

  /** unpack packed cG to fftbox
   * @param cG packed vector
   * @param gvecs g-coordinate for cG[i]
   * @param maxg  fft grid
   * @param fftbox unpacked data to be transformed
   */
  template<typename T> 
    inline void unpack4fftw(const Vector<std::complex<T> >& cG
        , const vector<TinyVector<int,3> >& gvecs
        , const TinyVector<int,3>& maxg
        , Array<std::complex<T>,3>& fftbox
        )
  {
    fftbox=std::complex<T>();
    //this is rather unsafe
#pragma omp parallel for
    for (int iG=0; iG<cG.size(); iG++) 
    {
      fftbox((gvecs[iG][0]+maxg[0])%maxg[0]
          ,(gvecs[iG][1]+maxg[1])%maxg[1]
          ,(gvecs[iG][2]+maxg[2])%maxg[2]) = cG[iG];
    }
  }

  /** fix phase
   * @param in input complex data
   * @param out output real data
   * @param twist vector to correct phase
   */
  template<typename T>
    inline void fix_phase_c2r(const Array<std::complex<T>,3>& in, Array<T,3>& out
        , const TinyVector<T,3>& twist
        )
  {
    const T two_pi=-2.0*M_PI;
    const int nx=in.size(0);
    const int ny=in.size(1);
    const int nz=in.size(2);
    T nx_i=1.0/static_cast<T>(nx);
    T ny_i=1.0/static_cast<T>(ny);
    T nz_i=1.0/static_cast<T>(nz);

#pragma omp parallel for firstprivate(nx_i,ny_i,nz_i)
    for (int ix=0; ix<nx; ix++) 
    {
      T s, c;
      const std::complex<T>* restrict in_ptr=in.data()+ix*ny*nz;
      T* restrict out_ptr=out.data()+ix*ny*nz;

      T rux=static_cast<T>(ix)*nx_i*twist[0];
      for (int iy=0; iy<ny; iy++) 
      {
        T ruy=static_cast<T>(iy)*ny_i*twist[1];
        for (int iz=0; iz<nz; iz++) 
        {
          T ruz=static_cast<T>(iz)*nz_i*twist[2];
          sincos(two_pi*(rux+ruy+ruz), &s, &c);
          *out_ptr = c*in_ptr->real()-s*in_ptr->imag();
          ++out_ptr;
          ++in_ptr;
        }
      }
    }
//#pragma omp parallel for
//          for (int ix=0; ix<nx; ix++) {
//            PosType ru;
//            ru[0] = (RealType)ix / (RealType)nx;
//            for (int iy=0; iy<ny; iy++) {
//              ru[1] = (RealType)iy / (RealType)ny;
//              for (int iz=0; iz<nz; iz++) {
//                ru[2] = (RealType)iz / (RealType)nz;
//                double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
//                double s, c;
//                sincos(phi, &s, &c);
//                complex<double> phase(c,s);
//                complex<double> z = phase*rawData(ix,iy,iz);
//                splineData(ix,iy,iz) = z.real();
//              }
//            }
//          }
  }

  /** rotate the state after 3dfft
   *
   * First, add the eikr phase factor.
   * Then, rotate the phase of the orbitals so that neither
   * the real part nor the imaginary part are very near
   * zero.  This sometimes happens in crystals with high
   * symmetry at special k-points.
   */
  template<typename T>
    inline void fix_phase_rotate(Array<std::complex<T>,3>& in
    , Array<T,3>& out, const TinyVector<T,3>& twist
    )
    {
      const T two_pi=-2.0*M_PI;
      const int nx=in.size(0);
      const int ny=in.size(1);
      const int nz=in.size(2);
      T nx_i=1.0/static_cast<T>(nx);
      T ny_i=1.0/static_cast<T>(ny);
      T nz_i=1.0/static_cast<T>(nz);

      T rNorm=0.0, iNorm=0.0;
#pragma omp parallel for reduction(+:rNorm,iNorm), firstprivate(nx_i,ny_i,nz_i)
      for (int ix=0; ix<nx; ix++) 
      {
        T s, c;
        std::complex<T>* restrict in_ptr=in.data()+ix*ny*nz;
        T rux=static_cast<T>(ix)*nx_i*twist[0];
        for (int iy=0; iy<ny; iy++)
        {
          T ruy=static_cast<T>(iy)*ny_i*twist[1];
          for (int iz=0; iz<nz; iz++) 
          {
            T ruz=static_cast<T>(iz)*nz_i*twist[2];
            sincos(two_pi*(rux+ruy+ruz), &s, &c);
            std::complex<T> eikr(c,s);
            *in_ptr *= eikr;
            rNorm += in_ptr->real()*in_ptr->real();
            iNorm += in_ptr->imag()*in_ptr->imag();
            ++in_ptr;
          }
        }
      }

      T arg = std::atan2(iNorm, rNorm);
      T phase_i,phase_r;
      sincos(0.125*M_PI-0.5*arg, &phase_i, &phase_r);
#pragma omp parallel for firstprivate(phase_r,phase_i)
      for (int ix=0; ix<nx; ix++)
      {
        const std::complex<T>* restrict in_ptr=in.data()+ix*ny*nz;
        T* restrict out_ptr=out.data()+ix*ny*nz;
        for (int iy=0; iy<ny; iy++)
          for (int iz=0; iz<nz; iz++) 
          {
            *out_ptr=phase_r*in_ptr->real()-phase_i*in_ptr->imag();
            ++in_ptr;
            ++out_ptr;
          }
      }
    }

  template<typename T>
  inline void fix_phase_rotate(const Array<std::complex<T>,3>& in
      , Array<std::complex<T>,3>& out, const TinyVector<T,3>& twist)
  {
    const T two_pi=-2.0*M_PI;
    const int nx=in.size(0);
    const int ny=in.size(1);
    const int nz=in.size(2);
    T nx_i=1.0/static_cast<T>(nx);
    T ny_i=1.0/static_cast<T>(ny);
    T nz_i=1.0/static_cast<T>(nz);

    T rNorm=0.0, iNorm=0.0;
//#pragma omp parallel for reduction(+:rNorm,iNorm), firstprivate(nx_i,ny_i,nz_i)
    for (int ix=0; ix<nx; ix++) 
    {
      T s, c, r, i;
      const std::complex<T>* restrict in_ptr=in.data()+ix*ny*nz;
      T rux=static_cast<T>(ix)*nx_i*twist[0];
      T rsum=0, isum=0;
      for (int iy=0; iy<ny; iy++)
      {
        T ruy=static_cast<T>(iy)*ny_i*twist[1];
        for (int iz=0; iz<nz; iz++) 
        {
          T ruz=static_cast<T>(iz)*nz_i*twist[2];
          sincos(two_pi*(rux+ruy+ruz), &s, &c);
          r = c*in_ptr->real()-s*in_ptr->imag();
          i = s*in_ptr->real()+c*in_ptr->imag();
          ++in_ptr;
          rsum += r*r;
          isum += i*i;
        }
      }
      rNorm += rsum;
      iNorm += isum;
    }

    T arg = std::atan2(iNorm, rNorm);
    //cerr << "Phase = " << arg/M_PI << " pi.\n";
    T phase_r, phase_i;
    sincos(0.5*(0.25*M_PI-arg), &phase_i, &phase_r);

//#pragma omp parallel for firstprivate(phase_r,phase_i)
    for (int ix=0; ix<nx; ix++)
    {
      const std::complex<T>* restrict in_ptr=in.data()+ix*ny*nz;
      std::complex<T>* restrict out_ptr=out.data()+ix*ny*nz;
      for (int iy=0; iy<ny; iy++)
        for (int iz=0; iz<nz; iz++) 
        {
          *out_ptr=std::complex<T>(
              phase_r*in_ptr->real()-phase_i*in_ptr->imag()
              , phase_i*in_ptr->real()+phase_r*in_ptr->imag());
          ++out_ptr;
          ++in_ptr;
        }
    }
  }

  inline bool EinsplineSetBuilder::bcastSortBands(int n, bool root)
  {
    PooledData<RealType> misc(n*3+1);
    bool isCore=false;
    if(root)
    {
      misc.rewind();
      for(int i=0; i<n; ++i)
      {
        int ti   = SortBands[i].TwistIndex;
        int bi   = SortBands[i].BandIndex;
        double e = SortBands[i].Energy;
        //PosType k= primLattice.k_cart(TwistAngles[ti]);
        PosType k= PrimCell.k_cart(TwistAngles[ti]);
        fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n", 
            ti, bi, e, k[0], k[1], k[2], myComm->rank());
        misc.put(ti);
	misc.put(bi);
	misc.put(e);
        isCore |= SortBands[i].IsCoreState;
      }
      misc.put(isCore);
    }
    myComm->bcast(misc);
    if(!root)
    {
      SortBands.resize(n);
      misc.rewind();
      for(int i=0; i<n; ++i)
      {
        misc.get(SortBands[i].TwistIndex);
	misc.get(SortBands[i].BandIndex);
	misc.get(SortBands[i].Energy);
      }
      misc.get(isCore);
    }

    return isCore;
  }

  bool EinsplineSetBuilder::ReadGvectors_ESHDF()
  {
    bool root=myComm->rank() ==0;

    //this is always ugly
    MeshSize = 0; 
    if(root)
    {
      HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize);
      h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
      h_mesh.read (H5FileID, "/electrons/mesh");
    }
    myComm->bcast(MeshSize);

    //check if psi_r exists
    bool hasPsig = (MeshSize[0] == 0);
    //if(root)
    //{
    //  hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_r");
    //  if(gid<0)
    //    hasPsig=true;
    //  else
    //    H5Dclose(gid);
    //}
    //myComm->bcast(hasPsig);

    if(hasPsig)
    {
      int numk=0;
      if(root)
      {
        HDFAttribIO<int> h_numk(numk);
        h_numk.read (H5FileID, "/electrons/number_of_kpoints");
      }
      myComm->bcast(numk);

      //    std::set<TinyVector<int,3> > Gset;
      // Read k-points for all G-vectors and take the union
      TinyVector<int,3> maxIndex(0,0,0);
      Gvecs.resize(numk);
      for (int ik=0; ik<numk; ik++) 
      {
        // ostringstream numGpath;
        // HDFAttribIO<int> h_numg(numG);
        // int numG=0;
        // numGpath << "/electrons/kpoint_" << ik << "/number_of_gvectors";
        int numg=0;
        if(root)
        {
          ostringstream Gpath;
          Gpath    << "/electrons/kpoint_" << ik << "/gvectors";
          HDFAttribIO<vector<TinyVector<int,3> > > h_Gvecs(Gvecs[ik]);
          h_Gvecs.read (H5FileID, Gpath.str().c_str());
          numg=Gvecs[ik].size();
        }
        myComm->bcast(numg);

        if(!root) Gvecs[ik].resize(numg);
        myComm->bcast(Gvecs[ik]);

        for (int ig=0; ig<Gvecs[ik].size(); ig++) 
        {
          maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[ik][ig][0]));
          maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[ik][ig][1]));
          maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[ik][ig][2]));
        }
        // for (int ig=0; ig<Gvecs.size(); ig++)
        // 	if (Gset.find(Gvecs[ig]) == Gset.end())
        // 	  Gset.insert(Gvecs[ig]);
      }
      MeshSize[0] = (int)std::ceil(4.0*MeshFactor*maxIndex[0]);
      MeshSize[1] = (int)std::ceil(4.0*MeshFactor*maxIndex[1]);
      MeshSize[2] = (int)std::ceil(4.0*MeshFactor*maxIndex[2]);
    }
    app_log() << "B-spline mesh factor is " << MeshFactor << endl;
    app_log() << "B-spline mesh size is (" << MeshSize[0] << ", "
      << MeshSize[1] << ", " << MeshSize[2] << ")\n";

    return hasPsig;
  }

}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5243 $   $Date: 2011-05-23 09:09:34 -0500 (Mon, 23 May 2011) $
 * $Id: einspine_helper.hpp 5243 2011-05-23 14:09:34Z jeongnim.kim $
 ***************************************************************************/
