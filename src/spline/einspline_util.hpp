//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file einspline_util.hpp
 * @brief utility functions for i/o and bcast of einspline objects
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_UTILITIES_H
#define QMCPLUSPLUS_EINSPLINE_UTILITIES_H

#include <Message/CommOperators.h>
#include <OhmmsData/FileUtility.h>
#include <io/hdf_archive.h>
#include <einspline/multi_bspline_copy.h>

namespace qmcplusplus
{
  ///handles i/o and bcast, testing for now
  template<typename T>
  inline void chunked_bcast(Communicate* comm, T* buffer, size_t ntot)
  {
    if(comm->size()==1) return;

    size_t chunk_size=(1<<30)/sizeof(T); //256 MB
    int n=static_cast<int>(ntot/chunk_size); 

    size_t offset=0;
    for(int i=0; i<n; ++i, offset+=chunk_size)
    {
      comm->bcast(buffer+offset,static_cast<int>(chunk_size));
    }

    if(offset<ntot)
    {
      comm->bcast(buffer+offset,static_cast<int>(ntot-offset));
    }
  }

  template<typename ENGT>
  inline void chunked_bcast(Communicate* comm, ENGT* buffer)
  {
    chunked_bcast(comm,buffer->coefs, buffer->coefs_size);
  }


  /** specialization of h5data_proxy for einspline_engine
   */
  template<typename ENGT>
    struct h5data_proxy<einspline_engine<ENGT> >
    : public h5_space_type<typename einspline_engine<ENGT>::value_type,4>
  {
    typedef typename einspline_engine<ENGT>::value_type value_type;
    using h5_space_type<value_type,4>::dims;
    using h5_space_type<value_type,4>::get_address;
    typedef einspline_engine<ENGT> data_type;

    data_type& ref_;

    inline h5data_proxy(data_type& a): ref_(a)
    {
      dims[0]=a.spliner->x_grid.num+3;
      dims[1]=a.spliner->y_grid.num+3;
      dims[2]=a.spliner->z_grid.num+3;
      dims[3]=a.spliner->z_stride;
    }

    inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      if(ref_.spliner) 
        return h5d_read(grp,aname,get_address(ref_.spliner->coefs),xfer_plist);
      else
        return false;
    }

    inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.spliner->coefs),xfer_plist);
    }
  };

  inline string make_spline_filename(const string& old, int spin, int twist, const TinyVector<int,3>& mesh)
  {
    string aname(old);
    if(getExtension(aname) == "h5")
    {
      aname.erase(aname.end()-3,aname.end());
    }
    ostringstream oo;
    oo<<".spin_"<< spin << ".tw" << twist <<".g"<<mesh[0]<<"x"<<mesh[1]<<"x"<<mesh[2]<<".h5";
    aname+=oo.str();
    return aname;
  }

  inline string make_spline_filename(const string& old, const Tensor<int,3>& tilematrix,int spin, int twist, const TinyVector<int,3>& mesh)
  {
    string aname(old);
    if(getExtension(aname) == "h5")
    {
      aname.erase(aname.end()-3,aname.end());
    }
    ostringstream oo;
    oo<<".tile_" 
      << tilematrix(0,0) <<tilematrix(0,1) <<tilematrix(0,2)
      << tilematrix(1,0) <<tilematrix(1,1) <<tilematrix(1,2)
      << tilematrix(2,0) <<tilematrix(2,1) <<tilematrix(2,2)
      << ".spin_"<< spin << ".tw" << twist <<".g"<<mesh[0]<<"x"<<mesh[1]<<"x"<<mesh[2]<<".h5";
    aname+=oo.str();
    return aname;
  }

  inline string make_spline_filename(const string& old, const Tensor<int,3>& tilematrix,int spin, int twist,
      int bg, const TinyVector<int,3>& mesh)
  {
    string aname(old);
    if(getExtension(aname) == "h5")
    {
      aname.erase(aname.end()-3,aname.end());
    }
    ostringstream oo;
    oo<<".tile_" 
      << tilematrix(0,0) <<tilematrix(0,1) <<tilematrix(0,2)
      << tilematrix(1,0) <<tilematrix(1,1) <<tilematrix(1,2)
      << tilematrix(2,0) <<tilematrix(2,1) <<tilematrix(2,2)
      << ".spin_"<< spin << ".tw" << twist 
      << ".b" << bg <<".g"<<mesh[0]<<"x"<<mesh[1]<<"x"<<mesh[2]<<".h5";
    aname+=oo.str();
    return aname;
  }

  template<typename T>
  string make_spline_filename(const string& old,int spin, const TinyVector<T,3>& twistangle, const TinyVector<int,3>& mesh)
  {
    string aname(old);
    if(getExtension(aname) == "h5")
    {
      aname.erase(aname.end()-3,aname.end());
    }

    char a[128];
    sprintf(a,".spin_%d.k%5.3f_%5.3f_%5.3f.g%dx%dx%d.h5",
        spin,twistangle[0],twistangle[1],twistangle[2],mesh[0],mesh[1],mesh[2]);
    aname+=a;
    return aname;
  }

  template<typename GT>
    void print_grid(GT& grid)
    {
      std::cout << grid.start << " " << grid.end << " " 
        << grid.num << " " 
        << grid.delta <<  " " << grid.delta_inv << std::endl; 
    }

  template<typename ENGT>
    void print_spliner(ENGT* spline)
    {
      cout << "xgrid       "; print_grid(spline->x_grid);
      cout << "ygrid       "; print_grid(spline->y_grid);
      cout << "zgrid       "; print_grid(spline->z_grid);
    }


  /** engine to create a einspline based on the input einspline
   */
  template <typename T>
  struct GridConvert
  {
    ///number of points in each direction including BC
    int N[3];
    ///offset
    int Offset[3];
    ///number of points of the original grid
    int BaseN[3];
    ///offset of the original grid, always 0
    int BaseOffset[3];
    ///spacing to be added or removed from the buffer
    TinyVector<T,3> Delta;

    template<typename ENGT1, typename ENGT2, typename PT>
      void create(ENGT1*& out, ENGT2* in, PT& lower, PT& upper, int num)
      {

        typedef typename bspline_engine_traits<ENGT1>::real_type real_type;

        Ugrid agrid[3];
        agrid[0]=in->x_grid;
        agrid[1]=in->y_grid;
        agrid[2]=in->z_grid;

        typename bspline_engine_traits<ENGT1>::BCType xyz_bc[3];
        xyz_bc[0].lCode=in->xBC.lCode;xyz_bc[0].rCode=in->xBC.rCode;
        xyz_bc[1].lCode=in->yBC.lCode;xyz_bc[1].rCode=in->yBC.rCode;
        xyz_bc[2].lCode=in->zBC.lCode;xyz_bc[2].rCode=in->zBC.rCode;
        
        for(int i=0; i<3; ++i)
        {
          int ngi=(int)(lower[i]*agrid[i].delta_inv);
          int ngf=(int)(upper[i]*agrid[i].delta_inv)+1;
          agrid[i].num  =std::min(agrid[i].num,ngf-ngi);
          ngf=agrid[i].num+ngi;
          agrid[i].start=static_cast<real_type>(ngi)*agrid[i].delta;
          agrid[i].end  =static_cast<real_type>(ngf)*agrid[i].delta;

          if (xyz_bc[i].lCode == PERIODIC || xyz_bc[i].lCode == ANTIPERIODIC)
            N[i] = agrid[i].num+3;
          else
            N[i] = agrid[i].num+2;

          Delta[i]=agrid[i].delta;
          Offset[i]=ngi;
        }
        out=einspline::create(out,agrid,xyz_bc,num);
      }
  };

  namespace einspline
  {
    template<typename IV>
      bool outOfBound(const IV& a)
      {
        for(int i=0; i<a.size(); ++i) if(a[i]<0.0 || a[i] >=1.0) return true;
        return false;
      }

    template<typename IV>
      bool validRange(const IV& low, const IV& up)
      {
        bool ok=low[0]<up[0];
        for(int i=1; i<low.size(); ++i) ok &= (low[i]<up[i]);
        return ok;
      }
  }

}
#endif
