namespace qmcplusplus
{
  /** class to use hyperslab with any container
  */
  template<typename CT, unsigned D> 
    struct HyperSlabProxy: public h5_space_type<typename CT::value_type,D>
  {
    ///typedef the pointer
    typedef h5_space_type<typename CT::value_type,D>::pointer pointer;
    ///reference to the physical container
    CT&  ref_;
    ///rank of the slab
    int slab_rank;
    ///local dimension of the hyperslab
    std::vector<hsize_t> slab_dims;
    ///offset of the hyperslab
    std::vector<hsize_t> slab_offset;

    Vector2HyperSlabProxy(CT& a): ref_(a) 
    {
      slab_rank= this->size();
      slab_dims.resize(slab_rank);
      slab_offset.resize(slab_rank);
    }

    inline int rank() const { return slab_rank;}
    inline hsize_t* gcount() { return this->dims;}
    inline hsize_t* count() { return &slab_dims[0];}
    inline hsize_t* offset() { return &slab_offset[0];}
    inline pointer data() { return get_address(&ref_[0]); }
    inline void resize() 
    { 
      hsize_t mreq=ref_.dims[0];
      for(int i=1;i<D; ++i) mreq *= ref_.dims[i];
      if(ref_.size() != mreq) ref_.resize(mreq); 
    }
    inline bool empty() const { return ref_.size() == 0;}
  };

  template<typename CT, unsigned D>
    struct HDFAttribIO<HyperSlabProxy<CT,D> >
    {
      typedef HyperSlabProxy<CT,D> data_type;
      data_type& ref_;

      HDFAttribIO(data_type& a):ref_(a) { }

      bool write(hid_t grp, const std::string& aname)
      {
        if(ref_.empty()) return; //cannot write an empty conainer
        hid_t dset_id=H5Dopen(grp,name);
        hid_t tid=get_h5_datatype(*ref_.data());
        if(dset_id<0)
        {
          hid_t sid1  = H5Screate_simple(ref_.rank(),ref_.gcount(),NULL);
          dset_id=H5Dcreate(grp,name,tid,sid1,H5P_DEFAULT);
          H5Sclose(sid1);
        }

        TinyVector<hsize_t,D> stride(1);
        hid_t memspace=H5Screate_simple(ref_.rank(),ref_.count(),NULL);
        hid_t filespace=H5Dget_space(dset_id);
        herr_t ret=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,ref_.offset(),stride.data(),ref_.count(),NULL);
        ret = H5Dwrite(dset_id,tid,memspace,filespace,H5P_DEFAULT,ref_.data());
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(dset_id);

        return ret != -1;
      }

      bool read(hid_t grp, const std::string& aname)
      {
        hid_t dataset = H5Dopen(grp,aname.c_str());
        if(dataset<0) return false;
        hid_t dataspace = H5Dget_space(dataset);
        int rank_in = H5Sget_simple_extent_ndims(dataspace);
        if(rank_in != ref_.rank()) return false;

        int status_n = H5Sget_simple_extent_dims(dataspace, ref_.count(), NULL);
        ref_.resize();

        hid_t memspace = H5Screate_simple(ref_.rank(), ref_.gcount(), NULL);
        herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, ref_.offset(),NULL,ref_.gcount(),NULL);
        hid_t tid=get_h5_datatype(*ref_.data());
        status = H5Dread(dataset, tid, memspace, dataspace, H5P_DEFAULT, ref_.data);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);
      }
    };
}

