
#ifndef QMCPLUSPLUS_AFQMC_SPARSEMATRIX_H
#define QMCPLUSPLUS_AFQMC_SPARSEMATRIX_H

#include<iostream>
#include<vector>
#include<tuple>
#include<assert.h>
#include<algorithm>
#include"AFQMC/config.0.h"

#define ASSERT_SPARSEMATRIX 

namespace qmcplusplus
{

// class that implements a sparse matrix in CSR format
template<class T>
class SparseMatrix
{
  public:

  typedef T            Type_t;
  typedef T            value_type;
  typedef T*           pointer;
  typedef const T*     const_pointer;
  typedef const int*   const_indxPtr;
  typedef int*           indxPtr;
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef SparseMatrix<T>  This_t;

  SparseMatrix<T>():vals(),colms(),myrows(),rowIndex(),nr(0),nc(0),compressed(false),zero_based(true),storage_format(0)
  {
  }

  SparseMatrix<T>(int n):vals(),colms(),myrows(),rowIndex(),nr(n),nc(n),compressed(false),zero_based(true),storage_format(0)
  {
  }

  SparseMatrix<T>(int n,int m):vals(),colms(),myrows(),rowIndex(),nr(n),nc(m),compressed(false),zero_based(true),storage_format(0)
  {
  }

  ~SparseMatrix<T>()
  {
  }

  SparseMatrix(const SparseMatrix<T> &rhs)
  {
    compressed=rhs.compressed;
    zero_based=true;
    nr=rhs.nr;
    nc=rhs.nc;
    vals=rhs.vals;
    myrows=rhs.myrows;
    colms=rhs.colms;
    rowIndex=rhs.rowIndex;
    storage_format=rhs.storage_format;
  }

  inline void reserve(int n)
  {
    vals.reserve(n);
    myrows.reserve(n);
    colms.reserve(n); 
    rowIndex.reserve(nr+1);
  }

  inline bool allocateMemoryAndReserve(int n)
  {
    reserve(n);
    return true;
  }

  inline bool initializeChildren()
  {
    return true;
  }

  inline void resize_arrays(int nnz)
  {
    vals.resize(nnz);
    myrows.resize(nnz);
    colms.resize(nnz);
    rowIndex.resize(nr+1);
  }

  inline void clear() { 
    vals.clear();
    colms.clear();
    myrows.clear();
    rowIndex.clear();
    compressed=false;
    zero_based=true;
  }

  inline void setDims(int n, int m)
  {
    nr=n;
    nc=m;
    compressed=false;
    zero_based=true;
    clear();
  }

  inline void setCompressed() 
  {
    compressed=true;
  }

  inline bool isCompressed() const
  {
    return compressed;
  }
  inline int size() const
  {
    return vals.size();
  }
  inline int rows() const
  {
    return nr;
  }
  inline int cols() const
  {
    return nc;
  }

  inline const_pointer values() const 
  {
    return vals.data();
  }

  inline pointer values() 
  {
    return vals.data();
  }

  inline const_indxPtr column_data() const 
  {
    return colms.data();
  }
  inline indxPtr column_data() 
  {
    return colms.data();
  }

  inline const_indxPtr row_data() const 
  {
    return myrows.data();
  }
  inline indxPtr row_data() 
  {
    return myrows.data();
  }

  inline const_indxPtr row_index() const 
  {
    return rowIndex.data();
  }
  inline indxPtr row_index() 
  {
    return rowIndex.data();
  }

  inline This_t& operator=(const SparseMatrix<T> &rhs) 
  { 
    compressed=rhs.compressed;
    zero_based=rhs.zero_based;
    nr=rhs.nr;
    nc=rhs.nc;
    vals=rhs.vals;
    myrows=rhs.myrows;
    colms=rhs.colms;
    rowIndex=rhs.rowIndex;
  }  

  inline int find_element(int i, int j) {
    for (int k = rowIndex[i]; k < rowIndex[i+1]; k++) {
      if (colms[k] == j) return k;
    }
    return -1;
  }

  // DANGER: This returns a reference, which could allow changes to the stored value.
  // If a zero element is changed, it will change zero everywhere in the matrix.
  // For now this method is only used for testing so it should not be a problem.
  inline Type_t& operator()(int i, int j)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc && compressed); 
#endif
    int idx = find_element(i,j);
    if (idx == -1) return zero;
    return vals[idx];
  }

  inline Type_t operator()( int i, int j) const
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc && compressed); 
#endif
    int idx = find_element(i,j);
    if (idx == -1) return 0;
    return vals[idx];
  }

  inline void add(const int i, const int j, const T& v, bool dummy=false) 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc);
#endif
    compressed=false;
    myrows.push_back(i);
    colms.push_back(j);
    vals.push_back(v);
  }

  inline bool remove_repeated() 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows.size() == colms.size() && myrows.size() == vals.size());
#endif

    compressed=false;
    for(std::vector<int>::iterator itri=myrows.begin(); itri<myrows.end(); itri++)
    {
      int ki = std::distance( myrows.begin(), itri ); 
      for(std::vector<int>::iterator itrj=itri+1; itrj<myrows.end(); itrj++)
      {
        int kj = std::distance( myrows.begin(), itrj ); 
        if( *itri == *itrj && colms[ki] == colms[kj] ) {
          if(vals[ki] != vals[kj]) {
            app_error()<<" Error in call to SparseMatrix::remove_repeated. Same indexes with different values. \n";
            app_error()<<"i: ri, ci, vi: " 
                       <<ki <<" "
                       <<*itri <<" "
                       <<colms[ki] <<" "
                       <<vals[ki] <<"\n"
                       <<"j: rj, cj, vj: " 
                       <<kj <<" "
                       <<*itrj <<" "
                       <<colms[kj] <<" "
                       <<vals[kj] <<std::endl;
            return false;
          }
          itrj = myrows.erase(itrj); 
          colms.erase( colms.begin()+kj );
          vals.erase( vals.begin()+kj ); 
        }
      }
    }
    return true;
  }

  inline void compress_old() 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows.size() == colms.size() && myrows.size() == vals.size());
#endif

   
    // This is not efficient. Write your own iterator to swap all arrays simultaneously during sort  
    // Simple options for now:
    // 1. use memory and efficient std::sort
    // 2. no memory but my inefficient algorithm???
    // Using #1 for now!!!
    
    // order along myrows
    int n=myrows.size(); 
    std::vector<std::tuple<int,int> > toSort;
    toSort.reserve(n);
    for(int i=0; i<n; i++) toSort.push_back(std::forward_as_tuple(myrows[i],i));  
    std::sort(toSort.begin(),toSort.end());
    std::vector<T> tmp;
    tmp=vals; 
    myrows=colms;
    for(int i=0; i<n; i++) {
      int k=std::get<1>(toSort[i]);
      colms[i] = myrows[k]; 
      vals[i] = tmp[k];
    }
    for(int i=0; i<n; i++) 
      myrows[i] = std::get<0>(toSort[i]); 

    if(!std::is_sorted(myrows.begin(),myrows.end())) 
      std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

    // define rowIndex

    rowIndex.resize(nr+1);
    int curr=-1;
    for(int n=0; n<myrows.size(); n++) {
      if( myrows[n] != curr ) {
        int old = curr;
        curr = myrows[n];
        for(int i=old+1; i<=curr; i++) rowIndex[i] = n;
      }
    }
    for(int i=myrows.back()+1; i<rowIndex.size(); i++)
      rowIndex[i] = vals.size();
   
    // order within each rowIndex block
    for(int k=0; k<nr; k++) {
      if(rowIndex[k] == rowIndex[k+1]) continue;       
      toSort.clear();
      tmp.clear();
      for(int i=rowIndex[k],p=0; i<rowIndex[k+1]; i++,p++) toSort.push_back(std::forward_as_tuple(colms[i],p));
      for(int i=rowIndex[k]; i<rowIndex[k+1]; i++) tmp.push_back(vals[i]);
      std::sort(toSort.begin(),toSort.end());
      for(int i=rowIndex[k],p=0; i<rowIndex[k+1]; i++,p++) {
        colms[i] = std::get<0>(toSort[p]);
        vals[i] = tmp[std::get<1>(toSort[p])];
      }      
    } 
   
    compressed=true;
  }

  inline void compress() 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows.size() == colms.size() && myrows.size() == vals.size());
#endif

    // order along myrows
    int n=myrows.size(); 
    sort_rows(0,n-1);
    if(!std::is_sorted(myrows.begin(),myrows.end())) 
      std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

    // define rowIndex
    rowIndex.resize(nr+1);
    int curr=-1;
    for(int n=0; n<myrows.size(); n++) {
      if( myrows[n] != curr ) {
        int old = curr;
        curr = myrows[n];
        for(int i=old+1; i<=curr; i++) rowIndex[i] = n;
      }
    }
    for(int i=myrows.back()+1; i<rowIndex.size(); i++)
      rowIndex[i] = vals.size();
   
    // order within each rowIndex block
    for(int k=0; k<nr; k++) {
      if(rowIndex[k] == rowIndex[k+1]) continue;       
      sort_colms(rowIndex[k],rowIndex[k+1]-1);
    } 
   
    compressed=true;
  }

  void sort_rows(int left, int right) {
    int i = left, j = right;
    auto pivot = myrows[(left + right) / 2];

    /* partition */
    while (i <= j) {
      while (myrows[i] < pivot)
        i++;
      while (myrows[j] > pivot)
        j--;
      if (i <= j) {
        std::swap(myrows[i],myrows[j]);
        std::swap(colms[i],colms[j]);
        std::swap(vals[i++],vals[j--]);
      }
    };

    /* recursion */
    if (left < j)
      sort_rows(left, j);
    if (i < right)
      sort_rows(i, right);
  }

  void sort_colms(int left, int right) {
    int i = left, j = right;
    auto pivot = colms[(left + right) / 2];

    /* partition */
    while (i <= j) {
      while (colms[i] < pivot)
        i++;
      while (colms[j] > pivot)
        j--;
      if (i <= j) {
        std::swap(colms[i],colms[j]);
        std::swap(vals[i++],vals[j--]);
      }
    };

    /* recursion */
    if (left < j)
      sort_colms(left, j);
    if (i < right)
      sort_colms(i, right);
  }

  inline void transpose() {
    assert(myrows.size() == colms.size() && myrows.size() == vals.size());
    for(std::vector<int>::iterator itR=myrows.begin(),itC=colms.begin(); itR!=myrows.end(); ++itR,++itC)
      std::swap(*itR,*itC);
    std::swap(nr,nc);
    compress();
  }

  inline void initFroms1D(std::vector<std::tuple<IndexType,RealType> >& V, bool sorted)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(nr==1);
#endif
    if(!sorted) 
      std::sort(V.begin(),V.end(),my_sort);
    myrows.clear();
    rowIndex.clear();
    vals.clear();
    colms.clear();
  
    int nnz=V.size();
    myrows.resize(nnz);
    vals.resize(nnz);
    colms.resize(nnz);
    rowIndex.resize(nr+1);
    rowIndex[0]=0;
    for(int i=0; i<V.size(); i++) {
      vals[i] = static_cast<T>(std::get<1>(V[i]));
      myrows[i] = 0; 
      colms[i] = std::get<0>(V[i]); 
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nc);
#endif
    }
    rowIndex[1]=V.size();
    compressed=true;
  }

  inline void initFroms1D(std::vector<s1D<std::complex<RealType> > >& V, bool sorted)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(nr==1);
#endif
    if(!sorted) 
      std::sort(V.begin(),V.end(),my_sort);
    myrows.clear();
    rowIndex.clear();
    vals.clear();
    colms.clear();
  
    int nnz=V.size();
    myrows.resize(nnz);
    vals.resize(nnz);
    colms.resize(nnz);
    rowIndex.resize(nr+1);
    rowIndex[0]=0;
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        vals[i] = std::get<1>(V[i]);
      } else {
       assert(false);
      }
      myrows[i] = 0; 
      colms[i] = std::get<0>(V[i]); 
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nc);
#endif
    }
    rowIndex[1]=V.size();
    compressed=true;
  }

  inline void initFroms2D(std::vector<s2D<std::complex<RealType> > >&  V, bool sorted)
  {
    if(!sorted) 
      std::sort(V.begin(),V.end(),my_sort);
    myrows.clear();
    rowIndex.clear();
    vals.clear();
    colms.clear();

    int nnz=V.size();
    myrows.resize(nnz);
    vals.resize(nnz);
    colms.resize(nnz);
    rowIndex.resize(nr+1);
    for(int i=0; i<V.size(); i++) {
      if( std::is_same<T,std::complex<double> >::value  ) {
        vals[i] = std::get<2>(V[i]);
      } else {
       assert(false);
      }
      myrows[i] = std::get<0>(V[i]);
      colms[i] = std::get<1>(V[i]);
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nr);
      assert(std::get<1>(V[i]) >= 0 && std::get<1>(V[i]) < nc);
#endif
    }
    int curr=-1;
    for(int n=0; n<myrows.size(); n++) {
      if( myrows[n] != curr ) {
        int old = curr;
        curr = myrows[n];
        for(int i=old+1; i<=curr; i++) rowIndex[i] = n;
      }
    }
    if (myrows.size() > 0) {
      for(int i=myrows.back()+1; i<rowIndex.size(); i++)
        rowIndex[i] = vals.size();
    }
    compressed=true;
  }

  inline void initFroms2D(std::vector<s2D<RealType> >&  V, bool sorted)
  {
    if(!sorted) 
      std::sort(V.begin(),V.end(),my_sort);
    myrows.clear();
    rowIndex.clear();
    vals.clear();
    colms.clear();

    int nnz=V.size();
    myrows.resize(nnz);
    vals.resize(nnz);
    colms.resize(nnz);
    rowIndex.resize(nr+1);
    for(int i=0; i<V.size(); i++) {
      vals[i] = static_cast<T>(std::get<2>(V[i]));
      myrows[i] = std::get<0>(V[i]);
      colms[i] = std::get<1>(V[i]);
#ifdef ASSERT_SPARSEMATRIX
      assert(std::get<0>(V[i]) >= 0 && std::get<0>(V[i]) < nr);
      assert(std::get<1>(V[i]) >= 0 && std::get<1>(V[i]) < nc);
#endif
    }
    int curr=-1;
    for(int n=0; n<myrows.size(); n++) {
      if( myrows[n] != curr ) {
        int old = curr;
        curr = myrows[n];
        for(int i=old+1; i<=curr; i++) rowIndex[i] = n;
      }
    }
    if (myrows.size() > 0) {
      for(int i=myrows.back()+1; i<rowIndex.size(); i++)
        rowIndex[i] = vals.size();
    }
    compressed=true;
  }

  inline void check()
  {
    for(int i=0; i<rowIndex.size()-1; i++)
    {
      if(rowIndex[i+1] < rowIndex[i]) std::cout<<"Error: SparseMatrix::check(): rowIndex. \n" <<std::endl; 
  
    }

  }

  inline SparseMatrix<T>& operator*=(const RealType rhs ) 
  {
    for(iterator it=vals.begin(); it!=vals.end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline SparseMatrix<T>& operator*=(const std::complex<RealType> rhs ) 
  {
    for(iterator it=vals.begin(); it!=vals.end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline void toZeroBase() {
    if(zero_based) return;
    zero_based=true;
    for (int& i : colms ) i--; 
    for (int& i : myrows ) i--; 
    for (int& i : rowIndex ) i--; 
  }

  inline void toOneBase() {
    if(!zero_based) return;
    zero_based=false;
    for (int& i : colms ) i++; 
    for (int& i : myrows ) i++; 
    for (int& i : rowIndex ) i++; 
  }

  friend std::ostream& operator<<(std::ostream& out, const SparseMatrix<T>& rhs)
  {
    for(int i=0; i<rhs.vals.size(); i++)
      out<<"(" <<rhs.myrows[i] <<"," <<rhs.colms[i] <<":" <<rhs.vals[i] <<")\n"; 
    return out;
  }

  friend std::istream& operator>>(std::istream& in, SparseMatrix<T>& rhs)
  {
    T v;
    int c,r;
    in>>r >>c >>v;
    rhs.vals.push_back(v); 
    rhs.myrows.push_back(r);
    rhs.colms.push_back(c);
    return in;
  }

  // this is ugly, but I need to code quickly 
  // so I'm doing this to avoid adding hdf5 support here 
  inline std::vector<T>* getVals() { return &vals; } 
  inline std::vector<int>* getRows() { return &myrows; }
  inline std::vector<int>* getCols() { return &colms; }
  inline std::vector<int>* getRowIndex() { return &rowIndex; }

  void setRowsFromRowIndex()
  {
    int shift = zero_based?0:1;
    myrows.resize(vals.size());
    for(int i=0; i<nr; i++)
     for(int j=rowIndex[i]; j<rowIndex[i+1]; j++)
      myrows[j]=i+shift;
  }
  bool zero_base() const { return zero_based; }
  int row_max() const { return max_in_row; }
  int format() const { return storage_format; }

  private:

  bool compressed;
  int nr,nc;
  std::vector<T> vals;
  std::vector<int> colms,myrows,rowIndex;
  bool zero_based;
  int storage_format; // 0: CSR, 1: Compressed Matrix (ESSL) 
  int max_in_row;
  Type_t zero; // zero for return value

  _mySort_snD_ my_sort;

};


}

#endif
