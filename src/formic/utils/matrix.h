///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/vector/matrix.h
///
/// \brief   header file for the formic::matrix class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_VECTOR_MATRIX_HEADER
#define FORMIC_VECTOR_MATRIX_HEADER

#include <string>
#include <utility>
#include <algorithm>
#include <complex>
#include <iterator>

#include <boost/format.hpp>

#include <formic/utils/reusable_array.h>
#include <formic/utils/numeric.h>
#include <formic/utils/lapack_interface.h>

namespace formic {

  // pre-declaration of non-const matrix and vectors
  template<class S> class Matrix;
  template<class S> class ColVec;
  template<class S> class RowVec;

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################  Row and Col Iterators  ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  An iterator used for iterating over the elements of a matrix row.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class RowConstIterator : public std::iterator<std::bidirectional_iterator_tag, const S> {
    private:
      const S * m_x;
      size_t m_n;
    public:
      RowConstIterator(const S * x, size_t n) : m_x(x), m_n(n) {}
      const S & operator*() const { return *m_x; }
      const S * operator->() const { return m_x; }
      RowConstIterator & operator++() { m_x += m_n; return *this; }                          // prefix ++
      RowConstIterator & operator--() { m_x -= m_n; return *this; }                          // prefix --
      RowConstIterator operator++(int) { RowConstIterator it(*this); ++(*this); return it; } // postfix ++
      RowConstIterator operator--(int) { RowConstIterator it(*this); --(*this); return it; } // postfix --
      bool operator!=(const RowConstIterator & other) const { return m_x != other.m_x; }
      bool operator==(const RowConstIterator & other) const { return m_x == other.m_x; }
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  An iterator used for iterating over the elements of a matrix row.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class RowIterator : public std::iterator<std::bidirectional_iterator_tag, S> {
    private:
      S * m_x;
      size_t m_n;
    public:
      RowIterator(S * x, size_t n) : m_x(x), m_n(n) {}
      S & operator*() { return *m_x; }
      S * operator->() { return m_x; }
      RowIterator & operator++() { m_x += m_n; return *this; }                     // prefix ++
      RowIterator & operator--() { m_x -= m_n; return *this; }                     // prefix --
      RowIterator operator++(int) { RowIterator it(*this); ++(*this); return it; } // postfix ++
      RowIterator operator--(int) { RowIterator it(*this); --(*this); return it; } // postfix --
      bool operator!=(const RowIterator & other) const { return m_x != other.m_x; }
      bool operator==(const RowIterator & other) const { return m_x == other.m_x; }
  };

//  ///////////////////////////////////////////////////////////////////////////////////////////////////
//  /// \brief  An iterator used for iterating over the elements of a matrix column.
//  ///
//  ///////////////////////////////////////////////////////////////////////////////////////////////////
//  template<class S> class ColConstIterator {
//    private:
//      const S * m_x;
//    public:
//      ColConstIterator(const S * x) : m_x(x) {}
//      const S & operator*() const { return *m_x; }
//      ColConstIterator & operator++() { m_x += 1; return *this; }                            // prefix ++
//      ColConstIterator & operator--() { m_x -= 1; return *this; }                            // prefix --
//      ColConstIterator operator++(int) { ColConstIterator it(*this); ++(*this); return it; } // postfix ++
//      ColConstIterator operator--(int) { ColConstIterator it(*this); --(*this); return it; } // postfix --
//  };
//
//  ///////////////////////////////////////////////////////////////////////////////////////////////////
//  /// \brief  An iterator used for iterating over the elements of a matrix column.
//  ///
//  ///////////////////////////////////////////////////////////////////////////////////////////////////
//  template<class S> class ColIterator {
//    private:
//      S * m_x;
//    public:
//      ColIterator(S * x) : m_x(x) {}
//      S & operator*() const { return *m_x; }
//      ColIterator & operator++() { m_x += 1; return *this; }                       // prefix ++
//      ColIterator & operator--() { m_x -= 1; return *this; }                       // prefix --
//      ColIterator operator++(int) { ColIterator it(*this); ++(*this); return it; } // postfix ++
//      ColIterator operator--(int) { ColIterator it(*this); --(*this); return it; } // postfix --
//  };

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################    ConstMatrix class    ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  The base class for matrices and vectors.  ConstMatrix cannot change the values of
  ///         the elements in the data array it points to, and can thus be safely used as a wrapper
  ///         around const data (see its constructor with external storage below).  The Matrix
  ///         class defined below cannot be used for this purpose as it is allowed to change the
  ///         values of the matrix elements, which is the reason for having ConstMatrix as its
  ///         own class.  Note that while the element values cannot be changed, ConstMatrix can
  ///         change the array it points to via the assignment operator, i.e.
  ///
  ///             ConstMatrix<double> cmat(2, 2, 2.0);
  ///             Matrix<double> mmat(3, 3, 3.0);
  ///             cmat = mat;
  ///
  ///         This works because Matrix is a child of ConstMatrix.  The result is that
  ///         both cmat and mat are now 3 by 3 and point to the same array of 3.0 values.
  ///         The array of 2.0 values is no longer used and is held in reserve by the
  ///         ReusableArray class until it is needed for something else or garbage collection
  ///         is asked for through the appropriate ReusableArray function.
  ///         From the user's perspective, the array of 2.0 values is gone.
  ///         Note that cmat can only perform operations that do not change the 3.0 values,
  ///         while mmat can edit them.  For example,
  ///
  ///             mmat(0,1) = 7.0; // allowed, now cmat(0,1) is also 7.0 due to shallow copy
  ///             cmat(0,1) = 5.0; // compile error due to attempt to write to const reference
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class ConstMatrix {

    //////////   data members   //////////
    protected:

      /// \brief  number of rows
      size_t m_n;

      /// \brief  number of columns
      size_t m_m;

      /// \brief  whether the matrix is actually a column vector
      bool m_is_col_vec;

      /// \brief  whether the matrix is actually a row vector
      bool m_is_row_vec;

      /// \brief  pointer to data array
      S * m_p;

      /// \brief  smart pointer to reusable array object
      boost::shared_ptr<formic::ReusableArray<S> > m_ra;

    //////////   initializers   //////////
    protected:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Initializes the matrix to use its internal data array.
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void initialize_with_internal_storage(const size_t n, const size_t m) {

        // initialize dimensions and data pointer
        m_n = n;
        m_m = m;
        m_p = 0;

        // drop any internally held data array
        m_ra.reset();

        // if necessary, get a new internal data array and point to it
        if ( n > 0 && m > 0 ) {
          boost::shared_ptr<formic::ReusableArray<S> > ra( new formic::ReusableArray<S>(n*m) );
          m_ra = ra;
          m_p = m_ra->cptr();
        }

      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Initializes the matrix to use an external data array.
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      p        pointer to the external array
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void initialize_with_external_storage(const size_t n, const size_t m, const S * const p) {

        // make sure the array is not null
        if ( p == 0 )
          throw formic::Exception("cannot initialize a matrix with a null array pointer");

        // error if we try to point back to our own internally held array
        this->ensure_p_and_ra_are_consistent();
        if ( m_ra && p >= this->begin() && p < this->end() )
          throw formic::Exception("illegal attempt to use the matrix's own internal array as the external storage array");

        // initialize dimensions
        m_n = n;
        m_m = m;

        // Point to the external data array.
        // Note that we cast away the const quality, but that if this is a ConstMatrix
        // the elements in the array will not be modified.
        m_p = (S *)p;

        // drop any internally held data array
        m_ra.reset();

      }

    //////////   constructors   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor without initialization of elements (also the defualt constructor)
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ConstMatrix(const size_t n = 0, const size_t m = 0) : m_is_col_vec(false), m_is_row_vec(false) {
        this->initialize_with_internal_storage(n, m);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor with initialization of elements
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      v        value to initialize all matrix elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ConstMatrix(const size_t n, const size_t m, const S v) : m_is_col_vec(false), m_is_row_vec(false) {
        this->initialize_with_internal_storage(n, m);
        std::fill(m_p, m_p + this->size(), v);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the matrix around the specified, externally held array
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      p        array holding the matrix elements
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ConstMatrix(const size_t n, const size_t m, const S * const p) : m_is_col_vec(false), m_is_row_vec(false) {
        this->initialize_with_external_storage(n, m, p);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copy constructor performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///         Because Matrix, ColVec, and RowVec are children of ConstMatrix, this constructor
      ///         can copy from any of them as well as from a ConstMatrix.
      ///
      /// \param[in]      other    the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ConstMatrix(const ConstMatrix & other) : m_is_col_vec(false), m_is_row_vec(false) {
        m_n = other.m_n;
        m_m = other.m_m;
        m_p = other.m_p;
        m_ra = other.m_ra;
      }

    //////////   assignment operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///         Because Matrix, ColVec, and RowVec are children of ConstMatrix, other
      ///         can be of those types as well as ConstMatrix.
      ///
      /// \param[in]      other    the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ConstMatrix & operator=(const ConstMatrix & other) {

        // assign
        if ( this != &other ) {
          m_n = other.m_n;
          m_m = other.m_m;
          m_p = other.m_p;
          m_ra = other.m_ra;
        }

        // check that the result makes sense
        this->check_vector_sanity("operator=");

        // return a reference to ourself
        return *this;

      }

    //////////   deep copy operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the matrix into the supplied array
      ///
      /// \param[in,out]  x      On input, a length this->size() array.
      ///                        On output, the matrix data has been copied into x.
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const ConstMatrix & operator>>=(S * const x) const {
        std::copy(this->begin(), this->end(), x);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the matrix into the supplied std::vector
      ///
      /// \param[out]     x      the vector to copy into, which is resized appropriately
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const ConstMatrix & operator>>=(std::vector<S> & x) const {
        x.resize(this->size());
        std::copy(this->begin(), this->end(), x.begin());
        return *this;
      }

    //////////   sanity checks   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  If the object is a RowVec or ColVec, checks it obeys their restrictions.
      ///
      /// \param[in]      calling_func   name of the function the check was called from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void check_vector_sanity(const std::string & calling_func) const {

        // we can't be both a row and a column vector
        if ( m_is_col_vec && m_is_row_vec )
          throw formic::Exception("object in %s claims to be both a RowVec and a ColVec, which is not allowed") % calling_func;

        // row vectors can only have one row
        if ( m_is_row_vec && this->size() > 0 && this->rows() != 1 )
          throw formic::Exception("RowVec with %i rows detected inside %s") % this->rows() % calling_func;

        // column vectors can only have one column
        if ( m_is_col_vec && this->size() > 0 && this->cols() != 1 )
          throw formic::Exception("ColVec with %i columns detected inside %s") % this->cols() % calling_func;

      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Check the matrix's dimensions
      ///
      /// \param[in]      n              how many rows we should have
      /// \param[in]      m              how many columns we should have
      /// \param[in]      calling_func   name of the function the check was called from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void check_dimensions(const size_t n, const size_t m, const std::string & calling_func) const {
        if ( n != m_n || m != m_m ) {
          const std::string s( m_is_row_vec || m_is_col_vec ? "vector" : "matrix" );
          throw formic::Exception("%s expected a %i by %i %s but found a %i by %i %s instead") % calling_func % n % m % s % m_n % m_m % s;
        }
      }

      /// \brief  checks that if the internally held array is in use, the data pointer points to it
      void ensure_p_and_ra_are_consistent() const {
        if ( m_ra )
          if ( m_p != m_ra->cptr() )
            throw formic::Exception("illegal situation in which m_ra is valid but m_p does not point to its array");
      }

    //////////   element access   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a const reference to the i,j element after checking that i and j are in bounds.
      ///
      /// \param[in]      i        the row index
      /// \param[in]      j        the column index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const S & at(const size_t i, const size_t j) const {
        if ( i >= m_n ) throw formic::Exception("out-of-range row index i = %i in formic::ConstMatrix::at for matrix of dimension %i by %i") % i % m_n % m_m;
        if ( j >= m_m ) throw formic::Exception("out-of-range col index j = %i in formic::ConstMatrix::at for matrix of dimension %i by %i") % j % m_n % m_m;
        return m_p[i+m_n*j];
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a const reference to the i,j element
      ///
      /// \param[in]      i        the row index
      /// \param[in]      j        the column index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const S & operator()(const size_t i, const size_t j) const {
        #ifndef NDEBUG
        return this->at(i,j); // do bounds checking if debugging
        #else
        return m_p[i+m_n*j];
        #endif
      }

    //////////   functions returning iterators   //////////
    public:

      /// \brief  return a const pointer to the beginning of the data array
      const S * begin() const { return m_p; }

      /// \brief  return a const pointer to the end of the data array
      const S * end() const { return m_p + this->size(); }

      /// \brief  return a const iterator to the beginning of the ith column
      const S * col_begin(const size_t i) const {
        if ( i >= m_m || this->size() == 0 )
          throw formic::Exception("out-of-bounds column index %i for a %i by %i matrix in formic::ConstMatrix::col_begin") % i % m_n % m_m;
        return m_p + i * m_n;
      }

      /// \brief  return a const iterator to the end of the ith column
      const S * col_end(const size_t i) const {
        if ( i >= m_m || this->size() == 0 )
          throw formic::Exception("out-of-bounds column index %i for a %i by %i matrix in formic::ConstMatrix::col_end") % i % m_n % m_m;
        return m_p + ( i + 1 ) * m_n;
      }

      /// \brief  return a const iterator to the beginning of the ith row
      formic::RowConstIterator<S> row_begin(const size_t i) const {
        if ( i >= m_n || this->size() == 0 )
          throw formic::Exception("out-of-bounds row index %i for a %i by %i matrix in formic::ConstMatrix::row_begin") % i % m_n % m_m;
        return formic::RowConstIterator<S>(m_p + i, m_n);
      }

      /// \brief  return a const iterator to the end of the ith row
      formic::RowConstIterator<S> row_end(const size_t i) const {
        if ( i >= m_n || this->size() == 0 )
          throw formic::Exception("out-of-bounds row index %i for a %i by %i matrix in formic::ConstMatrix::row_end") % i % m_n % m_m;
        return formic::RowConstIterator<S>(m_p + i + this->size(), m_n);
      }

    //////////   other public members   //////////
    public:

      /// \brief  return the number of elements
      size_t size() const { return m_n * m_m; }

      /// \brief  return the number of rows
      size_t rows() const { return m_n; }

      /// \brief  return the number of columns
      size_t cols() const { return m_m; }

      /// \brief  returns whether the object is actually a vector
      bool is_a_vector() const { return ( m_is_col_vec || m_is_row_vec ); }

      /// \brief  returns whether the object is actually a row vector
      bool is_row_vec() const { return m_is_row_vec; }

      /// \brief  returns whether the object is actually a row vector
      bool is_col_vec() const { return m_is_col_vec; }

      /// \brief  returns whether the matrix is using as its array something other than its own internal array
      bool storage_is_external() const {
        this->ensure_p_and_ra_are_consistent();
        return ( m_p && !m_ra );
      }

      /// \brief  returns whether this matrix's dimensions differ from that of other
      bool dimensions_differ_from(const ConstMatrix & other) const {
        return ( this->rows() != other.rows() || this->cols() != other.cols() );
      }

      /// \brief  returns the sum of square absolute values of the of the elements
      S norm2() const {
        S retval = formic::zero(S());
        for (const S * x = this->begin(); x != this->end(); x++)
          retval += formic::conj(*x) * (*x);
        return retval;
      }

      /// \brief  returns the compound index for element (i,j)
      size_t cmpd(const size_t i, const size_t j) const { return i + this->m_n * j; }

      /// \brief  gives the row and column indices for the supplied compound index
      void row_col_from_cmpd(const size_t c, size_t & i, size_t & j) const {
        if ( c >= this->size() )
          throw formic::Exception("out-of-bounds compound index %i for a %i by %i matrix in formic::ConstMatrix::row_col_from_cmpd") % c % m_n % m_m;
        i = c % m_n;
        j = c / m_n;
      }

      /// \brief  gives the row and column indices for the supplied iterator
      void row_col_from_iter(const S * const iter, size_t & i, size_t & j) const {
        if ( iter < this->begin() || iter >= this->end() )
          throw formic::Exception("out-of-bounds iterator for this %i by %i matrix in formic::ConstMatrix::row_col_from_iter") % m_n % m_m;
        this->row_col_from_cmpd(iter - this->begin(), i, j);
      }

      /// \brief  conversion operator to turn a 1 by 1 matrix into a scalar
      operator S() const {
        if ( this->size() != 1 )
          throw formic::Exception("illegal attempt to convert a %i by %i matrix to a scalar") % m_n % m_m;
        const S retval = m_p[0];
        return retval;
      }

      formic::Matrix<S> as_diag_mat() const;

      formic::Matrix<std::complex<double> > get_complex_form() const;

      std::string print(const std::string & fmt, const std::string & name = "") const;

      formic::Matrix<S> clone() const;

      formic::Matrix<S> operator-() const;

      formic::Matrix<S> operator+() const;

      formic::Matrix<S> t() const;

      formic::Matrix<S> c() const;

      formic::Matrix<S> conj() const;

      template<class V> formic::Matrix<S> slice(const V & row_inds, const V & col_inds) const;

      void svd(formic::Matrix<S> & u, formic::ColVec<S> & s, formic::Matrix<S> & vt) const;

      formic::Matrix<S> inv(const double thresh = -1.0) const;

      void nonsym_eig(formic::ColVec<std::complex<double> > & w, formic::Matrix<std::complex<double> > & v) const;
      
      void sym_eig(formic::ColVec<double> & w, formic::Matrix<S> & v) const;

      const S * max_abs_elem_iter() const;

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Returns the element with the largest absolute value and also gives its indices.
      ///        Note we return the element itself, not its absolute value.
      ///
      /// \param[out]     i        on exit, the row index of the desired element
      /// \param[out]     j        on exit, the column index of the desired element
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S max_abs_elem(size_t & i, size_t & j) const {
        const S * const best = this->max_abs_elem_iter();
        this->row_col_from_iter(best, i, j);
        return *best;
      }

      /// \brief Returns the element with the largest absolute value.  Note we return the element itself, not its absolute value.
      S max_abs_elem() const {
        return *this->max_abs_elem_iter();
      }

      const S * min_abs_elem_iter() const;

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Returns the element with the smallest absolute value and also gives its indices.
      ///        Note we return the element itself, not its absolute value.
      ///
      /// \param[out]     i        on exit, the row index of the desired element
      /// \param[out]     j        on exit, the column index of the desired element
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S min_abs_elem(size_t & i, size_t & j) const {
        const S * const best = this->min_abs_elem_iter();
        this->row_col_from_iter(best, i, j);
        return *best;
      }

      /// \brief Returns the element with the smallest absolute value.  Note we return the element itself, not its absolute value.
      S min_abs_elem() const {
        return *this->min_abs_elem_iter();
      }

      //S det() const {
      //  if ( m_n != m_m || this->size() == 0 )
      //    throw formic::Exception("cannot take determinant of a %i by %i matrix") % m_n % m_m;
      //  throw formic::Exception("ConstMatrix det not yet implemented");
      //}

      formic::ColVec<S> col_as_vec(const size_t j) const;

      formic::RowVec<S> row_as_vec(const size_t i) const;

  };

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################       Matrix class      ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  The base class for matrices and vectors that can modify their elements, which
  ///         inherits many features that do not modify elements from ConstMatrix.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class Matrix : public formic::ConstMatrix<S> {

    //////////   data members   //////////
    protected:

      // due to template base class and two phase lookup, we need to tell the compiler that these exist
      using formic::ConstMatrix<S>::m_n;
      using formic::ConstMatrix<S>::m_m;
      using formic::ConstMatrix<S>::m_is_col_vec;
      using formic::ConstMatrix<S>::m_is_row_vec;
      using formic::ConstMatrix<S>::m_p;
      using formic::ConstMatrix<S>::m_ra;

    //////////   constructors   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor without initialization of elements (also the defualt constructor)
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit Matrix(const size_t n = 0, const size_t m = 0) : formic::ConstMatrix<S>(n, m) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor with initialization of elements
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      v        value to initialize all matrix elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit Matrix(const size_t n, const size_t m, const S v) : formic::ConstMatrix<S>(n, m, v) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the matrix around the specified, externally held array
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      p        array holding the matrix elements
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit Matrix(const size_t n, const size_t m, S * const p) : formic::ConstMatrix<S>(n, m, p) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the matrix around the specified, externally held array, initializing
      ///         all elements to the specified value
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      m        the desired number of columns
      /// \param[in]      p        array holding the matrix elements
      /// \param[in]      v        value to initialize all matrix elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit Matrix(const size_t n, const size_t m, S * const p, const S v) : formic::ConstMatrix<S>(n, m, p) { (*this) = v; }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copy constructor performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///         Because ColVec and RowVec are children of Matrix, this constructor
      ///         can copy from either of them as well as from a Matrix.
      ///
      /// \param[in]      other    the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Matrix(const Matrix & other) : formic::ConstMatrix<S>(other) {}

    //////////   assignment operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///         Because ColVec and RowVec are children of Matrix, other
      ///         can be either of those types as well as type Matrix.
      ///
      /// \param[in]      other    the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Matrix & operator=(const Matrix & other) {
        ConstMatrix<S>::operator=(other);
        return *this;
      }

      /// \brief  sets all elements of the matrix equal to a
      formic::Matrix<S> & operator=(const S a) {
        std::fill(this->begin(), this->end(), a);
        return *this;
      }

    //////////   deep copy operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the matrix into the supplied array
      ///
      /// \param[in,out]  x      On input, a length this->size() array.
      ///                        On output, the matrix data has been copied into x.
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const Matrix & operator>>=(S * const x) const {
        std::copy(this->begin(), this->end(), x);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the supplied array into the matrix
      ///
      /// \param[in]      x      a length this->size() array to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Matrix & operator<<=(const S * const x) {
        std::copy(x, x + this->size(), this->begin());
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the matrix into the supplied std::vector
      ///
      /// \param[out]     x      the vector to copy into, which is resized appropriately
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const Matrix & operator>>=(std::vector<S> & x) const {
        x.resize(this->size());
        std::copy(this->begin(), this->end(), x.begin());
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copies the data in the supplied std::vector into the matrix
      ///
      /// \param[in]      x      a length this->size() std::vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Matrix & operator<<=(const std::vector<S> & x) {
        if ( this->size() != x.size() )
          throw formic::Exception("dimension mismatch in attempt to deep copy a length %i std::vector into a %i by %i formic::Matrix") % x.size() % m_n % m_m;
        std::copy(x.begin(), x.end(), this->begin());
        return *this;
      }

      /// \brief  deep copy from another matrix, resizing this matrix to match other's dimensions
      Matrix & operator<<=(const formic::ConstMatrix<S> & other) {
        if ( this->size() == other.size() ) {
          m_n = other.rows();
          m_m = other.cols();
          this->check_vector_sanity("operator<<=(ConstMatrix)");
        } else {
          // if the matrix has external storage, a size mismatch is an error
          if ( this->storage_is_external() )
            throw formic::Exception("illegal attempt to deep copy a %i by %i matrix into a %i by %i matrix that has external storage") % other.rows() % other.cols() % m_n % m_m;
          this->reset(other.rows(), other.cols());
        }
        std::copy(other.begin(), other.end(), this->begin());
        return *this;
      }

    //////////   element access   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a reference to the i,j element after checking that i and j are in bounds.
      ///
      /// \param[in]      i        the row index
      /// \param[in]      j        the column index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S & at(const size_t i, const size_t j) {
        if ( i >= m_n ) throw formic::Exception("out-of-range row index i = %i in formic::Matrix::at for matrix of dimension %i by %i") % i % m_n % m_m;
        if ( j >= m_m ) throw formic::Exception("out-of-range col index j = %i in formic::Matrix::at for matrix of dimension %i by %i") % j % m_n % m_m;
        return m_p[i+m_n*j];
      }

      // need the using keyword so that the const version of at in ConstMatrix isn't hidden by the one just above
      using ConstMatrix<S>::at;

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a reference to the i,j element
      ///
      /// \param[in]      i        the row index
      /// \param[in]      j        the column index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S & operator()(const size_t i, const size_t j) {
        #ifndef NDEBUG
        return this->at(i,j); // do bounds checking if debugging
        #else
        return m_p[i+m_n*j];
        #endif
      }

      // need the using keyword so that the const version of operator() in ConstMatrix isn't hidden by the one just above
      using ConstMatrix<S>::operator();

    //////////   functions returning iterators   //////////
    public:

      // need the using keyword so that the const version of iterator returning functions are not hidden
      using ConstMatrix<S>::begin;
      using ConstMatrix<S>::end;
      using ConstMatrix<S>::col_begin;
      using ConstMatrix<S>::col_end;
      using ConstMatrix<S>::row_begin;
      using ConstMatrix<S>::row_end;

      /// \brief  return a pointer to the beginning of the data array
      S * begin() { return m_p; }

      /// \brief  return a pointer to the end of the data array
      S * end()   { return m_p + this->size(); }

      /// \brief  return a const iterator to the beginning of the ith column
      S * col_begin(const size_t i) {
        if ( i >= m_m || this->size() == 0 )
          throw formic::Exception("out-of-bounds column index %i for a %i by %i matrix in formic::Matrix::col_begin") % i % m_n % m_m;
        return m_p + i * m_n;
      }

      /// \brief  return a const iterator to the end of the ith column
      S * col_end(const size_t i) {
        if ( i >= m_m || this->size() == 0 )
          throw formic::Exception("out-of-bounds column index %i for a %i by %i matrix in formic::Matrix::col_end") % i % m_n % m_m;
        return m_p + ( i + 1 ) * m_n;
      }

      /// \brief  return a const iterator to the beginning of the ith row
      formic::RowIterator<S> row_begin(const size_t i) {
        if ( i >= m_n || this->size() == 0 )
          throw formic::Exception("out-of-bounds row index %i for a %i by %i matrix in formic::Matrix::row_begin") % i % m_n % m_m;
        return formic::RowIterator<S>(m_p + i, m_n);
      }

      /// \brief  return a const iterator to the end of the ith row
      formic::RowIterator<S> row_end(const size_t i) {
        if ( i >= m_n || this->size() == 0 )
          throw formic::Exception("out-of-bounds row index %i for a %i by %i matrix in formic::Matrix::row_end") % i % m_n % m_m;
        return formic::RowIterator<S>(m_p + i + this->size(), m_n);
      }

    //////////   other public members   //////////
    public:

      /// \brief  clears the matrix, leaving it as a 0 by 0 empty matrix
      void clear() { *this = formic::Matrix<S>(); }

      /// \brief  clears the matrix and re-allocates it to be the specified size
      Matrix & reset(const int n, const int m) {
        this->clear();
        *this = formic::Matrix<S>(n,m);
        this->check_vector_sanity("reset(n,m)");
        return *this;
      }

      /// \brief  clears the matrix and re-allocates it to be the specified size, initializing its elements to v
      Matrix & reset(const int n, const int m, const S v) {
        this->clear();
        *this = formic::Matrix<S>(n,m,v);
        this->check_vector_sanity("reset(n,m,v)");
        return *this;
      }

      /// \brief  adds x to each element of the matrix
      formic::Matrix<S> & operator+=(const S x) {
        for (S * it = this->begin(); it != this->end(); it++)
          *it += x;
        return *this;
      }

      /// \brief  subtracts x from each element of the matrix
      formic::Matrix<S> & operator-=(const S x) {
        for (S * it = this->begin(); it != this->end(); it++)
          *it -= x;
        return *this;
      }

      /// \brief  multiplies each element of the matrix by x
      formic::Matrix<S> & operator*=(const S x) {
        for (S * it = this->begin(); it != this->end(); it++)
          *it *= x;
        return *this;
      }

      /// \brief  divides each element of the matrix by x
      formic::Matrix<S> & operator/=(const S x) {
        for (S * it = this->begin(); it != this->end(); it++)
          *it /= x;
        return *this;
      }

      /// \brief  adds to each element of this matrix the corresponding element of other
      formic::Matrix<S> & operator+=(const formic::ConstMatrix<S> & other) {
        if ( this->dimensions_differ_from(other) )
          throw formic::Exception("cannot += matrices of different dimensions");
        const S * jt = other.begin();
        for (S * it = this->begin(); it != this->end(); it++, jt++)
          *it += *jt;
        return *this;
      }

      /// \brief  subtracts from each element of this matrix the corresponding element of other
      formic::Matrix<S> & operator-=(const formic::ConstMatrix<S> & other) {
        if ( this->dimensions_differ_from(other) )
          throw formic::Exception("cannot -= matrices of different dimensions");
        const S * jt = other.begin();
        for (S * it = this->begin(); it != this->end(); it++, jt++)
          *it -= *jt;
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Adds the scaled elements of the supplied array to the matrix's elements
      ///
      /// \param[in]      a        scale factor to apply to the array elements
      /// \param[in]      x        pointer to the start of a length this->size() array
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      formic::Matrix<S> & axpy(const S a, const S * x) {
        for (S * it = this->begin(); it != this->end(); it++, x++)
          *it += a * (*x);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Adds the scaled elements of the supplied matrix to the matrix's elements
      ///
      /// \param[in]      a        factor to scale the elements of other by before adding them
      /// \param[in]      other    the other matrix
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      formic::Matrix<S> & axpy(const S a, const formic::ConstMatrix<S> & other) {
        if ( this->dimensions_differ_from(other) )
          throw formic::Exception("cannot axpy matrices of different dimensions");
        return this->axpy(a, other.begin());
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Scales the specified column by the specified constant
      ///
      /// \param[in]      i        which column to scale
      /// \param[in]      a        factor to scale the column by
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void scale_col_by(const size_t i, const S a) {
        if ( i >= this->cols() )
          throw formic::Exception("%i is out of bounds for Matrix::scale_col_by in a %i by %i matrix") % i % this->cols() % this->rows();
        S * start = this->begin() + i * this->rows();
        S * finish = start + this->rows();
        for ( ; start != finish; start++)
          *start *= a;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Scales the specified row by the specified constant
      ///
      /// \param[in]      i        which row to scale
      /// \param[in]      a        factor to scale the row by
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void scale_row_by(const size_t i, const S a) {
        if ( i >= this->rows() )
          throw formic::Exception("%i is out of bounds for Matrix::scale_row_by in a %i by %i matrix") % i % this->cols() % this->rows();
        S * start = this->begin() + i;
        S * finish = start + this->size();
        for ( ; start != finish; start += this->rows())
          *start *= a;
      }

      /// \brief  transpose the matrix in place
      formic::Matrix<S> & tip() {
        if ( this->is_a_vector() )
          throw formic::Exception("illegal attempt to transpose a vector in place");
        if ( m_n == m_m ) {
          for (size_t i = 0; i < m_n; i++)
          for (size_t j = 0; j < i; j++)
            std::swap((*this)(i,j), (*this)(j,i));
        } else {
          const formic::ConstMatrix<S> temp = this->clone();
          std::swap(m_n, m_m);
          for (size_t j = 0; j < m_m; j++)
          for (size_t i = 0; i < m_n; i++)
            (*this)(i,j) = temp(j,i);
        }
        return *this;
      }

      /// \brief  replace all elements with their complex conjugates (i.e. conjugate the matrix in-place)
      formic::Matrix<S> & conjip() {
        for (S * x = this->begin(); x != this->end(); x++)
          *x = formic::conj(*x);
        return *this;
      }

      /// \brief  conjugate transpose the matrix in place
      formic::Matrix<S> & cip() {
        if ( this->is_a_vector() )
          throw formic::Exception("illegal attempt to conjugate transpose a vector in place");
        this->conjip(); // conjugate
        this->tip();    // transpose
        return *this;
      }

      void put_vec_in_row(const size_t i, const formic::ConstMatrix<S> & v);

      void put_vec_in_col(const size_t j, const formic::ConstMatrix<S> & v);

      void random_fill(const double lb = -1.0, const double ub = 1.0);

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Extends the matrix by adding new rows and columns, after which the original values
      ///         will be in the upper left corner of the extended matrix
      ///
      /// \param[in]      p        number of rows to add
      /// \param[in]      q        number of columns to add
      /// \param[in]      val      value to initialize new elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void extend_by_p_rows_and_q_cols(const size_t p, const size_t q, const S val = formic::zero(S())) {
        formic::Matrix<S> old = *this;
        this->reset(this->rows() + p, this->cols() + q, val);
        for (size_t j = 0; j < old.cols(); j++)
          std::copy(old.col_begin(j), old.col_end(j), this->col_begin(j));
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Resize the matrix to new size, while keeping the old elements 
      ///        
      /// \param[in]      p       number of rows of the new matrix
      /// \param[in]      q       number of columns of the new matrix
      /// \param[in]      val     value to initialize new elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void conservativeResize(const size_t p, const size_t q, const S val = formic::zero(S())) {
        formic::Matrix<S> old = *this;

        // if the new matrix is of the same size of the old matrix, then do nothing 
        if ( p == this->rows() && q == this->cols() ) 
          return;

        // record old cols and rows
        size_t old_row = this->rows();
        size_t old_col = this->cols();

        // reset the matrix 
        this->reset(p, q, val);

        // copy data to the new matrix
        if ( old_col < q ) {
          for ( size_t j = 0; j < old_col; j++ ) {
            if ( old_row < p ) 
              std::copy(old.col_begin(j), old.col_end(j), this->col_begin(j));
            else
              std::copy(old.col_begin(j), old.col_begin(j)+p, this->col_begin(j));
          }
        }
        else {
          for ( size_t j = 0; j < q; j++) {
            if ( old_row < p )
              std::copy(old.col_begin(j), old.col_end(j), this->col_begin(j));
            else
              std::copy(old.col_begin(j), old.col_begin(j)+p, this->col_begin(j));
          }
        }
      }


      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Add the matrix to the specified part of this matrix beginning at element (i,j).
      ///
      /// \param[in]      i        row of position to add matrix to
      /// \param[in]      j        col of position to add matrix to
      /// \param[in]      a        factor to scale matrix by during addition
      /// \param[in]      mat      the matrix to add to a subsection of this matrix
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void axpy_submatrix(const size_t i, const size_t j, const S a, const formic::ConstMatrix<S> & mat) {
        if ( mat.size() == 0 )
          throw formic::Exception("illegal attempt to use a %i by %i matrix as the added matrix in axpy_submatrix") % mat.rows() % mat.cols();
        if ( this->size() == 0 )
          throw formic::Exception("illegal attempt to axpy_submatrix into a %i by %i matrix") % this->rows() % this->cols();
        if ( i + mat.rows() > this->rows() )
          throw formic::Exception("cannot axpy submatrix as row index %i is out of bounds when adding a %i by %i matrix into a subsection of a %i by %i matrix")
                % i % mat.rows() % mat.cols() % this->rows() % this->cols();
        if ( j + mat.cols() > this->cols() )
          throw formic::Exception("cannot axpy submatrix as col index %i is out of bounds when adding a %i by %i matrix into a subsection of a %i by %i matrix")
                % j % mat.rows() % mat.cols() % this->rows() % this->cols();
        // add one column at a time
        for (size_t k = j; k < j + mat.cols(); k++)
          formic::xaxpy(mat.rows(), a, mat.col_begin(k-j), 1, this->col_begin(k) + i, 1);
      }

      // !!!!!    DEPRECATED!  See inv function in ConstMatrix    !!!!!!
      ///// \brief  invert the matrix
      //formic::Matrix<S> & invert() {
      //  if ( m_n != m_m )
      //    throw formic::Exception("cannot invert a non-square matrix");
      //  if ( this->size() > 0 ) {
      //    std::vector<int> iwork(2 * m_n);
      //    formic::Matrix<S> work = this->clone();
      //    S det;
      //    formic::matrix_inverse_lu(m_n, det, this->begin(), work->begin(), &iwork(0));
      //  }
      //  return *this;
      //}

  };

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################    VectorBase class     ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  A class that will be the base for column and row vectors.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class VectorBase : public formic::Matrix<S> {

    //////////   data members   //////////
    protected:

      // due to template base class and two phase lookup, we need to tell the compiler that these exist
      using formic::ConstMatrix<S>::m_n;
      using formic::ConstMatrix<S>::m_m;
      using formic::ConstMatrix<S>::m_p;

    //////////   constructors   //////////
    protected:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  default constructor is protected because it should only be called by its children.
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      VectorBase() {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copy constructor is protected because it should only be called by its children.
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      VectorBase(const VectorBase & other) : formic::Matrix<S>(other) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Constructor from a matrix is protected because it should only be called by its children.
      ///
      /// \param[in]      mat      the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      VectorBase(const formic::Matrix<S> & mat) : formic::Matrix<S>(mat) {}

    //////////   assignment operators   //////////
    private:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment is private as it should not be called
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      VectorBase & operator=(const VectorBase & other) {
        throw formic::Exception("VectorBase assignment operator should not be called");
        return *this;
      }

    //////////   element access   //////////
    public:

      // note that these hide the matrix element access functions, as desired

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a reference to element i element after checking that i is in bounds.
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S & at(const size_t i) {
        if ( i >= this->size() ) throw formic::Exception("out-of-range row index i = %i in formic::VectorBase::at for vector of dimension %i by %i") % i % m_n % m_m;
        return m_p[i];
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a const reference to element i element after checking that i is in bounds.
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const S & at(const size_t i) const {
        if ( i >= this->size() ) throw formic::Exception("out-of-range row index i = %i in formic::VectorBase::at for vector of dimension %i by %i") % i % m_n % m_m;
        return m_p[i];
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a reference to element i
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S & operator()(const size_t i) {
        #ifndef NDEBUG
        return this->at(i); // do bounds checking if debugging
        #else
        return m_p[i];
        #endif
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a const reference to element i
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const S & operator()(const size_t i) const {
        #ifndef NDEBUG
        return this->at(i); // do bounds checking if debugging
        #else
        return m_p[i];
        #endif
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a reference to element i
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S & operator[](const size_t i) {
        #ifndef NDEBUG
        return this->at(i); // do bounds checking if debugging
        #else
        return m_p[i];
        #endif
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Returns a const reference to element i
      ///
      /// \param[in]      i        the row index
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      const S & operator[](const size_t i) const {
        #ifndef NDEBUG
        return this->at(i); // do bounds checking if debugging
        #else
        return m_p[i];
        #endif
      }

    //////////   other public members   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Returns the element with the largest absolute value and also gives its index.
      ///        Note we return the element itself, not its absolute value.
      ///
      /// \param[out]     i        on exit, the index of the desired element
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S max_abs_elem(size_t & i) const {
        const S * const best = this->max_abs_elem_iter();
        size_t j, k;
        this->row_col_from_iter(best, j, k);
        i = std::max(j, k);
        return *best;
      }

      /// \brief Returns the element with the largest absolute value.  Note we return the element itself, not its absolute value.
      S max_abs_elem() const {
        return *this->max_abs_elem_iter();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Returns the element with the smallest absolute value and also gives its index.
      ///        Note we return the element itself, not its absolute value.
      ///
      /// \param[out]     i        on exit, the index of the desired element
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      S min_abs_elem(size_t & i) const {
        const S * const best = this->min_abs_elem_iter();
        size_t j, k;
        this->row_col_from_iter(best, j, k);
        i = std::max(j, k);
        return *best;
      }

      /// \brief Returns the element with the smallest absolute value.  Note we return the element itself, not its absolute value.
      S min_abs_elem() const {
        return *this->min_abs_elem_iter();
      }

  };

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################       ColVec class      ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  A class for column vectors, which are n by 1 matrices that inherit most of their
  ///         features from VectorBase.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class ColVec : public formic::VectorBase<S> {

    //////////   data members   //////////
    protected:

      // due to template base class and two phase lookup, we need to tell the compiler that these exist
      using formic::ConstMatrix<S>::m_is_col_vec;
      using formic::ConstMatrix<S>::m_is_row_vec;

    //////////   initializers   //////////
    private:

      /// \brief  sets the boolean flags appropriate for ColVec
      void init_vector_flags() {
        m_is_col_vec = true;
        m_is_row_vec = false;
      }

    //////////   constructors   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor without initialization of elements (also the defualt constructor)
      ///
      /// \param[in]      n        the desired number of rows
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ColVec(const size_t n = 0) {
        this->initialize_with_internal_storage(n, 1);
        this->init_vector_flags();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor with initialization of elements
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      v        value to initialize all elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ColVec(const size_t n, const S v) {
        this->initialize_with_internal_storage(n, 1);
        this->init_vector_flags();
        *this = v;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the vector around the specified, externally held array
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      p        array holding the matrix elements
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ColVec(const size_t n, S * const p) {
        this->initialize_with_external_storage(n, 1, p);
        this->init_vector_flags();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the matrix around the specified, externally held array, initializing
      ///         all elements to the specified value
      ///
      /// \param[in]      n        the desired number of rows
      /// \param[in]      p        array holding the matrix elements
      /// \param[in]      v        value to initialize all matrix elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit ColVec(const size_t n, S * const p, const S v) {
        this->initialize_with_external_storage(n, 1, p);
        this->init_vector_flags();
        *this = v;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copy constructor performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ColVec(const ColVec & other) : formic::VectorBase<S>(other) {
        this->init_vector_flags();
        this->check_vector_sanity("ColVec(ColVec)");
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Constructor from a matrix, which allows implicit conversion from Matrix to ColVec.
      ///         Error if mat is not a single column.
      ///
      /// \param[in]      mat      the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ColVec(const formic::Matrix<S> & mat) : formic::VectorBase<S>(mat) {
        if ( mat.is_row_vec() )
          throw formic::Exception("illegal attempt to construct ColVec from RowVec");
        this->init_vector_flags();
        this->check_vector_sanity("ColVec(Matrix)");
      }

    //////////   assignment operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ColVec & operator=(const ColVec & other) {
        ConstMatrix<S>::operator=(other);
        this->check_vector_sanity("ColVec::operator=(ColVec)");
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and mat will have the
      ///         same dimensions and will point to the same data array.
      ///         Error if mat is not a single column.
      ///
      /// \param[in]      mat      the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ColVec & operator=(const formic::Matrix<S> & mat) {
        if ( mat.is_row_vec() )
          throw formic::Exception("illegal attempt to assign RowVec to ColVec");
        ConstMatrix<S>::operator=(mat);
        this->check_vector_sanity("ColVec::operator=(Matrix)");
        return *this;
      }

      /// \brief  sets all elements equal to a
      ColVec & operator=(const S a) {
        std::fill(this->begin(), this->end(), a);
        return *this;
      }

    //////////   other public members   //////////
    public:

      // note that the Matrix::reset functions are hidden by those below, as desired

      /// \brief  clears the vector and re-allocates it to be the specified size
      ColVec & reset(const int n) {
        this->clear();
        *this = formic::ColVec<S>(n);
        this->check_vector_sanity("reset(n)");
        return *this;
      }

      /// \brief  clears the vector and re-allocates it to be the specified size, initializing its elements to v
      ColVec & reset(const int n, const S v) {
        this->clear();
        *this = formic::ColVec<S>(n, v);
        this->check_vector_sanity("reset(n,v)");
        return *this;
      }

      /// \brief  returns a new clone of the vector
      ColVec clone() const { return formic::ConstMatrix<S>::clone(); }

      /// \brief  resize the vector while keeping old elements
      void conservativeResize(const size_t p, const S val = formic::zero(S())) {
        formic::ColVec<S> old = *this;

        // if the new matrix is of the same size of the old matrix, then do nothing
        if ( p == this->size() )
          return;

        // record old size
        size_t old_size = this->size();

        // reset the vector
        this->reset(p, val);

        // copy data to the new vector
        for ( size_t j = 0; j < (old_size < p ? old_size : p); j++) 
          this->at(j) = old.at(j);
      }

  };

  //###############################################################################################\\
  //###############################################################################################\\
  //##############################                         ########################################\\
  //##############################       RowVec class      ########################################\\
  //##############################                         ########################################\\
  //###############################################################################################\\
  //###############################################################################################\\

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  A class for column vectors, which are n by 1 matrices that inherit most of their
  ///         features from VectorBase.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> class RowVec : public formic::VectorBase<S> {

    //////////   data members   //////////
    protected:

      // due to template base class and two phase lookup, we need to tell the compiler that these exist
      using formic::ConstMatrix<S>::m_is_col_vec;
      using formic::ConstMatrix<S>::m_is_row_vec;

    //////////   initializers   //////////
    private:

      /// \brief  sets the boolean flags appropriate for RowVec
      void init_vector_flags() {
        m_is_col_vec = false;
        m_is_row_vec = true;
      }

    //////////   constructors   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor without initialization of elements (also the defualt constructor)
      ///
      /// \param[in]      n        the desired number of columns
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit RowVec(const size_t n = 0) {
        this->initialize_with_internal_storage(1, n);
        this->init_vector_flags();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  constructor with initialization of elements
      ///
      /// \param[in]      n        the desired number of columns
      /// \param[in]      v        value to initialize all elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit RowVec(const size_t n, const S v) {
        this->initialize_with_internal_storage(1, n);
        this->init_vector_flags();
        *this = v;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the vector around the specified, externally held array
      ///
      /// \param[in]      n        the desired number of columns
      /// \param[in]      p        array holding the matrix elements
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit RowVec(const size_t n, S * const p) {
        this->initialize_with_external_storage(1, n, p);
        this->init_vector_flags();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  construct the matrix around the specified, externally held array, initializing
      ///         all elements to the specified value
      ///
      /// \param[in]      n        the desired number of columns
      /// \param[in]      p        array holding the matrix elements
      /// \param[in]      v        value to initialize all matrix elements to
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit RowVec(const size_t n, S * const p, const S v) {
        this->initialize_with_external_storage(1, n, p);
        this->init_vector_flags();
        *this = v;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Copy constructor performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      RowVec(const RowVec & other) : formic::VectorBase<S>(other) {
        this->init_vector_flags();
        this->check_vector_sanity("RowVec(RowVec)");
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Constructor from a matrix, which allows implicit conversion from Matrix to RowVec.
      ///         Error if mat is not a single row.
      ///
      /// \param[in]      mat      the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      RowVec(const formic::Matrix<S> & mat) : formic::VectorBase<S>(mat) {
        if ( mat.is_col_vec() )
          throw formic::Exception("illegal attempt to construct RowVec from ColVec");
        this->init_vector_flags();
        this->check_vector_sanity("RowVec(Matrix)");
      }

    //////////   assignment operators   //////////
    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and other will have the
      ///         same dimensions and will point to the same data array.
      ///
      /// \param[in]      other    the vector to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      RowVec & operator=(const RowVec & other) {
        ConstMatrix<S>::operator=(other);
        this->check_vector_sanity("RowVec::operator=(RowVec)");
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief  Assignment performs a shallow copy, after which this and mat will have the
      ///         same dimensions and will point to the same data array.
      ///         Error if mat is not a single column.
      ///
      /// \param[in]      mat      the matrix to copy from
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      RowVec & operator=(const formic::Matrix<S> & mat) {
        if ( mat.is_col_vec() )
          throw formic::Exception("illegal attempt to assign ColVec to RowVec");
        ConstMatrix<S>::operator=(mat);
        this->check_vector_sanity("RowVec::operator=(Matrix)");
        return *this;
      }

      /// \brief  sets all elements equal to a
      RowVec & operator=(const S a) {
        std::fill(this->begin(), this->end(), a);
        return *this;
      }

    //////////   other public members   //////////
    public:

      // note that the Matrix::reset functions are hidden by those below, as desired

      /// \brief  clears the vector and re-allocates it to be the specified size
      RowVec & reset(const int n) {
        this->clear();
        *this = formic::RowVec<S>(n);
        this->check_vector_sanity("reset(n)");
        return *this;
      }

      /// \brief  clears the vector and re-allocates it to be the specified size, initializing its elements to v
      RowVec & reset(const int n, const S v) {
        this->clear();
        *this = formic::RowVec<S>(n, v);
        this->check_vector_sanity("reset(n,v)");
        return *this;
      }

      /// \brief  returns a new clone of the vector
      RowVec clone() const { return formic::ConstMatrix<S>::clone(); }

  };

  //###############################################################################################\\
  //###############################################################################################\\
  //###########################                                      ##############################\\
  //###########################       Some Functions on Matrices     ##############################\\
  //###########################                                      ##############################\\
  //###############################################################################################\\
  //###############################################################################################\\

  /// \brief  check that the two matrices' dimensions allow them to be used together in matrix matrix multiplication
  template<class S> inline void check_dimensions_for_matrix_multiply(const formic::ConstMatrix<S> & m1, const formic::ConstMatrix<S> & m2) {
    if ( m1.size() == 0 )
      throw formic::Exception("cannot multiply matrices: found a zero dimension in the first matrix, which is %i by %i") % m1.rows() % m1.cols();
    if ( m2.size() == 0 )
      throw formic::Exception("cannot multiply matrices: found a zero dimension in the second matrix, which is %i by %i") % m2.rows() % m2.cols();
    if ( m1.cols() != m2.rows() )
      throw formic::Exception("cannot multiply matrices: inner dimensions of (%i by %i) and (%i by %i) matrices do not match.") % m1.rows() % m1.cols() % m2.rows() % m2.cols();
  }

  /// \brief  return the n by n identity matrix
  template<class S> inline formic::Matrix<S> identity_matrix(const size_t n) {
    formic::Matrix<S> retval(n, n, formic::zero(S()));
    for (size_t i = 0; i < n; i++)
      retval(i,i) = formic::unity(S());
    return retval;
  }

  /// \brief  return the square root of the sum of the square moduli of the matrx's elements
  template<class S> inline double abs(const formic::ConstMatrix<S> & mat) {
    return std::abs( std::sqrt( double( std::abs( mat.norm2() ) ) ) );
  }

  template<class S> formic::Matrix<S> matrix_exponent(const formic::ConstMatrix<S> & m, const double tol);

  template<class S> formic::Matrix<S> matrix_exponent_der_adj(const formic::ConstMatrix<S> & m, const formic::ConstMatrix<S> & a, const double tol);

  template<class S> S dotc(const formic::ConstMatrix<S> & mat1, const formic::ConstMatrix<S> & mat2);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  creates and returns a matrix with random, uniformly distributed elements
  ///
  /// \param[in]      n      number of rows for the matrix
  /// \param[in]      m      number of cols for the matrix
  /// \param[in]      lb     lower bound for random number range (default is -1.0)
  /// \param[in]      ub     upper bound for random number range (default if  1.0)
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> formic::Matrix<S> randMatrix(const int n, const int m, const double lb = -1.0, const double ub = 1.0) {
    formic::Matrix<S> retval(n, m);
    retval.random_fill(lb, ub);
    return retval;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  creates and returns a column vector with random, uniformly distributed elements
  ///
  /// \param[in]      n      number of rows for the vector
  /// \param[in]      lb     lower bound for random number range (default is -1.0)
  /// \param[in]      ub     upper bound for random number range (default if  1.0)
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> formic::Matrix<S> randColVec(const int n, const double lb = -1.0, const double ub = 1.0) {
    formic::ColVec<S> retval(n);
    retval.random_fill(lb, ub);
    return retval;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  creates and returns a row vector with random, uniformly distributed elements
  ///
  /// \param[in]      n      number of rows for the vector
  /// \param[in]      lb     lower bound for random number range (default is -1.0)
  /// \param[in]      ub     upper bound for random number range (default if  1.0)
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class S> formic::Matrix<S> randRowVec(const int n, const double lb = -1.0, const double ub = 1.0) {
    formic::RowVec<S> retval(n);
    retval.random_fill(lb, ub);
    return retval;
  }

} // end of namespace formic

//###############################################################################################\\
//###############################################################################################\\
//###########################                                    ################################\\
//###########################     Matrix Arithmetic Operators    ################################\\
//###########################                                    ################################\\
//###############################################################################################\\
//###############################################################################################\\

/// \brief  the unary + operator returns a deep copy of the matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::operator+() const {
  return this->clone();
}

/// \brief  the unary - operator returns -1 times a deep copy of the matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::operator-() const {
  formic::Matrix<S> retval(m_n,m_m);
  const S * jt = this->begin();
  for (S * it = retval.begin(); it != retval.end(); it++, jt++)
    *it = -(*jt);
  return retval;
}

/// \brief  this binary + operator returns a new matrix that is the sum of the two arguments
template<class S> inline formic::Matrix<S> operator+(const formic::ConstMatrix<S> & m1, const formic::ConstMatrix<S> & m2) {
  return m1.clone() += m2;
}

/// \brief  this binary - operator returns a new matrix that is the difference of the two arguments
template<class S> inline formic::Matrix<S> operator-(const formic::ConstMatrix<S> & m1, const formic::ConstMatrix<S> & m2) {
  return m1.clone() -= m2;
}

/// \brief  this binary + operator returns a new matrix that is equal to the input matrix with x added to each of its elements
template<class S> inline formic::Matrix<S> operator+(const formic::ConstMatrix<S> & m1, const S x) {
  return m1.clone() += x;
}

/// \brief  this binary - operator returns a new matrix that is equal to the input matrix with x subtracted from each of its elements
template<class S> inline formic::Matrix<S> operator-(const formic::ConstMatrix<S> & m1, const S x) {
  return m1.clone() -= x;
}

/// \brief  this binary * operator returns a new matrix that is equal to the input matrix with each element multiplied by x
template<class S> inline formic::Matrix<S> operator*(const formic::ConstMatrix<S> & m1, const S x) {
  return m1.clone() *= x;
}

/// \brief  this binary / operator returns a new matrix that is equal to the input matrix with each element divided by x
template<class S> inline formic::Matrix<S> operator/(const formic::ConstMatrix<S> & m1, const S x) {
  return m1.clone() /= x;
}

/// \brief  this binary + operator returns a new matrix that is equal to the input matrix with x added to each element
template<class S> inline formic::Matrix<S> operator+(const S x, const formic::ConstMatrix<S> & m1) {
  return m1.clone() += x;
}

/// \brief  this binary - operator returns a new matrix, each of whose elements is equal to x minus the corresponding element in m1
template<class S> inline formic::Matrix<S> operator-(const S x, const formic::ConstMatrix<S> & m1) {
  return formic::Matrix<S>(m1.rows(), m1.cols(), x) -= m1;
}

/// \brief  this binary * operator returns a new matrix that is equal to the input matrix with each element multiplied by x
template<class S> inline formic::Matrix<S> operator*(const S x, const formic::ConstMatrix<S> & m1) {
  return m1.clone() *= x;
}

/// \brief  this binary / operator returns a new matrix, each of whose elements is equal to x divided by the corresponding element in m1
template<class S> inline formic::Matrix<S> operator/(const S x, const formic::ConstMatrix<S> & m1) {
  formic::Matrix<S> retval(m1.rows(), m1.cols());
  const S * jt = m1.begin();
  for (S * it = retval.begin(); it != retval.end(); it++, jt++)
    *it = x / (*jt);
  return retval;
}

/// \brief  this binary * operator returns a new matrix that is equal to the matrix product between the two input matrices
template<class S> inline formic::Matrix<S> operator*(const formic::ConstMatrix<S> & m1, const formic::ConstMatrix<S> & m2) {
  formic::check_dimensions_for_matrix_multiply(m1, m2);
  formic::Matrix<S> retval(m1.rows(), m2.cols(), formic::zero(S()));
  for ( int j = 0; j < m2.cols(); j++)
  for ( int k = 0; k < m1.cols(); k++)
  for ( int i = 0; i < m1.rows(); i++)
    retval(i,j) += m1(i,k) * m2(k,j);
  return retval;
}

/// \brief  this binary * operator returns a new matrix that is equal to the matrix product between the two input matrices
template<> inline formic::Matrix<double> operator*(const formic::ConstMatrix<double> & m1, const formic::ConstMatrix<double> & m2) {
  formic::check_dimensions_for_matrix_multiply(m1, m2);
  formic::Matrix<double> retval(m1.rows(), m2.cols(), formic::zero(double()));
  formic::xgemm('N', 'N', m1.rows(), m2.cols(), m1.cols(), formic::unity(double()),
                m1.begin(), m1.rows(), m2.begin(), m2.rows(), formic::zero(double()), retval.begin(), retval.rows());
  return retval;
}

/// \brief  this binary * operator returns a new matrix that is equal to the matrix product between the two input matrices
template<> inline formic::Matrix<std::complex<double> > operator*(const formic::ConstMatrix<std::complex<double> > & m1, const formic::ConstMatrix<std::complex<double> > & m2) {
  formic::check_dimensions_for_matrix_multiply(m1, m2);
  formic::Matrix<std::complex<double> > retval(m1.rows(), m2.cols(), formic::zero(std::complex<double>()));
  formic::xgemm('N', 'N', m1.rows(), m2.cols(), m1.cols(), formic::unity(std::complex<double>()),
                m1.begin(), m1.rows(), m2.begin(), m2.rows(), formic::zero(std::complex<double>()), retval.begin(), retval.rows());
  return retval;
}

//template<class S> inline
//formic::Vector<S> operator*(const formic::ConstMatrix<S> & m1, const formic::Vector<S> & v) {
//  if ( m1.cols() != v.size() )
//    throw formic::Exception("cannot multiply: dimensions of (%i by %i) matrix and (%i by 1) vector do not match.")
//          % m1.rows() % m1.cols() % v.size();
//  if ( v.size() == 0 )
//    throw formic::Exception("cannot multiply: vector is zero length");
//  if ( m1.rows() == 0 )
//    throw formic::Exception("cannot multiply: matrix has no rows");
//  return formic::Vector<S>(m1.rows()) <<= &( m1 * formic::ConstMatrix<S>(v.size(), 1, &v.at(0)) ).at(0,0);
//}

//template<class S> inline
//formic::Matrix<S> operator*(const formic::Vector<S> & v, const formic::ConstMatrix<S> & m1) {
//  if ( m1.rows() != v.size() )
//    throw formic::Exception("cannot multiply: dimensions of (1 by %i) vector and (%i by %i) matrix do not match.")
//          % v.size() % m1.rows() % m1.cols();
//  if ( v.size() == 0 )
//    throw formic::Exception("cannot multiply: vector is zero length");
//  if ( m1.cols() == 0 )
//    throw formic::Exception("cannot multiply: matrix has no columns");
//  return formic::Vector<S>(m1.cols()) <<= &( formic::ConstMatrix<S>(1, v.size(), &v.at(0)) * m1 ).at(0,0);
//}

//###############################################################################################\\
//###############################################################################################\\
//###########################                                    ################################\\
//###########################     Member Function Definitions    ################################\\
//###########################                                    ################################\\
//###############################################################################################\\
//###############################################################################################\\

/// \brief  return a deep copy of the matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::clone() const {
  return formic::Matrix<S>() <<= (*this);
}

/// \brief  return a copy of the matrix in which all elements have been complex conjugated
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::conj() const {
  return this->clone().conjip();
}

/// \brief  return the matrix transpose in a new matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::t() const {
  formic::Matrix<S> retval(m_m, m_n);
  for (size_t i = 0; i < m_n; i++)
  for (size_t j = 0; j < m_m; j++)
    retval(j,i) = (*this)(i,j);
  return retval;
}

/// \brief  return the matrix conjugate transpose in a new matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::c() const {
  formic::Matrix<S> retval(m_m, m_n);
  for (size_t i = 0; i < m_n; i++)
  for (size_t j = 0; j < m_m; j++)
    retval(j,i) = formic::conj((*this)(i,j));
  return retval;
}

/// \brief  returns a new diagonal matrix whose diagonal elements are equal to those in this vector-like matrix
template<class S> formic::Matrix<S> formic::ConstMatrix<S>::as_diag_mat() const {
  if ( this->size() == 0 || ( m_n != 1 && m_m != 1 ) )
    throw formic::Exception("illegal attempt to convert a %i by %i matrix to a diagonal matrix") % m_n % m_m;
  formic::Matrix<S> retval(this->size(), this->size(), formic::zero(S()));
  for (size_t i = 0; i < this->size(); i++)
    retval.at(i,i) = m_p[i];
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a string holding the matrix formatted for printing
///
/// \param[in]      fmt      format string for how to format numbers, e.g. "%12.4f"
/// \param[in]      name     a name to call the matrix by
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> std::string formic::ConstMatrix<S>::print(const std::string & fmt, const std::string & name) const {
  std::string retval;
  if ( name.size() > 0 )
    retval += "printing matrix \"" + name + "\":\n";
  for (size_t i = 0; i < m_n; i++) {
    for (size_t j = 0; j < m_m; j++)
      retval += "  " + formic::format_number(fmt, (*this)(i,j));
    retval += "\n";
  }
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a matrix resulting from slicing out certain rows and columns from this matrix
///
///         Note that V must support the size() function and element access by operator[]
///
/// \param[in]      row_inds   a vector telling which rows to select
/// \param[in]      col_inds   a vector telling which columns to select
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> template<class V> formic::Matrix<S> formic::ConstMatrix<S>::slice(const V & row_inds, const V & col_inds) const {

  // get new matrix dimensions and check them for sanity
  const int n = row_inds.size();
  const int m = col_inds.size();
  if ( n > m_n )
    throw formic::Exception("too many row indices (%i) in formic::ConstMatrix::slice for matrix of dimension %i by %i") % n % m_n % m_m;
  if ( m > m_m )
    throw formic::Exception("too many col indices (%i) in formic::ConstMatrix::slice for matrix of dimension %i by %i") % m % m_n % m_m;

  // check row indices for sanity
  for (int i = 0; i < n; i++)
    if ( size_t(row_inds[i]) >= this->rows() )
      throw formic::Exception("row index %i is out of bounds for this %i by %i matrix in formic::ConstMatrix::slice") % row_inds[i] % m_n % m_m;

  // check column indices for sanity
  for (int i = 0; i < m; i++)
    if ( size_t(col_inds[i]) >= this->cols() )
      throw formic::Exception("column index %i is out of bounds for this %i by %i matrix in formic::ConstMatrix::slice") % col_inds[i] % m_n % m_m;

  // build new matrix
  formic::Matrix<S> retval(n, m);
  for (size_t j = 0; j < m; j++)
  for (size_t i = 0; i < n; i++)
    retval(i,j) = (*this)(row_inds[i], col_inds[j]);

  // return new matrix
  return retval;

}

//template<class S>
//formic::Matrix<S> formic::Vector<S>::t() const { return formic::Matrix<S>(m_n, 1, _p).t(); }

//template<class S>
//formic::Matrix<S> formic::Vector<S>::m() const { return formic::Matrix<S>(1, m_n, _p).t(); }

#endif
