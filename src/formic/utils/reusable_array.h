///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/reusable_array.h
///
/// \brief   header file for the formic::ReusableArray class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_VECTOR_REUSABLE_ARRAY_HEADER
#define FORMIC_VECTOR_REUSABLE_ARRAY_HEADER

#include <map>
#include <stack>

#include <boost/shared_array.hpp>

namespace formic {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  A class to provide reusable arrays that, upon destruction, have their underlying
  ///         data arrays kept available for future use to help avoid unecessary reallocations
  ///         of memory.  A garbage collection function is available to force all currently
  ///         unused arrays to be deallocated.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class DataType> class ReusableArray {

    private:

      /// \brief  length of the array
      size_t m_length;

      /// \brief  smart pointer to the underlying data array
      boost::shared_array<DataType> m_array;

      /// \brief  a container for data arrays that are currently unused
      static std::map<size_t, std::stack<boost::shared_array<DataType> > > available_arrays;

    private:

      // disable assignment operator
      ReusableArray & operator=(const ReusableArray &);

      // disable copy construction
      ReusableArray(const ReusableArray &);

      static void thread_lock() {}

      static void thread_unlock() {}

    public:

      /// \brief  return the length
      size_t length() { return m_length; }

      /// \brief  return a pointer to the underlying data array
      DataType * cptr() { return m_array.get(); }

      ReusableArray(const size_t n);

      ~ReusableArray();

      static void garbage_collect();

  };

  void reusable_array_garbage_collect();

}

#endif
