///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/reusable_array.cpp
///
/// \brief   implementation file for the formic::ReusableArray class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<complex>

#include<formic/utils/reusable_array.h>

template<class S> std::map<size_t, std::stack<boost::shared_array<S> > > formic::ReusableArray<S>::available_arrays;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Constructor first attempts to find an unused data array of the correct size.
///         Otherwise, a new data array is allocated.
///
/// \param[in]      n        desired length of the array
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::ReusableArray<S>::ReusableArray(const size_t n)
  : m_length(n),
    m_array()
{

  if ( n > 0 ) {

    this->thread_lock();

    // get the array stack for this length
    std::stack<boost::shared_array<S> > & st = available_arrays[n];

    // if there is an array of the correct length available, use it
    if ( st.size() > 0 ) {

      //if ( formic::mpi::rank() == 0 )
      //  formic::of << boost::format("using a stacked array of size %i") % m_length << std::endl;
      m_array = st.top();
      st.pop();

    // otherwise, allocate a new array
    } else {

      //if ( formic::mpi::rank() == 0 )
      //  formic::of << boost::format("allocating new array of size %i") % n << std::endl;
      boost::shared_array<S> temp_array( new S[n] );
      m_array = temp_array;

    }

    this->thread_unlock();

  }

}

template formic::ReusableArray<int>::ReusableArray(const size_t n);
template formic::ReusableArray<double>::ReusableArray(const size_t n);
template formic::ReusableArray<std::complex<double> >::ReusableArray(const size_t n);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Destructor saves the underlying data array for later (if it's length is not zero).
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> formic::ReusableArray<S>::~ReusableArray() {
  this->thread_lock();
  if (m_length > 0) {
    //if ( formic::mpi::rank() == 0 )
    //  formic::of << boost::format("stacking an array of size %i") % m_length << std::endl;
    available_arrays[m_length].push(m_array);
  }
  m_array.reset();
  this->thread_unlock();
}

template formic::ReusableArray<int>::~ReusableArray();
template formic::ReusableArray<double>::~ReusableArray();
template formic::ReusableArray<std::complex<double> >::~ReusableArray();

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Deallocates unused data arrays of this type.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::ReusableArray<S>::garbage_collect() {
  ReusableArray<S>::thread_lock();
  //if ( formic::mpi::rank() == 0 )
  //  for (typename std::map<size_t, std::stack<boost::shared_array<S> > >::iterator mi = available_arrays.begin(); mi != available_arrays.end(); mi++)
  //    formic::of << boost::format("deallocating %i arrays of size %i") % mi->second.size() % mi->first << std::endl;
  available_arrays.clear();
  ReusableArray<S>::thread_unlock();
}

template void formic::ReusableArray<int>::garbage_collect();
template void formic::ReusableArray<double>::garbage_collect();
template void formic::ReusableArray<std::complex<double> >::garbage_collect();

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Deallocates unused data arrays of all types.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::reusable_array_garbage_collect() {
  formic::ReusableArray<int>::garbage_collect();
  formic::ReusableArray<double>::garbage_collect();
  formic::ReusableArray<std::complex<double> >::garbage_collect();
}
