///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/mpi/interface.cpp
///
/// \brief   definitions of mpi interface functions
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "formic/utils/mpi_interface.h"
#include "formic/utils/archive.h"

//#include <unistd.h> // needed for gethostname

//MPI_Comm MPI_COMM_WORLD;
//MPI_Op MPI_SUM;
//MPI_Datatype MPI_CHAR;
//MPI_Datatype MPI_SHORT;
//MPI_Datatype MPI_INT;
//MPI_Datatype MPI_LONG;
//MPI_Datatype MPI_SIGNED_CHAR;
//MPI_Datatype MPI_UNSIGNED_CHAR;
//MPI_Datatype MPI_UNSIGNED_SHORT;
//MPI_Datatype MPI_UNSIGNED;
//MPI_Datatype MPI_UNSIGNED_LONG;
//MPI_Datatype MPI_FLOAT;
//MPI_Datatype MPI_DOUBLE;
//MPI_Datatype MPI_LONG_DOUBLE;
//MPI_Datatype MPI_BOOL;
//MPI_Datatype MPI_COMPLEX;
//MPI_Datatype MPI_DOUBLE_COMPLEX;
//MPI_Datatype MPI_LONG_DOUBLE_COMPLEX;


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   initializes MPI
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::mpi::init(int argc, char **argv) {

//  MPI_Init(argc, argv);

  const int requested = MPI_THREAD_SERIALIZED;
  int provided = 0;
  MPI_Init_thread(&argc, &argv, requested, &provided);
  if ( provided != requested) {
    //std::printf("MPI_THREAD_SINGLE     = %10i\n", MPI_THREAD_SINGLE);
    //std::printf("MPI_THREAD_FUNNELED   = %10i\n", MPI_THREAD_FUNNELED);
    //std::printf("MPI_THREAD_SERIALIZED = %10i\n", MPI_THREAD_SERIALIZED);
    //std::printf("MPI_THREAD_MULTIPLE   = %10i\n", MPI_THREAD_MULTIPLE);
    //std::printf("requested             = %10i\n", requested);
    //std::printf("provided              = %10i\n", provided);
    throw formic::Exception("requested mpi level of thread support not provided");
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   finalizes MPI
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::mpi::finalize() {
  MPI_Finalize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///// \brief   function to return the hostname for this process
/////
///// \return the hostname as given by gethostname
/////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//std::string formic::mpi::get_hostname() {
//
//  // allocate vector to hold the character string
//  std::vector<char> hostname_vec(1024);
//
//  // load the hostname into the character string
//  const int rc = gethostname(&hostname_vec.at(0), 1000);
//  if (rc != 0)
//    throw formic::Exception("gethostname failed in formic::mpi::get_hostname with error code %i") % rc;
//
//  // return the hostname
//  return std::string(&hostname_vec.at(0));
//
//}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   broadcasts data from one process to others
///
/// \param[in,out] data     a pointer to the data to be broadcast
/// \param[in]     n        the number of data elements
/// \param[in]     root     which processor to broadcast from
/// \param[in]     comm     the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::bcast(T * data, size_t n, int root, const MPI_Comm & comm) {
  const int block_size = 100000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Bcast((void *)data, block_size, formic::mpi::datatype<T>(), root, comm);
    data += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Bcast((void *)data, int(n), formic::mpi::datatype<T>(), root, comm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   reduces data on to one process
///
/// \param[in]     send_buf   a pointer to the data to be reduced
/// \param[out]    recv_buf   a pointer to where the reduced result will be stored
/// \param[in]     n          the number of data elements
/// \param[in]     op         the reduce operation to use
/// \param[in]     root       which processor to reduce to
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::reduce(const T * send_buf, T * recv_buf, size_t n,
                                         const MPI_Op & op, int root, const MPI_Comm & comm) {

  const int block_size = 100000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Reduce(const_cast<T*>(send_buf), recv_buf, block_size, formic::mpi::datatype<T>(), op, root, comm);
    send_buf += block_size;
    recv_buf += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Reduce(const_cast<T*>(send_buf), recv_buf, int(n), formic::mpi::datatype<T>(), op, root, comm);

  //std::copy(send_buf, send_buf+n, recv_buf);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   reduces data on to all processes
///
/// \param[in]     send_buf   a pointer to the data to be reduced
/// \param[out]    recv_buf   a pointer to where the reduced result will be stored
/// \param[in]     n          the number of data elements
/// \param[in]     op         the reduce operation to use
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::allreduce(const T * send_buf, T * recv_buf, size_t n,
                                            const MPI_Op & op, const MPI_Comm & comm) {

  const int block_size = 100000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Allreduce(const_cast<T*>(send_buf), recv_buf, block_size, formic::mpi::datatype<T>(), op, comm);
    send_buf += block_size;
    recv_buf += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Allreduce(const_cast<T*>(send_buf), recv_buf, int(n), formic::mpi::datatype<T>(), op, comm);

  //std::copy(send_buf, send_buf+n, recv_buf);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   sends data to another process
///
/// \param[in]     buf        a pointer to the data to be sent
/// \param[in]     n          the number of data elements
/// \param[in]     dest       rank of the process to send to
/// \param[in]     tag        a tag to label the transmission
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::send(const T * buf, size_t n, const int dest, const int tag, const MPI_Comm & comm) {

  const int block_size = 100000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Send(const_cast<T*>(buf), block_size, formic::mpi::datatype<T>(), dest, tag, comm);
    buf += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Send(const_cast<T*>(buf), int(n), formic::mpi::datatype<T>(), dest, tag, comm);

  //throw formic::Exception("cannot use formic::mpi::send without MPI");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   receive data from another process
///
/// \param[out]    buf        a pointer to where the received data will be placed
/// \param[in]     n          the number of data elements
/// \param[in]     source     rank of the sending process
/// \param[in]     tag        a tag to label the transmission
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::recv(T * buf, size_t n, const int source, const int tag, const MPI_Comm & comm) {

  const int block_size = 100000000/sizeof(T);
  MPI_Status status;
  while ( n > size_t(block_size) ) {
    MPI_Recv(buf, block_size, formic::mpi::datatype<T>(), source, tag, comm, &status);
    buf += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Recv(buf, int(n), formic::mpi::datatype<T>(), source, tag, comm, &status);

  //throw formic::Exception("cannot use formic::mpi::recv without MPI");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   scatter data from one process to others
///
/// \param[in]     send_buf   a pointer to the data to be scattered
/// \param[out]    recv_buf   a pointer to where the scattered data will be placed
/// \param[in]     n          the number of data elements each process should receive
/// \param[in]     root       which processor to scatter from
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::scatter(const T * send_buf, T * recv_buf, size_t n, int root, const MPI_Comm & comm) {

  const int block_size = 1000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Scatter(const_cast<T*>(send_buf), block_size, formic::mpi::datatype<T>(), recv_buf, block_size, formic::mpi::datatype<T>(), root, comm);
    int size;
    MPI_Comm_size(comm, &size);
    send_buf += size_t(size) * block_size;
    recv_buf += block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Scatter(const_cast<T*>(send_buf), int(n), formic::mpi::datatype<T>(), recv_buf, int(n), formic::mpi::datatype<T>(), root, comm);

  //std::copy(send_buf, send_buf+n, recv_buf);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   gather data from all processes on to one
///
/// \param[in]     send_buf   a pointer to the data to be gathered
/// \param[out]    recv_buf   a pointer to where the gathered data will be placed
/// \param[in]     n          the number of data elements each process will send
/// \param[in]     root       which processor to gather to
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::gather(const T * send_buf, T * recv_buf, size_t n, int root, const MPI_Comm & comm) {

  const int block_size = 1000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Gather(const_cast<T*>(send_buf), block_size, formic::mpi::datatype<T>(), recv_buf, block_size, formic::mpi::datatype<T>(), root, comm);
    send_buf += block_size;
    int size;
    MPI_Comm_size(comm, &size);
    recv_buf += size_t(size) * block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Gather(const_cast<T*>(send_buf), int(n), formic::mpi::datatype<T>(), recv_buf, int(n), formic::mpi::datatype<T>(), root, comm);

  //std::copy(send_buf, send_buf+n, recv_buf);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   gather data from all processes and place the result on all processes
///
/// \param[in]     send_buf   a pointer to the data to be gathered
/// \param[out]    recv_buf   a pointer to where the gathered data will be placed
/// \param[in]     n          the number of data elements each process will send
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::allgather(const T * send_buf, T * recv_buf, size_t n, const MPI_Comm & comm) {

  const int block_size = 1000000/sizeof(T);
  while ( n > size_t(block_size) ) {
    MPI_Allgather(const_cast<T*>(send_buf), block_size, formic::mpi::datatype<T>(), recv_buf, block_size, formic::mpi::datatype<T>(), comm);
    send_buf += block_size;
    int size;
    MPI_Comm_size(comm, &size);
    recv_buf += size_t(size) * block_size;
    n -= size_t(block_size);
  }
  if (n > 0)
    MPI_Allgather(const_cast<T*>(send_buf), int(n), formic::mpi::datatype<T>(), recv_buf, int(n), formic::mpi::datatype<T>(), comm);

  //std::copy(send_buf, send_buf+n, recv_buf);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   broadcasts a vector from one process to others
///
/// \param[in,out] v        the vector to be broadcast
/// \param[in]     root     which processor to broadcast from
/// \param[in]     comm     the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::bcast(std::vector<T> & v, int root, const MPI_Comm & comm) {
  size_t n = v.size();
  formic::mpi::bcast(&n, 1, root, comm);
  v.resize(n);
  if (n > 0)
    formic::mpi::bcast(&v.at(0), n, root, comm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   reduces a vector on to one process
///
/// \param[in]     send_v     the vector to reduce
/// \param[out]    recv_v     the vector where the result of the reduction will be stored
/// \param[in]     op         the reduce operation to use
/// \param[in]     root       which processor to reduce to
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::reduce(const std::vector<T> & send_v, std::vector<T> & recv_v,
                                         const MPI_Op & op, int root, const MPI_Comm & comm) {
  const size_t n = send_v.size();
  if (formic::mpi::rank() == root)
    recv_v.resize(n);
  if (n > 0)
    formic::mpi::reduce(&send_v.at(0), &recv_v.at(0), n, op, root, comm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   reduces a vector on to all processes
///
/// \param[in]     send_v     the vector to reduce
/// \param[out]    recv_v     the vector where the result of the reduction will be stored
/// \param[in]     op         the reduce operation to use
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::allreduce(const std::vector<T> & send_v, std::vector<T> & recv_v,
                                            const MPI_Op & op, const MPI_Comm & comm) {
  const size_t n = send_v.size();
  recv_v.resize(n);
  if (n > 0)
    formic::mpi::allreduce(&send_v.at(0), &recv_v.at(0), n, op, comm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   sends a vector to another process
///
/// \param[in]     v          the vector to send
/// \param[in]     dest       rank of the process to send to
/// \param[in]     tag        a tag to label the transmission
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::send(const std::vector<T> & v, const int dest, const int tag,
                                       const MPI_Comm & comm) {
  const size_t n = v.size();
  formic::mpi::send(&n, 1, dest, tag, comm);
  if (n > 0)
    formic::mpi::send(&v.at(0), n, dest, tag, comm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   receive a vector from another process
///
/// \param[out]    v          the vector that the received vector will overwrite
/// \param[in]     source     rank of the sending process
/// \param[in]     tag        a tag to label the transmission
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::recv(std::vector<T> & v, const int source, const int tag,
                                       const MPI_Comm & comm) {
  size_t n;
  formic::mpi::recv(&n, 1, source, tag, comm);
  v.resize(n);
  if (n > 0)
    formic::mpi::recv(&v.at(0), n, source, tag, comm);
}

template void formic::mpi::bcast(bool                  *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(char                  *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(signed short          *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(signed int            *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(signed long           *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(signed char           *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(unsigned char         *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(unsigned short        *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(unsigned int          *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(unsigned long int     *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(float                 *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(double                *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(long double           *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(std::complex<float>   *, size_t, int, const MPI_Comm &);
template void formic::mpi::bcast(std::complex<double>  *, size_t, int, const MPI_Comm &);

template void formic::mpi::reduce(const char                *, char                *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const signed short        *, signed short        *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const signed int          *, signed int          *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const signed long         *, signed long         *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const signed char         *, signed char         *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const unsigned char       *, unsigned char       *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const unsigned short      *, unsigned short      *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const unsigned int        *, unsigned int        *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const unsigned long int   *, unsigned long int   *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const float               *, float               *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const double              *, double              *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const long double         *, long double         *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::complex<float> *, std::complex<float> *, size_t, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::complex<double>*, std::complex<double>*, size_t, const MPI_Op&, int, const MPI_Comm&);

template void formic::mpi::allreduce(const char                *, char                *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const signed short        *, signed short        *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const signed int          *, signed int          *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const signed long         *, signed long         *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const signed char         *, signed char         *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const unsigned char       *, unsigned char       *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const unsigned short      *, unsigned short      *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const unsigned int        *, unsigned int        *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const unsigned long int   *, unsigned long int   *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const float               *, float               *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const double              *, double              *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const long double         *, long double         *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::complex<float> *, std::complex<float> *, size_t, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::complex<double>*, std::complex<double>*, size_t, const MPI_Op&, const MPI_Comm&);

template void formic::mpi::send(const bool                *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const char                *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const signed short        *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const signed int          *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const signed long         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const signed char         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const unsigned char       *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const unsigned short      *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const unsigned int        *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const unsigned long int   *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const float               *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const double              *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const long double         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const std::complex<float> *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::send(const std::complex<double>*, size_t, const int, const int, const MPI_Comm&);

template void formic::mpi::recv(bool                *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(char                *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(signed short        *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(signed int          *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(signed long         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(signed char         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(unsigned char       *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(unsigned short      *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(unsigned int        *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(unsigned long int   *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(float               *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(double              *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(long double         *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(std::complex<float> *, size_t, const int, const int, const MPI_Comm&);
template void formic::mpi::recv(std::complex<double>*, size_t, const int, const int, const MPI_Comm&);

template void formic::mpi::scatter(const bool                *, bool                *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const char                *, char                *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const signed short        *, signed short        *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const signed int          *, signed int          *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const signed long         *, signed long         *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const signed char         *, signed char         *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const unsigned char       *, unsigned char       *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const unsigned short      *, unsigned short      *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const unsigned int        *, unsigned int        *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const unsigned long int   *, unsigned long int   *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const float               *, float               *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const double              *, double              *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const long double         *, long double         *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const std::complex<float> *, std::complex<float> *, size_t, int, const MPI_Comm&);
template void formic::mpi::scatter(const std::complex<double>*, std::complex<double>*, size_t, int, const MPI_Comm&);

template void formic::mpi::gather(const bool                *, bool                *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const char                *, char                *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const signed short        *, signed short        *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const signed int          *, signed int          *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const signed long         *, signed long         *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const signed char         *, signed char         *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const unsigned char       *, unsigned char       *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const unsigned short      *, unsigned short      *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const unsigned int        *, unsigned int        *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const unsigned long int   *, unsigned long int   *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const float               *, float               *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const double              *, double              *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const long double         *, long double         *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const std::complex<float> *, std::complex<float> *, size_t, int, const MPI_Comm&);
template void formic::mpi::gather(const std::complex<double>*, std::complex<double>*, size_t, int, const MPI_Comm&);

template void formic::mpi::allgather(const bool                *, bool                *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const char                *, char                *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const signed short        *, signed short        *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const signed int          *, signed int          *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const signed long         *, signed long         *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const signed char         *, signed char         *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const unsigned char       *, unsigned char       *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const unsigned short      *, unsigned short      *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const unsigned int        *, unsigned int        *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const unsigned long int   *, unsigned long int   *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const float               *, float               *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const double              *, double              *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const long double         *, long double         *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const std::complex<float> *, std::complex<float> *, size_t, const MPI_Comm&);
template void formic::mpi::allgather(const std::complex<double>*, std::complex<double>*, size_t, const MPI_Comm&);

template void formic::mpi::bcast(std::vector<char                  > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<signed short          > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<signed int            > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<signed long           > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<signed char           > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<unsigned char         > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<unsigned short        > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<unsigned int          > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<unsigned long int     > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<float                 > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<double                > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<long double           > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<std::complex<float>   > &, int, const MPI_Comm &);
template void formic::mpi::bcast(std::vector<std::complex<double>  > &, int, const MPI_Comm &);

template void formic::mpi::reduce(const std::vector<char                 >&, std::vector<char                 >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<signed short         >&, std::vector<signed short         >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<signed int           >&, std::vector<signed int           >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<signed long          >&, std::vector<signed long          >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<signed char          >&, std::vector<signed char          >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<unsigned char        >&, std::vector<unsigned char        >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<unsigned short       >&, std::vector<unsigned short       >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<unsigned int         >&, std::vector<unsigned int         >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<unsigned long int    >&, std::vector<unsigned long int    >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<float                >&, std::vector<float                >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<double               >&, std::vector<double               >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<long double          >&, std::vector<long double          >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<std::complex<float>  >&, std::vector<std::complex<float>  >&, const MPI_Op&, int, const MPI_Comm&);
template void formic::mpi::reduce(const std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, const MPI_Op&, int, const MPI_Comm&);

template void formic::mpi::allreduce(const std::vector<char                 >&, std::vector<char                 >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<signed short         >&, std::vector<signed short         >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<signed int           >&, std::vector<signed int           >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<signed long          >&, std::vector<signed long          >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<signed char          >&, std::vector<signed char          >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<unsigned char        >&, std::vector<unsigned char        >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<unsigned short       >&, std::vector<unsigned short       >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<unsigned int         >&, std::vector<unsigned int         >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<unsigned long int    >&, std::vector<unsigned long int    >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<float                >&, std::vector<float                >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<double               >&, std::vector<double               >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<long double          >&, std::vector<long double          >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<std::complex<float>  >&, std::vector<std::complex<float>  >&, const MPI_Op&, const MPI_Comm&);
template void formic::mpi::allreduce(const std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, const MPI_Op&, const MPI_Comm&);

template void formic::mpi::send(const std::vector<char                  > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<signed short          > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<signed int            > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<signed long           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<signed char           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<unsigned char         > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<unsigned short        > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<unsigned int          > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<unsigned long int     > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<float                 > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<double                > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<long double           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<std::complex<float>   > &, const int, const int, const MPI_Comm &);
template void formic::mpi::send(const std::vector<std::complex<double>  > &, const int, const int, const MPI_Comm &);

template void formic::mpi::recv(std::vector<char                  > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<signed short          > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<signed int            > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<signed long           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<signed char           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<unsigned char         > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<unsigned short        > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<unsigned int          > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<unsigned long int     > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<float                 > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<double                > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<long double           > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<std::complex<float>   > &, const int, const int, const MPI_Comm &);
template void formic::mpi::recv(std::vector<std::complex<double>  > &, const int, const int, const MPI_Comm &);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   reads an object from an archive an broadcasts it
///
/// \param[in,out] arch      the archive to read from
/// \param[out]    val       the object that is read
/// \param[in]     error_msg  an error message to display if the read fails
/// \param[in]     root       which processor to broadcast from
/// \param[in]     comm       the communicator to use (defaults to MPI_COMM_WORLD)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> void formic::mpi::read_and_bcast(formic::Archive & arch,
                                                   T & val,
                                                   const std::string & error_msg,
                                                   int root,
                                                   const MPI_Comm & comm) {

  // get MPI info
  const int nproc  = formic::mpi::size(comm);
  const int myrank = formic::mpi::rank(comm);

  // read the object on the root process
  if (myrank == root)
    if ( !(arch >> val) )
      throw formic::Exception("%s") % error_msg;

  // broadcast the object
  formic::mpi::bcast(val, root, comm);

}

template void formic::mpi::read_and_bcast(formic::Archive &, bool &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, int &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, long int &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, double &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::complex<double> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::string &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::vector<int> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::vector<long int> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::vector<double> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::vector<std::complex<double> > &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, std::vector<std::string> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, formic::Matrix<double> &, const std::string &, int, const MPI_Comm &);
template void formic::mpi::read_and_bcast(formic::Archive &, formic::Matrix<std::complex<double> > &, const std::string &, int, const MPI_Comm &);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   broadcast a string
///
/// \param[inout]  datum      the item to be broadcast
/// \param[in]     root       the process from which to broadcast
/// \param[in]     comm       the communicator to use
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::mpi::bcast(std::string & s, int root, const MPI_Comm & comm) {
  std::vector<char> v(s.begin(), s.end());
  formic::mpi::bcast(v, root, comm);
  s = std::string(v.begin(), v.end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   broadcast a vector of strings
///
/// \param[inout]  datum      the item to be broadcast
/// \param[in]     root       the process from which to broadcast
/// \param[in]     comm       the communicator to use
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::mpi::bcast(std::vector<std::string> & v, int root, const MPI_Comm & comm) {
  size_t n = v.size();
  formic::mpi::bcast(n, root, comm);
  v.resize(n);
  for (size_t i = 0; i < n; i++)
    formic::mpi::bcast(v[i], root, comm);
}
