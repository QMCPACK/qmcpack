/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#ifndef MPI_PROTO_H_INCLUDED
#define MPI_PROTO_H_INCLUDED

/* Begin Prototypes */
int MPI_DUP_FN(MPI_Comm oldcomm, int keyval, void *extra_state, void *attribute_val_in,
               void *attribute_val_out, int *flag) MPICH_API_PUBLIC;
int MPI_Status_f082c(const MPI_F08_status *f08_status, MPI_Status *c_status) MPICH_API_PUBLIC;
int MPI_Status_c2f08(const MPI_Status *c_status, MPI_F08_status *f08_status) MPICH_API_PUBLIC;
int MPI_Status_f082f(const MPI_F08_status *f08_status, MPI_Fint *f_status) MPICH_API_PUBLIC;
int MPI_Status_f2f08(const MPI_Fint *f_status, MPI_F08_status *f08_status) MPICH_API_PUBLIC;
int MPI_Type_create_f90_integer(int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Attr_delete(MPI_Comm comm, int keyval) MPICH_API_PUBLIC;
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) MPICH_API_PUBLIC;
int MPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) MPICH_API_PUBLIC;
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                           MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                           void *extra_state) MPICH_API_PUBLIC;
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) MPICH_API_PUBLIC;
int MPI_Comm_free_keyval(int *comm_keyval) MPICH_API_PUBLIC;
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag)
    MPICH_API_PUBLIC;
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) MPICH_API_PUBLIC;
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                      void *extra_state) MPICH_API_PUBLIC;
int MPI_Keyval_free(int *keyval) MPICH_API_PUBLIC;
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                           MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                           void *extra_state) MPICH_API_PUBLIC;
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) MPICH_API_PUBLIC;
int MPI_Type_free_keyval(int *type_keyval) MPICH_API_PUBLIC;
int MPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag)
    MPICH_API_PUBLIC;
int MPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val)
    MPICH_API_PUBLIC;
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                          MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                          void *extra_state) MPICH_API_PUBLIC;
int MPI_Win_delete_attr(MPI_Win win, int win_keyval) MPICH_API_PUBLIC;
int MPI_Win_free_keyval(int *win_keyval) MPICH_API_PUBLIC;
int MPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag) MPICH_API_PUBLIC;
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) MPICH_API_PUBLIC;
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                       int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                   MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                        const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                  MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Allreduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[],
                  MPI_Datatype recvtype, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                       MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                       const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                  const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                  MPICH_API_PUBLIC;
int MPI_Alltoallw_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                       const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                       const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                       MPI_Info info, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Barrier(MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Barrier_init(MPI_Comm comm, MPI_Info info, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Bcast_init(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                   MPI_Info info, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Exscan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    MPI_Comm comm, MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Gather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Gatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                     MPI_Comm comm, MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                   MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                           MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                            MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                int root, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                              MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm, MPI_Request *request)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int displs[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Neighbor_allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, const int recvcounts[], const int displs[],
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                          int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Info info, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                           MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                           const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                           const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                           const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                           MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallw_init(const void *sendbuf, const int sendcounts[],
                                const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                                const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               int root, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                     MPI_Op op)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_block_init(const void *sendbuf, void *recvbuf, int recvcount,
                                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_init(const void *sendbuf, void *recvbuf, const int recvcounts[],
                            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                            MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm)
             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Scatter_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Scatterv_init(const void *sendbuf, const int sendcounts[], const int displs[],
                      MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) MPICH_API_PUBLIC;
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_free(MPI_Comm *comm) MPICH_API_PUBLIC;
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info *info_used) MPICH_API_PUBLIC;
int MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) MPICH_API_PUBLIC;
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group) MPICH_API_PUBLIC;
int MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Comm_idup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm, MPI_Request *request)
    MPICH_API_PUBLIC;
int MPI_Comm_rank(MPI_Comm comm, int *rank) MPICH_API_PUBLIC;
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) MPICH_API_PUBLIC;
int MPI_Comm_remote_size(MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info) MPICH_API_PUBLIC;
int MPI_Comm_set_name(MPI_Comm comm, const char *comm_name) MPICH_API_PUBLIC;
int MPI_Comm_size(MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int MPI_Comm_test_inter(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                         int remote_leader, int tag, MPI_Comm *newintercomm) MPICH_API_PUBLIC;
int MPI_Intercomm_create_from_groups(MPI_Group local_group, int local_leader,
                                     MPI_Group remote_group, int remote_leader,
                                     const char *stringtag, MPI_Info info,
                                     MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                     MPICH_API_PUBLIC;
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) MPICH_API_PUBLIC;
int MPIX_Comm_test_threadcomm(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int MPIX_Comm_revoke(MPI_Comm comm) MPICH_API_PUBLIC;
int MPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPIX_Comm_failure_ack(MPI_Comm comm) MPICH_API_PUBLIC;
int MPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) MPICH_API_PUBLIC;
int MPIX_Comm_agree(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int MPIX_Comm_get_failed(MPI_Comm comm, MPI_Group *failedgrp) MPICH_API_PUBLIC;
int MPI_Get_address(const void *location, MPI_Aint *address) MPICH_API_PUBLIC;
int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) MPICH_API_PUBLIC;
int MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count) MPICH_API_PUBLIC;
int MPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
             int *position, MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Pack_external(const char *datarep, const void *inbuf, int incount, MPI_Datatype datatype,
                      void *outbuf, MPI_Aint outsize, MPI_Aint *position) MPICH_API_PUBLIC;
int MPI_Pack_external_size(const char *datarep, int incount, MPI_Datatype datatype, MPI_Aint *size)
    MPICH_API_PUBLIC;
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count) MPICH_API_PUBLIC;
int MPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
    MPICH_API_PUBLIC;
int MPI_Type_commit(MPI_Datatype *datatype) MPICH_API_PUBLIC;
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                           const int array_of_distribs[], const int array_of_dargs[],
                           const int array_of_psizes[], int order, MPI_Datatype oldtype,
                           MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                             const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                             MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hindexed_block(int count, int blocklength,
                                   const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                   MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                            MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                           MPICH_API_PUBLIC;
int MPI_Type_create_subarray(int ndims, const int array_of_sizes[], const int array_of_subsizes[],
                             const int array_of_starts[], int order, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_free(MPI_Datatype *datatype) MPICH_API_PUBLIC;
int MPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                          int max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
                          MPI_Datatype array_of_datatypes[]) MPICH_API_PUBLIC;
int MPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                          int *num_datatypes, int *combiner) MPICH_API_PUBLIC;
int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) MPICH_API_PUBLIC;
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
    MPICH_API_PUBLIC;
int MPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) MPICH_API_PUBLIC;
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent)
    MPICH_API_PUBLIC;
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
    MPICH_API_PUBLIC;
int MPI_Type_get_value_index(MPI_Datatype value_type, MPI_Datatype index_type,
                             MPI_Datatype *pair_type) MPICH_API_PUBLIC;
int MPI_Type_indexed(int count, const int array_of_blocklengths[],
                     const int array_of_displacements[], MPI_Datatype oldtype,
                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) MPICH_API_PUBLIC;
int MPI_Type_set_name(MPI_Datatype datatype, const char *type_name) MPICH_API_PUBLIC;
int MPI_Type_size(MPI_Datatype datatype, int *size) MPICH_API_PUBLIC;
int MPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) MPICH_API_PUBLIC;
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                    MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
               MPI_Datatype datatype, MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                        MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                        MPICH_API_PUBLIC;
int MPI_Address(void *location, MPI_Aint *address) MPICH_API_PUBLIC;
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) MPICH_API_PUBLIC;
int MPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) MPICH_API_PUBLIC;
int MPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) MPICH_API_PUBLIC;
int MPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                      MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                    MPI_Datatype array_of_types[], MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPIX_Type_iov_len(MPI_Datatype datatype, MPI_Count max_iov_bytes, MPI_Count *iov_len,
                      MPI_Count *actual_iov_bytes) MPICH_API_PUBLIC;
int MPIX_Type_iov(MPI_Datatype datatype, MPI_Count iov_offset, MPIX_Iov *iov, MPI_Count max_iov_len,
                  MPI_Count *actual_iov_len) MPICH_API_PUBLIC;
int MPI_Add_error_class(int *errorclass) MPICH_API_PUBLIC;
int MPI_Add_error_code(int errorclass, int *errorcode) MPICH_API_PUBLIC;
int MPI_Add_error_string(int errorcode, const char *string) MPICH_API_PUBLIC;
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) MPICH_API_PUBLIC;
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                               MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int MPI_Errhandler_free(MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Error_class(int errorcode, int *errorclass) MPICH_API_PUBLIC;
int MPI_Error_string(int errorcode, char *string, int *resultlen) MPICH_API_PUBLIC;
int MPI_File_call_errhandler(MPI_File fh, int errorcode) MPICH_API_PUBLIC;
int MPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                               MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int MPI_Remove_error_class(int errorclass) MPICH_API_PUBLIC;
int MPI_Remove_error_code(int errorcode) MPICH_API_PUBLIC;
int MPI_Remove_error_string(int errorcode) MPICH_API_PUBLIC;
int MPI_Session_call_errhandler(MPI_Session session, int errorcode) MPICH_API_PUBLIC;
int MPI_Session_create_errhandler(MPI_Session_errhandler_function *session_errhandler_fn,
                                  MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Session_get_errhandler(MPI_Session session, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Session_set_errhandler(MPI_Session session, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int MPI_Win_call_errhandler(MPI_Win win, int errorcode) MPICH_API_PUBLIC;
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                              MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int MPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                          MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
MPI_Fint MPI_Comm_c2f(MPI_Comm comm) MPICH_API_PUBLIC;
MPI_Comm MPI_Comm_f2c(MPI_Fint comm) MPICH_API_PUBLIC;
MPI_Fint MPI_Errhandler_c2f(MPI_Errhandler errhandler) MPICH_API_PUBLIC;
MPI_Errhandler MPI_Errhandler_f2c(MPI_Fint errhandler) MPICH_API_PUBLIC;
MPI_Fint MPI_Group_c2f(MPI_Group group) MPICH_API_PUBLIC;
MPI_Group MPI_Group_f2c(MPI_Fint group) MPICH_API_PUBLIC;
MPI_Fint MPI_Info_c2f(MPI_Info info) MPICH_API_PUBLIC;
MPI_Info MPI_Info_f2c(MPI_Fint info) MPICH_API_PUBLIC;
MPI_Fint MPI_Message_c2f(MPI_Message message) MPICH_API_PUBLIC;
MPI_Message MPI_Message_f2c(MPI_Fint message) MPICH_API_PUBLIC;
MPI_Fint MPI_Op_c2f(MPI_Op op) MPICH_API_PUBLIC;
MPI_Op MPI_Op_f2c(MPI_Fint op) MPICH_API_PUBLIC;
MPI_Fint MPI_Request_c2f(MPI_Request request) MPICH_API_PUBLIC;
MPI_Request MPI_Request_f2c(MPI_Fint request) MPICH_API_PUBLIC;
MPI_Fint MPI_Session_c2f(MPI_Session session) MPICH_API_PUBLIC;
MPI_Session MPI_Session_f2c(MPI_Fint session) MPICH_API_PUBLIC;
MPI_Fint MPI_Type_c2f(MPI_Datatype datatype) MPICH_API_PUBLIC;
MPI_Datatype MPI_Type_f2c(MPI_Fint datatype) MPICH_API_PUBLIC;
MPI_Fint MPI_Win_c2f(MPI_Win win) MPICH_API_PUBLIC;
MPI_Win MPI_Win_f2c(MPI_Fint win) MPICH_API_PUBLIC;
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) MPICH_API_PUBLIC;
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) MPICH_API_PUBLIC;
int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Group_free(MPI_Group *group) MPICH_API_PUBLIC;
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Group_rank(MPI_Group group, int *rank) MPICH_API_PUBLIC;
int MPI_Group_size(MPI_Group group, int *size) MPICH_API_PUBLIC;
int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                              int ranks2[]) MPICH_API_PUBLIC;
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) MPICH_API_PUBLIC;
int MPI_Info_create(MPI_Info *info) MPICH_API_PUBLIC;
int MPI_Info_create_env(int argc, char *argv[], MPI_Info *info) MPICH_API_PUBLIC;
int MPI_Info_delete(MPI_Info info, const char *key) MPICH_API_PUBLIC;
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo) MPICH_API_PUBLIC;
int MPI_Info_free(MPI_Info *info) MPICH_API_PUBLIC;
int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag)
    MPICH_API_PUBLIC;
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys) MPICH_API_PUBLIC;
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key) MPICH_API_PUBLIC;
int MPI_Info_get_string(MPI_Info info, const char *key, int *buflen, char *value, int *flag)
    MPICH_API_PUBLIC;
int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag)
    MPICH_API_PUBLIC;
int MPI_Info_set(MPI_Info info, const char *key, const char *value) MPICH_API_PUBLIC;
int MPIX_Info_set_hex(MPI_Info info, const char *key, const void *value, int value_size)
    MPICH_API_PUBLIC;
int MPI_Abort(MPI_Comm comm, int errorcode) MPICH_API_PUBLIC;
int MPI_Comm_create_from_group(MPI_Group group, const char *stringtag, MPI_Info info,
                               MPI_Errhandler errhandler, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Finalize(void) MPICH_API_PUBLIC;
int MPI_Finalized(int *flag) MPICH_API_PUBLIC;
int MPI_Group_from_session_pset(MPI_Session session, const char *pset_name, MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int MPI_Init(int *argc, char ***argv) MPICH_API_PUBLIC;
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided) MPICH_API_PUBLIC;
int MPI_Initialized(int *flag) MPICH_API_PUBLIC;
int MPI_Is_thread_main(int *flag) MPICH_API_PUBLIC;
int MPI_Query_thread(int *provided) MPICH_API_PUBLIC;
int MPI_Session_finalize(MPI_Session *session) MPICH_API_PUBLIC;
int MPI_Session_get_info(MPI_Session session, MPI_Info *info_used) MPICH_API_PUBLIC;
int MPI_Session_get_nth_pset(MPI_Session session, MPI_Info info, int n, int *pset_len,
                             char *pset_name) MPICH_API_PUBLIC;
int MPI_Session_get_num_psets(MPI_Session session, MPI_Info info, int *npset_names)
    MPICH_API_PUBLIC;
int MPI_Session_get_pset_info(MPI_Session session, const char *pset_name, MPI_Info *info)
    MPICH_API_PUBLIC;
int MPI_Session_init(MPI_Info info, MPI_Errhandler errhandler, MPI_Session *session)
    MPICH_API_PUBLIC;
MPI_Aint MPI_Aint_add(MPI_Aint base, MPI_Aint disp) MPICH_API_PUBLIC;
MPI_Aint MPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) MPICH_API_PUBLIC;
int MPI_Get_library_version(char *version, int *resultlen) MPICH_API_PUBLIC;
int MPI_Get_processor_name(char *name, int *resultlen) MPICH_API_PUBLIC;
int MPI_Get_version(int *version, int *subversion) MPICH_API_PUBLIC;
int MPI_Pcontrol(const int level, ...) MPICH_API_PUBLIC;
int MPIX_GPU_query_support(int gpu_type, int *is_supported) MPICH_API_PUBLIC;
int MPIX_Query_cuda_support(void) MPICH_API_PUBLIC;
int MPIX_Query_ze_support(void) MPICH_API_PUBLIC;
int MPIX_Query_hip_support(void) MPICH_API_PUBLIC;
int MPI_Op_commutative(MPI_Op op, int *commute) MPICH_API_PUBLIC;
int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) MPICH_API_PUBLIC;
int MPI_Op_free(MPI_Op *op) MPICH_API_PUBLIC;
int MPI_Parrived(MPI_Request request, int partition, int *flag) MPICH_API_PUBLIC;
int MPI_Pready(int partition, MPI_Request request) MPICH_API_PUBLIC;
int MPI_Pready_list(int length, const int array_of_partitions[], MPI_Request request)
    MPICH_API_PUBLIC;
int MPI_Pready_range(int partition_low, int partition_high, MPI_Request request) MPICH_API_PUBLIC;
int MPI_Precv_init(void *buf, int partitions, MPI_Count count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_API_PUBLIC;
int MPI_Psend_init(const void *buf, int partitions, MPI_Count count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_API_PUBLIC;
int MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Buffer_attach(void *buffer, int size) MPICH_API_PUBLIC;
int MPI_Buffer_detach(void *buffer_addr, int *size) MPICH_API_PUBLIC;
int MPI_Buffer_flush(void) MPICH_API_PUBLIC;
int MPI_Buffer_iflush(MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Comm_attach_buffer(MPI_Comm comm, void *buffer, int size) MPICH_API_PUBLIC;
int MPI_Comm_detach_buffer(MPI_Comm comm, void *buffer_addr, int *size) MPICH_API_PUBLIC;
int MPI_Comm_flush_buffer(MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Comm_iflush_buffer(MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
              MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Isendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int MPI_Isendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                          int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status)
    MPICH_API_PUBLIC;
int MPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
              MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
             MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                  MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                 MPI_Comm comm, MPI_Status *status)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                         int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Session_attach_buffer(MPI_Session session, void *buffer, int size) MPICH_API_PUBLIC;
int MPI_Session_detach_buffer(MPI_Session session, void *buffer_addr, int *size) MPICH_API_PUBLIC;
int MPI_Session_flush_buffer(MPI_Session session) MPICH_API_PUBLIC;
int MPI_Session_iflush_buffer(MPI_Session session, MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Cancel(MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Grequest_complete(MPI_Request request) MPICH_API_PUBLIC;
int MPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                       MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                       MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Request_free(MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Request_get_status_all(int count, MPI_Request array_of_requests[], int *flag,
                               MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int MPI_Request_get_status_any(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                               MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Request_get_status_some(int incount, MPI_Request array_of_requests[], int *outcount,
                                int array_of_indices[], MPI_Status *array_of_statuses)
                                MPICH_API_PUBLIC;
int MPI_Start(MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Startall(int count, MPI_Request array_of_requests[]) MPICH_API_PUBLIC;
int MPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) MPICH_API_PUBLIC;
int MPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) MPICH_API_PUBLIC;
int MPI_Status_get_error(MPI_Status *status, int *error) MPICH_API_PUBLIC;
int MPI_Status_get_source(MPI_Status *status, int *source) MPICH_API_PUBLIC;
int MPI_Status_get_tag(MPI_Status *status, int *tag) MPICH_API_PUBLIC;
int MPI_Status_set_error(MPI_Status *status, int error) MPICH_API_PUBLIC;
int MPI_Status_set_source(MPI_Status *status, int source) MPICH_API_PUBLIC;
int MPI_Status_set_tag(MPI_Status *status, int tag) MPICH_API_PUBLIC;
int MPI_Status_set_cancelled(MPI_Status *status, int flag) MPICH_API_PUBLIC;
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Test_cancelled(const MPI_Status *status, int *flag) MPICH_API_PUBLIC;
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int MPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int MPI_Wait(MPI_Request *request, MPI_Status *status) MPICH_API_PUBLIC;
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
    MPICH_API_PUBLIC;
int MPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status)
    MPICH_API_PUBLIC;
int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                   int target_rank, MPI_Aint target_disp, int target_count,
                   MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) MPICH_API_PUBLIC;
int MPI_Compare_and_swap(const void *origin_addr, const void *compare_addr, void *result_addr,
                         MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Win win)
                         MPICH_API_PUBLIC;
int MPI_Fetch_and_op(const void *origin_addr, void *result_addr, MPI_Datatype datatype,
                     int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win)
                     MPICH_API_PUBLIC;
int MPI_Free_mem(void *base) MPICH_API_PUBLIC;
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
            MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win)
            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                       void *result_addr, int result_count, MPI_Datatype result_datatype,
                       int target_rank, MPI_Aint target_disp, int target_count,
                       MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
            MPI_Win win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Raccumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win,
             MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                        void *result_addr, int result_count, MPI_Datatype result_datatype,
                        int target_rank, MPI_Aint target_disp, int target_count,
                        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request)
             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                     MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                            void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) MPICH_API_PUBLIC;
int MPI_Win_complete(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                   MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_detach(MPI_Win win, const void *base) MPICH_API_PUBLIC;
int MPI_Win_fence(int assert, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_flush(int rank, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_flush_all(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_flush_local(int rank, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_flush_local_all(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_free(MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_get_group(MPI_Win win, MPI_Group *group) MPICH_API_PUBLIC;
int MPI_Win_get_info(MPI_Win win, MPI_Info *info_used) MPICH_API_PUBLIC;
int MPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) MPICH_API_PUBLIC;
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_lock_all(int assert, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_post(MPI_Group group, int assert, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_set_info(MPI_Win win, MPI_Info info) MPICH_API_PUBLIC;
int MPI_Win_set_name(MPI_Win win, const char *win_name) MPICH_API_PUBLIC;
int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr)
    MPICH_API_PUBLIC;
int MPI_Win_start(MPI_Group group, int assert, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_sync(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_test(MPI_Win win, int *flag) MPICH_API_PUBLIC;
int MPI_Win_unlock(int rank, MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_unlock_all(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Win_wait(MPI_Win win) MPICH_API_PUBLIC;
int MPI_Close_port(const char *port_name) MPICH_API_PUBLIC;
int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                    MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Comm_disconnect(MPI_Comm *comm) MPICH_API_PUBLIC;
int MPI_Comm_get_parent(MPI_Comm *parent) MPICH_API_PUBLIC;
int MPI_Comm_join(int fd, MPI_Comm *intercomm) MPICH_API_PUBLIC;
int MPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                   MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) MPICH_API_PUBLIC;
int MPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                            const int array_of_maxprocs[], const MPI_Info array_of_info[], int root,
                            MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                            MPICH_API_PUBLIC;
int MPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) MPICH_API_PUBLIC;
int MPI_Open_port(MPI_Info info, char *port_name) MPICH_API_PUBLIC;
int MPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
    MPICH_API_PUBLIC;
int MPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
    MPICH_API_PUBLIC;
int MPIX_Stream_create(MPI_Info info, MPIX_Stream *stream) MPICH_API_PUBLIC;
int MPIX_Stream_free(MPIX_Stream *stream) MPICH_API_PUBLIC;
int MPIX_Stream_comm_create(MPI_Comm comm, MPIX_Stream stream, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPIX_Stream_comm_create_multiplex(MPI_Comm comm, int count, MPIX_Stream array_of_streams[],
                                      MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPIX_Comm_get_stream(MPI_Comm comm, int idx, MPIX_Stream *stream) MPICH_API_PUBLIC;
int MPIX_Stream_progress(MPIX_Stream stream) MPICH_API_PUBLIC;
int MPIX_Start_progress_thread(MPIX_Stream stream) MPICH_API_PUBLIC;
int MPIX_Stop_progress_thread(MPIX_Stream stream) MPICH_API_PUBLIC;
int MPIX_Stream_send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, int source_stream_index, int dest_stream_index)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                     MPI_Comm comm, int source_stream_index, int dest_stream_index,
                     MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Send_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Recv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, MPI_Status *status)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Isend_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Irecv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Wait_enqueue(MPI_Request *request, MPI_Status *status) MPICH_API_PUBLIC;
int MPIX_Waitall_enqueue(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
    MPICH_API_PUBLIC;
int MPIX_Allreduce_enqueue(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                           MPI_Op op, MPI_Comm comm)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPIX_Threadcomm_init(MPI_Comm comm, int num_threads, MPI_Comm *newthreadcomm) MPICH_API_PUBLIC;
int MPIX_Threadcomm_free(MPI_Comm *threadcomm) MPICH_API_PUBLIC;
int MPIX_Threadcomm_start(MPI_Comm threadcomm) MPICH_API_PUBLIC;
int MPIX_Threadcomm_finish(MPI_Comm threadcomm) MPICH_API_PUBLIC;
double MPI_Wtick(void) MPICH_API_PUBLIC;
double MPI_Wtime(void) MPICH_API_PUBLIC;
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) MPICH_API_PUBLIC;
int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                    int reorder, MPI_Comm *comm_cart) MPICH_API_PUBLIC;
int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[])
    MPICH_API_PUBLIC;
int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank)
    MPICH_API_PUBLIC;
int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) MPICH_API_PUBLIC;
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)
    MPICH_API_PUBLIC;
int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) MPICH_API_PUBLIC;
int MPI_Cartdim_get(MPI_Comm comm, int *ndims) MPICH_API_PUBLIC;
int MPI_Dims_create(int nnodes, int ndims, int dims[]) MPICH_API_PUBLIC;
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                          const int destinations[], const int weights[], MPI_Info info, int reorder,
                          MPI_Comm *comm_dist_graph) MPICH_API_PUBLIC;
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                   const int sourceweights[], int outdegree,
                                   const int destinations[], const int destweights[], MPI_Info info,
                                   int reorder, MPI_Comm *comm_dist_graph) MPICH_API_PUBLIC;
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                             int maxoutdegree, int destinations[], int destweights[])
                             MPICH_API_PUBLIC;
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted)
    MPICH_API_PUBLIC;
int MPI_Get_hw_resource_info(MPI_Info *hw_info) MPICH_API_PUBLIC;
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                     int reorder, MPI_Comm *comm_graph) MPICH_API_PUBLIC;
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[])
    MPICH_API_PUBLIC;
int MPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank)
    MPICH_API_PUBLIC;
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[])
    MPICH_API_PUBLIC;
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) MPICH_API_PUBLIC;
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) MPICH_API_PUBLIC;
int MPI_Topo_test(MPI_Comm comm, int *status) MPICH_API_PUBLIC;
/* End Prototypes */

int MPI_T_category_changed(int *update_number) MPICH_API_PUBLIC;
int MPI_T_category_get_categories(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int MPI_T_category_get_cvars(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int MPI_T_category_get_events(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int MPI_T_category_get_index(const char *name, int *cat_index) MPICH_API_PUBLIC;
int MPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                            int *num_cvars, int *num_pvars, int *num_categories) MPICH_API_PUBLIC;
int MPI_T_category_get_num(int *num_cat) MPICH_API_PUBLIC;
int MPI_T_category_get_num_events(int cat_index, int *num_events) MPICH_API_PUBLIC;
int MPI_T_category_get_pvars(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int MPI_T_cvar_get_index(const char *name, int *cvar_index) MPICH_API_PUBLIC;
int MPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *bind, int *scope) MPICH_API_PUBLIC;
int MPI_T_cvar_get_num(int *num_cvar) MPICH_API_PUBLIC;
int MPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                            int *count) MPICH_API_PUBLIC;
int MPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) MPICH_API_PUBLIC;
int MPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) MPICH_API_PUBLIC;
int MPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) MPICH_API_PUBLIC;
int MPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len) MPICH_API_PUBLIC;
int MPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len)
    MPICH_API_PUBLIC;
int MPI_T_event_callback_get_info(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info *info_used) MPICH_API_PUBLIC;
int MPI_T_event_callback_set_info(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info info) MPICH_API_PUBLIC;
int MPI_T_event_copy(MPI_T_event_instance event_instance, void *buffer) MPICH_API_PUBLIC;
int MPI_T_event_get_index(const char *name, int *event_index) MPICH_API_PUBLIC;
int MPI_T_event_get_info(int event_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype array_of_datatypes[], MPI_Aint array_of_displacements[],
                         int *num_elements, MPI_T_enum *enumtype, MPI_Info *info, char *desc,
                         int *desc_len, int *bind) MPICH_API_PUBLIC;
int MPI_T_event_get_num(int *num_events) MPICH_API_PUBLIC;
int MPI_T_event_get_source(MPI_T_event_instance event_instance, int *source_index)
    MPICH_API_PUBLIC;
int MPI_T_event_get_timestamp(MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
    MPICH_API_PUBLIC;
int MPI_T_event_handle_alloc(int event_index, void *obj_handle, MPI_Info info,
                             MPI_T_event_registration *event_registration) MPICH_API_PUBLIC;
int MPI_T_event_handle_free(MPI_T_event_registration event_registration, void *user_data,
                            MPI_T_event_free_cb_function free_cb_function) MPICH_API_PUBLIC;
int MPI_T_event_handle_get_info(MPI_T_event_registration event_registration, MPI_Info *info_used)
    MPICH_API_PUBLIC;
int MPI_T_event_handle_set_info(MPI_T_event_registration event_registration, MPI_Info info)
    MPICH_API_PUBLIC;
int MPI_T_event_read(MPI_T_event_instance event_instance, int element_index, void *buffer)
    MPICH_API_PUBLIC;
int MPI_T_event_register_callback(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                  MPI_T_event_cb_function event_cb_function) MPICH_API_PUBLIC;
int MPI_T_event_set_dropped_handler(MPI_T_event_registration event_registration,
                                    MPI_T_event_dropped_cb_function dropped_cb_function)
                                    MPICH_API_PUBLIC;
int MPI_T_finalize(void) MPICH_API_PUBLIC;
int MPI_T_init_thread(int required, int *provided) MPICH_API_PUBLIC;
int MPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) MPICH_API_PUBLIC;
int MPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *bind, int *readonly, int *continuous, int *atomic) MPICH_API_PUBLIC;
int MPI_T_pvar_get_num(int *num_pvar) MPICH_API_PUBLIC;
int MPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                            MPI_T_pvar_handle *handle, int *count) MPICH_API_PUBLIC;
int MPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle) MPICH_API_PUBLIC;
int MPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
    MPICH_API_PUBLIC;
int MPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
    MPICH_API_PUBLIC;
int MPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int MPI_T_pvar_session_create(MPI_T_pvar_session *session) MPICH_API_PUBLIC;
int MPI_T_pvar_session_free(MPI_T_pvar_session *session) MPICH_API_PUBLIC;
int MPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int MPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int MPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf)
    MPICH_API_PUBLIC;
int MPI_T_source_get_info(int source_index, char *name, int *name_len, char *desc, int *desc_len,
                          MPI_T_source_order *ordering, MPI_Count *ticks_per_second,
                          MPI_Count *max_ticks, MPI_Info *info) MPICH_API_PUBLIC;
int MPI_T_source_get_num(int *num_sources) MPICH_API_PUBLIC;
int MPI_T_source_get_timestamp(int source_index, MPI_Count *timestamp) MPICH_API_PUBLIC;
int MPIX_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn,
                        MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                        void *extra_state, MPI_Request *request) MPICH_API_PUBLIC;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                               MPI_Grequest_free_function *free_fn,
                               MPI_Grequest_cancel_function *cancel_fn,
                               MPIX_Grequest_poll_function *poll_fn,
                               MPIX_Grequest_wait_function *wait_fn,
                               MPIX_Grequest_class *greq_class) MPICH_API_PUBLIC;
int MPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                 MPI_Request *request) MPICH_API_PUBLIC;

int MPI_Allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                         void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                         MPI_Info info, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                     MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                          void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                          MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Allreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Allreduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                        MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                         const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                         const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                         MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                    MPICH_API_PUBLIC;
int MPI_Alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                         const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                         const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                         const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                         MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Bcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Bcast_init_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Exscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Exscan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Gather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                 MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Gather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                      MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Gatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Gatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                       MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Iallgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Iallgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Iallreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ialltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ialltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Ialltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                     MPI_Request *request) MPICH_API_PUBLIC;
int MPI_Ibcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                 MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Iexscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Igather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Igatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Ineighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm, MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ineighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Ineighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request)
                              MPICH_API_PUBLIC;
int MPI_Ireduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ireduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Ireduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                                MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Iscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Iscatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Iscatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Neighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                              MPI_Datatype recvtype, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Neighbor_allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                   void *recvbuf, const MPI_Count recvcounts[],
                                   const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                   MPI_Info info, MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                            void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                 MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                             const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                             const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                  const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                                  const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                  MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                             const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                             const MPI_Datatype recvtypes[], MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Neighbor_alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                  const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                  void *recvbuf, const MPI_Count recvcounts[],
                                  const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                  MPICH_API_PUBLIC;
int MPI_Reduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_local_c(const void *inbuf, void *inoutbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_block_init_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                    MPI_Request *request)
                                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Reduce_scatter_init_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                              MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int MPI_Scatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Scatter_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Scatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Scatterv_init_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                        MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int MPI_Get_count_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int MPI_Get_elements_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int MPI_Pack_c(const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
               MPI_Count outsize, MPI_Count *position, MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Pack_external_c(const char *datarep, const void *inbuf, MPI_Count incount,
                        MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                        MPI_Count *position) MPICH_API_PUBLIC;
int MPI_Pack_external_size_c(const char *datarep, MPI_Count incount, MPI_Datatype datatype,
                             MPI_Count *size) MPICH_API_PUBLIC;
int MPI_Pack_size_c(MPI_Count incount, MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size)
    MPICH_API_PUBLIC;
int MPI_Status_set_elements_c(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
    MPICH_API_PUBLIC;
int MPI_Type_contiguous_c(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype *newtype)
    MPICH_API_PUBLIC;
int MPI_Type_create_darray_c(int size, int rank, int ndims, const MPI_Count array_of_gsizes[],
                             const int array_of_distribs[], const int array_of_dargs[],
                             const int array_of_psizes[], int order, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hindexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                               const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                               MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hindexed_block_c(MPI_Count count, MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_hvector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                              MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_indexed_block_c(MPI_Count count, MPI_Count blocklength,
                                    const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_resized_c(MPI_Datatype oldtype, MPI_Count lb, MPI_Count extent,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_create_struct_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                             const MPI_Count array_of_displacements[],
                             const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                             MPICH_API_PUBLIC;
int MPI_Type_create_subarray_c(int ndims, const MPI_Count array_of_sizes[],
                               const MPI_Count array_of_subsizes[],
                               const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                               MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_get_contents_c(MPI_Datatype datatype, MPI_Count max_integers, MPI_Count max_addresses,
                            MPI_Count max_large_counts, MPI_Count max_datatypes,
                            int array_of_integers[], MPI_Aint array_of_addresses[],
                            MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                            MPICH_API_PUBLIC;
int MPI_Type_get_envelope_c(MPI_Datatype datatype, MPI_Count *num_integers,
                            MPI_Count *num_addresses, MPI_Count *num_large_counts,
                            MPI_Count *num_datatypes, int *combiner) MPICH_API_PUBLIC;
int MPI_Type_get_extent_c(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
    MPICH_API_PUBLIC;
int MPI_Type_get_true_extent_c(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
    MPICH_API_PUBLIC;
int MPI_Type_indexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                       const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                       MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Type_size_c(MPI_Datatype datatype, MPI_Count *size) MPICH_API_PUBLIC;
int MPI_Type_vector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                      MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int MPI_Unpack_c(const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                 MPI_Count outcount, MPI_Datatype datatype, MPI_Comm comm) MPICH_API_PUBLIC;
int MPI_Unpack_external_c(const char datarep[], const void *inbuf, MPI_Count insize,
                          MPI_Count *position, void *outbuf, MPI_Count outcount,
                          MPI_Datatype datatype) MPICH_API_PUBLIC;
int MPI_Op_create_c(MPI_User_function_c *user_fn, int commute, MPI_Op *op) MPICH_API_PUBLIC;
int MPI_Bsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Bsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Buffer_attach_c(void *buffer, MPI_Count size) MPICH_API_PUBLIC;
int MPI_Buffer_detach_c(void *buffer_addr, MPI_Count *size) MPICH_API_PUBLIC;
int MPI_Comm_attach_buffer_c(MPI_Comm comm, void *buffer, MPI_Count size) MPICH_API_PUBLIC;
int MPI_Comm_detach_buffer_c(MPI_Comm comm, void *buffer_addr, MPI_Count *size) MPICH_API_PUBLIC;
int MPI_Ibsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Imrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                 MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Irsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Isendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                    int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                    int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int MPI_Isendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                            int sendtag, int source, int recvtag, MPI_Comm comm,
                            MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Issend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Mrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
               MPI_Comm comm, MPI_Status *status)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Recv_init_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Send_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Sendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                   int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int MPI_Sendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int sendtag,
                           int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Session_attach_buffer_c(MPI_Session session, void *buffer, MPI_Count size)
    MPICH_API_PUBLIC;
int MPI_Session_detach_buffer_c(MPI_Session session, void *buffer_addr, MPI_Count *size)
    MPICH_API_PUBLIC;
int MPI_Ssend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Ssend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Accumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                     int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                     MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Get_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, MPI_Count target_count,
              MPI_Datatype target_datatype, MPI_Win win)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Get_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                         MPI_Win win)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Put_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, MPI_Count target_count,
              MPI_Datatype target_datatype, MPI_Win win)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Raccumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                      MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rget_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Rget_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                          MPI_Win win, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int MPI_Rput_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPI_Win_allocate_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                       void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_allocate_shared_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                              void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_create_c(void *base, MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                     MPI_Win *win) MPICH_API_PUBLIC;
int MPI_Win_shared_query_c(MPI_Win win, int rank, MPI_Aint *size, MPI_Aint *disp_unit,
                           void *baseptr) MPICH_API_PUBLIC;
int MPIX_Stream_send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Stream_irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Send_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Recv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Status *status)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Isend_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Irecv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int MPIX_Allreduce_enqueue_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;

int PMPI_Status_f082c(const MPI_F08_status *f08_status, MPI_Status *c_status) MPICH_API_PUBLIC;
int PMPI_Status_c2f08(const MPI_Status *c_status, MPI_F08_status *f08_status) MPICH_API_PUBLIC;
int PMPI_Status_f082f(const MPI_F08_status *f08_status, MPI_Fint *f_status) MPICH_API_PUBLIC;
int PMPI_Status_f2f08(const MPI_Fint *f_status, MPI_F08_status *f08_status) MPICH_API_PUBLIC;
int PMPI_Type_create_f90_integer(int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Attr_delete(MPI_Comm comm, int keyval) MPICH_API_PUBLIC;
int PMPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) MPICH_API_PUBLIC;
int PMPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) MPICH_API_PUBLIC;
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state) MPICH_API_PUBLIC;
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) MPICH_API_PUBLIC;
int PMPI_Comm_free_keyval(int *comm_keyval) MPICH_API_PUBLIC;
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) MPICH_API_PUBLIC;
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                       void *extra_state) MPICH_API_PUBLIC;
int PMPI_Keyval_free(int *keyval) MPICH_API_PUBLIC;
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                            void *extra_state) MPICH_API_PUBLIC;
int PMPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) MPICH_API_PUBLIC;
int PMPI_Type_free_keyval(int *type_keyval) MPICH_API_PUBLIC;
int PMPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val)
    MPICH_API_PUBLIC;
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state) MPICH_API_PUBLIC;
int PMPI_Win_delete_attr(MPI_Win win, int win_keyval) MPICH_API_PUBLIC;
int PMPI_Win_free_keyval(int *win_keyval) MPICH_API_PUBLIC;
int PMPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) MPICH_API_PUBLIC;
int PMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                        int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                          void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                      MPI_Datatype recvtype, MPI_Comm comm)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                         const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                         MPI_Comm comm, MPI_Info info, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                           void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                           MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                           MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Allreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Allreduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                        MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Allreduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                       int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                         void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                         MPI_Info info, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                        MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                        const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                          const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                   MPICH_API_PUBLIC;
int PMPI_Alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                     MPICH_API_PUBLIC;
int PMPI_Alltoallw_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                        const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                        const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                        MPI_Info info, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                          const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                          MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Barrier(MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Barrier_init(MPI_Comm comm, MPI_Info info, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bcast_init(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                    MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bcast_init_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                      MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Exscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Exscan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Exscan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Gather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Gather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Gather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Gatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Gatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                      MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Gatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Iallgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                      MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Iallgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                       MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Iallreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ialltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Ialltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                      MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                    const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                    MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ialltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                      MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ibarrier(MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Ibcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                  MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Iexscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Igather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Igatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                    int root, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ineighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                               MPI_Comm comm, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const int recvcounts[], const int displs[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Ineighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                void *recvbuf, const MPI_Count recvcounts[],
                                const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm, MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                               const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                               MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                             MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ineighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                               const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                               void *recvbuf, const MPI_Count recvcounts[],
                               const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                               MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                 int root, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ireduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                   MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ireduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                           MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                               MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Ireduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                                 MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Iscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Iscatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Iscatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                 MPI_Info info, MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                   void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, const int recvcounts[], const int displs[],
                                  MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Neighbor_allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                    void *recvbuf, const MPI_Count recvcounts[],
                                    const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                    MPI_Info info, MPI_Request *request)
                                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              MPI_Datatype recvtype, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                 const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                                 MPI_Info info, MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                   const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallw_init(const void *sendbuf, const int sendcounts[],
                                 const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                 void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                                 const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Neighbor_alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                   const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                   void *recvbuf, const MPI_Count recvcounts[],
                                   const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                   MPICH_API_PUBLIC;
int PMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                int root, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                     MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                      MPI_Op op)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_local_c(const void *inbuf, void *inoutbuf, MPI_Count count, MPI_Datatype datatype,
                        MPI_Op op)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_block_init(const void *sendbuf, void *recvbuf, int recvcount,
                                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_block_init_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                     MPI_Request *request)
                                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_init(const void *sendbuf, void *recvbuf, const int recvcounts[],
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                             MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Reduce_scatter_init_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                               MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Scan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Scan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Scan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Scatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Scatter_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Scatter_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Scatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Scatterv_init(const void *sendbuf, const int sendcounts[], const int displs[],
                       MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Scatterv_init_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                         MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                         MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7) MPICH_API_PUBLIC;
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) MPICH_API_PUBLIC;
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_free(MPI_Comm *comm) MPICH_API_PUBLIC;
int PMPI_Comm_get_info(MPI_Comm comm, MPI_Info *info_used) MPICH_API_PUBLIC;
int PMPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) MPICH_API_PUBLIC;
int PMPI_Comm_group(MPI_Comm comm, MPI_Group *group) MPICH_API_PUBLIC;
int PMPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Comm_idup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm, MPI_Request *request)
    MPICH_API_PUBLIC;
int PMPI_Comm_rank(MPI_Comm comm, int *rank) MPICH_API_PUBLIC;
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) MPICH_API_PUBLIC;
int PMPI_Comm_remote_size(MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int PMPI_Comm_set_info(MPI_Comm comm, MPI_Info info) MPICH_API_PUBLIC;
int PMPI_Comm_set_name(MPI_Comm comm, const char *comm_name) MPICH_API_PUBLIC;
int PMPI_Comm_size(MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int PMPI_Comm_test_inter(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                          int remote_leader, int tag, MPI_Comm *newintercomm) MPICH_API_PUBLIC;
int PMPI_Intercomm_create_from_groups(MPI_Group local_group, int local_leader,
                                      MPI_Group remote_group, int remote_leader,
                                      const char *stringtag, MPI_Info info,
                                      MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                      MPICH_API_PUBLIC;
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) MPICH_API_PUBLIC;
int PMPIX_Comm_test_threadcomm(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int PMPIX_Comm_revoke(MPI_Comm comm) MPICH_API_PUBLIC;
int PMPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPIX_Comm_failure_ack(MPI_Comm comm) MPICH_API_PUBLIC;
int PMPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) MPICH_API_PUBLIC;
int PMPIX_Comm_agree(MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int PMPIX_Comm_get_failed(MPI_Comm comm, MPI_Group *failedgrp) MPICH_API_PUBLIC;
int PMPI_Get_address(const void *location, MPI_Aint *address) MPICH_API_PUBLIC;
int PMPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) MPICH_API_PUBLIC;
int PMPI_Get_count_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int PMPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count)
    MPICH_API_PUBLIC;
int PMPI_Get_elements_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int PMPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
    MPICH_API_PUBLIC;
int PMPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
              int *position, MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Pack_c(const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
                MPI_Count outsize, MPI_Count *position, MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Pack_external(const char *datarep, const void *inbuf, int incount, MPI_Datatype datatype,
                       void *outbuf, MPI_Aint outsize, MPI_Aint *position) MPICH_API_PUBLIC;
int PMPI_Pack_external_c(const char *datarep, const void *inbuf, MPI_Count incount,
                         MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                         MPI_Count *position) MPICH_API_PUBLIC;
int PMPI_Pack_external_size(const char *datarep, int incount, MPI_Datatype datatype,
                            MPI_Aint *size) MPICH_API_PUBLIC;
int PMPI_Pack_external_size_c(const char *datarep, MPI_Count incount, MPI_Datatype datatype,
                              MPI_Count *size) MPICH_API_PUBLIC;
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int PMPI_Pack_size_c(MPI_Count incount, MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size)
    MPICH_API_PUBLIC;
int PMPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count)
    MPICH_API_PUBLIC;
int PMPI_Status_set_elements_c(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
    MPICH_API_PUBLIC;
int PMPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
    MPICH_API_PUBLIC;
int PMPI_Type_commit(MPI_Datatype *datatype) MPICH_API_PUBLIC;
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_contiguous_c(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype *newtype)
    MPICH_API_PUBLIC;
int PMPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                            const int array_of_distribs[], const int array_of_dargs[],
                            const int array_of_psizes[], int order, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_darray_c(int size, int rank, int ndims, const MPI_Count array_of_gsizes[],
                              const int array_of_distribs[], const int array_of_dargs[],
                              const int array_of_psizes[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_hindexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                                const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_hindexed_block(int count, int blocklength,
                                    const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_hindexed_block_c(MPI_Count count, MPI_Count blocklength,
                                      const MPI_Count array_of_displacements[],
                                      MPI_Datatype oldtype, MPI_Datatype *newtype)
                                      MPICH_API_PUBLIC;
int PMPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_hvector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                               MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_indexed_block_c(MPI_Count count, MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                             MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_resized_c(MPI_Datatype oldtype, MPI_Count lb, MPI_Count extent,
                               MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_struct(int count, const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                            MPICH_API_PUBLIC;
int PMPI_Type_create_struct_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                              const MPI_Count array_of_displacements[],
                              const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                              MPICH_API_PUBLIC;
int PMPI_Type_create_subarray(int ndims, const int array_of_sizes[], const int array_of_subsizes[],
                              const int array_of_starts[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_create_subarray_c(int ndims, const MPI_Count array_of_sizes[],
                                const MPI_Count array_of_subsizes[],
                                const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                                MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_free(MPI_Datatype *datatype) MPICH_API_PUBLIC;
int PMPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                           int max_datatypes, int array_of_integers[],
                           MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[])
                           MPICH_API_PUBLIC;
int PMPI_Type_get_contents_c(MPI_Datatype datatype, MPI_Count max_integers, MPI_Count max_addresses,
                             MPI_Count max_large_counts, MPI_Count max_datatypes,
                             int array_of_integers[], MPI_Aint array_of_addresses[],
                             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                             MPICH_API_PUBLIC;
int PMPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                           int *num_datatypes, int *combiner) MPICH_API_PUBLIC;
int PMPI_Type_get_envelope_c(MPI_Datatype datatype, MPI_Count *num_integers,
                             MPI_Count *num_addresses, MPI_Count *num_large_counts,
                             MPI_Count *num_datatypes, int *combiner) MPICH_API_PUBLIC;
int PMPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) MPICH_API_PUBLIC;
int PMPI_Type_get_extent_c(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
    MPICH_API_PUBLIC;
int PMPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
    MPICH_API_PUBLIC;
int PMPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) MPICH_API_PUBLIC;
int PMPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent)
    MPICH_API_PUBLIC;
int PMPI_Type_get_true_extent_c(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
    MPICH_API_PUBLIC;
int PMPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
    MPICH_API_PUBLIC;
int PMPI_Type_get_value_index(MPI_Datatype value_type, MPI_Datatype index_type,
                              MPI_Datatype *pair_type) MPICH_API_PUBLIC;
int PMPI_Type_indexed(int count, const int array_of_blocklengths[],
                      const int array_of_displacements[], MPI_Datatype oldtype,
                      MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_indexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                        const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                        MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) MPICH_API_PUBLIC;
int PMPI_Type_set_name(MPI_Datatype datatype, const char *type_name) MPICH_API_PUBLIC;
int PMPI_Type_size(MPI_Datatype datatype, int *size) MPICH_API_PUBLIC;
int PMPI_Type_size_c(MPI_Datatype datatype, MPI_Count *size) MPICH_API_PUBLIC;
int PMPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) MPICH_API_PUBLIC;
int PMPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_vector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                       MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
                MPI_Datatype datatype, MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Unpack_c(const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                  MPI_Count outcount, MPI_Datatype datatype, MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                         MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                         MPICH_API_PUBLIC;
int PMPI_Unpack_external_c(const char datarep[], const void *inbuf, MPI_Count insize,
                           MPI_Count *position, void *outbuf, MPI_Count outcount,
                           MPI_Datatype datatype) MPICH_API_PUBLIC;
int PMPI_Address(void *location, MPI_Aint *address) MPICH_API_PUBLIC;
int PMPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) MPICH_API_PUBLIC;
int PMPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) MPICH_API_PUBLIC;
int PMPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) MPICH_API_PUBLIC;
int PMPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                       MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                      MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                     MPI_Datatype array_of_types[], MPI_Datatype *newtype) MPICH_API_PUBLIC;
int PMPIX_Type_iov_len(MPI_Datatype datatype, MPI_Count max_iov_bytes, MPI_Count *iov_len,
                       MPI_Count *actual_iov_bytes) MPICH_API_PUBLIC;
int PMPIX_Type_iov(MPI_Datatype datatype, MPI_Count iov_offset, MPIX_Iov *iov,
                   MPI_Count max_iov_len, MPI_Count *actual_iov_len) MPICH_API_PUBLIC;
int PMPI_Add_error_class(int *errorclass) MPICH_API_PUBLIC;
int PMPI_Add_error_code(int errorclass, int *errorcode) MPICH_API_PUBLIC;
int PMPI_Add_error_string(int errorcode, const char *string) MPICH_API_PUBLIC;
int PMPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) MPICH_API_PUBLIC;
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int PMPI_Errhandler_free(MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Error_class(int errorcode, int *errorclass) MPICH_API_PUBLIC;
int PMPI_Error_string(int errorcode, char *string, int *resultlen) MPICH_API_PUBLIC;
int PMPI_File_call_errhandler(MPI_File fh, int errorcode) MPICH_API_PUBLIC;
int PMPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int PMPI_Remove_error_class(int errorclass) MPICH_API_PUBLIC;
int PMPI_Remove_error_code(int errorcode) MPICH_API_PUBLIC;
int PMPI_Remove_error_string(int errorcode) MPICH_API_PUBLIC;
int PMPI_Session_call_errhandler(MPI_Session session, int errorcode) MPICH_API_PUBLIC;
int PMPI_Session_create_errhandler(MPI_Session_errhandler_function *session_errhandler_fn,
                                   MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Session_get_errhandler(MPI_Session session, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Session_set_errhandler(MPI_Session session, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int PMPI_Win_call_errhandler(MPI_Win win, int errorcode) MPICH_API_PUBLIC;
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int PMPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                           MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) MPICH_API_PUBLIC;
MPI_Fint PMPI_Comm_c2f(MPI_Comm comm) MPICH_API_PUBLIC;
MPI_Comm PMPI_Comm_f2c(MPI_Fint comm) MPICH_API_PUBLIC;
MPI_Fint PMPI_Errhandler_c2f(MPI_Errhandler errhandler) MPICH_API_PUBLIC;
MPI_Errhandler PMPI_Errhandler_f2c(MPI_Fint errhandler) MPICH_API_PUBLIC;
MPI_Fint PMPI_Group_c2f(MPI_Group group) MPICH_API_PUBLIC;
MPI_Group PMPI_Group_f2c(MPI_Fint group) MPICH_API_PUBLIC;
MPI_Fint PMPI_Info_c2f(MPI_Info info) MPICH_API_PUBLIC;
MPI_Info PMPI_Info_f2c(MPI_Fint info) MPICH_API_PUBLIC;
MPI_Fint PMPI_Message_c2f(MPI_Message message) MPICH_API_PUBLIC;
MPI_Message PMPI_Message_f2c(MPI_Fint message) MPICH_API_PUBLIC;
MPI_Fint PMPI_Op_c2f(MPI_Op op) MPICH_API_PUBLIC;
MPI_Op PMPI_Op_f2c(MPI_Fint op) MPICH_API_PUBLIC;
MPI_Fint PMPI_Request_c2f(MPI_Request request) MPICH_API_PUBLIC;
MPI_Request PMPI_Request_f2c(MPI_Fint request) MPICH_API_PUBLIC;
MPI_Fint PMPI_Session_c2f(MPI_Session session) MPICH_API_PUBLIC;
MPI_Session PMPI_Session_f2c(MPI_Fint session) MPICH_API_PUBLIC;
MPI_Fint PMPI_Type_c2f(MPI_Datatype datatype) MPICH_API_PUBLIC;
MPI_Datatype PMPI_Type_f2c(MPI_Fint datatype) MPICH_API_PUBLIC;
MPI_Fint PMPI_Win_c2f(MPI_Win win) MPICH_API_PUBLIC;
MPI_Win PMPI_Win_f2c(MPI_Fint win) MPICH_API_PUBLIC;
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) MPICH_API_PUBLIC;
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_free(MPI_Group *group) MPICH_API_PUBLIC;
int PMPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Group_rank(MPI_Group group, int *rank) MPICH_API_PUBLIC;
int PMPI_Group_size(MPI_Group group, int *size) MPICH_API_PUBLIC;
int PMPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                               int ranks2[]) MPICH_API_PUBLIC;
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) MPICH_API_PUBLIC;
int PMPI_Info_create(MPI_Info *info) MPICH_API_PUBLIC;
int PMPI_Info_create_env(int argc, char *argv[], MPI_Info *info) MPICH_API_PUBLIC;
int PMPI_Info_delete(MPI_Info info, const char *key) MPICH_API_PUBLIC;
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo) MPICH_API_PUBLIC;
int PMPI_Info_free(MPI_Info *info) MPICH_API_PUBLIC;
int PMPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Info_get_nkeys(MPI_Info info, int *nkeys) MPICH_API_PUBLIC;
int PMPI_Info_get_nthkey(MPI_Info info, int n, char *key) MPICH_API_PUBLIC;
int PMPI_Info_get_string(MPI_Info info, const char *key, int *buflen, char *value, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag)
    MPICH_API_PUBLIC;
int PMPI_Info_set(MPI_Info info, const char *key, const char *value) MPICH_API_PUBLIC;
int PMPIX_Info_set_hex(MPI_Info info, const char *key, const void *value, int value_size)
    MPICH_API_PUBLIC;
int PMPI_Abort(MPI_Comm comm, int errorcode) MPICH_API_PUBLIC;
int PMPI_Comm_create_from_group(MPI_Group group, const char *stringtag, MPI_Info info,
                                MPI_Errhandler errhandler, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Finalize(void) MPICH_API_PUBLIC;
int PMPI_Finalized(int *flag) MPICH_API_PUBLIC;
int PMPI_Group_from_session_pset(MPI_Session session, const char *pset_name, MPI_Group *newgroup)
    MPICH_API_PUBLIC;
int PMPI_Init(int *argc, char ***argv) MPICH_API_PUBLIC;
int PMPI_Init_thread(int *argc, char ***argv, int required, int *provided) MPICH_API_PUBLIC;
int PMPI_Initialized(int *flag) MPICH_API_PUBLIC;
int PMPI_Is_thread_main(int *flag) MPICH_API_PUBLIC;
int PMPI_Query_thread(int *provided) MPICH_API_PUBLIC;
int PMPI_Session_finalize(MPI_Session *session) MPICH_API_PUBLIC;
int PMPI_Session_get_info(MPI_Session session, MPI_Info *info_used) MPICH_API_PUBLIC;
int PMPI_Session_get_nth_pset(MPI_Session session, MPI_Info info, int n, int *pset_len,
                              char *pset_name) MPICH_API_PUBLIC;
int PMPI_Session_get_num_psets(MPI_Session session, MPI_Info info, int *npset_names)
    MPICH_API_PUBLIC;
int PMPI_Session_get_pset_info(MPI_Session session, const char *pset_name, MPI_Info *info)
    MPICH_API_PUBLIC;
int PMPI_Session_init(MPI_Info info, MPI_Errhandler errhandler, MPI_Session *session)
    MPICH_API_PUBLIC;
MPI_Aint PMPI_Aint_add(MPI_Aint base, MPI_Aint disp) MPICH_API_PUBLIC;
MPI_Aint PMPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) MPICH_API_PUBLIC;
int PMPI_Get_library_version(char *version, int *resultlen) MPICH_API_PUBLIC;
int PMPI_Get_processor_name(char *name, int *resultlen) MPICH_API_PUBLIC;
int PMPI_Get_version(int *version, int *subversion) MPICH_API_PUBLIC;
int PMPI_Pcontrol(const int level, ...) MPICH_API_PUBLIC;
int PMPIX_GPU_query_support(int gpu_type, int *is_supported) MPICH_API_PUBLIC;
int PMPIX_Query_cuda_support(void) MPICH_API_PUBLIC;
int PMPIX_Query_ze_support(void) MPICH_API_PUBLIC;
int PMPIX_Query_hip_support(void) MPICH_API_PUBLIC;
int PMPI_T_category_changed(int *update_number) MPICH_API_PUBLIC;
int PMPI_T_category_get_categories(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int PMPI_T_category_get_cvars(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int PMPI_T_category_get_events(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int PMPI_T_category_get_index(const char *name, int *cat_index) MPICH_API_PUBLIC;
int PMPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                             int *num_cvars, int *num_pvars, int *num_categories) MPICH_API_PUBLIC;
int PMPI_T_category_get_num(int *num_cat) MPICH_API_PUBLIC;
int PMPI_T_category_get_num_events(int cat_index, int *num_events) MPICH_API_PUBLIC;
int PMPI_T_category_get_pvars(int cat_index, int len, int indices[]) MPICH_API_PUBLIC;
int PMPI_T_cvar_get_index(const char *name, int *cvar_index) MPICH_API_PUBLIC;
int PMPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *bind, int *scope) MPICH_API_PUBLIC;
int PMPI_T_cvar_get_num(int *num_cvar) MPICH_API_PUBLIC;
int PMPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                             int *count) MPICH_API_PUBLIC;
int PMPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) MPICH_API_PUBLIC;
int PMPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) MPICH_API_PUBLIC;
int PMPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) MPICH_API_PUBLIC;
int PMPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len)
    MPICH_API_PUBLIC;
int PMPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len)
    MPICH_API_PUBLIC;
int PMPI_T_event_callback_get_info(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info *info_used)
                                   MPICH_API_PUBLIC;
int PMPI_T_event_callback_set_info(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info) MPICH_API_PUBLIC;
int PMPI_T_event_copy(MPI_T_event_instance event_instance, void *buffer) MPICH_API_PUBLIC;
int PMPI_T_event_get_index(const char *name, int *event_index) MPICH_API_PUBLIC;
int PMPI_T_event_get_info(int event_index, char *name, int *name_len, int *verbosity,
                          MPI_Datatype array_of_datatypes[], MPI_Aint array_of_displacements[],
                          int *num_elements, MPI_T_enum *enumtype, MPI_Info *info, char *desc,
                          int *desc_len, int *bind) MPICH_API_PUBLIC;
int PMPI_T_event_get_num(int *num_events) MPICH_API_PUBLIC;
int PMPI_T_event_get_source(MPI_T_event_instance event_instance, int *source_index)
    MPICH_API_PUBLIC;
int PMPI_T_event_get_timestamp(MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
    MPICH_API_PUBLIC;
int PMPI_T_event_handle_alloc(int event_index, void *obj_handle, MPI_Info info,
                              MPI_T_event_registration *event_registration) MPICH_API_PUBLIC;
int PMPI_T_event_handle_free(MPI_T_event_registration event_registration, void *user_data,
                             MPI_T_event_free_cb_function free_cb_function) MPICH_API_PUBLIC;
int PMPI_T_event_handle_get_info(MPI_T_event_registration event_registration, MPI_Info *info_used)
    MPICH_API_PUBLIC;
int PMPI_T_event_handle_set_info(MPI_T_event_registration event_registration, MPI_Info info)
    MPICH_API_PUBLIC;
int PMPI_T_event_read(MPI_T_event_instance event_instance, int element_index, void *buffer)
    MPICH_API_PUBLIC;
int PMPI_T_event_register_callback(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                   MPI_T_event_cb_function event_cb_function) MPICH_API_PUBLIC;
int PMPI_T_event_set_dropped_handler(MPI_T_event_registration event_registration,
                                     MPI_T_event_dropped_cb_function dropped_cb_function)
                                     MPICH_API_PUBLIC;
int PMPI_T_finalize(void) MPICH_API_PUBLIC;
int PMPI_T_init_thread(int required, int *provided) MPICH_API_PUBLIC;
int PMPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) MPICH_API_PUBLIC;
int PMPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *bind, int *readonly, int *continuous, int *atomic) MPICH_API_PUBLIC;
int PMPI_T_pvar_get_num(int *num_pvar) MPICH_API_PUBLIC;
int PMPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                             MPI_T_pvar_handle *handle, int *count) MPICH_API_PUBLIC;
int PMPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle)
    MPICH_API_PUBLIC;
int PMPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
    MPICH_API_PUBLIC;
int PMPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
    MPICH_API_PUBLIC;
int PMPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int PMPI_T_pvar_session_create(MPI_T_pvar_session *session) MPICH_API_PUBLIC;
int PMPI_T_pvar_session_free(MPI_T_pvar_session *session) MPICH_API_PUBLIC;
int PMPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int PMPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int PMPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf)
    MPICH_API_PUBLIC;
int PMPI_T_source_get_info(int source_index, char *name, int *name_len, char *desc, int *desc_len,
                           MPI_T_source_order *ordering, MPI_Count *ticks_per_second,
                           MPI_Count *max_ticks, MPI_Info *info) MPICH_API_PUBLIC;
int PMPI_T_source_get_num(int *num_sources) MPICH_API_PUBLIC;
int PMPI_T_source_get_timestamp(int source_index, MPI_Count *timestamp) MPICH_API_PUBLIC;
int PMPI_Op_commutative(MPI_Op op, int *commute) MPICH_API_PUBLIC;
int PMPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) MPICH_API_PUBLIC;
int PMPI_Op_create_c(MPI_User_function_c *user_fn, int commute, MPI_Op *op) MPICH_API_PUBLIC;
int PMPI_Op_free(MPI_Op *op) MPICH_API_PUBLIC;
int PMPI_Parrived(MPI_Request request, int partition, int *flag) MPICH_API_PUBLIC;
int PMPI_Pready(int partition, MPI_Request request) MPICH_API_PUBLIC;
int PMPI_Pready_list(int length, const int array_of_partitions[], MPI_Request request)
    MPICH_API_PUBLIC;
int PMPI_Pready_range(int partition_low, int partition_high, MPI_Request request) MPICH_API_PUBLIC;
int PMPI_Precv_init(void *buf, int partitions, MPI_Count count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_API_PUBLIC;
int PMPI_Psend_init(const void *buf, int partitions, MPI_Count count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_API_PUBLIC;
int PMPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Bsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Buffer_attach(void *buffer, int size) MPICH_API_PUBLIC;
int PMPI_Buffer_attach_c(void *buffer, MPI_Count size) MPICH_API_PUBLIC;
int PMPI_Buffer_detach(void *buffer_addr, int *size) MPICH_API_PUBLIC;
int PMPI_Buffer_detach_c(void *buffer_addr, MPI_Count *size) MPICH_API_PUBLIC;
int PMPI_Buffer_flush(void) MPICH_API_PUBLIC;
int PMPI_Buffer_iflush(MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Comm_attach_buffer(MPI_Comm comm, void *buffer, int size) MPICH_API_PUBLIC;
int PMPI_Comm_attach_buffer_c(MPI_Comm comm, void *buffer, MPI_Count size) MPICH_API_PUBLIC;
int PMPI_Comm_detach_buffer(MPI_Comm comm, void *buffer_addr, int *size) MPICH_API_PUBLIC;
int PMPI_Comm_detach_buffer_c(MPI_Comm comm, void *buffer_addr, MPI_Count *size) MPICH_API_PUBLIC;
int PMPI_Comm_flush_buffer(MPI_Comm comm) MPICH_API_PUBLIC;
int PMPI_Comm_iflush_buffer(MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Ibsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                 MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Imrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                  MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status)
    MPICH_API_PUBLIC;
int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Irsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Isendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int PMPI_Isendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                     int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                     int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int PMPI_Isendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                           int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Isendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                             int sendtag, int source, int recvtag, MPI_Comm comm,
                             MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Issend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status)
    MPICH_API_PUBLIC;
int PMPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Mrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                 MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                MPI_Comm comm, MPI_Status *status)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                   MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Recv_init_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Send_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Status *status)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int PMPI_Sendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                    int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                    int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int PMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                          int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Sendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                            int sendtag, int source, int recvtag, MPI_Comm comm,
                            MPI_Status *status)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Session_attach_buffer(MPI_Session session, void *buffer, int size) MPICH_API_PUBLIC;
int PMPI_Session_attach_buffer_c(MPI_Session session, void *buffer, MPI_Count size)
    MPICH_API_PUBLIC;
int PMPI_Session_detach_buffer(MPI_Session session, void *buffer_addr, int *size) MPICH_API_PUBLIC;
int PMPI_Session_detach_buffer_c(MPI_Session session, void *buffer_addr, MPI_Count *size)
    MPICH_API_PUBLIC;
int PMPI_Session_flush_buffer(MPI_Session session) MPICH_API_PUBLIC;
int PMPI_Session_iflush_buffer(MPI_Session session, MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Ssend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Ssend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Cancel(MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Grequest_complete(MPI_Request request) MPICH_API_PUBLIC;
int PMPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Request_free(MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Request_get_status_all(int count, MPI_Request array_of_requests[], int *flag,
                                MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int PMPI_Request_get_status_any(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                                MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Request_get_status_some(int incount, MPI_Request array_of_requests[], int *outcount,
                                 int array_of_indices[], MPI_Status *array_of_statuses)
                                 MPICH_API_PUBLIC;
int PMPI_Start(MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Startall(int count, MPI_Request array_of_requests[]) MPICH_API_PUBLIC;
int PMPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) MPICH_API_PUBLIC;
int PMPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) MPICH_API_PUBLIC;
int PMPI_Status_get_error(MPI_Status *status, int *error) MPICH_API_PUBLIC;
int PMPI_Status_get_source(MPI_Status *status, int *source) MPICH_API_PUBLIC;
int PMPI_Status_get_tag(MPI_Status *status, int *tag) MPICH_API_PUBLIC;
int PMPI_Status_set_error(MPI_Status *status, int error) MPICH_API_PUBLIC;
int PMPI_Status_set_source(MPI_Status *status, int source) MPICH_API_PUBLIC;
int PMPI_Status_set_tag(MPI_Status *status, int tag) MPICH_API_PUBLIC;
int PMPI_Status_set_cancelled(MPI_Status *status, int flag) MPICH_API_PUBLIC;
int PMPI_Test(MPI_Request *request, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Test_cancelled(const MPI_Status *status, int *flag) MPICH_API_PUBLIC;
int PMPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                 MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int PMPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                 MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int PMPI_Wait(MPI_Request *request, MPI_Status *status) MPICH_API_PUBLIC;
int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
    MPICH_API_PUBLIC;
int PMPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status)
    MPICH_API_PUBLIC;
int PMPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int PMPIX_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                         void *extra_state, MPI_Request *request) MPICH_API_PUBLIC;
int PMPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class) MPICH_API_PUBLIC;
int PMPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                  MPI_Request *request) MPICH_API_PUBLIC;
int PMPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Accumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                      MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) MPICH_API_PUBLIC;
int PMPI_Compare_and_swap(const void *origin_addr, const void *compare_addr, void *result_addr,
                          MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Fetch_and_op(const void *origin_addr, void *result_addr, MPI_Datatype datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win)
                      MPICH_API_PUBLIC;
int PMPI_Free_mem(void *base) MPICH_API_PUBLIC;
int PMPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win)
             MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Get_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                        void *result_addr, int result_count, MPI_Datatype result_datatype,
                        int target_rank, MPI_Aint target_disp, int target_count,
                        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Get_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                          MPI_Win win)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
             MPI_Win win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Put_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Raccumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                     int target_rank, MPI_Aint target_disp, int target_count,
                     MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Raccumulate_c(const void *origin_addr, MPI_Count origin_count,
                       MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                       MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
              MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rget_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                         void *result_addr, int result_count, MPI_Datatype result_datatype,
                         int target_rank, MPI_Aint target_disp, int target_count,
                         MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Rget_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                           MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                           MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                           MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                           MPI_Win win, MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int PMPI_Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
              MPI_Win win, MPI_Request *request)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Rput_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                      MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_allocate_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                        void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_allocate_shared_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                               void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) MPICH_API_PUBLIC;
int PMPI_Win_complete(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                    MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_create_c(void *base, MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                      MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_detach(MPI_Win win, const void *base) MPICH_API_PUBLIC;
int PMPI_Win_fence(int assert, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_flush(int rank, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_flush_all(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_flush_local(int rank, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_flush_local_all(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_free(MPI_Win *win) MPICH_API_PUBLIC;
int PMPI_Win_get_group(MPI_Win win, MPI_Group *group) MPICH_API_PUBLIC;
int PMPI_Win_get_info(MPI_Win win, MPI_Info *info_used) MPICH_API_PUBLIC;
int PMPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) MPICH_API_PUBLIC;
int PMPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_lock_all(int assert, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_post(MPI_Group group, int assert, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_set_info(MPI_Win win, MPI_Info info) MPICH_API_PUBLIC;
int PMPI_Win_set_name(MPI_Win win, const char *win_name) MPICH_API_PUBLIC;
int PMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr)
    MPICH_API_PUBLIC;
int PMPI_Win_shared_query_c(MPI_Win win, int rank, MPI_Aint *size, MPI_Aint *disp_unit,
                            void *baseptr) MPICH_API_PUBLIC;
int PMPI_Win_start(MPI_Group group, int assert, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_sync(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_test(MPI_Win win, int *flag) MPICH_API_PUBLIC;
int PMPI_Win_unlock(int rank, MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_unlock_all(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Win_wait(MPI_Win win) MPICH_API_PUBLIC;
int PMPI_Close_port(const char *port_name) MPICH_API_PUBLIC;
int PMPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                      MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Comm_disconnect(MPI_Comm *comm) MPICH_API_PUBLIC;
int PMPI_Comm_get_parent(MPI_Comm *parent) MPICH_API_PUBLIC;
int PMPI_Comm_join(int fd, MPI_Comm *intercomm) MPICH_API_PUBLIC;
int PMPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                    MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) MPICH_API_PUBLIC;
int PMPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                             MPICH_API_PUBLIC;
int PMPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) MPICH_API_PUBLIC;
int PMPI_Open_port(MPI_Info info, char *port_name) MPICH_API_PUBLIC;
int PMPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
    MPICH_API_PUBLIC;
int PMPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
    MPICH_API_PUBLIC;
int PMPIX_Stream_create(MPI_Info info, MPIX_Stream *stream) MPICH_API_PUBLIC;
int PMPIX_Stream_free(MPIX_Stream *stream) MPICH_API_PUBLIC;
int PMPIX_Stream_comm_create(MPI_Comm comm, MPIX_Stream stream, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int PMPIX_Stream_comm_create_multiplex(MPI_Comm comm, int count, MPIX_Stream array_of_streams[],
                                       MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPIX_Comm_get_stream(MPI_Comm comm, int idx, MPIX_Stream *stream) MPICH_API_PUBLIC;
int PMPIX_Stream_progress(MPIX_Stream stream) MPICH_API_PUBLIC;
int PMPIX_Start_progress_thread(MPIX_Stream stream) MPICH_API_PUBLIC;
int PMPIX_Stop_progress_thread(MPIX_Stream stream) MPICH_API_PUBLIC;
int PMPIX_Stream_send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm, int source_stream_index, int dest_stream_index,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Stream_irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, int source_stream_index, int dest_stream_index,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Send_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Send_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Recv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Status *status)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Recv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, MPI_Status *status)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Isend_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Isend_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                          int tag, MPI_Comm comm, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Irecv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Irecv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                          MPI_Comm comm, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_API_PUBLIC;
int PMPIX_Wait_enqueue(MPI_Request *request, MPI_Status *status) MPICH_API_PUBLIC;
int PMPIX_Waitall_enqueue(int count, MPI_Request array_of_requests[],
                          MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int PMPIX_Allreduce_enqueue(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                            MPI_Op op, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPIX_Allreduce_enqueue_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4) MPICH_API_PUBLIC;
int PMPIX_Threadcomm_init(MPI_Comm comm, int num_threads, MPI_Comm *newthreadcomm)
    MPICH_API_PUBLIC;
int PMPIX_Threadcomm_free(MPI_Comm *threadcomm) MPICH_API_PUBLIC;
int PMPIX_Threadcomm_start(MPI_Comm threadcomm) MPICH_API_PUBLIC;
int PMPIX_Threadcomm_finish(MPI_Comm threadcomm) MPICH_API_PUBLIC;
double PMPI_Wtick(void) MPICH_API_PUBLIC;
double PMPI_Wtime(void) MPICH_API_PUBLIC;
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) MPICH_API_PUBLIC;
int PMPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                     int reorder, MPI_Comm *comm_cart) MPICH_API_PUBLIC;
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[])
    MPICH_API_PUBLIC;
int PMPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank)
    MPICH_API_PUBLIC;
int PMPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) MPICH_API_PUBLIC;
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)
    MPICH_API_PUBLIC;
int PMPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) MPICH_API_PUBLIC;
int PMPI_Cartdim_get(MPI_Comm comm, int *ndims) MPICH_API_PUBLIC;
int PMPI_Dims_create(int nnodes, int ndims, int dims[]) MPICH_API_PUBLIC;
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                           const int destinations[], const int weights[], MPI_Info info,
                           int reorder, MPI_Comm *comm_dist_graph) MPICH_API_PUBLIC;
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                    const int sourceweights[], int outdegree,
                                    const int destinations[], const int destweights[],
                                    MPI_Info info, int reorder, MPI_Comm *comm_dist_graph)
                                    MPICH_API_PUBLIC;
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                              int maxoutdegree, int destinations[], int destweights[])
                              MPICH_API_PUBLIC;
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted)
    MPICH_API_PUBLIC;
int PMPI_Get_hw_resource_info(MPI_Info *hw_info) MPICH_API_PUBLIC;
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                      int reorder, MPI_Comm *comm_graph) MPICH_API_PUBLIC;
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[])
    MPICH_API_PUBLIC;
int PMPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank)
    MPICH_API_PUBLIC;
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[])
    MPICH_API_PUBLIC;
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) MPICH_API_PUBLIC;
int PMPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) MPICH_API_PUBLIC;
int PMPI_Topo_test(MPI_Comm comm, int *status) MPICH_API_PUBLIC;

enum QMPI_Functions_enum {
    MPI_ATTR_DELETE_T,
    MPI_ATTR_GET_T,
    MPI_ATTR_PUT_T,
    MPI_COMM_CREATE_KEYVAL_T,
    MPI_COMM_DELETE_ATTR_T,
    MPI_COMM_FREE_KEYVAL_T,
    MPI_COMM_GET_ATTR_T,
    MPI_COMM_SET_ATTR_T,
    MPI_KEYVAL_CREATE_T,
    MPI_KEYVAL_FREE_T,
    MPI_TYPE_CREATE_KEYVAL_T,
    MPI_TYPE_DELETE_ATTR_T,
    MPI_TYPE_FREE_KEYVAL_T,
    MPI_TYPE_GET_ATTR_T,
    MPI_TYPE_SET_ATTR_T,
    MPI_WIN_CREATE_KEYVAL_T,
    MPI_WIN_DELETE_ATTR_T,
    MPI_WIN_FREE_KEYVAL_T,
    MPI_WIN_GET_ATTR_T,
    MPI_WIN_SET_ATTR_T,
    MPI_ALLGATHER_T,
    MPI_ALLGATHER_C_T,
    MPI_ALLGATHER_INIT_T,
    MPI_ALLGATHER_INIT_C_T,
    MPI_ALLGATHERV_T,
    MPI_ALLGATHERV_C_T,
    MPI_ALLGATHERV_INIT_T,
    MPI_ALLGATHERV_INIT_C_T,
    MPI_ALLREDUCE_T,
    MPI_ALLREDUCE_C_T,
    MPI_ALLREDUCE_INIT_T,
    MPI_ALLREDUCE_INIT_C_T,
    MPI_ALLTOALL_T,
    MPI_ALLTOALL_C_T,
    MPI_ALLTOALL_INIT_T,
    MPI_ALLTOALL_INIT_C_T,
    MPI_ALLTOALLV_T,
    MPI_ALLTOALLV_C_T,
    MPI_ALLTOALLV_INIT_T,
    MPI_ALLTOALLV_INIT_C_T,
    MPI_ALLTOALLW_T,
    MPI_ALLTOALLW_C_T,
    MPI_ALLTOALLW_INIT_T,
    MPI_ALLTOALLW_INIT_C_T,
    MPI_BARRIER_T,
    MPI_BARRIER_INIT_T,
    MPI_BCAST_T,
    MPI_BCAST_C_T,
    MPI_BCAST_INIT_T,
    MPI_BCAST_INIT_C_T,
    MPI_EXSCAN_T,
    MPI_EXSCAN_C_T,
    MPI_EXSCAN_INIT_T,
    MPI_EXSCAN_INIT_C_T,
    MPI_GATHER_T,
    MPI_GATHER_C_T,
    MPI_GATHER_INIT_T,
    MPI_GATHER_INIT_C_T,
    MPI_GATHERV_T,
    MPI_GATHERV_C_T,
    MPI_GATHERV_INIT_T,
    MPI_GATHERV_INIT_C_T,
    MPI_IALLGATHER_T,
    MPI_IALLGATHER_C_T,
    MPI_IALLGATHERV_T,
    MPI_IALLGATHERV_C_T,
    MPI_IALLREDUCE_T,
    MPI_IALLREDUCE_C_T,
    MPI_IALLTOALL_T,
    MPI_IALLTOALL_C_T,
    MPI_IALLTOALLV_T,
    MPI_IALLTOALLV_C_T,
    MPI_IALLTOALLW_T,
    MPI_IALLTOALLW_C_T,
    MPI_IBARRIER_T,
    MPI_IBCAST_T,
    MPI_IBCAST_C_T,
    MPI_IEXSCAN_T,
    MPI_IEXSCAN_C_T,
    MPI_IGATHER_T,
    MPI_IGATHER_C_T,
    MPI_IGATHERV_T,
    MPI_IGATHERV_C_T,
    MPI_INEIGHBOR_ALLGATHER_T,
    MPI_INEIGHBOR_ALLGATHER_C_T,
    MPI_INEIGHBOR_ALLGATHERV_T,
    MPI_INEIGHBOR_ALLGATHERV_C_T,
    MPI_INEIGHBOR_ALLTOALL_T,
    MPI_INEIGHBOR_ALLTOALL_C_T,
    MPI_INEIGHBOR_ALLTOALLV_T,
    MPI_INEIGHBOR_ALLTOALLV_C_T,
    MPI_INEIGHBOR_ALLTOALLW_T,
    MPI_INEIGHBOR_ALLTOALLW_C_T,
    MPI_IREDUCE_T,
    MPI_IREDUCE_C_T,
    MPI_IREDUCE_SCATTER_T,
    MPI_IREDUCE_SCATTER_C_T,
    MPI_IREDUCE_SCATTER_BLOCK_T,
    MPI_IREDUCE_SCATTER_BLOCK_C_T,
    MPI_ISCAN_T,
    MPI_ISCAN_C_T,
    MPI_ISCATTER_T,
    MPI_ISCATTER_C_T,
    MPI_ISCATTERV_T,
    MPI_ISCATTERV_C_T,
    MPI_NEIGHBOR_ALLGATHER_T,
    MPI_NEIGHBOR_ALLGATHER_C_T,
    MPI_NEIGHBOR_ALLGATHER_INIT_T,
    MPI_NEIGHBOR_ALLGATHER_INIT_C_T,
    MPI_NEIGHBOR_ALLGATHERV_T,
    MPI_NEIGHBOR_ALLGATHERV_C_T,
    MPI_NEIGHBOR_ALLGATHERV_INIT_T,
    MPI_NEIGHBOR_ALLGATHERV_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALL_T,
    MPI_NEIGHBOR_ALLTOALL_C_T,
    MPI_NEIGHBOR_ALLTOALL_INIT_T,
    MPI_NEIGHBOR_ALLTOALL_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALLV_T,
    MPI_NEIGHBOR_ALLTOALLV_C_T,
    MPI_NEIGHBOR_ALLTOALLV_INIT_T,
    MPI_NEIGHBOR_ALLTOALLV_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALLW_T,
    MPI_NEIGHBOR_ALLTOALLW_C_T,
    MPI_NEIGHBOR_ALLTOALLW_INIT_T,
    MPI_NEIGHBOR_ALLTOALLW_INIT_C_T,
    MPI_REDUCE_T,
    MPI_REDUCE_C_T,
    MPI_REDUCE_INIT_T,
    MPI_REDUCE_INIT_C_T,
    MPI_REDUCE_LOCAL_T,
    MPI_REDUCE_LOCAL_C_T,
    MPI_REDUCE_SCATTER_T,
    MPI_REDUCE_SCATTER_C_T,
    MPI_REDUCE_SCATTER_BLOCK_T,
    MPI_REDUCE_SCATTER_BLOCK_C_T,
    MPI_REDUCE_SCATTER_BLOCK_INIT_T,
    MPI_REDUCE_SCATTER_BLOCK_INIT_C_T,
    MPI_REDUCE_SCATTER_INIT_T,
    MPI_REDUCE_SCATTER_INIT_C_T,
    MPI_SCAN_T,
    MPI_SCAN_C_T,
    MPI_SCAN_INIT_T,
    MPI_SCAN_INIT_C_T,
    MPI_SCATTER_T,
    MPI_SCATTER_C_T,
    MPI_SCATTER_INIT_T,
    MPI_SCATTER_INIT_C_T,
    MPI_SCATTERV_T,
    MPI_SCATTERV_C_T,
    MPI_SCATTERV_INIT_T,
    MPI_SCATTERV_INIT_C_T,
    MPI_COMM_COMPARE_T,
    MPI_COMM_CREATE_T,
    MPI_COMM_CREATE_GROUP_T,
    MPI_COMM_DUP_T,
    MPI_COMM_DUP_WITH_INFO_T,
    MPI_COMM_FREE_T,
    MPI_COMM_GET_INFO_T,
    MPI_COMM_GET_NAME_T,
    MPI_COMM_GROUP_T,
    MPI_COMM_IDUP_T,
    MPI_COMM_IDUP_WITH_INFO_T,
    MPI_COMM_RANK_T,
    MPI_COMM_REMOTE_GROUP_T,
    MPI_COMM_REMOTE_SIZE_T,
    MPI_COMM_SET_INFO_T,
    MPI_COMM_SET_NAME_T,
    MPI_COMM_SIZE_T,
    MPI_COMM_SPLIT_T,
    MPI_COMM_SPLIT_TYPE_T,
    MPI_COMM_TEST_INTER_T,
    MPI_INTERCOMM_CREATE_T,
    MPI_INTERCOMM_CREATE_FROM_GROUPS_T,
    MPI_INTERCOMM_MERGE_T,
    MPIX_COMM_TEST_THREADCOMM_T,
    MPIX_COMM_REVOKE_T,
    MPIX_COMM_SHRINK_T,
    MPIX_COMM_FAILURE_ACK_T,
    MPIX_COMM_FAILURE_GET_ACKED_T,
    MPIX_COMM_AGREE_T,
    MPIX_COMM_GET_FAILED_T,
    MPI_GET_ADDRESS_T,
    MPI_GET_COUNT_T,
    MPI_GET_COUNT_C_T,
    MPI_GET_ELEMENTS_T,
    MPI_GET_ELEMENTS_C_T,
    MPI_GET_ELEMENTS_X_T,
    MPI_PACK_T,
    MPI_PACK_C_T,
    MPI_PACK_EXTERNAL_T,
    MPI_PACK_EXTERNAL_C_T,
    MPI_PACK_EXTERNAL_SIZE_T,
    MPI_PACK_EXTERNAL_SIZE_C_T,
    MPI_PACK_SIZE_T,
    MPI_PACK_SIZE_C_T,
    MPI_STATUS_SET_ELEMENTS_T,
    MPI_STATUS_SET_ELEMENTS_C_T,
    MPI_STATUS_SET_ELEMENTS_X_T,
    MPI_TYPE_COMMIT_T,
    MPI_TYPE_CONTIGUOUS_T,
    MPI_TYPE_CONTIGUOUS_C_T,
    MPI_TYPE_CREATE_DARRAY_T,
    MPI_TYPE_CREATE_DARRAY_C_T,
    MPI_TYPE_CREATE_HINDEXED_T,
    MPI_TYPE_CREATE_HINDEXED_C_T,
    MPI_TYPE_CREATE_HINDEXED_BLOCK_T,
    MPI_TYPE_CREATE_HINDEXED_BLOCK_C_T,
    MPI_TYPE_CREATE_HVECTOR_T,
    MPI_TYPE_CREATE_HVECTOR_C_T,
    MPI_TYPE_CREATE_INDEXED_BLOCK_T,
    MPI_TYPE_CREATE_INDEXED_BLOCK_C_T,
    MPI_TYPE_CREATE_RESIZED_T,
    MPI_TYPE_CREATE_RESIZED_C_T,
    MPI_TYPE_CREATE_STRUCT_T,
    MPI_TYPE_CREATE_STRUCT_C_T,
    MPI_TYPE_CREATE_SUBARRAY_T,
    MPI_TYPE_CREATE_SUBARRAY_C_T,
    MPI_TYPE_DUP_T,
    MPI_TYPE_FREE_T,
    MPI_TYPE_GET_CONTENTS_T,
    MPI_TYPE_GET_CONTENTS_C_T,
    MPI_TYPE_GET_ENVELOPE_T,
    MPI_TYPE_GET_ENVELOPE_C_T,
    MPI_TYPE_GET_EXTENT_T,
    MPI_TYPE_GET_EXTENT_C_T,
    MPI_TYPE_GET_EXTENT_X_T,
    MPI_TYPE_GET_NAME_T,
    MPI_TYPE_GET_TRUE_EXTENT_T,
    MPI_TYPE_GET_TRUE_EXTENT_C_T,
    MPI_TYPE_GET_TRUE_EXTENT_X_T,
    MPI_TYPE_GET_VALUE_INDEX_T,
    MPI_TYPE_INDEXED_T,
    MPI_TYPE_INDEXED_C_T,
    MPI_TYPE_MATCH_SIZE_T,
    MPI_TYPE_SET_NAME_T,
    MPI_TYPE_SIZE_T,
    MPI_TYPE_SIZE_C_T,
    MPI_TYPE_SIZE_X_T,
    MPI_TYPE_VECTOR_T,
    MPI_TYPE_VECTOR_C_T,
    MPI_UNPACK_T,
    MPI_UNPACK_C_T,
    MPI_UNPACK_EXTERNAL_T,
    MPI_UNPACK_EXTERNAL_C_T,
    MPI_ADDRESS_T,
    MPI_TYPE_EXTENT_T,
    MPI_TYPE_LB_T,
    MPI_TYPE_UB_T,
    MPI_TYPE_HINDEXED_T,
    MPI_TYPE_HVECTOR_T,
    MPI_TYPE_STRUCT_T,
    MPIX_TYPE_IOV_LEN_T,
    MPIX_TYPE_IOV_T,
    MPI_ADD_ERROR_CLASS_T,
    MPI_ADD_ERROR_CODE_T,
    MPI_ADD_ERROR_STRING_T,
    MPI_COMM_CALL_ERRHANDLER_T,
    MPI_COMM_CREATE_ERRHANDLER_T,
    MPI_COMM_GET_ERRHANDLER_T,
    MPI_COMM_SET_ERRHANDLER_T,
    MPI_ERRHANDLER_FREE_T,
    MPI_ERROR_CLASS_T,
    MPI_ERROR_STRING_T,
    MPI_FILE_CALL_ERRHANDLER_T,
    MPI_FILE_CREATE_ERRHANDLER_T,
    MPI_FILE_GET_ERRHANDLER_T,
    MPI_FILE_SET_ERRHANDLER_T,
    MPI_REMOVE_ERROR_CLASS_T,
    MPI_REMOVE_ERROR_CODE_T,
    MPI_REMOVE_ERROR_STRING_T,
    MPI_SESSION_CALL_ERRHANDLER_T,
    MPI_SESSION_CREATE_ERRHANDLER_T,
    MPI_SESSION_GET_ERRHANDLER_T,
    MPI_SESSION_SET_ERRHANDLER_T,
    MPI_WIN_CALL_ERRHANDLER_T,
    MPI_WIN_CREATE_ERRHANDLER_T,
    MPI_WIN_GET_ERRHANDLER_T,
    MPI_WIN_SET_ERRHANDLER_T,
    MPI_ERRHANDLER_CREATE_T,
    MPI_ERRHANDLER_GET_T,
    MPI_ERRHANDLER_SET_T,
    MPI_COMM_C2F_T,
    MPI_COMM_F2C_T,
    MPI_ERRHANDLER_C2F_T,
    MPI_ERRHANDLER_F2C_T,
    MPI_GROUP_C2F_T,
    MPI_GROUP_F2C_T,
    MPI_INFO_C2F_T,
    MPI_INFO_F2C_T,
    MPI_MESSAGE_C2F_T,
    MPI_MESSAGE_F2C_T,
    MPI_OP_C2F_T,
    MPI_OP_F2C_T,
    MPI_REQUEST_C2F_T,
    MPI_REQUEST_F2C_T,
    MPI_SESSION_C2F_T,
    MPI_SESSION_F2C_T,
    MPI_TYPE_C2F_T,
    MPI_TYPE_F2C_T,
    MPI_WIN_C2F_T,
    MPI_WIN_F2C_T,
    MPI_GROUP_COMPARE_T,
    MPI_GROUP_DIFFERENCE_T,
    MPI_GROUP_EXCL_T,
    MPI_GROUP_FREE_T,
    MPI_GROUP_INCL_T,
    MPI_GROUP_INTERSECTION_T,
    MPI_GROUP_RANGE_EXCL_T,
    MPI_GROUP_RANGE_INCL_T,
    MPI_GROUP_RANK_T,
    MPI_GROUP_SIZE_T,
    MPI_GROUP_TRANSLATE_RANKS_T,
    MPI_GROUP_UNION_T,
    MPI_INFO_CREATE_T,
    MPI_INFO_CREATE_ENV_T,
    MPI_INFO_DELETE_T,
    MPI_INFO_DUP_T,
    MPI_INFO_FREE_T,
    MPI_INFO_GET_T,
    MPI_INFO_GET_NKEYS_T,
    MPI_INFO_GET_NTHKEY_T,
    MPI_INFO_GET_STRING_T,
    MPI_INFO_GET_VALUELEN_T,
    MPI_INFO_SET_T,
    MPIX_INFO_SET_HEX_T,
    MPI_ABORT_T,
    MPI_COMM_CREATE_FROM_GROUP_T,
    MPI_FINALIZE_T,
    MPI_FINALIZED_T,
    MPI_GROUP_FROM_SESSION_PSET_T,
    MPI_INIT_T,
    MPI_INIT_THREAD_T,
    MPI_INITIALIZED_T,
    MPI_IS_THREAD_MAIN_T,
    MPI_QUERY_THREAD_T,
    MPI_SESSION_FINALIZE_T,
    MPI_SESSION_GET_INFO_T,
    MPI_SESSION_GET_NTH_PSET_T,
    MPI_SESSION_GET_NUM_PSETS_T,
    MPI_SESSION_GET_PSET_INFO_T,
    MPI_SESSION_INIT_T,
    MPI_AINT_ADD_T,
    MPI_AINT_DIFF_T,
    MPI_GET_LIBRARY_VERSION_T,
    MPI_GET_PROCESSOR_NAME_T,
    MPI_GET_VERSION_T,
    MPI_PCONTROL_T,
    MPIX_GPU_QUERY_SUPPORT_T,
    MPIX_QUERY_CUDA_SUPPORT_T,
    MPIX_QUERY_ZE_SUPPORT_T,
    MPIX_QUERY_HIP_SUPPORT_T,
    MPI_T_CATEGORY_CHANGED_T,
    MPI_T_CATEGORY_GET_CATEGORIES_T,
    MPI_T_CATEGORY_GET_CVARS_T,
    MPI_T_CATEGORY_GET_EVENTS_T,
    MPI_T_CATEGORY_GET_INDEX_T,
    MPI_T_CATEGORY_GET_INFO_T,
    MPI_T_CATEGORY_GET_NUM_T,
    MPI_T_CATEGORY_GET_NUM_EVENTS_T,
    MPI_T_CATEGORY_GET_PVARS_T,
    MPI_T_CVAR_GET_INDEX_T,
    MPI_T_CVAR_GET_INFO_T,
    MPI_T_CVAR_GET_NUM_T,
    MPI_T_CVAR_HANDLE_ALLOC_T,
    MPI_T_CVAR_HANDLE_FREE_T,
    MPI_T_CVAR_READ_T,
    MPI_T_CVAR_WRITE_T,
    MPI_T_ENUM_GET_INFO_T,
    MPI_T_ENUM_GET_ITEM_T,
    MPI_T_EVENT_CALLBACK_GET_INFO_T,
    MPI_T_EVENT_CALLBACK_SET_INFO_T,
    MPI_T_EVENT_COPY_T,
    MPI_T_EVENT_GET_INDEX_T,
    MPI_T_EVENT_GET_INFO_T,
    MPI_T_EVENT_GET_NUM_T,
    MPI_T_EVENT_GET_SOURCE_T,
    MPI_T_EVENT_GET_TIMESTAMP_T,
    MPI_T_EVENT_HANDLE_ALLOC_T,
    MPI_T_EVENT_HANDLE_FREE_T,
    MPI_T_EVENT_HANDLE_GET_INFO_T,
    MPI_T_EVENT_HANDLE_SET_INFO_T,
    MPI_T_EVENT_READ_T,
    MPI_T_EVENT_REGISTER_CALLBACK_T,
    MPI_T_EVENT_SET_DROPPED_HANDLER_T,
    MPI_T_FINALIZE_T,
    MPI_T_INIT_THREAD_T,
    MPI_T_PVAR_GET_INDEX_T,
    MPI_T_PVAR_GET_INFO_T,
    MPI_T_PVAR_GET_NUM_T,
    MPI_T_PVAR_HANDLE_ALLOC_T,
    MPI_T_PVAR_HANDLE_FREE_T,
    MPI_T_PVAR_READ_T,
    MPI_T_PVAR_READRESET_T,
    MPI_T_PVAR_RESET_T,
    MPI_T_PVAR_SESSION_CREATE_T,
    MPI_T_PVAR_SESSION_FREE_T,
    MPI_T_PVAR_START_T,
    MPI_T_PVAR_STOP_T,
    MPI_T_PVAR_WRITE_T,
    MPI_T_SOURCE_GET_INFO_T,
    MPI_T_SOURCE_GET_NUM_T,
    MPI_T_SOURCE_GET_TIMESTAMP_T,
    MPI_OP_COMMUTATIVE_T,
    MPI_OP_CREATE_T,
    MPI_OP_CREATE_C_T,
    MPI_OP_FREE_T,
    MPI_PARRIVED_T,
    MPI_PREADY_T,
    MPI_PREADY_LIST_T,
    MPI_PREADY_RANGE_T,
    MPI_PRECV_INIT_T,
    MPI_PSEND_INIT_T,
    MPI_BSEND_T,
    MPI_BSEND_C_T,
    MPI_BSEND_INIT_T,
    MPI_BSEND_INIT_C_T,
    MPI_BUFFER_ATTACH_T,
    MPI_BUFFER_ATTACH_C_T,
    MPI_BUFFER_DETACH_T,
    MPI_BUFFER_DETACH_C_T,
    MPI_BUFFER_FLUSH_T,
    MPI_BUFFER_IFLUSH_T,
    MPI_COMM_ATTACH_BUFFER_T,
    MPI_COMM_ATTACH_BUFFER_C_T,
    MPI_COMM_DETACH_BUFFER_T,
    MPI_COMM_DETACH_BUFFER_C_T,
    MPI_COMM_FLUSH_BUFFER_T,
    MPI_COMM_IFLUSH_BUFFER_T,
    MPI_IBSEND_T,
    MPI_IBSEND_C_T,
    MPI_IMPROBE_T,
    MPI_IMRECV_T,
    MPI_IMRECV_C_T,
    MPI_IPROBE_T,
    MPI_IRECV_T,
    MPI_IRECV_C_T,
    MPI_IRSEND_T,
    MPI_IRSEND_C_T,
    MPI_ISEND_T,
    MPI_ISEND_C_T,
    MPI_ISENDRECV_T,
    MPI_ISENDRECV_C_T,
    MPI_ISENDRECV_REPLACE_T,
    MPI_ISENDRECV_REPLACE_C_T,
    MPI_ISSEND_T,
    MPI_ISSEND_C_T,
    MPI_MPROBE_T,
    MPI_MRECV_T,
    MPI_MRECV_C_T,
    MPI_PROBE_T,
    MPI_RECV_T,
    MPI_RECV_C_T,
    MPI_RECV_INIT_T,
    MPI_RECV_INIT_C_T,
    MPI_RSEND_T,
    MPI_RSEND_C_T,
    MPI_RSEND_INIT_T,
    MPI_RSEND_INIT_C_T,
    MPI_SEND_T,
    MPI_SEND_C_T,
    MPI_SEND_INIT_T,
    MPI_SEND_INIT_C_T,
    MPI_SENDRECV_T,
    MPI_SENDRECV_C_T,
    MPI_SENDRECV_REPLACE_T,
    MPI_SENDRECV_REPLACE_C_T,
    MPI_SESSION_ATTACH_BUFFER_T,
    MPI_SESSION_ATTACH_BUFFER_C_T,
    MPI_SESSION_DETACH_BUFFER_T,
    MPI_SESSION_DETACH_BUFFER_C_T,
    MPI_SESSION_FLUSH_BUFFER_T,
    MPI_SESSION_IFLUSH_BUFFER_T,
    MPI_SSEND_T,
    MPI_SSEND_C_T,
    MPI_SSEND_INIT_T,
    MPI_SSEND_INIT_C_T,
    MPI_CANCEL_T,
    MPI_GREQUEST_COMPLETE_T,
    MPI_GREQUEST_START_T,
    MPI_REQUEST_FREE_T,
    MPI_REQUEST_GET_STATUS_T,
    MPI_REQUEST_GET_STATUS_ALL_T,
    MPI_REQUEST_GET_STATUS_ANY_T,
    MPI_REQUEST_GET_STATUS_SOME_T,
    MPI_START_T,
    MPI_STARTALL_T,
    MPI_STATUS_C2F_T,
    MPI_STATUS_F2C_T,
    MPI_STATUS_GET_ERROR_T,
    MPI_STATUS_GET_SOURCE_T,
    MPI_STATUS_GET_TAG_T,
    MPI_STATUS_SET_ERROR_T,
    MPI_STATUS_SET_SOURCE_T,
    MPI_STATUS_SET_TAG_T,
    MPI_STATUS_SET_CANCELLED_T,
    MPI_TEST_T,
    MPI_TEST_CANCELLED_T,
    MPI_TESTALL_T,
    MPI_TESTANY_T,
    MPI_TESTSOME_T,
    MPI_WAIT_T,
    MPI_WAITALL_T,
    MPI_WAITANY_T,
    MPI_WAITSOME_T,
    MPIX_GREQUEST_START_T,
    MPIX_GREQUEST_CLASS_CREATE_T,
    MPIX_GREQUEST_CLASS_ALLOCATE_T,
    MPI_ACCUMULATE_T,
    MPI_ACCUMULATE_C_T,
    MPI_ALLOC_MEM_T,
    MPI_COMPARE_AND_SWAP_T,
    MPI_FETCH_AND_OP_T,
    MPI_FREE_MEM_T,
    MPI_GET_T,
    MPI_GET_C_T,
    MPI_GET_ACCUMULATE_T,
    MPI_GET_ACCUMULATE_C_T,
    MPI_PUT_T,
    MPI_PUT_C_T,
    MPI_RACCUMULATE_T,
    MPI_RACCUMULATE_C_T,
    MPI_RGET_T,
    MPI_RGET_C_T,
    MPI_RGET_ACCUMULATE_T,
    MPI_RGET_ACCUMULATE_C_T,
    MPI_RPUT_T,
    MPI_RPUT_C_T,
    MPI_WIN_ALLOCATE_T,
    MPI_WIN_ALLOCATE_C_T,
    MPI_WIN_ALLOCATE_SHARED_T,
    MPI_WIN_ALLOCATE_SHARED_C_T,
    MPI_WIN_ATTACH_T,
    MPI_WIN_COMPLETE_T,
    MPI_WIN_CREATE_T,
    MPI_WIN_CREATE_C_T,
    MPI_WIN_CREATE_DYNAMIC_T,
    MPI_WIN_DETACH_T,
    MPI_WIN_FENCE_T,
    MPI_WIN_FLUSH_T,
    MPI_WIN_FLUSH_ALL_T,
    MPI_WIN_FLUSH_LOCAL_T,
    MPI_WIN_FLUSH_LOCAL_ALL_T,
    MPI_WIN_FREE_T,
    MPI_WIN_GET_GROUP_T,
    MPI_WIN_GET_INFO_T,
    MPI_WIN_GET_NAME_T,
    MPI_WIN_LOCK_T,
    MPI_WIN_LOCK_ALL_T,
    MPI_WIN_POST_T,
    MPI_WIN_SET_INFO_T,
    MPI_WIN_SET_NAME_T,
    MPI_WIN_SHARED_QUERY_T,
    MPI_WIN_SHARED_QUERY_C_T,
    MPI_WIN_START_T,
    MPI_WIN_SYNC_T,
    MPI_WIN_TEST_T,
    MPI_WIN_UNLOCK_T,
    MPI_WIN_UNLOCK_ALL_T,
    MPI_WIN_WAIT_T,
    MPI_CLOSE_PORT_T,
    MPI_COMM_ACCEPT_T,
    MPI_COMM_CONNECT_T,
    MPI_COMM_DISCONNECT_T,
    MPI_COMM_GET_PARENT_T,
    MPI_COMM_JOIN_T,
    MPI_COMM_SPAWN_T,
    MPI_COMM_SPAWN_MULTIPLE_T,
    MPI_LOOKUP_NAME_T,
    MPI_OPEN_PORT_T,
    MPI_PUBLISH_NAME_T,
    MPI_UNPUBLISH_NAME_T,
    MPIX_STREAM_CREATE_T,
    MPIX_STREAM_FREE_T,
    MPIX_STREAM_COMM_CREATE_T,
    MPIX_STREAM_COMM_CREATE_MULTIPLEX_T,
    MPIX_COMM_GET_STREAM_T,
    MPIX_STREAM_PROGRESS_T,
    MPIX_START_PROGRESS_THREAD_T,
    MPIX_STOP_PROGRESS_THREAD_T,
    MPIX_STREAM_SEND_T,
    MPIX_STREAM_SEND_C_T,
    MPIX_STREAM_ISEND_T,
    MPIX_STREAM_ISEND_C_T,
    MPIX_STREAM_RECV_T,
    MPIX_STREAM_RECV_C_T,
    MPIX_STREAM_IRECV_T,
    MPIX_STREAM_IRECV_C_T,
    MPIX_SEND_ENQUEUE_T,
    MPIX_SEND_ENQUEUE_C_T,
    MPIX_RECV_ENQUEUE_T,
    MPIX_RECV_ENQUEUE_C_T,
    MPIX_ISEND_ENQUEUE_T,
    MPIX_ISEND_ENQUEUE_C_T,
    MPIX_IRECV_ENQUEUE_T,
    MPIX_IRECV_ENQUEUE_C_T,
    MPIX_WAIT_ENQUEUE_T,
    MPIX_WAITALL_ENQUEUE_T,
    MPIX_ALLREDUCE_ENQUEUE_T,
    MPIX_ALLREDUCE_ENQUEUE_C_T,
    MPIX_THREADCOMM_INIT_T,
    MPIX_THREADCOMM_FREE_T,
    MPIX_THREADCOMM_START_T,
    MPIX_THREADCOMM_FINISH_T,
    MPI_WTICK_T,
    MPI_WTIME_T,
    MPI_CART_COORDS_T,
    MPI_CART_CREATE_T,
    MPI_CART_GET_T,
    MPI_CART_MAP_T,
    MPI_CART_RANK_T,
    MPI_CART_SHIFT_T,
    MPI_CART_SUB_T,
    MPI_CARTDIM_GET_T,
    MPI_DIMS_CREATE_T,
    MPI_DIST_GRAPH_CREATE_T,
    MPI_DIST_GRAPH_CREATE_ADJACENT_T,
    MPI_DIST_GRAPH_NEIGHBORS_T,
    MPI_DIST_GRAPH_NEIGHBORS_COUNT_T,
    MPI_GET_HW_RESOURCE_INFO_T,
    MPI_GRAPH_CREATE_T,
    MPI_GRAPH_GET_T,
    MPI_GRAPH_MAP_T,
    MPI_GRAPH_NEIGHBORS_T,
    MPI_GRAPH_NEIGHBORS_COUNT_T,
    MPI_GRAPHDIMS_GET_T,
    MPI_TOPO_TEST_T,
    MPI_LAST_FUNC_T
};

int QMPI_Attr_delete(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval)
    MPICH_API_PUBLIC;
int QMPI_Attr_get(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval, void *attribute_val,
                  int *flag) MPICH_API_PUBLIC;
int QMPI_Attr_put(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
                  void *attribute_val) MPICH_API_PUBLIC;
int QMPI_Comm_create_keyval(QMPI_Context context, int tool_id,
                            MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state) MPICH_API_PUBLIC;
int QMPI_Comm_delete_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval)
    MPICH_API_PUBLIC;
int QMPI_Comm_free_keyval(QMPI_Context context, int tool_id, int *comm_keyval) MPICH_API_PUBLIC;
int QMPI_Comm_get_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval,
                       void *attribute_val, int *flag) MPICH_API_PUBLIC;
int QMPI_Comm_set_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval,
                       void *attribute_val) MPICH_API_PUBLIC;
int QMPI_Keyval_create(QMPI_Context context, int tool_id, MPI_Copy_function *copy_fn,
                       MPI_Delete_function *delete_fn, int *keyval, void *extra_state)
                       MPICH_API_PUBLIC;
int QMPI_Keyval_free(QMPI_Context context, int tool_id, int *keyval) MPICH_API_PUBLIC;
int QMPI_Type_create_keyval(QMPI_Context context, int tool_id,
                            MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                            void *extra_state) MPICH_API_PUBLIC;
int QMPI_Type_delete_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                          int type_keyval) MPICH_API_PUBLIC;
int QMPI_Type_free_keyval(QMPI_Context context, int tool_id, int *type_keyval) MPICH_API_PUBLIC;
int QMPI_Type_get_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype, int type_keyval,
                       void *attribute_val, int *flag) MPICH_API_PUBLIC;
int QMPI_Type_set_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype, int type_keyval,
                       void *attribute_val) MPICH_API_PUBLIC;
int QMPI_Win_create_keyval(QMPI_Context context, int tool_id,
                           MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state) MPICH_API_PUBLIC;
int QMPI_Win_delete_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval)
    MPICH_API_PUBLIC;
int QMPI_Win_free_keyval(QMPI_Context context, int tool_id, int *win_keyval) MPICH_API_PUBLIC;
int QMPI_Win_get_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
                      void *attribute_val, int *flag) MPICH_API_PUBLIC;
int QMPI_Win_set_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
                      void *attribute_val) MPICH_API_PUBLIC;
int QMPI_Allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Allgather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Allgather_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                        MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Allgather_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                          MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int displs[], MPI_Datatype recvtype, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                      MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Allgatherv_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                         MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                         const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Allgatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                           MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                           const MPI_Count recvcounts[], const MPI_Aint displs[],
                           MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                           MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Allreduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Allreduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Allreduce_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                        int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Allreduce_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                          MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Alltoall_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                       MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Alltoall_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                         MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                         MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Alltoallv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
                   const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                   MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                     MPI_Datatype recvtype, MPI_Comm comm)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Alltoallv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                        const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                        void *recvbuf, const int recvcounts[], const int rdispls[],
                        MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Alltoallv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                          MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                          const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Alltoallw(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                   const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                   MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                     MPICH_API_PUBLIC;
int QMPI_Alltoallw_init(QMPI_Context context, int tool_id, const void *sendbuf,
                        const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
                        void *recvbuf, const int recvcounts[], const int rdispls[],
                        const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                        MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Alltoallw_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                          const MPI_Datatype sendtypes[], void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                          MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Barrier(QMPI_Context context, int tool_id, MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Barrier_init(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Bcast(QMPI_Context context, int tool_id, void *buffer, int count, MPI_Datatype datatype,
               int root, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bcast_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                 MPI_Datatype datatype, int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bcast_init(QMPI_Context context, int tool_id, void *buffer, int count,
                    MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bcast_init_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                      MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Exscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Exscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                  MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Exscan_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Exscan_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                       MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Gather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Gather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                  MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Gather_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                     MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                     int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Gather_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                       MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                       MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Gatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Gatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                   const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Gatherv_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                      MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                      const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                      MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Gatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                        MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                        const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                        MPI_Info info, MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Iallgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                    MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Iallgather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                      MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Iallgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                     MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                     const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Iallgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                       MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                       const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                       MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Iallreduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                    int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Iallreduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                      MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                      MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ialltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ialltoall_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ialltoallv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                    const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Ialltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                      const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Ialltoallw(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                    const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                    const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                    MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Ialltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                      const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                      MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Ibarrier(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Request *request)
    MPICH_API_PUBLIC;
int QMPI_Ibcast(QMPI_Context context, int tool_id, void *buffer, int count, MPI_Datatype datatype,
                int root, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Ibcast_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                  MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Iexscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Iexscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                   MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Igather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Igather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Igatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
                  MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Igatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                    MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Ineighbor_allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                             MPI_Datatype sendtype, void *recvbuf, int recvcount,
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ineighbor_allgather_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                               MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ineighbor_allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                              MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                              const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
                              MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Ineighbor_allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                const MPI_Count recvcounts[], const MPI_Aint displs[],
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                              MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                              MPI_Request *request)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoallv(QMPI_Context context, int tool_id, const void *sendbuf,
                             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int rdispls[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                               MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                               const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoallw(QMPI_Context context, int tool_id, const void *sendbuf,
                             const int sendcounts[], const MPI_Aint sdispls[],
                             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                             MPI_Comm comm, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Ineighbor_alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                               const MPI_Datatype sendtypes[], void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                               const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request)
                               MPICH_API_PUBLIC;
int QMPI_Ireduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                 MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ireduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                   MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ireduce_scatter(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                         const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ireduce_scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                           const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
                           MPI_Comm comm, MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ireduce_scatter_block(QMPI_Context context, int tool_id, const void *sendbuf,
                               void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
                               MPI_Comm comm, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Ireduce_scatter_block_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Iscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Iscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                 MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                 MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Iscatter(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Iscatter_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Iscatterv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Iscatterv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                     void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgather_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                              MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgather_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgather_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                   MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                   MPI_Info info, MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int displs[], MPI_Datatype recvtype, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgatherv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                  int sendcount, MPI_Datatype sendtype, void *recvbuf,
                                  const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Neighbor_allgatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                    MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                    const MPI_Count recvcounts[], const MPI_Aint displs[],
                                    MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                    MPI_Request *request)
                                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,9) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                           MPI_Datatype sendtype, void *recvbuf, int recvcount,
                           MPI_Datatype recvtype, MPI_Comm comm)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf,
                             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoall_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                MPI_Request *request)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoall_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                  MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                  MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                  MPI_Info info, MPI_Request *request)
                                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallv(QMPI_Context context, int tool_id, const void *sendbuf,
                            const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int rdispls[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                              const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                                 void *recvbuf, const int recvcounts[], const int rdispls[],
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                                   MPI_Datatype sendtype, void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,10) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallw(QMPI_Context context, int tool_id, const void *sendbuf,
                            const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              const MPI_Datatype sendtypes[], void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallw_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 const int sendcounts[], const MPI_Aint sdispls[],
                                 const MPI_Datatype sendtypes[], void *recvbuf,
                                 const int recvcounts[], const MPI_Aint rdispls[],
                                 const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Neighbor_alltoallw_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                                   const MPI_Datatype sendtypes[], void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Reduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                  MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                       MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_local(QMPI_Context context, int tool_id, const void *inbuf, void *inoutbuf,
                      int count, MPI_Datatype datatype, MPI_Op op)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_local_c(QMPI_Context context, int tool_id, const void *inbuf, void *inoutbuf,
                        MPI_Count count, MPI_Datatype datatype, MPI_Op op)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                        const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                          const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
                          MPI_Comm comm)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_block(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                              int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_block_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                MPI_Op op, MPI_Comm comm)
                                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_block_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                   void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_block_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                     void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                             const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                             MPI_Comm comm, MPI_Info info, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Reduce_scatter_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype,
                               MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Scan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Scan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Scan_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Scan_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPI_Scatter(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Scatter_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                      MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Scatter_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                        MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Scatterv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                  const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, int root, MPI_Comm comm)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Scatterv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                    const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                    MPI_Comm comm)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Scatterv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                       const int sendcounts[], const int displs[], MPI_Datatype sendtype,
                       void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Scatterv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                         const MPI_Count sendcounts[], const MPI_Aint displs[],
                         MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                         MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_ATTR_POINTER_WITH_TYPE_TAG(7,9) MPICH_API_PUBLIC;
int QMPI_Comm_compare(QMPI_Context context, int tool_id, MPI_Comm comm1, MPI_Comm comm2,
                      int *result) MPICH_API_PUBLIC;
int QMPI_Comm_create(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
                     MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_create_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
                           int tag, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_dup(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int QMPI_Comm_dup_with_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                            MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_free(QMPI_Context context, int tool_id, MPI_Comm *comm) MPICH_API_PUBLIC;
int QMPI_Comm_get_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info *info_used)
    MPICH_API_PUBLIC;
int QMPI_Comm_get_name(QMPI_Context context, int tool_id, MPI_Comm comm, char *comm_name,
                       int *resultlen) MPICH_API_PUBLIC;
int QMPI_Comm_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *group)
    MPICH_API_PUBLIC;
int QMPI_Comm_idup(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm,
                   MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Comm_idup_with_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                             MPI_Comm *newcomm, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Comm_rank(QMPI_Context context, int tool_id, MPI_Comm comm, int *rank) MPICH_API_PUBLIC;
int QMPI_Comm_remote_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *group)
    MPICH_API_PUBLIC;
int QMPI_Comm_remote_size(QMPI_Context context, int tool_id, MPI_Comm comm, int *size)
    MPICH_API_PUBLIC;
int QMPI_Comm_set_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info)
    MPICH_API_PUBLIC;
int QMPI_Comm_set_name(QMPI_Context context, int tool_id, MPI_Comm comm, const char *comm_name)
    MPICH_API_PUBLIC;
int QMPI_Comm_size(QMPI_Context context, int tool_id, MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int QMPI_Comm_split(QMPI_Context context, int tool_id, MPI_Comm comm, int color, int key,
                    MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_split_type(QMPI_Context context, int tool_id, MPI_Comm comm, int split_type, int key,
                         MPI_Info info, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_test_inter(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag)
    MPICH_API_PUBLIC;
int QMPI_Intercomm_create(QMPI_Context context, int tool_id, MPI_Comm local_comm, int local_leader,
                          MPI_Comm peer_comm, int remote_leader, int tag, MPI_Comm *newintercomm)
                          MPICH_API_PUBLIC;
int QMPI_Intercomm_create_from_groups(QMPI_Context context, int tool_id, MPI_Group local_group,
                                      int local_leader, MPI_Group remote_group, int remote_leader,
                                      const char *stringtag, MPI_Info info,
                                      MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                      MPICH_API_PUBLIC;
int QMPI_Intercomm_merge(QMPI_Context context, int tool_id, MPI_Comm intercomm, int high,
                         MPI_Comm *newintracomm) MPICH_API_PUBLIC;
int QMPIX_Comm_test_threadcomm(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag)
    MPICH_API_PUBLIC;
int QMPIX_Comm_revoke(QMPI_Context context, int tool_id, MPI_Comm comm) MPICH_API_PUBLIC;
int QMPIX_Comm_shrink(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm)
    MPICH_API_PUBLIC;
int QMPIX_Comm_failure_ack(QMPI_Context context, int tool_id, MPI_Comm comm) MPICH_API_PUBLIC;
int QMPIX_Comm_failure_get_acked(QMPI_Context context, int tool_id, MPI_Comm comm,
                                 MPI_Group *failedgrp) MPICH_API_PUBLIC;
int QMPIX_Comm_agree(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag) MPICH_API_PUBLIC;
int QMPIX_Comm_get_failed(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *failedgrp)
    MPICH_API_PUBLIC;
int QMPI_Get_address(QMPI_Context context, int tool_id, const void *location, MPI_Aint *address)
    MPICH_API_PUBLIC;
int QMPI_Get_count(QMPI_Context context, int tool_id, const MPI_Status *status,
                   MPI_Datatype datatype, int *count) MPICH_API_PUBLIC;
int QMPI_Get_count_c(QMPI_Context context, int tool_id, const MPI_Status *status,
                     MPI_Datatype datatype, MPI_Count *count) MPICH_API_PUBLIC;
int QMPI_Get_elements(QMPI_Context context, int tool_id, const MPI_Status *status,
                      MPI_Datatype datatype, int *count) MPICH_API_PUBLIC;
int QMPI_Get_elements_c(QMPI_Context context, int tool_id, const MPI_Status *status,
                        MPI_Datatype datatype, MPI_Count *count) MPICH_API_PUBLIC;
int QMPI_Get_elements_x(QMPI_Context context, int tool_id, const MPI_Status *status,
                        MPI_Datatype datatype, MPI_Count *count) MPICH_API_PUBLIC;
int QMPI_Pack(QMPI_Context context, int tool_id, const void *inbuf, int incount,
              MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm)
              MPICH_API_PUBLIC;
int QMPI_Pack_c(QMPI_Context context, int tool_id, const void *inbuf, MPI_Count incount,
                MPI_Datatype datatype, void *outbuf, MPI_Count outsize, MPI_Count *position,
                MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Pack_external(QMPI_Context context, int tool_id, const char *datarep, const void *inbuf,
                       int incount, MPI_Datatype datatype, void *outbuf, MPI_Aint outsize,
                       MPI_Aint *position) MPICH_API_PUBLIC;
int QMPI_Pack_external_c(QMPI_Context context, int tool_id, const char *datarep, const void *inbuf,
                         MPI_Count incount, MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                         MPI_Count *position) MPICH_API_PUBLIC;
int QMPI_Pack_external_size(QMPI_Context context, int tool_id, const char *datarep, int incount,
                            MPI_Datatype datatype, MPI_Aint *size) MPICH_API_PUBLIC;
int QMPI_Pack_external_size_c(QMPI_Context context, int tool_id, const char *datarep,
                              MPI_Count incount, MPI_Datatype datatype, MPI_Count *size)
                              MPICH_API_PUBLIC;
int QMPI_Pack_size(QMPI_Context context, int tool_id, int incount, MPI_Datatype datatype,
                   MPI_Comm comm, int *size) MPICH_API_PUBLIC;
int QMPI_Pack_size_c(QMPI_Context context, int tool_id, MPI_Count incount, MPI_Datatype datatype,
                     MPI_Comm comm, MPI_Count *size) MPICH_API_PUBLIC;
int QMPI_Status_set_elements(QMPI_Context context, int tool_id, MPI_Status *status,
                             MPI_Datatype datatype, int count) MPICH_API_PUBLIC;
int QMPI_Status_set_elements_c(QMPI_Context context, int tool_id, MPI_Status *status,
                               MPI_Datatype datatype, MPI_Count count) MPICH_API_PUBLIC;
int QMPI_Status_set_elements_x(QMPI_Context context, int tool_id, MPI_Status *status,
                               MPI_Datatype datatype, MPI_Count count) MPICH_API_PUBLIC;
int QMPI_Type_commit(QMPI_Context context, int tool_id, MPI_Datatype *datatype) MPICH_API_PUBLIC;
int QMPI_Type_contiguous(QMPI_Context context, int tool_id, int count, MPI_Datatype oldtype,
                         MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_contiguous_c(QMPI_Context context, int tool_id, MPI_Count count, MPI_Datatype oldtype,
                           MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_darray(QMPI_Context context, int tool_id, int size, int rank, int ndims,
                            const int array_of_gsizes[], const int array_of_distribs[],
                            const int array_of_dargs[], const int array_of_psizes[], int order,
                            MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_darray_c(QMPI_Context context, int tool_id, int size, int rank, int ndims,
                              const MPI_Count array_of_gsizes[], const int array_of_distribs[],
                              const int array_of_dargs[], const int array_of_psizes[], int order,
                              MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_hindexed(QMPI_Context context, int tool_id, int count,
                              const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_hindexed_c(QMPI_Context context, int tool_id, MPI_Count count,
                                const MPI_Count array_of_blocklengths[],
                                const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_hindexed_block(QMPI_Context context, int tool_id, int count, int blocklength,
                                    const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_hindexed_block_c(QMPI_Context context, int tool_id, MPI_Count count,
                                      MPI_Count blocklength,
                                      const MPI_Count array_of_displacements[],
                                      MPI_Datatype oldtype, MPI_Datatype *newtype)
                                      MPICH_API_PUBLIC;
int QMPI_Type_create_hvector(QMPI_Context context, int tool_id, int count, int blocklength,
                             MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                             MPICH_API_PUBLIC;
int QMPI_Type_create_hvector_c(QMPI_Context context, int tool_id, MPI_Count count,
                               MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype,
                               MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_indexed_block(QMPI_Context context, int tool_id, int count, int blocklength,
                                   const int array_of_displacements[], MPI_Datatype oldtype,
                                   MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_indexed_block_c(QMPI_Context context, int tool_id, MPI_Count count,
                                     MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_resized(QMPI_Context context, int tool_id, MPI_Datatype oldtype, MPI_Aint lb,
                             MPI_Aint extent, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_resized_c(QMPI_Context context, int tool_id, MPI_Datatype oldtype,
                               MPI_Count lb, MPI_Count extent, MPI_Datatype *newtype)
                               MPICH_API_PUBLIC;
int QMPI_Type_create_struct(QMPI_Context context, int tool_id, int count,
                            const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                            MPICH_API_PUBLIC;
int QMPI_Type_create_struct_c(QMPI_Context context, int tool_id, MPI_Count count,
                              const MPI_Count array_of_blocklengths[],
                              const MPI_Count array_of_displacements[],
                              const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                              MPICH_API_PUBLIC;
int QMPI_Type_create_subarray(QMPI_Context context, int tool_id, int ndims,
                              const int array_of_sizes[], const int array_of_subsizes[],
                              const int array_of_starts[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_create_subarray_c(QMPI_Context context, int tool_id, int ndims,
                                const MPI_Count array_of_sizes[],
                                const MPI_Count array_of_subsizes[],
                                const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                                MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_dup(QMPI_Context context, int tool_id, MPI_Datatype oldtype, MPI_Datatype *newtype)
    MPICH_API_PUBLIC;
int QMPI_Type_free(QMPI_Context context, int tool_id, MPI_Datatype *datatype) MPICH_API_PUBLIC;
int QMPI_Type_get_contents(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                           int max_integers, int max_addresses, int max_datatypes,
                           int array_of_integers[], MPI_Aint array_of_addresses[],
                           MPI_Datatype array_of_datatypes[]) MPICH_API_PUBLIC;
int QMPI_Type_get_contents_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                             MPI_Count max_integers, MPI_Count max_addresses,
                             MPI_Count max_large_counts, MPI_Count max_datatypes,
                             int array_of_integers[], MPI_Aint array_of_addresses[],
                             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                             MPICH_API_PUBLIC;
int QMPI_Type_get_envelope(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                           int *num_integers, int *num_addresses, int *num_datatypes,
                           int *combiner) MPICH_API_PUBLIC;
int QMPI_Type_get_envelope_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                             MPI_Count *num_integers, MPI_Count *num_addresses,
                             MPI_Count *num_large_counts, MPI_Count *num_datatypes, int *combiner)
                             MPICH_API_PUBLIC;
int QMPI_Type_get_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *lb,
                         MPI_Aint *extent) MPICH_API_PUBLIC;
int QMPI_Type_get_extent_c(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *lb,
                           MPI_Count *extent) MPICH_API_PUBLIC;
int QMPI_Type_get_extent_x(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *lb,
                           MPI_Count *extent) MPICH_API_PUBLIC;
int QMPI_Type_get_name(QMPI_Context context, int tool_id, MPI_Datatype datatype, char *type_name,
                       int *resultlen) MPICH_API_PUBLIC;
int QMPI_Type_get_true_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                              MPI_Aint *true_lb, MPI_Aint *true_extent) MPICH_API_PUBLIC;
int QMPI_Type_get_true_extent_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                                MPI_Count *true_lb, MPI_Count *true_extent) MPICH_API_PUBLIC;
int QMPI_Type_get_true_extent_x(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                                MPI_Count *true_lb, MPI_Count *true_extent) MPICH_API_PUBLIC;
int QMPI_Type_get_value_index(QMPI_Context context, int tool_id, MPI_Datatype value_type,
                              MPI_Datatype index_type, MPI_Datatype *pair_type) MPICH_API_PUBLIC;
int QMPI_Type_indexed(QMPI_Context context, int tool_id, int count,
                      const int array_of_blocklengths[], const int array_of_displacements[],
                      MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_indexed_c(QMPI_Context context, int tool_id, MPI_Count count,
                        const MPI_Count array_of_blocklengths[],
                        const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                        MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_match_size(QMPI_Context context, int tool_id, int typeclass, int size,
                         MPI_Datatype *datatype) MPICH_API_PUBLIC;
int QMPI_Type_set_name(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                       const char *type_name) MPICH_API_PUBLIC;
int QMPI_Type_size(QMPI_Context context, int tool_id, MPI_Datatype datatype, int *size)
    MPICH_API_PUBLIC;
int QMPI_Type_size_c(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *size)
    MPICH_API_PUBLIC;
int QMPI_Type_size_x(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *size)
    MPICH_API_PUBLIC;
int QMPI_Type_vector(QMPI_Context context, int tool_id, int count, int blocklength, int stride,
                     MPI_Datatype oldtype, MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_vector_c(QMPI_Context context, int tool_id, MPI_Count count, MPI_Count blocklength,
                       MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                       MPICH_API_PUBLIC;
int QMPI_Unpack(QMPI_Context context, int tool_id, const void *inbuf, int insize, int *position,
                void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Unpack_c(QMPI_Context context, int tool_id, const void *inbuf, MPI_Count insize,
                  MPI_Count *position, void *outbuf, MPI_Count outcount, MPI_Datatype datatype,
                  MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Unpack_external(QMPI_Context context, int tool_id, const char datarep[], const void *inbuf,
                         MPI_Aint insize, MPI_Aint *position, void *outbuf, int outcount,
                         MPI_Datatype datatype) MPICH_API_PUBLIC;
int QMPI_Unpack_external_c(QMPI_Context context, int tool_id, const char datarep[],
                           const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                           MPI_Count outcount, MPI_Datatype datatype) MPICH_API_PUBLIC;
int QMPI_Address(QMPI_Context context, int tool_id, void *location, MPI_Aint *address)
    MPICH_API_PUBLIC;
int QMPI_Type_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *extent)
    MPICH_API_PUBLIC;
int QMPI_Type_lb(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *displacement)
    MPICH_API_PUBLIC;
int QMPI_Type_ub(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *displacement)
    MPICH_API_PUBLIC;
int QMPI_Type_hindexed(QMPI_Context context, int tool_id, int count, int array_of_blocklengths[],
                       MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                       MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPI_Type_hvector(QMPI_Context context, int tool_id, int count, int blocklength,
                      MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                      MPICH_API_PUBLIC;
int QMPI_Type_struct(QMPI_Context context, int tool_id, int count, int array_of_blocklengths[],
                     MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[],
                     MPI_Datatype *newtype) MPICH_API_PUBLIC;
int QMPIX_Type_iov_len(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                       MPI_Count max_iov_bytes, MPI_Count *iov_len, MPI_Count *actual_iov_bytes)
                       MPICH_API_PUBLIC;
int QMPIX_Type_iov(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count iov_offset,
                   MPIX_Iov *iov, MPI_Count max_iov_len, MPI_Count *actual_iov_len)
                   MPICH_API_PUBLIC;
int QMPI_Add_error_class(QMPI_Context context, int tool_id, int *errorclass) MPICH_API_PUBLIC;
int QMPI_Add_error_code(QMPI_Context context, int tool_id, int errorclass, int *errorcode)
    MPICH_API_PUBLIC;
int QMPI_Add_error_string(QMPI_Context context, int tool_id, int errorcode, const char *string)
    MPICH_API_PUBLIC;
int QMPI_Comm_call_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode)
    MPICH_API_PUBLIC;
int QMPI_Comm_create_errhandler(QMPI_Context context, int tool_id,
                                MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Comm_get_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm,
                             MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Comm_set_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm,
                             MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int QMPI_Errhandler_free(QMPI_Context context, int tool_id, MPI_Errhandler *errhandler)
    MPICH_API_PUBLIC;
int QMPI_Error_class(QMPI_Context context, int tool_id, int errorcode, int *errorclass)
    MPICH_API_PUBLIC;
int QMPI_Error_string(QMPI_Context context, int tool_id, int errorcode, char *string,
                      int *resultlen) MPICH_API_PUBLIC;
int QMPI_File_call_errhandler(QMPI_Context context, int tool_id, MPI_File fh, int errorcode)
    MPICH_API_PUBLIC;
int QMPI_File_create_errhandler(QMPI_Context context, int tool_id,
                                MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_File_get_errhandler(QMPI_Context context, int tool_id, MPI_File file,
                             MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_File_set_errhandler(QMPI_Context context, int tool_id, MPI_File file,
                             MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int QMPI_Remove_error_class(QMPI_Context context, int tool_id, int errorclass) MPICH_API_PUBLIC;
int QMPI_Remove_error_code(QMPI_Context context, int tool_id, int errorcode) MPICH_API_PUBLIC;
int QMPI_Remove_error_string(QMPI_Context context, int tool_id, int errorcode) MPICH_API_PUBLIC;
int QMPI_Session_call_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                 int errorcode) MPICH_API_PUBLIC;
int QMPI_Session_create_errhandler(QMPI_Context context, int tool_id,
                                   MPI_Session_errhandler_function *session_errhandler_fn,
                                   MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Session_get_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Session_set_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int QMPI_Win_call_errhandler(QMPI_Context context, int tool_id, MPI_Win win, int errorcode)
    MPICH_API_PUBLIC;
int QMPI_Win_create_errhandler(QMPI_Context context, int tool_id,
                               MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Win_get_errhandler(QMPI_Context context, int tool_id, MPI_Win win,
                            MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Win_set_errhandler(QMPI_Context context, int tool_id, MPI_Win win,
                            MPI_Errhandler errhandler) MPICH_API_PUBLIC;
int QMPI_Errhandler_create(QMPI_Context context, int tool_id,
                           MPI_Comm_errhandler_function *comm_errhandler_fn,
                           MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Errhandler_get(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler *errhandler) MPICH_API_PUBLIC;
int QMPI_Errhandler_set(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler errhandler) MPICH_API_PUBLIC;
MPI_Fint QMPI_Comm_c2f(QMPI_Context context, int tool_id, MPI_Comm comm) MPICH_API_PUBLIC;
MPI_Comm QMPI_Comm_f2c(QMPI_Context context, int tool_id, MPI_Fint comm) MPICH_API_PUBLIC;
MPI_Fint QMPI_Errhandler_c2f(QMPI_Context context, int tool_id, MPI_Errhandler errhandler)
    MPICH_API_PUBLIC;
MPI_Errhandler QMPI_Errhandler_f2c(QMPI_Context context, int tool_id, MPI_Fint errhandler)
    MPICH_API_PUBLIC;
MPI_Fint QMPI_Group_c2f(QMPI_Context context, int tool_id, MPI_Group group) MPICH_API_PUBLIC;
MPI_Group QMPI_Group_f2c(QMPI_Context context, int tool_id, MPI_Fint group) MPICH_API_PUBLIC;
MPI_Fint QMPI_Info_c2f(QMPI_Context context, int tool_id, MPI_Info info) MPICH_API_PUBLIC;
MPI_Info QMPI_Info_f2c(QMPI_Context context, int tool_id, MPI_Fint info) MPICH_API_PUBLIC;
MPI_Fint QMPI_Message_c2f(QMPI_Context context, int tool_id, MPI_Message message) MPICH_API_PUBLIC;
MPI_Message QMPI_Message_f2c(QMPI_Context context, int tool_id, MPI_Fint message) MPICH_API_PUBLIC;
MPI_Fint QMPI_Op_c2f(QMPI_Context context, int tool_id, MPI_Op op) MPICH_API_PUBLIC;
MPI_Op QMPI_Op_f2c(QMPI_Context context, int tool_id, MPI_Fint op) MPICH_API_PUBLIC;
MPI_Fint QMPI_Request_c2f(QMPI_Context context, int tool_id, MPI_Request request) MPICH_API_PUBLIC;
MPI_Request QMPI_Request_f2c(QMPI_Context context, int tool_id, MPI_Fint request) MPICH_API_PUBLIC;
MPI_Fint QMPI_Session_c2f(QMPI_Context context, int tool_id, MPI_Session session) MPICH_API_PUBLIC;
MPI_Session QMPI_Session_f2c(QMPI_Context context, int tool_id, MPI_Fint session) MPICH_API_PUBLIC;
MPI_Fint QMPI_Type_c2f(QMPI_Context context, int tool_id, MPI_Datatype datatype) MPICH_API_PUBLIC;
MPI_Datatype QMPI_Type_f2c(QMPI_Context context, int tool_id, MPI_Fint datatype) MPICH_API_PUBLIC;
MPI_Fint QMPI_Win_c2f(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
MPI_Win QMPI_Win_f2c(QMPI_Context context, int tool_id, MPI_Fint win) MPICH_API_PUBLIC;
int QMPI_Group_compare(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                       int *result) MPICH_API_PUBLIC;
int QMPI_Group_difference(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                          MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_excl(QMPI_Context context, int tool_id, MPI_Group group, int n, const int ranks[],
                    MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_free(QMPI_Context context, int tool_id, MPI_Group *group) MPICH_API_PUBLIC;
int QMPI_Group_incl(QMPI_Context context, int tool_id, MPI_Group group, int n, const int ranks[],
                    MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_intersection(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                            MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_range_excl(QMPI_Context context, int tool_id, MPI_Group group, int n,
                          int ranges[][3], MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_range_incl(QMPI_Context context, int tool_id, MPI_Group group, int n,
                          int ranges[][3], MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Group_rank(QMPI_Context context, int tool_id, MPI_Group group, int *rank)
    MPICH_API_PUBLIC;
int QMPI_Group_size(QMPI_Context context, int tool_id, MPI_Group group, int *size)
    MPICH_API_PUBLIC;
int QMPI_Group_translate_ranks(QMPI_Context context, int tool_id, MPI_Group group1, int n,
                               const int ranks1[], MPI_Group group2, int ranks2[])
                               MPICH_API_PUBLIC;
int QMPI_Group_union(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                     MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Info_create(QMPI_Context context, int tool_id, MPI_Info *info) MPICH_API_PUBLIC;
int QMPI_Info_create_env(QMPI_Context context, int tool_id, int argc, char *argv[], MPI_Info *info)
    MPICH_API_PUBLIC;
int QMPI_Info_delete(QMPI_Context context, int tool_id, MPI_Info info, const char *key)
    MPICH_API_PUBLIC;
int QMPI_Info_dup(QMPI_Context context, int tool_id, MPI_Info info, MPI_Info *newinfo)
    MPICH_API_PUBLIC;
int QMPI_Info_free(QMPI_Context context, int tool_id, MPI_Info *info) MPICH_API_PUBLIC;
int QMPI_Info_get(QMPI_Context context, int tool_id, MPI_Info info, const char *key, int valuelen,
                  char *value, int *flag) MPICH_API_PUBLIC;
int QMPI_Info_get_nkeys(QMPI_Context context, int tool_id, MPI_Info info, int *nkeys)
    MPICH_API_PUBLIC;
int QMPI_Info_get_nthkey(QMPI_Context context, int tool_id, MPI_Info info, int n, char *key)
    MPICH_API_PUBLIC;
int QMPI_Info_get_string(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                         int *buflen, char *value, int *flag) MPICH_API_PUBLIC;
int QMPI_Info_get_valuelen(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                           int *valuelen, int *flag) MPICH_API_PUBLIC;
int QMPI_Info_set(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                  const char *value) MPICH_API_PUBLIC;
int QMPIX_Info_set_hex(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                       const void *value, int value_size) MPICH_API_PUBLIC;
int QMPI_Abort(QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode) MPICH_API_PUBLIC;
int QMPI_Comm_create_from_group(QMPI_Context context, int tool_id, MPI_Group group,
                                const char *stringtag, MPI_Info info, MPI_Errhandler errhandler,
                                MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Finalize(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPI_Finalized(QMPI_Context context, int tool_id, int *flag) MPICH_API_PUBLIC;
int QMPI_Group_from_session_pset(QMPI_Context context, int tool_id, MPI_Session session,
                                 const char *pset_name, MPI_Group *newgroup) MPICH_API_PUBLIC;
int QMPI_Init(QMPI_Context context, int tool_id, int *argc, char ***argv) MPICH_API_PUBLIC;
int QMPI_Init_thread(QMPI_Context context, int tool_id, int *argc, char ***argv, int required,
                     int *provided) MPICH_API_PUBLIC;
int QMPI_Initialized(QMPI_Context context, int tool_id, int *flag) MPICH_API_PUBLIC;
int QMPI_Is_thread_main(QMPI_Context context, int tool_id, int *flag) MPICH_API_PUBLIC;
int QMPI_Query_thread(QMPI_Context context, int tool_id, int *provided) MPICH_API_PUBLIC;
int QMPI_Session_finalize(QMPI_Context context, int tool_id, MPI_Session *session)
    MPICH_API_PUBLIC;
int QMPI_Session_get_info(QMPI_Context context, int tool_id, MPI_Session session,
                          MPI_Info *info_used) MPICH_API_PUBLIC;
int QMPI_Session_get_nth_pset(QMPI_Context context, int tool_id, MPI_Session session, MPI_Info info,
                              int n, int *pset_len, char *pset_name) MPICH_API_PUBLIC;
int QMPI_Session_get_num_psets(QMPI_Context context, int tool_id, MPI_Session session,
                               MPI_Info info, int *npset_names) MPICH_API_PUBLIC;
int QMPI_Session_get_pset_info(QMPI_Context context, int tool_id, MPI_Session session,
                               const char *pset_name, MPI_Info *info) MPICH_API_PUBLIC;
int QMPI_Session_init(QMPI_Context context, int tool_id, MPI_Info info, MPI_Errhandler errhandler,
                      MPI_Session *session) MPICH_API_PUBLIC;
MPI_Aint QMPI_Aint_add(QMPI_Context context, int tool_id, MPI_Aint base, MPI_Aint disp)
    MPICH_API_PUBLIC;
MPI_Aint QMPI_Aint_diff(QMPI_Context context, int tool_id, MPI_Aint addr1, MPI_Aint addr2)
    MPICH_API_PUBLIC;
int QMPI_Get_library_version(QMPI_Context context, int tool_id, char *version, int *resultlen)
    MPICH_API_PUBLIC;
int QMPI_Get_processor_name(QMPI_Context context, int tool_id, char *name, int *resultlen)
    MPICH_API_PUBLIC;
int QMPI_Get_version(QMPI_Context context, int tool_id, int *version, int *subversion)
    MPICH_API_PUBLIC;
int QMPI_Pcontrol(QMPI_Context context, int tool_id, const int level, ...) MPICH_API_PUBLIC;
int QMPIX_GPU_query_support(QMPI_Context context, int tool_id, int gpu_type, int *is_supported)
    MPICH_API_PUBLIC;
int QMPIX_Query_cuda_support(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPIX_Query_ze_support(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPIX_Query_hip_support(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPI_T_category_changed(QMPI_Context context, int tool_id, int *update_number)
    MPICH_API_PUBLIC;
int QMPI_T_category_get_categories(QMPI_Context context, int tool_id, int cat_index, int len,
                                   int indices[]) MPICH_API_PUBLIC;
int QMPI_T_category_get_cvars(QMPI_Context context, int tool_id, int cat_index, int len,
                              int indices[]) MPICH_API_PUBLIC;
int QMPI_T_category_get_events(QMPI_Context context, int tool_id, int cat_index, int len,
                               int indices[]) MPICH_API_PUBLIC;
int QMPI_T_category_get_index(QMPI_Context context, int tool_id, const char *name, int *cat_index)
    MPICH_API_PUBLIC;
int QMPI_T_category_get_info(QMPI_Context context, int tool_id, int cat_index, char *name,
                             int *name_len, char *desc, int *desc_len, int *num_cvars,
                             int *num_pvars, int *num_categories) MPICH_API_PUBLIC;
int QMPI_T_category_get_num(QMPI_Context context, int tool_id, int *num_cat) MPICH_API_PUBLIC;
int QMPI_T_category_get_num_events(QMPI_Context context, int tool_id, int cat_index,
                                   int *num_events) MPICH_API_PUBLIC;
int QMPI_T_category_get_pvars(QMPI_Context context, int tool_id, int cat_index, int len,
                              int indices[]) MPICH_API_PUBLIC;
int QMPI_T_cvar_get_index(QMPI_Context context, int tool_id, const char *name, int *cvar_index)
    MPICH_API_PUBLIC;
int QMPI_T_cvar_get_info(QMPI_Context context, int tool_id, int cvar_index, char *name,
                         int *name_len, int *verbosity, MPI_Datatype *datatype,
                         MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *scope)
                         MPICH_API_PUBLIC;
int QMPI_T_cvar_get_num(QMPI_Context context, int tool_id, int *num_cvar) MPICH_API_PUBLIC;
int QMPI_T_cvar_handle_alloc(QMPI_Context context, int tool_id, int cvar_index, void *obj_handle,
                             MPI_T_cvar_handle *handle, int *count) MPICH_API_PUBLIC;
int QMPI_T_cvar_handle_free(QMPI_Context context, int tool_id, MPI_T_cvar_handle *handle)
    MPICH_API_PUBLIC;
int QMPI_T_cvar_read(QMPI_Context context, int tool_id, MPI_T_cvar_handle handle, void *buf)
    MPICH_API_PUBLIC;
int QMPI_T_cvar_write(QMPI_Context context, int tool_id, MPI_T_cvar_handle handle, const void *buf)
    MPICH_API_PUBLIC;
int QMPI_T_enum_get_info(QMPI_Context context, int tool_id, MPI_T_enum enumtype, int *num,
                         char *name, int *name_len) MPICH_API_PUBLIC;
int QMPI_T_enum_get_item(QMPI_Context context, int tool_id, MPI_T_enum enumtype, int indx,
                         int *value, char *name, int *name_len) MPICH_API_PUBLIC;
int QMPI_T_event_callback_get_info(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info *info_used)
                                   MPICH_API_PUBLIC;
int QMPI_T_event_callback_set_info(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info) MPICH_API_PUBLIC;
int QMPI_T_event_copy(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                      void *buffer) MPICH_API_PUBLIC;
int QMPI_T_event_get_index(QMPI_Context context, int tool_id, const char *name, int *event_index)
    MPICH_API_PUBLIC;
int QMPI_T_event_get_info(QMPI_Context context, int tool_id, int event_index, char *name,
                          int *name_len, int *verbosity, MPI_Datatype array_of_datatypes[],
                          MPI_Aint array_of_displacements[], int *num_elements,
                          MPI_T_enum *enumtype, MPI_Info *info, char *desc, int *desc_len,
                          int *bind) MPICH_API_PUBLIC;
int QMPI_T_event_get_num(QMPI_Context context, int tool_id, int *num_events) MPICH_API_PUBLIC;
int QMPI_T_event_get_source(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                            int *source_index) MPICH_API_PUBLIC;
int QMPI_T_event_get_timestamp(QMPI_Context context, int tool_id,
                               MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
                               MPICH_API_PUBLIC;
int QMPI_T_event_handle_alloc(QMPI_Context context, int tool_id, int event_index, void *obj_handle,
                              MPI_Info info, MPI_T_event_registration *event_registration)
                              MPICH_API_PUBLIC;
int QMPI_T_event_handle_free(QMPI_Context context, int tool_id,
                             MPI_T_event_registration event_registration, void *user_data,
                             MPI_T_event_free_cb_function free_cb_function) MPICH_API_PUBLIC;
int QMPI_T_event_handle_get_info(QMPI_Context context, int tool_id,
                                 MPI_T_event_registration event_registration, MPI_Info *info_used)
                                 MPICH_API_PUBLIC;
int QMPI_T_event_handle_set_info(QMPI_Context context, int tool_id,
                                 MPI_T_event_registration event_registration, MPI_Info info)
                                 MPICH_API_PUBLIC;
int QMPI_T_event_read(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                      int element_index, void *buffer) MPICH_API_PUBLIC;
int QMPI_T_event_register_callback(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                   MPI_T_event_cb_function event_cb_function) MPICH_API_PUBLIC;
int QMPI_T_event_set_dropped_handler(QMPI_Context context, int tool_id,
                                     MPI_T_event_registration event_registration,
                                     MPI_T_event_dropped_cb_function dropped_cb_function)
                                     MPICH_API_PUBLIC;
int QMPI_T_finalize(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPI_T_init_thread(QMPI_Context context, int tool_id, int required, int *provided)
    MPICH_API_PUBLIC;
int QMPI_T_pvar_get_index(QMPI_Context context, int tool_id, const char *name, int var_class,
                          int *pvar_index) MPICH_API_PUBLIC;
int QMPI_T_pvar_get_info(QMPI_Context context, int tool_id, int pvar_index, char *name,
                         int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype,
                         MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly,
                         int *continuous, int *atomic) MPICH_API_PUBLIC;
int QMPI_T_pvar_get_num(QMPI_Context context, int tool_id, int *num_pvar) MPICH_API_PUBLIC;
int QMPI_T_pvar_handle_alloc(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                             int pvar_index, void *obj_handle, MPI_T_pvar_handle *handle,
                             int *count) MPICH_API_PUBLIC;
int QMPI_T_pvar_handle_free(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                            MPI_T_pvar_handle *handle) MPICH_API_PUBLIC;
int QMPI_T_pvar_read(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                     MPI_T_pvar_handle handle, void *buf) MPICH_API_PUBLIC;
int QMPI_T_pvar_readreset(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                          MPI_T_pvar_handle handle, void *buf) MPICH_API_PUBLIC;
int QMPI_T_pvar_reset(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int QMPI_T_pvar_session_create(QMPI_Context context, int tool_id, MPI_T_pvar_session *session)
    MPICH_API_PUBLIC;
int QMPI_T_pvar_session_free(QMPI_Context context, int tool_id, MPI_T_pvar_session *session)
    MPICH_API_PUBLIC;
int QMPI_T_pvar_start(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int QMPI_T_pvar_stop(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                     MPI_T_pvar_handle handle) MPICH_API_PUBLIC;
int QMPI_T_pvar_write(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle, const void *buf) MPICH_API_PUBLIC;
int QMPI_T_source_get_info(QMPI_Context context, int tool_id, int source_index, char *name,
                           int *name_len, char *desc, int *desc_len, MPI_T_source_order *ordering,
                           MPI_Count *ticks_per_second, MPI_Count *max_ticks, MPI_Info *info)
                           MPICH_API_PUBLIC;
int QMPI_T_source_get_num(QMPI_Context context, int tool_id, int *num_sources) MPICH_API_PUBLIC;
int QMPI_T_source_get_timestamp(QMPI_Context context, int tool_id, int source_index,
                                MPI_Count *timestamp) MPICH_API_PUBLIC;
int QMPI_Op_commutative(QMPI_Context context, int tool_id, MPI_Op op, int *commute)
    MPICH_API_PUBLIC;
int QMPI_Op_create(QMPI_Context context, int tool_id, MPI_User_function *user_fn, int commute,
                   MPI_Op *op) MPICH_API_PUBLIC;
int QMPI_Op_create_c(QMPI_Context context, int tool_id, MPI_User_function_c *user_fn, int commute,
                     MPI_Op *op) MPICH_API_PUBLIC;
int QMPI_Op_free(QMPI_Context context, int tool_id, MPI_Op *op) MPICH_API_PUBLIC;
int QMPI_Parrived(QMPI_Context context, int tool_id, MPI_Request request, int partition, int *flag)
    MPICH_API_PUBLIC;
int QMPI_Pready(QMPI_Context context, int tool_id, int partition, MPI_Request request)
    MPICH_API_PUBLIC;
int QMPI_Pready_list(QMPI_Context context, int tool_id, int length, const int array_of_partitions[],
                     MPI_Request request) MPICH_API_PUBLIC;
int QMPI_Pready_range(QMPI_Context context, int tool_id, int partition_low, int partition_high,
                      MPI_Request request) MPICH_API_PUBLIC;
int QMPI_Precv_init(QMPI_Context context, int tool_id, void *buf, int partitions, MPI_Count count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_API_PUBLIC;
int QMPI_Psend_init(QMPI_Context context, int tool_id, const void *buf, int partitions,
                    MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                    MPI_Info info, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,6) MPICH_API_PUBLIC;
int QMPI_Bsend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bsend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Bsend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Buffer_attach(QMPI_Context context, int tool_id, void *buffer, int size) MPICH_API_PUBLIC;
int QMPI_Buffer_attach_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count size)
    MPICH_API_PUBLIC;
int QMPI_Buffer_detach(QMPI_Context context, int tool_id, void *buffer_addr, int *size)
    MPICH_API_PUBLIC;
int QMPI_Buffer_detach_c(QMPI_Context context, int tool_id, void *buffer_addr, MPI_Count *size)
    MPICH_API_PUBLIC;
int QMPI_Buffer_flush(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPI_Buffer_iflush(QMPI_Context context, int tool_id, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Comm_attach_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer,
                            int size) MPICH_API_PUBLIC;
int QMPI_Comm_attach_buffer_c(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer,
                              MPI_Count size) MPICH_API_PUBLIC;
int QMPI_Comm_detach_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer_addr,
                            int *size) MPICH_API_PUBLIC;
int QMPI_Comm_detach_buffer_c(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer_addr,
                              MPI_Count *size) MPICH_API_PUBLIC;
int QMPI_Comm_flush_buffer(QMPI_Context context, int tool_id, MPI_Comm comm) MPICH_API_PUBLIC;
int QMPI_Comm_iflush_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Request *request)
    MPICH_API_PUBLIC;
int QMPI_Ibsend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Ibsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Improbe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm, int *flag,
                 MPI_Message *message, MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Imrecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
                MPI_Message *message, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Imrecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                  MPI_Datatype datatype, MPI_Message *message, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Iprobe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm, int *flag,
                MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Irecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Irecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                 MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Irsend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Irsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Isend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm, MPI_Request *request)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Isend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Isendrecv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                   MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(8,10) MPICH_API_PUBLIC;
int QMPI_Isendrecv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                     MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(8,10) MPICH_API_PUBLIC;
int QMPI_Isendrecv_replace(QMPI_Context context, int tool_id, void *buf, int count,
                           MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                           MPI_Comm comm, MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Isendrecv_replace_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                             MPI_Comm comm, MPI_Request *request)
                             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Issend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Issend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Mprobe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
                MPI_Message *message, MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Mrecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
               MPI_Message *message, MPI_Status *status)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Mrecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                 MPI_Datatype datatype, MPI_Message *message, MPI_Status *status)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Probe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
               MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Recv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
              int source, int tag, MPI_Comm comm, MPI_Status *status)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Recv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Recv_init(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Recv_init_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                     MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                     MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rsend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rsend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rsend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Send(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
              int dest, int tag, MPI_Comm comm)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Send_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Send_init(QMPI_Context context, int tool_id, const void *buf, int count,
                   MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                   MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Send_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                     MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                     MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Sendrecv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                  MPI_Status *status)
                  MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(8,10) MPICH_API_PUBLIC;
int QMPI_Sendrecv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                    MPI_Comm comm, MPI_Status *status)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(8,10) MPICH_API_PUBLIC;
int QMPI_Sendrecv_replace(QMPI_Context context, int tool_id, void *buf, int count,
                          MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                          MPI_Comm comm, MPI_Status *status)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Sendrecv_replace_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                            MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                            MPI_Comm comm, MPI_Status *status)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Session_attach_buffer(QMPI_Context context, int tool_id, MPI_Session session, void *buffer,
                               int size) MPICH_API_PUBLIC;
int QMPI_Session_attach_buffer_c(QMPI_Context context, int tool_id, MPI_Session session,
                                 void *buffer, MPI_Count size) MPICH_API_PUBLIC;
int QMPI_Session_detach_buffer(QMPI_Context context, int tool_id, MPI_Session session,
                               void *buffer_addr, int *size) MPICH_API_PUBLIC;
int QMPI_Session_detach_buffer_c(QMPI_Context context, int tool_id, MPI_Session session,
                                 void *buffer_addr, MPI_Count *size) MPICH_API_PUBLIC;
int QMPI_Session_flush_buffer(QMPI_Context context, int tool_id, MPI_Session session)
    MPICH_API_PUBLIC;
int QMPI_Session_iflush_buffer(QMPI_Context context, int tool_id, MPI_Session session,
                               MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Ssend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Ssend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                 MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Ssend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Ssend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Cancel(QMPI_Context context, int tool_id, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Grequest_complete(QMPI_Context context, int tool_id, MPI_Request request)
    MPICH_API_PUBLIC;
int QMPI_Grequest_start(QMPI_Context context, int tool_id, MPI_Grequest_query_function *query_fn,
                        MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Request_free(QMPI_Context context, int tool_id, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Request_get_status(QMPI_Context context, int tool_id, MPI_Request request, int *flag,
                            MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Request_get_status_all(QMPI_Context context, int tool_id, int count,
                                MPI_Request array_of_requests[], int *flag,
                                MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int QMPI_Request_get_status_any(QMPI_Context context, int tool_id, int count,
                                MPI_Request array_of_requests[], int *indx, int *flag,
                                MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Request_get_status_some(QMPI_Context context, int tool_id, int incount,
                                 MPI_Request array_of_requests[], int *outcount,
                                 int array_of_indices[], MPI_Status *array_of_statuses)
                                 MPICH_API_PUBLIC;
int QMPI_Start(QMPI_Context context, int tool_id, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Startall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[])
    MPICH_API_PUBLIC;
int QMPI_Status_c2f(QMPI_Context context, int tool_id, const MPI_Status *c_status,
                    MPI_Fint *f_status) MPICH_API_PUBLIC;
int QMPI_Status_f2c(QMPI_Context context, int tool_id, const MPI_Fint *f_status,
                    MPI_Status *c_status) MPICH_API_PUBLIC;
int QMPI_Status_get_error(QMPI_Context context, int tool_id, MPI_Status *status, int *error)
    MPICH_API_PUBLIC;
int QMPI_Status_get_source(QMPI_Context context, int tool_id, MPI_Status *status, int *source)
    MPICH_API_PUBLIC;
int QMPI_Status_get_tag(QMPI_Context context, int tool_id, MPI_Status *status, int *tag)
    MPICH_API_PUBLIC;
int QMPI_Status_set_error(QMPI_Context context, int tool_id, MPI_Status *status, int error)
    MPICH_API_PUBLIC;
int QMPI_Status_set_source(QMPI_Context context, int tool_id, MPI_Status *status, int source)
    MPICH_API_PUBLIC;
int QMPI_Status_set_tag(QMPI_Context context, int tool_id, MPI_Status *status, int tag)
    MPICH_API_PUBLIC;
int QMPI_Status_set_cancelled(QMPI_Context context, int tool_id, MPI_Status *status, int flag)
    MPICH_API_PUBLIC;
int QMPI_Test(QMPI_Context context, int tool_id, MPI_Request *request, int *flag,
              MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Test_cancelled(QMPI_Context context, int tool_id, const MPI_Status *status, int *flag)
    MPICH_API_PUBLIC;
int QMPI_Testall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *flag, MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int QMPI_Testany(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *indx, int *flag, MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Testsome(QMPI_Context context, int tool_id, int incount, MPI_Request array_of_requests[],
                  int *outcount, int array_of_indices[], MPI_Status *array_of_statuses)
                  MPICH_API_PUBLIC;
int QMPI_Wait(QMPI_Context context, int tool_id, MPI_Request *request, MPI_Status *status)
    MPICH_API_PUBLIC;
int QMPI_Waitall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 MPI_Status *array_of_statuses) MPICH_API_PUBLIC;
int QMPI_Waitany(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *indx, MPI_Status *status) MPICH_API_PUBLIC;
int QMPI_Waitsome(QMPI_Context context, int tool_id, int incount, MPI_Request array_of_requests[],
                  int *outcount, int array_of_indices[], MPI_Status *array_of_statuses)
                  MPICH_API_PUBLIC;
int QMPIX_Grequest_start(QMPI_Context context, int tool_id, MPI_Grequest_query_function *query_fn,
                         MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                         void *extra_state, MPI_Request *request) MPICH_API_PUBLIC;
int QMPIX_Grequest_class_create(QMPI_Context context, int tool_id,
                                MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class) MPICH_API_PUBLIC;
int QMPIX_Grequest_class_allocate(QMPI_Context context, int tool_id, MPIX_Grequest_class greq_class,
                                  void *extra_state, MPI_Request *request) MPICH_API_PUBLIC;
int QMPI_Accumulate(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
                    MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                    int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                    MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                      MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
                      MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                      MPI_Op op, MPI_Win win)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Alloc_mem(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Info info, void *baseptr)
    MPICH_API_PUBLIC;
int QMPI_Compare_and_swap(QMPI_Context context, int tool_id, const void *origin_addr,
                          const void *compare_addr, void *result_addr, MPI_Datatype datatype,
                          int target_rank, MPI_Aint target_disp, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Fetch_and_op(QMPI_Context context, int tool_id, const void *origin_addr, void *result_addr,
                      MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Op op,
                      MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Free_mem(QMPI_Context context, int tool_id, void *base) MPICH_API_PUBLIC;
int QMPI_Get(QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)
             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Get_c(QMPI_Context context, int tool_id, void *origin_addr, MPI_Count origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Get_accumulate(QMPI_Context context, int tool_id, const void *origin_addr,
                        int origin_count, MPI_Datatype origin_datatype, void *result_addr,
                        int result_count, MPI_Datatype result_datatype, int target_rank,
                        MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
                        MPI_Op op, MPI_Win win)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Get_accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                          MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
                          MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
                          MPI_Aint target_disp, MPI_Count target_count,
                          MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Put(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)
             MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Put_c(QMPI_Context context, int tool_id, const void *origin_addr, MPI_Count origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
               MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Raccumulate(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
                     MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                     int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                     MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Raccumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                       MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
                       MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                       MPI_Op op, MPI_Win win, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rget(QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
              MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rget_c(QMPI_Context context, int tool_id, void *origin_addr, MPI_Count origin_count,
                MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rget_accumulate(QMPI_Context context, int tool_id, const void *origin_addr,
                         int origin_count, MPI_Datatype origin_datatype, void *result_addr,
                         int result_count, MPI_Datatype result_datatype, int target_rank,
                         MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
                         MPI_Op op, MPI_Win win, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Rget_accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                           MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
                           MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
                           MPI_Aint target_disp, MPI_Count target_count,
                           MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                           MPI_Request *request)
                           MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8) MPICH_API_PUBLIC;
int QMPI_Rput(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
              MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
              MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Rput_c(QMPI_Context context, int tool_id, const void *origin_addr, MPI_Count origin_count,
                MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win,
                MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPI_Win_allocate(QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
                      MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win) MPICH_API_PUBLIC;
int QMPI_Win_allocate_c(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Aint disp_unit,
                        MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                        MPICH_API_PUBLIC;
int QMPI_Win_allocate_shared(QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
                             MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                             MPICH_API_PUBLIC;
int QMPI_Win_allocate_shared_c(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Aint disp_unit,
                               MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                               MPICH_API_PUBLIC;
int QMPI_Win_attach(QMPI_Context context, int tool_id, MPI_Win win, void *base, MPI_Aint size)
    MPICH_API_PUBLIC;
int QMPI_Win_complete(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_create(QMPI_Context context, int tool_id, void *base, MPI_Aint size, int disp_unit,
                    MPI_Info info, MPI_Comm comm, MPI_Win *win) MPICH_API_PUBLIC;
int QMPI_Win_create_c(QMPI_Context context, int tool_id, void *base, MPI_Aint size,
                      MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win)
                      MPICH_API_PUBLIC;
int QMPI_Win_create_dynamic(QMPI_Context context, int tool_id, MPI_Info info, MPI_Comm comm,
                            MPI_Win *win) MPICH_API_PUBLIC;
int QMPI_Win_detach(QMPI_Context context, int tool_id, MPI_Win win, const void *base)
    MPICH_API_PUBLIC;
int QMPI_Win_fence(QMPI_Context context, int tool_id, int assert, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_flush(QMPI_Context context, int tool_id, int rank, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_flush_all(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_flush_local(QMPI_Context context, int tool_id, int rank, MPI_Win win)
    MPICH_API_PUBLIC;
int QMPI_Win_flush_local_all(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_free(QMPI_Context context, int tool_id, MPI_Win *win) MPICH_API_PUBLIC;
int QMPI_Win_get_group(QMPI_Context context, int tool_id, MPI_Win win, MPI_Group *group)
    MPICH_API_PUBLIC;
int QMPI_Win_get_info(QMPI_Context context, int tool_id, MPI_Win win, MPI_Info *info_used)
    MPICH_API_PUBLIC;
int QMPI_Win_get_name(QMPI_Context context, int tool_id, MPI_Win win, char *win_name,
                      int *resultlen) MPICH_API_PUBLIC;
int QMPI_Win_lock(QMPI_Context context, int tool_id, int lock_type, int rank, int assert,
                  MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_lock_all(QMPI_Context context, int tool_id, int assert, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_post(QMPI_Context context, int tool_id, MPI_Group group, int assert, MPI_Win win)
    MPICH_API_PUBLIC;
int QMPI_Win_set_info(QMPI_Context context, int tool_id, MPI_Win win, MPI_Info info)
    MPICH_API_PUBLIC;
int QMPI_Win_set_name(QMPI_Context context, int tool_id, MPI_Win win, const char *win_name)
    MPICH_API_PUBLIC;
int QMPI_Win_shared_query(QMPI_Context context, int tool_id, MPI_Win win, int rank, MPI_Aint *size,
                          int *disp_unit, void *baseptr) MPICH_API_PUBLIC;
int QMPI_Win_shared_query_c(QMPI_Context context, int tool_id, MPI_Win win, int rank,
                            MPI_Aint *size, MPI_Aint *disp_unit, void *baseptr) MPICH_API_PUBLIC;
int QMPI_Win_start(QMPI_Context context, int tool_id, MPI_Group group, int assert, MPI_Win win)
    MPICH_API_PUBLIC;
int QMPI_Win_sync(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_test(QMPI_Context context, int tool_id, MPI_Win win, int *flag) MPICH_API_PUBLIC;
int QMPI_Win_unlock(QMPI_Context context, int tool_id, int rank, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_unlock_all(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Win_wait(QMPI_Context context, int tool_id, MPI_Win win) MPICH_API_PUBLIC;
int QMPI_Close_port(QMPI_Context context, int tool_id, const char *port_name) MPICH_API_PUBLIC;
int QMPI_Comm_accept(QMPI_Context context, int tool_id, const char *port_name, MPI_Info info,
                     int root, MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_connect(QMPI_Context context, int tool_id, const char *port_name, MPI_Info info,
                      int root, MPI_Comm comm, MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Comm_disconnect(QMPI_Context context, int tool_id, MPI_Comm *comm) MPICH_API_PUBLIC;
int QMPI_Comm_get_parent(QMPI_Context context, int tool_id, MPI_Comm *parent) MPICH_API_PUBLIC;
int QMPI_Comm_join(QMPI_Context context, int tool_id, int fd, MPI_Comm *intercomm)
    MPICH_API_PUBLIC;
int QMPI_Comm_spawn(QMPI_Context context, int tool_id, const char *command, char *argv[],
                    int maxprocs, MPI_Info info, int root, MPI_Comm comm, MPI_Comm *intercomm,
                    int array_of_errcodes[]) MPICH_API_PUBLIC;
int QMPI_Comm_spawn_multiple(QMPI_Context context, int tool_id, int count,
                             char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                             MPICH_API_PUBLIC;
int QMPI_Lookup_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                     char *port_name) MPICH_API_PUBLIC;
int QMPI_Open_port(QMPI_Context context, int tool_id, MPI_Info info, char *port_name)
    MPICH_API_PUBLIC;
int QMPI_Publish_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                      const char *port_name) MPICH_API_PUBLIC;
int QMPI_Unpublish_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                        const char *port_name) MPICH_API_PUBLIC;
int QMPIX_Stream_create(QMPI_Context context, int tool_id, MPI_Info info, MPIX_Stream *stream)
    MPICH_API_PUBLIC;
int QMPIX_Stream_free(QMPI_Context context, int tool_id, MPIX_Stream *stream) MPICH_API_PUBLIC;
int QMPIX_Stream_comm_create(QMPI_Context context, int tool_id, MPI_Comm comm, MPIX_Stream stream,
                             MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPIX_Stream_comm_create_multiplex(QMPI_Context context, int tool_id, MPI_Comm comm, int count,
                                       MPIX_Stream array_of_streams[], MPI_Comm *newcomm)
                                       MPICH_API_PUBLIC;
int QMPIX_Comm_get_stream(QMPI_Context context, int tool_id, MPI_Comm comm, int idx,
                          MPIX_Stream *stream) MPICH_API_PUBLIC;
int QMPIX_Stream_progress(QMPI_Context context, int tool_id, MPIX_Stream stream) MPICH_API_PUBLIC;
int QMPIX_Start_progress_thread(QMPI_Context context, int tool_id, MPIX_Stream stream)
    MPICH_API_PUBLIC;
int QMPIX_Stop_progress_thread(QMPI_Context context, int tool_id, MPIX_Stream stream)
    MPICH_API_PUBLIC;
int QMPIX_Stream_send(QMPI_Context context, int tool_id, const void *buf, int count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      int source_stream_index, int dest_stream_index)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_send_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                        MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                        int source_stream_index, int dest_stream_index)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_isend(QMPI_Context context, int tool_id, const void *buf, int count,
                       MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                       int source_stream_index, int dest_stream_index, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_isend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                         int source_stream_index, int dest_stream_index, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_recv(QMPI_Context context, int tool_id, void *buf, int count,
                      MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                      int source_stream_index, int dest_stream_index, MPI_Status *status)
                      MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_recv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                        MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                        int source_stream_index, int dest_stream_index, MPI_Status *status)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_irecv(QMPI_Context context, int tool_id, void *buf, int count,
                       MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                       int source_stream_index, int dest_stream_index, MPI_Request *request)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Stream_irecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                         MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                         int source_stream_index, int dest_stream_index, MPI_Request *request)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Send_enqueue(QMPI_Context context, int tool_id, const void *buf, int count,
                       MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                       MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Send_enqueue_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Recv_enqueue(QMPI_Context context, int tool_id, void *buf, int count,
                       MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                       MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Recv_enqueue_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                         MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                         MPI_Status *status)
                         MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Isend_enqueue(QMPI_Context context, int tool_id, const void *buf, int count,
                        MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Isend_enqueue_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                          MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Irecv_enqueue(QMPI_Context context, int tool_id, void *buf, int count,
                        MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                        MPI_Request *request)
                        MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Irecv_enqueue_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                          MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                          MPI_Request *request)
                          MPICH_ATTR_POINTER_WITH_TYPE_TAG(3,5) MPICH_API_PUBLIC;
int QMPIX_Wait_enqueue(QMPI_Context context, int tool_id, MPI_Request *request, MPI_Status *status)
    MPICH_API_PUBLIC;
int QMPIX_Waitall_enqueue(QMPI_Context context, int tool_id, int count,
                          MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
                          MPICH_API_PUBLIC;
int QMPIX_Allreduce_enqueue(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                            int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                            MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPIX_Allreduce_enqueue_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                              MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                              MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6) MPICH_API_PUBLIC;
int QMPIX_Threadcomm_init(QMPI_Context context, int tool_id, MPI_Comm comm, int num_threads,
                          MPI_Comm *newthreadcomm) MPICH_API_PUBLIC;
int QMPIX_Threadcomm_free(QMPI_Context context, int tool_id, MPI_Comm *threadcomm)
    MPICH_API_PUBLIC;
int QMPIX_Threadcomm_start(QMPI_Context context, int tool_id, MPI_Comm threadcomm)
    MPICH_API_PUBLIC;
int QMPIX_Threadcomm_finish(QMPI_Context context, int tool_id, MPI_Comm threadcomm)
    MPICH_API_PUBLIC;
double QMPI_Wtick(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
double QMPI_Wtime(QMPI_Context context, int tool_id) MPICH_API_PUBLIC;
int QMPI_Cart_coords(QMPI_Context context, int tool_id, MPI_Comm comm, int rank, int maxdims,
                     int coords[]) MPICH_API_PUBLIC;
int QMPI_Cart_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int ndims,
                     const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart)
                     MPICH_API_PUBLIC;
int QMPI_Cart_get(QMPI_Context context, int tool_id, MPI_Comm comm, int maxdims, int dims[],
                  int periods[], int coords[]) MPICH_API_PUBLIC;
int QMPI_Cart_map(QMPI_Context context, int tool_id, MPI_Comm comm, int ndims, const int dims[],
                  const int periods[], int *newrank) MPICH_API_PUBLIC;
int QMPI_Cart_rank(QMPI_Context context, int tool_id, MPI_Comm comm, const int coords[], int *rank)
    MPICH_API_PUBLIC;
int QMPI_Cart_shift(QMPI_Context context, int tool_id, MPI_Comm comm, int direction, int disp,
                    int *rank_source, int *rank_dest) MPICH_API_PUBLIC;
int QMPI_Cart_sub(QMPI_Context context, int tool_id, MPI_Comm comm, const int remain_dims[],
                  MPI_Comm *newcomm) MPICH_API_PUBLIC;
int QMPI_Cartdim_get(QMPI_Context context, int tool_id, MPI_Comm comm, int *ndims)
    MPICH_API_PUBLIC;
int QMPI_Dims_create(QMPI_Context context, int tool_id, int nnodes, int ndims, int dims[])
    MPICH_API_PUBLIC;
int QMPI_Dist_graph_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int n,
                           const int sources[], const int degrees[], const int destinations[],
                           const int weights[], MPI_Info info, int reorder,
                           MPI_Comm *comm_dist_graph) MPICH_API_PUBLIC;
int QMPI_Dist_graph_create_adjacent(QMPI_Context context, int tool_id, MPI_Comm comm_old,
                                    int indegree, const int sources[], const int sourceweights[],
                                    int outdegree, const int destinations[],
                                    const int destweights[], MPI_Info info, int reorder,
                                    MPI_Comm *comm_dist_graph) MPICH_API_PUBLIC;
int QMPI_Dist_graph_neighbors(QMPI_Context context, int tool_id, MPI_Comm comm, int maxindegree,
                              int sources[], int sourceweights[], int maxoutdegree,
                              int destinations[], int destweights[]) MPICH_API_PUBLIC;
int QMPI_Dist_graph_neighbors_count(QMPI_Context context, int tool_id, MPI_Comm comm, int *indegree,
                                    int *outdegree, int *weighted) MPICH_API_PUBLIC;
int QMPI_Get_hw_resource_info(QMPI_Context context, int tool_id, MPI_Info *hw_info)
    MPICH_API_PUBLIC;
int QMPI_Graph_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int nnodes,
                      const int indx[], const int edges[], int reorder, MPI_Comm *comm_graph)
                      MPICH_API_PUBLIC;
int QMPI_Graph_get(QMPI_Context context, int tool_id, MPI_Comm comm, int maxindex, int maxedges,
                   int indx[], int edges[]) MPICH_API_PUBLIC;
int QMPI_Graph_map(QMPI_Context context, int tool_id, MPI_Comm comm, int nnodes, const int indx[],
                   const int edges[], int *newrank) MPICH_API_PUBLIC;
int QMPI_Graph_neighbors(QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
                         int maxneighbors, int neighbors[]) MPICH_API_PUBLIC;
int QMPI_Graph_neighbors_count(QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
                               int *nneighbors) MPICH_API_PUBLIC;
int QMPI_Graphdims_get(QMPI_Context context, int tool_id, MPI_Comm comm, int *nnodes, int *nedges)
    MPICH_API_PUBLIC;
int QMPI_Topo_test(QMPI_Context context, int tool_id, MPI_Comm comm, int *status) MPICH_API_PUBLIC;

typedef int (QMPI_Attr_delete_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval);
typedef int (QMPI_Attr_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
             void *attribute_val, int *flag);
typedef int (QMPI_Attr_put_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
             void *attribute_val);
typedef int (QMPI_Comm_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Comm_copy_attr_function *comm_copy_attr_fn,
             MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
             void *extra_state);
typedef int (QMPI_Comm_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval);
typedef int (QMPI_Comm_free_keyval_t) (QMPI_Context context, int tool_id, int *comm_keyval);
typedef int (QMPI_Comm_get_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval, void *attribute_val, int *flag);
typedef int (QMPI_Comm_set_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval, void *attribute_val);
typedef int (QMPI_Keyval_create_t) (QMPI_Context context, int tool_id, MPI_Copy_function *copy_fn,
             MPI_Delete_function *delete_fn, int *keyval, void *extra_state);
typedef int (QMPI_Keyval_free_t) (QMPI_Context context, int tool_id, int *keyval);
typedef int (QMPI_Type_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Type_copy_attr_function *type_copy_attr_fn,
             MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
             void *extra_state);
typedef int (QMPI_Type_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval);
typedef int (QMPI_Type_free_keyval_t) (QMPI_Context context, int tool_id, int *type_keyval);
typedef int (QMPI_Type_get_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval, void *attribute_val, int *flag);
typedef int (QMPI_Type_set_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval, void *attribute_val);
typedef int (QMPI_Win_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Win_copy_attr_function *win_copy_attr_fn,
             MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval, void *extra_state);
typedef int (QMPI_Win_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Win win,
             int win_keyval);
typedef int (QMPI_Win_free_keyval_t) (QMPI_Context context, int tool_id, int *win_keyval);
typedef int (QMPI_Win_get_attr_t) (QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
             void *attribute_val, int *flag);
typedef int (QMPI_Win_set_attr_t) (QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
             void *attribute_val);
typedef int (QMPI_Allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgather_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allgather_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm);
typedef int (QMPI_Allgatherv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Allgatherv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allreduce_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Allreduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Allreduce_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allreduce_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoall_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoall_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoallv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Alltoallw_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallw_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Barrier_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPI_Barrier_init_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Bcast_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm);
typedef int (QMPI_Bcast_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm);
typedef int (QMPI_Bcast_init_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Bcast_init_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Exscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Exscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Exscan_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Exscan_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gather_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Gather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Gather_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gather_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Gatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Gatherv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Gatherv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Iallgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallreduce_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Iallreduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ialltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ialltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibarrier_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ibcast_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibcast_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iexscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iexscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Igather_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Igather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Igatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Igatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_block_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_block_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
             MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Iscatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatterv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatterv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgather_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgather_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm);
typedef int (QMPI_Neighbor_allgatherv_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_allgatherv_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoall_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoall_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallv_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const int sendcounts[], const int sdispls[],
             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallv_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
             MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
             const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallw_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallw_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
             const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Reduce_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
typedef int (QMPI_Reduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm);
typedef int (QMPI_Reduce_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_local_t) (QMPI_Context context, int tool_id, const void *inbuf,
             void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op);
typedef int (QMPI_Reduce_local_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             void *inoutbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op);
typedef int (QMPI_Reduce_scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_block_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
             MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Scan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Scan_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scan_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatter_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatter_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatterv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatterv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatterv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Scatterv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Comm_compare_t) (QMPI_Context context, int tool_id, MPI_Comm comm1,
             MPI_Comm comm2, int *result);
typedef int (QMPI_Comm_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
             MPI_Comm *newcomm);
typedef int (QMPI_Comm_create_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group group, int tag, MPI_Comm *newcomm);
typedef int (QMPI_Comm_dup_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Comm *newcomm);
typedef int (QMPI_Comm_dup_with_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info, MPI_Comm *newcomm);
typedef int (QMPI_Comm_free_t) (QMPI_Context context, int tool_id, MPI_Comm *comm);
typedef int (QMPI_Comm_get_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info *info_used);
typedef int (QMPI_Comm_get_name_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             char *comm_name, int *resultlen);
typedef int (QMPI_Comm_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *group);
typedef int (QMPI_Comm_idup_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm,
             MPI_Request *request);
typedef int (QMPI_Comm_idup_with_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info, MPI_Comm *newcomm, MPI_Request *request);
typedef int (QMPI_Comm_rank_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *rank);
typedef int (QMPI_Comm_remote_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *group);
typedef int (QMPI_Comm_remote_size_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *size);
typedef int (QMPI_Comm_set_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info);
typedef int (QMPI_Comm_set_name_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const char *comm_name);
typedef int (QMPI_Comm_size_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *size);
typedef int (QMPI_Comm_split_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int color,
             int key, MPI_Comm *newcomm);
typedef int (QMPI_Comm_split_type_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
typedef int (QMPI_Comm_test_inter_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *flag);
typedef int (QMPI_Intercomm_create_t) (QMPI_Context context, int tool_id, MPI_Comm local_comm,
             int local_leader, MPI_Comm peer_comm, int remote_leader, int tag,
             MPI_Comm *newintercomm);
typedef int (QMPI_Intercomm_create_from_groups_t) (QMPI_Context context, int tool_id,
             MPI_Group local_group, int local_leader, MPI_Group remote_group, int remote_leader,
             const char *stringtag, MPI_Info info, MPI_Errhandler errhandler,
             MPI_Comm *newintercomm);
typedef int (QMPI_Intercomm_merge_t) (QMPI_Context context, int tool_id, MPI_Comm intercomm,
             int high, MPI_Comm *newintracomm);
typedef int (QMPIX_Comm_test_threadcomm_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *flag);
typedef int (QMPIX_Comm_revoke_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPIX_Comm_shrink_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Comm *newcomm);
typedef int (QMPIX_Comm_failure_ack_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPIX_Comm_failure_get_acked_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *failedgrp);
typedef int (QMPIX_Comm_agree_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *flag);
typedef int (QMPIX_Comm_get_failed_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *failedgrp);
typedef int (QMPI_Get_address_t) (QMPI_Context context, int tool_id, const void *location,
             MPI_Aint *address);
typedef int (QMPI_Get_count_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, int *count);
typedef int (QMPI_Get_count_c_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Get_elements_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, int *count);
typedef int (QMPI_Get_elements_c_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Get_elements_x_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Pack_t) (QMPI_Context context, int tool_id, const void *inbuf, int incount,
             MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm);
typedef int (QMPI_Pack_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             MPI_Count incount, MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
             MPI_Count *position, MPI_Comm comm);
typedef int (QMPI_Pack_external_t) (QMPI_Context context, int tool_id, const char *datarep,
             const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, MPI_Aint outsize,
             MPI_Aint *position);
typedef int (QMPI_Pack_external_c_t) (QMPI_Context context, int tool_id, const char *datarep,
             const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
             MPI_Count outsize, MPI_Count *position);
typedef int (QMPI_Pack_external_size_t) (QMPI_Context context, int tool_id, const char *datarep,
             int incount, MPI_Datatype datatype, MPI_Aint *size);
typedef int (QMPI_Pack_external_size_c_t) (QMPI_Context context, int tool_id, const char *datarep,
             MPI_Count incount, MPI_Datatype datatype, MPI_Count *size);
typedef int (QMPI_Pack_size_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Datatype datatype, MPI_Comm comm, int *size);
typedef int (QMPI_Pack_size_c_t) (QMPI_Context context, int tool_id, MPI_Count incount,
             MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size);
typedef int (QMPI_Status_set_elements_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, int count);
typedef int (QMPI_Status_set_elements_c_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, MPI_Count count);
typedef int (QMPI_Status_set_elements_x_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, MPI_Count count);
typedef int (QMPI_Type_commit_t) (QMPI_Context context, int tool_id, MPI_Datatype *datatype);
typedef int (QMPI_Type_contiguous_t) (QMPI_Context context, int tool_id, int count,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_contiguous_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_darray_t) (QMPI_Context context, int tool_id, int size, int rank,
             int ndims, const int array_of_gsizes[], const int array_of_distribs[],
             const int array_of_dargs[], const int array_of_psizes[], int order,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_darray_c_t) (QMPI_Context context, int tool_id, int size, int rank,
             int ndims, const MPI_Count array_of_gsizes[], const int array_of_distribs[],
             const int array_of_dargs[], const int array_of_psizes[], int order,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_block_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_block_c_t) (QMPI_Context context, int tool_id,
             MPI_Count count, MPI_Count blocklength, const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hvector_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hvector_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_indexed_block_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, const int array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_indexed_block_c_t) (QMPI_Context context, int tool_id,
             MPI_Count count, MPI_Count blocklength, const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_resized_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Aint lb, MPI_Aint extent, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_resized_c_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Count lb, MPI_Count extent, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_struct_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
             const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPI_Type_create_struct_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPI_Type_create_subarray_t) (QMPI_Context context, int tool_id, int ndims,
             const int array_of_sizes[], const int array_of_subsizes[], const int array_of_starts[],
             int order, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_subarray_c_t) (QMPI_Context context, int tool_id, int ndims,
             const MPI_Count array_of_sizes[], const MPI_Count array_of_subsizes[],
             const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_dup_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_free_t) (QMPI_Context context, int tool_id, MPI_Datatype *datatype);
typedef int (QMPI_Type_get_contents_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int max_integers, int max_addresses, int max_datatypes, int array_of_integers[],
             MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]);
typedef int (QMPI_Type_get_contents_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count max_integers, MPI_Count max_addresses, MPI_Count max_large_counts,
             MPI_Count max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[]);
typedef int (QMPI_Type_get_envelope_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int *num_integers, int *num_addresses, int *num_datatypes, int *combiner);
typedef int (QMPI_Type_get_envelope_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *num_integers, MPI_Count *num_addresses, MPI_Count *num_large_counts,
             MPI_Count *num_datatypes, int *combiner);
typedef int (QMPI_Type_get_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *lb, MPI_Aint *extent);
typedef int (QMPI_Type_get_extent_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *lb, MPI_Count *extent);
typedef int (QMPI_Type_get_extent_x_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *lb, MPI_Count *extent);
typedef int (QMPI_Type_get_name_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             char *type_name, int *resultlen);
typedef int (QMPI_Type_get_true_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *true_lb, MPI_Aint *true_extent);
typedef int (QMPI_Type_get_true_extent_c_t) (QMPI_Context context, int tool_id,
             MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent);
typedef int (QMPI_Type_get_true_extent_x_t) (QMPI_Context context, int tool_id,
             MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent);
typedef int (QMPI_Type_get_value_index_t) (QMPI_Context context, int tool_id,
             MPI_Datatype value_type, MPI_Datatype index_type, MPI_Datatype *pair_type);
typedef int (QMPI_Type_indexed_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const int array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_indexed_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_match_size_t) (QMPI_Context context, int tool_id, int typeclass, int size,
             MPI_Datatype *datatype);
typedef int (QMPI_Type_set_name_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             const char *type_name);
typedef int (QMPI_Type_size_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int *size);
typedef int (QMPI_Type_size_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *size);
typedef int (QMPI_Type_size_x_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *size);
typedef int (QMPI_Type_vector_t) (QMPI_Context context, int tool_id, int count, int blocklength,
             int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_vector_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Unpack_t) (QMPI_Context context, int tool_id, const void *inbuf, int insize,
             int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm);
typedef int (QMPI_Unpack_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             MPI_Count insize, MPI_Count *position, void *outbuf, MPI_Count outcount,
             MPI_Datatype datatype, MPI_Comm comm);
typedef int (QMPI_Unpack_external_t) (QMPI_Context context, int tool_id, const char datarep[],
             const void *inbuf, MPI_Aint insize, MPI_Aint *position, void *outbuf, int outcount,
             MPI_Datatype datatype);
typedef int (QMPI_Unpack_external_c_t) (QMPI_Context context, int tool_id, const char datarep[],
             const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
             MPI_Count outcount, MPI_Datatype datatype);
typedef int (QMPI_Address_t) (QMPI_Context context, int tool_id, void *location,
             MPI_Aint *address);
typedef int (QMPI_Type_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *extent);
typedef int (QMPI_Type_lb_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *displacement);
typedef int (QMPI_Type_ub_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *displacement);
typedef int (QMPI_Type_hindexed_t) (QMPI_Context context, int tool_id, int count,
             int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_hvector_t) (QMPI_Context context, int tool_id, int count, int blocklength,
             MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_struct_t) (QMPI_Context context, int tool_id, int count,
             int array_of_blocklengths[], MPI_Aint array_of_displacements[],
             MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPIX_Type_iov_len_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count max_iov_bytes, MPI_Count *iov_len, MPI_Count *actual_iov_bytes);
typedef int (QMPIX_Type_iov_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count iov_offset, MPIX_Iov *iov, MPI_Count max_iov_len,
             MPI_Count *actual_iov_len);
typedef int (QMPI_Add_error_class_t) (QMPI_Context context, int tool_id, int *errorclass);
typedef int (QMPI_Add_error_code_t) (QMPI_Context context, int tool_id, int errorclass,
             int *errorcode);
typedef int (QMPI_Add_error_string_t) (QMPI_Context context, int tool_id, int errorcode,
             const char *string);
typedef int (QMPI_Comm_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int errorcode);
typedef int (QMPI_Comm_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Comm_errhandler_function *comm_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Comm_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Comm_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler errhandler);
typedef int (QMPI_Errhandler_free_t) (QMPI_Context context, int tool_id,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Error_class_t) (QMPI_Context context, int tool_id, int errorcode,
             int *errorclass);
typedef int (QMPI_Error_string_t) (QMPI_Context context, int tool_id, int errorcode, char *string,
             int *resultlen);
typedef int (QMPI_File_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_File fh,
             int errorcode);
typedef int (QMPI_File_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_File_errhandler_function *file_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_File_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_File file,
             MPI_Errhandler *errhandler);
typedef int (QMPI_File_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_File file,
             MPI_Errhandler errhandler);
typedef int (QMPI_Remove_error_class_t) (QMPI_Context context, int tool_id, int errorclass);
typedef int (QMPI_Remove_error_code_t) (QMPI_Context context, int tool_id, int errorcode);
typedef int (QMPI_Remove_error_string_t) (QMPI_Context context, int tool_id, int errorcode);
typedef int (QMPI_Session_call_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Session session, int errorcode);
typedef int (QMPI_Session_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Session_errhandler_function *session_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Session_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Session_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Errhandler errhandler);
typedef int (QMPI_Win_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             int errorcode);
typedef int (QMPI_Win_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Win_errhandler_function *win_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Win_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Win_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Errhandler errhandler);
typedef int (QMPI_Errhandler_create_t) (QMPI_Context context, int tool_id,
             MPI_Comm_errhandler_function *comm_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Errhandler_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Errhandler_set_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler errhandler);
typedef MPI_Fint (QMPI_Comm_c2f_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef MPI_Comm (QMPI_Comm_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint comm);
typedef MPI_Fint (QMPI_Errhandler_c2f_t) (QMPI_Context context, int tool_id,
                  MPI_Errhandler errhandler);
typedef MPI_Errhandler (QMPI_Errhandler_f2c_t) (QMPI_Context context, int tool_id,
                        MPI_Fint errhandler);
typedef MPI_Fint (QMPI_Group_c2f_t) (QMPI_Context context, int tool_id, MPI_Group group);
typedef MPI_Group (QMPI_Group_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint group);
typedef MPI_Fint (QMPI_Info_c2f_t) (QMPI_Context context, int tool_id, MPI_Info info);
typedef MPI_Info (QMPI_Info_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint info);
typedef MPI_Fint (QMPI_Message_c2f_t) (QMPI_Context context, int tool_id, MPI_Message message);
typedef MPI_Message (QMPI_Message_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint message);
typedef MPI_Fint (QMPI_Op_c2f_t) (QMPI_Context context, int tool_id, MPI_Op op);
typedef MPI_Op (QMPI_Op_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint op);
typedef MPI_Fint (QMPI_Request_c2f_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef MPI_Request (QMPI_Request_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint request);
typedef MPI_Fint (QMPI_Session_c2f_t) (QMPI_Context context, int tool_id, MPI_Session session);
typedef MPI_Session (QMPI_Session_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint session);
typedef MPI_Fint (QMPI_Type_c2f_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype);
typedef MPI_Datatype (QMPI_Type_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint datatype);
typedef MPI_Fint (QMPI_Win_c2f_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef MPI_Win (QMPI_Win_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint win);
typedef int (QMPI_Group_compare_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, int *result);
typedef int (QMPI_Group_difference_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Group_excl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             const int ranks[], MPI_Group *newgroup);
typedef int (QMPI_Group_free_t) (QMPI_Context context, int tool_id, MPI_Group *group);
typedef int (QMPI_Group_incl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             const int ranks[], MPI_Group *newgroup);
typedef int (QMPI_Group_intersection_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Group_range_excl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             int ranges[][3], MPI_Group *newgroup);
typedef int (QMPI_Group_range_incl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             int ranges[][3], MPI_Group *newgroup);
typedef int (QMPI_Group_rank_t) (QMPI_Context context, int tool_id, MPI_Group group, int *rank);
typedef int (QMPI_Group_size_t) (QMPI_Context context, int tool_id, MPI_Group group, int *size);
typedef int (QMPI_Group_translate_ranks_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             int n, const int ranks1[], MPI_Group group2, int ranks2[]);
typedef int (QMPI_Group_union_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Info_create_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Info_create_env_t) (QMPI_Context context, int tool_id, int argc, char *argv[],
             MPI_Info *info);
typedef int (QMPI_Info_delete_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key);
typedef int (QMPI_Info_dup_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Info *newinfo);
typedef int (QMPI_Info_free_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Info_get_t) (QMPI_Context context, int tool_id, MPI_Info info, const char *key,
             int valuelen, char *value, int *flag);
typedef int (QMPI_Info_get_nkeys_t) (QMPI_Context context, int tool_id, MPI_Info info, int *nkeys);
typedef int (QMPI_Info_get_nthkey_t) (QMPI_Context context, int tool_id, MPI_Info info, int n,
             char *key);
typedef int (QMPI_Info_get_string_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, int *buflen, char *value, int *flag);
typedef int (QMPI_Info_get_valuelen_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, int *valuelen, int *flag);
typedef int (QMPI_Info_set_t) (QMPI_Context context, int tool_id, MPI_Info info, const char *key,
             const char *value);
typedef int (QMPIX_Info_set_hex_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, const void *value, int value_size);
typedef int (QMPI_Abort_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode);
typedef int (QMPI_Comm_create_from_group_t) (QMPI_Context context, int tool_id, MPI_Group group,
             const char *stringtag, MPI_Info info, MPI_Errhandler errhandler, MPI_Comm *newcomm);
typedef int (QMPI_Finalize_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Finalized_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Group_from_session_pset_t) (QMPI_Context context, int tool_id,
             MPI_Session session, const char *pset_name, MPI_Group *newgroup);
typedef int (QMPI_Init_t) (QMPI_Context context, int tool_id, int *argc, char ***argv);
typedef int (QMPI_Init_thread_t) (QMPI_Context context, int tool_id, int *argc, char ***argv,
             int required, int *provided);
typedef int (QMPI_Initialized_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Is_thread_main_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Query_thread_t) (QMPI_Context context, int tool_id, int *provided);
typedef int (QMPI_Session_finalize_t) (QMPI_Context context, int tool_id, MPI_Session *session);
typedef int (QMPI_Session_get_info_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info *info_used);
typedef int (QMPI_Session_get_nth_pset_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info info, int n, int *pset_len, char *pset_name);
typedef int (QMPI_Session_get_num_psets_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info info, int *npset_names);
typedef int (QMPI_Session_get_pset_info_t) (QMPI_Context context, int tool_id, MPI_Session session,
             const char *pset_name, MPI_Info *info);
typedef int (QMPI_Session_init_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Errhandler errhandler, MPI_Session *session);
typedef MPI_Aint (QMPI_Aint_add_t) (QMPI_Context context, int tool_id, MPI_Aint base,
                  MPI_Aint disp);
typedef MPI_Aint (QMPI_Aint_diff_t) (QMPI_Context context, int tool_id, MPI_Aint addr1,
                  MPI_Aint addr2);
typedef int (QMPI_Get_library_version_t) (QMPI_Context context, int tool_id, char *version,
             int *resultlen);
typedef int (QMPI_Get_processor_name_t) (QMPI_Context context, int tool_id, char *name,
             int *resultlen);
typedef int (QMPI_Get_version_t) (QMPI_Context context, int tool_id, int *version,
             int *subversion);
typedef int (QMPI_Pcontrol_t) (QMPI_Context context, int tool_id, const int level, ...);
typedef int (QMPIX_GPU_query_support_t) (QMPI_Context context, int tool_id, int gpu_type,
             int *is_supported);
typedef int (QMPIX_Query_cuda_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPIX_Query_ze_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPIX_Query_hip_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_T_category_changed_t) (QMPI_Context context, int tool_id, int *update_number);
typedef int (QMPI_T_category_get_categories_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_cvars_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_events_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *cat_index);
typedef int (QMPI_T_category_get_info_t) (QMPI_Context context, int tool_id, int cat_index,
             char *name, int *name_len, char *desc, int *desc_len, int *num_cvars, int *num_pvars,
             int *num_categories);
typedef int (QMPI_T_category_get_num_t) (QMPI_Context context, int tool_id, int *num_cat);
typedef int (QMPI_T_category_get_num_events_t) (QMPI_Context context, int tool_id, int cat_index,
             int *num_events);
typedef int (QMPI_T_category_get_pvars_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_cvar_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *cvar_index);
typedef int (QMPI_T_cvar_get_info_t) (QMPI_Context context, int tool_id, int cvar_index, char *name,
             int *name_len, int *verbosity, MPI_Datatype *datatype, MPI_T_enum *enumtype,
             char *desc, int *desc_len, int *bind, int *scope);
typedef int (QMPI_T_cvar_get_num_t) (QMPI_Context context, int tool_id, int *num_cvar);
typedef int (QMPI_T_cvar_handle_alloc_t) (QMPI_Context context, int tool_id, int cvar_index,
             void *obj_handle, MPI_T_cvar_handle *handle, int *count);
typedef int (QMPI_T_cvar_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_cvar_handle *handle);
typedef int (QMPI_T_cvar_read_t) (QMPI_Context context, int tool_id, MPI_T_cvar_handle handle,
             void *buf);
typedef int (QMPI_T_cvar_write_t) (QMPI_Context context, int tool_id, MPI_T_cvar_handle handle,
             const void *buf);
typedef int (QMPI_T_enum_get_info_t) (QMPI_Context context, int tool_id, MPI_T_enum enumtype,
             int *num, char *name, int *name_len);
typedef int (QMPI_T_enum_get_item_t) (QMPI_Context context, int tool_id, MPI_T_enum enumtype,
             int indx, int *value, char *name, int *name_len);
typedef int (QMPI_T_event_callback_get_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety,
             MPI_Info *info_used);
typedef int (QMPI_T_event_callback_set_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety,
             MPI_Info info);
typedef int (QMPI_T_event_copy_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, void *buffer);
typedef int (QMPI_T_event_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *event_index);
typedef int (QMPI_T_event_get_info_t) (QMPI_Context context, int tool_id, int event_index,
             char *name, int *name_len, int *verbosity, MPI_Datatype array_of_datatypes[],
             MPI_Aint array_of_displacements[], int *num_elements, MPI_T_enum *enumtype,
             MPI_Info *info, char *desc, int *desc_len, int *bind);
typedef int (QMPI_T_event_get_num_t) (QMPI_Context context, int tool_id, int *num_events);
typedef int (QMPI_T_event_get_source_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, int *source_index);
typedef int (QMPI_T_event_get_timestamp_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, MPI_Count *event_timestamp);
typedef int (QMPI_T_event_handle_alloc_t) (QMPI_Context context, int tool_id, int event_index,
             void *obj_handle, MPI_Info info, MPI_T_event_registration *event_registration);
typedef int (QMPI_T_event_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, void *user_data,
             MPI_T_event_free_cb_function free_cb_function);
typedef int (QMPI_T_event_handle_get_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_Info *info_used);
typedef int (QMPI_T_event_handle_set_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_Info info);
typedef int (QMPI_T_event_read_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, int element_index, void *buffer);
typedef int (QMPI_T_event_register_callback_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety, MPI_Info info,
             void *user_data, MPI_T_event_cb_function event_cb_function);
typedef int (QMPI_T_event_set_dropped_handler_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration,
             MPI_T_event_dropped_cb_function dropped_cb_function);
typedef int (QMPI_T_finalize_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_T_init_thread_t) (QMPI_Context context, int tool_id, int required,
             int *provided);
typedef int (QMPI_T_pvar_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int var_class, int *pvar_index);
typedef int (QMPI_T_pvar_get_info_t) (QMPI_Context context, int tool_id, int pvar_index, char *name,
             int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype,
             MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly,
             int *continuous, int *atomic);
typedef int (QMPI_T_pvar_get_num_t) (QMPI_Context context, int tool_id, int *num_pvar);
typedef int (QMPI_T_pvar_handle_alloc_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, int pvar_index, void *obj_handle,
             MPI_T_pvar_handle *handle, int *count);
typedef int (QMPI_T_pvar_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, MPI_T_pvar_handle *handle);
typedef int (QMPI_T_pvar_read_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle, void *buf);
typedef int (QMPI_T_pvar_readreset_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
typedef int (QMPI_T_pvar_reset_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_session_create_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session *session);
typedef int (QMPI_T_pvar_session_free_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session *session);
typedef int (QMPI_T_pvar_start_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_stop_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_write_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle, const void *buf);
typedef int (QMPI_T_source_get_info_t) (QMPI_Context context, int tool_id, int source_index,
             char *name, int *name_len, char *desc, int *desc_len, MPI_T_source_order *ordering,
             MPI_Count *ticks_per_second, MPI_Count *max_ticks, MPI_Info *info);
typedef int (QMPI_T_source_get_num_t) (QMPI_Context context, int tool_id, int *num_sources);
typedef int (QMPI_T_source_get_timestamp_t) (QMPI_Context context, int tool_id, int source_index,
             MPI_Count *timestamp);
typedef int (QMPI_Op_commutative_t) (QMPI_Context context, int tool_id, MPI_Op op, int *commute);
typedef int (QMPI_Op_create_t) (QMPI_Context context, int tool_id, MPI_User_function *user_fn,
             int commute, MPI_Op *op);
typedef int (QMPI_Op_create_c_t) (QMPI_Context context, int tool_id, MPI_User_function_c *user_fn,
             int commute, MPI_Op *op);
typedef int (QMPI_Op_free_t) (QMPI_Context context, int tool_id, MPI_Op *op);
typedef int (QMPI_Parrived_t) (QMPI_Context context, int tool_id, MPI_Request request,
             int partition, int *flag);
typedef int (QMPI_Pready_t) (QMPI_Context context, int tool_id, int partition,
             MPI_Request request);
typedef int (QMPI_Pready_list_t) (QMPI_Context context, int tool_id, int length,
             const int array_of_partitions[], MPI_Request request);
typedef int (QMPI_Pready_range_t) (QMPI_Context context, int tool_id, int partition_low,
             int partition_high, MPI_Request request);
typedef int (QMPI_Precv_init_t) (QMPI_Context context, int tool_id, void *buf, int partitions,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Psend_init_t) (QMPI_Context context, int tool_id, const void *buf, int partitions,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Bsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Bsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Bsend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Bsend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Buffer_attach_t) (QMPI_Context context, int tool_id, void *buffer, int size);
typedef int (QMPI_Buffer_attach_c_t) (QMPI_Context context, int tool_id, void *buffer,
             MPI_Count size);
typedef int (QMPI_Buffer_detach_t) (QMPI_Context context, int tool_id, void *buffer_addr,
             int *size);
typedef int (QMPI_Buffer_detach_c_t) (QMPI_Context context, int tool_id, void *buffer_addr,
             MPI_Count *size);
typedef int (QMPI_Buffer_flush_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Buffer_iflush_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Comm_attach_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer, int size);
typedef int (QMPI_Comm_attach_buffer_c_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer, MPI_Count size);
typedef int (QMPI_Comm_detach_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer_addr, int *size);
typedef int (QMPI_Comm_detach_buffer_c_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer_addr, MPI_Count *size);
typedef int (QMPI_Comm_flush_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPI_Comm_iflush_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ibsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Improbe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             int *flag, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Imrecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Request *request);
typedef int (QMPI_Imrecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Request *request);
typedef int (QMPI_Iprobe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             int *flag, MPI_Status *status);
typedef int (QMPI_Irecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isendrecv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_replace_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_replace_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Issend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Issend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Mprobe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Mrecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Mrecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Probe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Recv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Recv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Recv_init_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Recv_init_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Rsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Rsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Rsend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Rsend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Send_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Send_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Send_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Send_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Sendrecv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_replace_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_replace_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
             MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Session_attach_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             void *buffer, int size);
typedef int (QMPI_Session_attach_buffer_c_t) (QMPI_Context context, int tool_id,
             MPI_Session session, void *buffer, MPI_Count size);
typedef int (QMPI_Session_detach_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             void *buffer_addr, int *size);
typedef int (QMPI_Session_detach_buffer_c_t) (QMPI_Context context, int tool_id,
             MPI_Session session, void *buffer_addr, MPI_Count *size);
typedef int (QMPI_Session_flush_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session);
typedef int (QMPI_Session_iflush_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Request *request);
typedef int (QMPI_Ssend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Ssend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Ssend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ssend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Cancel_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Grequest_complete_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef int (QMPI_Grequest_start_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, void *extra_state, MPI_Request *request);
typedef int (QMPI_Request_free_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Request_get_status_t) (QMPI_Context context, int tool_id, MPI_Request request,
             int *flag, MPI_Status *status);
typedef int (QMPI_Request_get_status_all_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *flag, MPI_Status *array_of_statuses);
typedef int (QMPI_Request_get_status_any_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *indx, int *flag, MPI_Status *status);
typedef int (QMPI_Request_get_status_some_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPI_Start_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Startall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[]);
typedef int (QMPI_Status_c2f_t) (QMPI_Context context, int tool_id, const MPI_Status *c_status,
             MPI_Fint *f_status);
typedef int (QMPI_Status_f2c_t) (QMPI_Context context, int tool_id, const MPI_Fint *f_status,
             MPI_Status *c_status);
typedef int (QMPI_Status_get_error_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int *error);
typedef int (QMPI_Status_get_source_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int *source);
typedef int (QMPI_Status_get_tag_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int *tag);
typedef int (QMPI_Status_set_error_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int error);
typedef int (QMPI_Status_set_source_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int source);
typedef int (QMPI_Status_set_tag_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int tag);
typedef int (QMPI_Status_set_cancelled_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int flag);
typedef int (QMPI_Test_t) (QMPI_Context context, int tool_id, MPI_Request *request, int *flag,
             MPI_Status *status);
typedef int (QMPI_Test_cancelled_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             int *flag);
typedef int (QMPI_Testall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *flag, MPI_Status *array_of_statuses);
typedef int (QMPI_Testany_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *indx, int *flag, MPI_Status *status);
typedef int (QMPI_Testsome_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPI_Wait_t) (QMPI_Context context, int tool_id, MPI_Request *request,
             MPI_Status *status);
typedef int (QMPI_Waitall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], MPI_Status *array_of_statuses);
typedef int (QMPI_Waitany_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *indx, MPI_Status *status);
typedef int (QMPI_Waitsome_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPIX_Grequest_start_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, MPIX_Grequest_poll_function *poll_fn,
             MPIX_Grequest_wait_function *wait_fn, void *extra_state, MPI_Request *request);
typedef int (QMPIX_Grequest_class_create_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, MPIX_Grequest_poll_function *poll_fn,
             MPIX_Grequest_wait_function *wait_fn, MPIX_Grequest_class *greq_class);
typedef int (QMPIX_Grequest_class_allocate_t) (QMPI_Context context, int tool_id,
             MPIX_Grequest_class greq_class, void *extra_state, MPI_Request *request);
typedef int (QMPI_Accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
typedef int (QMPI_Accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win);
typedef int (QMPI_Alloc_mem_t) (QMPI_Context context, int tool_id, MPI_Aint size, MPI_Info info,
             void *baseptr);
typedef int (QMPI_Compare_and_swap_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             const void *compare_addr, void *result_addr, MPI_Datatype datatype, int target_rank,
             MPI_Aint target_disp, MPI_Win win);
typedef int (QMPI_Fetch_and_op_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             void *result_addr, MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
             MPI_Op op, MPI_Win win);
typedef int (QMPI_Free_mem_t) (QMPI_Context context, int tool_id, void *base);
typedef int (QMPI_Get_t) (QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win);
typedef int (QMPI_Get_c_t) (QMPI_Context context, int tool_id, void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win);
typedef int (QMPI_Get_accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count,
             MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
typedef int (QMPI_Get_accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
             MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win);
typedef int (QMPI_Put_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Win win);
typedef int (QMPI_Put_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win);
typedef int (QMPI_Raccumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
             MPI_Request *request);
typedef int (QMPI_Raccumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_t) (QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_c_t) (QMPI_Context context, int tool_id, void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count,
             MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
             MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rput_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rput_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Win_allocate_t) (QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
             MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_c_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_shared_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_shared_c_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_attach_t) (QMPI_Context context, int tool_id, MPI_Win win, void *base,
             MPI_Aint size);
typedef int (QMPI_Win_complete_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_create_t) (QMPI_Context context, int tool_id, void *base, MPI_Aint size,
             int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_create_c_t) (QMPI_Context context, int tool_id, void *base, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_create_dynamic_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_detach_t) (QMPI_Context context, int tool_id, MPI_Win win, const void *base);
typedef int (QMPI_Win_fence_t) (QMPI_Context context, int tool_id, int assert, MPI_Win win);
typedef int (QMPI_Win_flush_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_flush_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_flush_local_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_flush_local_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_free_t) (QMPI_Context context, int tool_id, MPI_Win *win);
typedef int (QMPI_Win_get_group_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Group *group);
typedef int (QMPI_Win_get_info_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Info *info_used);
typedef int (QMPI_Win_get_name_t) (QMPI_Context context, int tool_id, MPI_Win win, char *win_name,
             int *resultlen);
typedef int (QMPI_Win_lock_t) (QMPI_Context context, int tool_id, int lock_type, int rank,
             int assert, MPI_Win win);
typedef int (QMPI_Win_lock_all_t) (QMPI_Context context, int tool_id, int assert, MPI_Win win);
typedef int (QMPI_Win_post_t) (QMPI_Context context, int tool_id, MPI_Group group, int assert,
             MPI_Win win);
typedef int (QMPI_Win_set_info_t) (QMPI_Context context, int tool_id, MPI_Win win, MPI_Info info);
typedef int (QMPI_Win_set_name_t) (QMPI_Context context, int tool_id, MPI_Win win,
             const char *win_name);
typedef int (QMPI_Win_shared_query_t) (QMPI_Context context, int tool_id, MPI_Win win, int rank,
             MPI_Aint *size, int *disp_unit, void *baseptr);
typedef int (QMPI_Win_shared_query_c_t) (QMPI_Context context, int tool_id, MPI_Win win, int rank,
             MPI_Aint *size, MPI_Aint *disp_unit, void *baseptr);
typedef int (QMPI_Win_start_t) (QMPI_Context context, int tool_id, MPI_Group group, int assert,
             MPI_Win win);
typedef int (QMPI_Win_sync_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_test_t) (QMPI_Context context, int tool_id, MPI_Win win, int *flag);
typedef int (QMPI_Win_unlock_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_unlock_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_wait_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Close_port_t) (QMPI_Context context, int tool_id, const char *port_name);
typedef int (QMPI_Comm_accept_t) (QMPI_Context context, int tool_id, const char *port_name,
             MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm);
typedef int (QMPI_Comm_connect_t) (QMPI_Context context, int tool_id, const char *port_name,
             MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm);
typedef int (QMPI_Comm_disconnect_t) (QMPI_Context context, int tool_id, MPI_Comm *comm);
typedef int (QMPI_Comm_get_parent_t) (QMPI_Context context, int tool_id, MPI_Comm *parent);
typedef int (QMPI_Comm_join_t) (QMPI_Context context, int tool_id, int fd, MPI_Comm *intercomm);
typedef int (QMPI_Comm_spawn_t) (QMPI_Context context, int tool_id, const char *command,
             char *argv[], int maxprocs, MPI_Info info, int root, MPI_Comm comm,
             MPI_Comm *intercomm, int array_of_errcodes[]);
typedef int (QMPI_Comm_spawn_multiple_t) (QMPI_Context context, int tool_id, int count,
             char *array_of_commands[], char **array_of_argv[], const int array_of_maxprocs[],
             const MPI_Info array_of_info[], int root, MPI_Comm comm, MPI_Comm *intercomm,
             int array_of_errcodes[]);
typedef int (QMPI_Lookup_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, char *port_name);
typedef int (QMPI_Open_port_t) (QMPI_Context context, int tool_id, MPI_Info info, char *port_name);
typedef int (QMPI_Publish_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, const char *port_name);
typedef int (QMPI_Unpublish_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, const char *port_name);
typedef int (QMPIX_Stream_create_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPIX_Stream *stream);
typedef int (QMPIX_Stream_free_t) (QMPI_Context context, int tool_id, MPIX_Stream *stream);
typedef int (QMPIX_Stream_comm_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPIX_Stream stream, MPI_Comm *newcomm);
typedef int (QMPIX_Stream_comm_create_multiplex_t) (QMPI_Context context, int tool_id,
             MPI_Comm comm, int count, MPIX_Stream array_of_streams[], MPI_Comm *newcomm);
typedef int (QMPIX_Comm_get_stream_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int idx,
             MPIX_Stream *stream);
typedef int (QMPIX_Stream_progress_t) (QMPI_Context context, int tool_id, MPIX_Stream stream);
typedef int (QMPIX_Start_progress_thread_t) (QMPI_Context context, int tool_id,
             MPIX_Stream stream);
typedef int (QMPIX_Stop_progress_thread_t) (QMPI_Context context, int tool_id, MPIX_Stream stream);
typedef int (QMPIX_Stream_send_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index);
typedef int (QMPIX_Stream_send_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             int source_stream_index, int dest_stream_index);
typedef int (QMPIX_Stream_isend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_isend_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             int source_stream_index, int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_recv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Status *status);
typedef int (QMPIX_Stream_recv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Status *status);
typedef int (QMPIX_Stream_irecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_irecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Send_enqueue_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPIX_Send_enqueue_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPIX_Recv_enqueue_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPIX_Recv_enqueue_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPIX_Isend_enqueue_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPIX_Isend_enqueue_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPIX_Irecv_enqueue_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPIX_Irecv_enqueue_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPIX_Wait_enqueue_t) (QMPI_Context context, int tool_id, MPI_Request *request,
             MPI_Status *status);
typedef int (QMPIX_Waitall_enqueue_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], MPI_Status *array_of_statuses);
typedef int (QMPIX_Allreduce_enqueue_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPIX_Allreduce_enqueue_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPIX_Threadcomm_init_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int num_threads, MPI_Comm *newthreadcomm);
typedef int (QMPIX_Threadcomm_free_t) (QMPI_Context context, int tool_id, MPI_Comm *threadcomm);
typedef int (QMPIX_Threadcomm_start_t) (QMPI_Context context, int tool_id, MPI_Comm threadcomm);
typedef int (QMPIX_Threadcomm_finish_t) (QMPI_Context context, int tool_id, MPI_Comm threadcomm);
typedef double (QMPI_Wtick_t) (QMPI_Context context, int tool_id);
typedef double (QMPI_Wtime_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Cart_coords_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
             int maxdims, int coords[]);
typedef int (QMPI_Cart_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int ndims,
             const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart);
typedef int (QMPI_Cart_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int maxdims,
             int dims[], int periods[], int coords[]);
typedef int (QMPI_Cart_map_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int ndims,
             const int dims[], const int periods[], int *newrank);
typedef int (QMPI_Cart_rank_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const int coords[], int *rank);
typedef int (QMPI_Cart_shift_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int direction,
             int disp, int *rank_source, int *rank_dest);
typedef int (QMPI_Cart_sub_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const int remain_dims[], MPI_Comm *newcomm);
typedef int (QMPI_Cartdim_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *ndims);
typedef int (QMPI_Dims_create_t) (QMPI_Context context, int tool_id, int nnodes, int ndims,
             int dims[]);
typedef int (QMPI_Dist_graph_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int n,
             const int sources[], const int degrees[], const int destinations[],
             const int weights[], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
typedef int (QMPI_Dist_graph_create_adjacent_t) (QMPI_Context context, int tool_id,
             MPI_Comm comm_old, int indegree, const int sources[], const int sourceweights[],
             int outdegree, const int destinations[], const int destweights[], MPI_Info info,
             int reorder, MPI_Comm *comm_dist_graph);
typedef int (QMPI_Dist_graph_neighbors_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int maxindegree, int sources[], int sourceweights[], int maxoutdegree,
             int destinations[], int destweights[]);
typedef int (QMPI_Dist_graph_neighbors_count_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *indegree, int *outdegree, int *weighted);
typedef int (QMPI_Get_hw_resource_info_t) (QMPI_Context context, int tool_id, MPI_Info *hw_info);
typedef int (QMPI_Graph_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int nnodes,
             const int indx[], const int edges[], int reorder, MPI_Comm *comm_graph);
typedef int (QMPI_Graph_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int maxindex,
             int maxedges, int indx[], int edges[]);
typedef int (QMPI_Graph_map_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int nnodes,
             const int indx[], const int edges[], int *newrank);
typedef int (QMPI_Graph_neighbors_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
             int maxneighbors, int neighbors[]);
typedef int (QMPI_Graph_neighbors_count_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int rank, int *nneighbors);
typedef int (QMPI_Graphdims_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *nnodes,
             int *nedges);
typedef int (QMPI_Topo_test_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *status);

#endif /* MPI_PROTO_H_INCLUDED */
