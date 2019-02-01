#ifndef FAKE_MPI_H
#define FAKE_MPI_H
// Fake MPI for serial execution

#include <limits.h> // HOST_NAME_MAX
#include <unistd.h> // gethostname
#include <string.h> // strlen
#include <assert.h> // assert
#include <stdlib.h> // exit
#include <time.h> // clock_gettime
#include <stdio.h> // puts for debug

const int MPI_MAX_PROCESSOR_NAME = HOST_NAME_MAX;
const int MPI_MAX_INFO_KEY = 128;
const int MPI_MAX_INFO_VAL = 128;
const int MPI_MAX_ERROR_STRING = 128;

enum {//: int { // error classes
	MPI_SUCCESS,          // No error
	MPI_ERR_BUFFER,       // Invalid buffer pointer
	MPI_ERR_COUNT,        // Invalid count argument
	MPI_ERR_TYPE,         // Invalid datatype argument
	MPI_ERR_TAG,          // Invalid tag argument
	MPI_ERR_COMM,         // Invalid communicator
	MPI_ERR_RANK,         // Invalid rank
	MPI_ERR_REQUEST,      // Invalid request (handle)
	MPI_ERR_ROOT,         // Invalid root
	MPI_ERR_GROUP,        // Invalid group
	MPI_ERR_OP,           // Invalid operation
	MPI_ERR_TOPOLOGY,     // Invalid topology
	MPI_ERR_DIMS,         // Invalid dimension argument
	MPI_ERR_ARG,          // Invalid argument of some other kind
	MPI_ERR_UNKNOWN,      // Unknown error
	MPI_ERR_TRUNCATE,     // Message truncated on receive
	MPI_ERR_OTHER,        // Known error not in this list
	MPI_ERR_INTERN,       // Internal MPI (implementation) error
	MPI_ERR_IN_STATUS,    // Error code is in status
	MPI_ERR_PENDING,      // Pending request
	MPI_ERR_KEYVAL,       // Invalid keyval has been passed
	MPI_ERR_NO_MEM,       // MPI_ALLOC_MEM failed because memory is exhausted
	MPI_ERR_BASE,         // Invalid base passed to MPI_FREE_MEM
	MPI_ERR_INFO_KEY,     // Key longer than MPI_MAX_INFO_KEY
	MPI_ERR_INFO_VALUE,   // Value longer than MPI_MAX_INFO_VAL
	MPI_ERR_INFO_NOKEY,   // Invalid key passed to MPI_INFO_DELETE
	MPI_ERR_SPAWN,        // Error in spawning processes
	MPI_ERR_PORT,         // Invalid port name passed to MPI_COMM_CONNECT
	MPI_ERR_SERVICE,      // Invalid service name passed to MPI_UNPUBLISH_NAME
	MPI_ERR_NAME,         // Invalid service name passed to MPI_LOOKUP_NAME
	MPI_ERR_WIN,          // Invalid win argument
	MPI_ERR_SIZE,         // Invalid size argument
	MPI_ERR_DISP,         // Invalid disp argument
	MPI_ERR_INFO,         // Invalid info argument
	MPI_ERR_LOCKTYPE,     // Invalid locktype argument
	MPI_ERR_ASSERT,       // Invalid assert argument
	MPI_ERR_RMA_CONFLICT, // Conflicting accesses to window
	MPI_ERR_RMA_SYNC,     // Wrong synchronization of RMA calls 
	MPI_ERR_RMA_RANGE,    // Target memory is not part of the window (in the case of a window created with MPI_WIN_CREATE_DYNAMIC, target memory is not attached)
	MPI_ERR_RMA_ATTACH,   // Memory cannot be attached (e.g., because of resource exhaustion)
	MPI_ERR_RMA_SHARED,   // Memory cannot be shared (e.g., some process in the group of the specified communicator cannot expose shared memory)
	MPI_ERR_RMA_FLAVOR,   // Passed window has the wrong flavor for the called function
	MPI_ERR_FILE,         // Invalid file handle
	MPI_ERR_NOT_SAME,     // Collective argument not identical on all processes, or collective routines called in a different order by different processes
	MPI_ERR_AMODE,        // Error related to the amode passed to MPI_FILE_OPEN
	MPI_ERR_UNSUPPORTED_DATAREP, // Unsupported datarep passed to MPI_FILE_SET_VIEW
	MPI_ERR_UNSUPPORTED_OPERATION, // Unsupported operation, such as seeking on a file which supports sequential access only
	MPI_ERR_NO_SUCH_FILE, // File does not exist
	MPI_ERR_FILE_EXISTS,  // File exists
	MPI_ERR_BAD_FILE,     // Invalid file name (e.g., path name too long)
	MPI_ERR_ACCESS,       // Permission denied
	MPI_ERR_NO_SPACE,     // Not enough space
	MPI_ERR_QUOTA,        // Quota exceeded
	MPI_ERR_READ_ONLY,    // Read-only file or file system
	MPI_ERR_FILE_IN_USE,  // File operation could not be completed, as the file is currently open by some process
	MPI_ERR_DUP_DATAREP,  // Conversion functions could not be registered because a data representation identifier that was already defined was passed to MPI_REGISTER_DATAREP
	MPI_ERR_CONVERSION,   // An error occurred in a user supplied data conversion function.
	MPI_ERR_IO,           // Other I/O error
	MPI_ERR_LASTCODE      // Last error code
}; // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node222.htm

//struct MPI_Group{};

typedef enum {
	MPI_GROUP_NULL = 0,
	MPI_GROUP_EMPTY
} MPI_Group;

// MPI Info handling.
//  Would not be too hard to create an std::map implementation
//struct MPI_Info{};
//const struct MPI_Info MPI_INFO_NULL;

typedef enum { //MPI_Info {
	MPI_INFO_NULL
} MPI_Info ;
//typedef struct {} MPI_Info;

int MPI_Info_create(MPI_Info *info);

int MPI_Info_free(MPI_Info *info);

int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo);

int MPI_Info_delete(MPI_Info info, const char *key);

int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag);

int MPI_Info_get_nkeys(MPI_Info info, int *nkeys);

int MPI_Info_get_nthkey(MPI_Info info, int n, char *key);

int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag);

int MPI_Info_set(MPI_Info info, const char *key, const char *value);

const int MPI_MAX_PORT_NAME = 128;

int MPI_Open_port(MPI_Info info, char *port_name);

int MPI_Close_port(const char *port_name);

enum { // level of thread support
	MPI_THREAD_SINGLE,
	MPI_THREAD_FUNNELED,
	MPI_THREAD_SERIALIZED,
	MPI_THREAD_MULTIPLE
}; // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm

typedef unsigned long long MPI_Aint;

struct MPI_Datatype_impl_;
typedef struct MPI_Datatype_impl_* MPI_Datatype;

struct MPI_Datatype_impl_{
//	MPI_Aint lb;
//	MPI_Aint extent;
	int count; // [in] number of blocks (integer) -- also number of entries in arrays array_of_types , array_of_displacements and array_of_blocklengths 
	int* blocklens; // [in] number of elements in each block (array) 
	MPI_Aint* indices; // [in] byte displacement of each block (array) 
	MPI_Datatype* old_types; // [in] type of elements in each block (array of handles to datatype objects)
};

static MPI_Datatype MPI_CHAR = NULL;
static MPI_Datatype MPI_INT = NULL;
static MPI_Datatype MPI_DOUBLE = NULL;
static MPI_Datatype MPI_BYTE = NULL;
static MPI_Datatype MPI_FLOAT_INT = NULL;

typedef enum {
	MPI_DATATYPE_NULL = 0,
//	MPI_CHAR,
	MPI_SHORT,
//	MPI_INT,
	MPI_LONG,
	MPI_LONG_LONG_INT,
	MPI_LONG_LONG,
	MPI_SIGNED_CHAR,
	MPI_UNSIGNED_CHAR,
	MPI_UNSIGNED_SHORT,
	MPI_UNSIGNED,
	MPI_UNSIGNED_LONG,
	MPI_UNSIGNED_LONG_LONG,
	MPI_FLOAT,
//	MPI_DOUBLE,
	MPI_LONG_DOUBLE,
	MPI_WCHAR,
	MPI_C_BOOL,
	MPI_INT8_T,
	MPI_INT16_T,
	MPI_INT32_T,
	MPI_INT64_T,
	MPI_UINT8_T,
	MPI_UINT16_T,
	MPI_UINT32_T,
	MPI_UINT64_T,
	MPI_AINT,
	MPI_COUNT,
	MPI_OFFSET,
	MPI_C_COMPLEX,
	MPI_C_FLOAT_COMPLEX,
	MPI_C_DOUBLE_COMPLEX,
//	MPI_BYTE,
	MPI_PACKED,
	MPI_CXX_BOOL,
	MPI_CXX_FLOAT_COMPLEX,
	MPI_CXX_DOUBLE_COMPLEX,
	MPI_CXX_LONG_DOUBLE_COMPLEX,
//	MPI_FLOAT_INT,
	MPI_DOUBLE_INT,
	MPI_LONG_INT,
	MPI_2INT,
	MPI_SHORT_INT,
	MPI_LONG_DOUBLE_INT,
} MPI_Datatype_;


inline int MPI_Type_create_struct(int count, const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
{
    return MPI_SUCCESS;
}

//inline int MPI_Type_size(MPI_Datatype datatype, int *size)
//{
//    if (size) *size = 0;
//    return MPI_SUCCESS;
//}

const int MPI_MAX_OBJECT_NAME = 128;

inline int MPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen)
{
    if (resultlen) *resultlen = 0;
    return MPI_SUCCESS;
}

inline int MPI_Type_set_name(MPI_Datatype datatype, const char *type_name)
{
    return MPI_SUCCESS;
}

enum {
	MPI_IDENT,
	MPI_CONGRUENT,
	MPI_SIMILAR,
	MPI_UNEQUAL
};

//struct MPI_Comm_impl_{
//	errorhandler_;
//};
struct MPI_Comm_impl_;
typedef struct MPI_Comm_impl_* MPI_Comm;

// http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node218.htm
typedef void MPI_Comm_errhandler_function(MPI_Comm *, int *, ...); 

// typedef enum {
//	MPI_COMM_NULL = 0, 
//	MPI_COMM_WORLD = 1, 
//	MPI_COMM_SELF = 2,
//} MPI_Comm;

struct MPI_Errhandler_impl_{
	MPI_Comm_errhandler_function* func_;
};

typedef struct MPI_Errhandler_impl_* MPI_Errhandler;

void fatal_error_(MPI_Comm * comm, int * errcode, ...);
void no_op_error_(MPI_Comm * comm, int * errcode, ...){}


static MPI_Errhandler MPI_ERRORS_ARE_FATAL;
static MPI_Errhandler MPI_ERRORS_RETURN;
//static struct MPI_Errhandler_impl_ MPI_ERRORS_ARE_FATAL = {&fatal_error_};
//static struct MPI_Errhandler_impl_ MPI_ERRORS_RETURN    = {&no_op_error_};

static MPI_Errhandler* MPI_ERRHANDLER_NULL = NULL;

struct MPI_Comm_impl_{
	MPI_Errhandler errhandler_;
};

extern MPI_Comm MPI_COMM_NULL = NULL;//{{&fatal_error_}};
extern MPI_Comm MPI_COMM_WORLD = NULL;
extern MPI_Comm MPI_COMM_SELF = NULL;

//typedef MPI_Comm_impl* MPI_Comm;

//typedef enum { // Predefined error handlers
//	MPI_ERRORS_ARE_FATAL, // The handler, when called, causes the program to abort on all executing processes. This has the same effect as if MPI_ABORT was called by the process that invoked the handler. 
//	MPI_ERRORS_RETURN     // The handler has no effect other than returning the error code to the user. 
//} MPI_Errhandler;


static const void* MPI_IN_PLACE = (const void*)(-1);

typedef enum { // redefined operations are supplied for MPI_REDUCE 
	MPI_OP_NULL, // TODO: is this an operator?
	MPI_MAX, // maximum
	MPI_MIN, // minimum
	MPI_SUM, // sum
	MPI_PROD, // product
	MPI_LAND, // logical and
	MPI_BAND, // bitwise and
	MPI_LOR, // logical or
	MPI_BOR, // bitwise or 
	MPI_LXOR, // logical exclusive or (xor)
	MPI_BXOR, // bitwise excluse or (xor)
	MPI_MAXLOC, // max value and location
	MPI_MINLOC, // min value and location
//	MPI_REPLACE, // TODO: is this an operator?
	MPI_NO_OP, // TODO: is this an operator?
	MPI_OP_LASTCODE // not standard?
} MPI_Op; // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node112.htm

//  MPI Requests
typedef enum {
	MPI_REQUEST_NULL
} MPI_Request;

typedef struct {
    int MPI_SOURCE;
    int MPI_TAG;
} MPI_Status;

//static MPI_Status *MPI_STATUS_IGNORE = NULL;

//static char *MPI_ARGV_NULL[] = {};
//static int MPI_ERRCODES_IGNORE[] = {};

typedef struct {} MPI_Message;
//[[noreturn]]
static inline 
int MPI_Abort( // Terminates MPI execution environment
	MPI_Comm comm, // [in] communicator of tasks to abort
	int errorcode  // [in] error code to return to invoking environment
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node225.htm
	exit(errorcode);
//	return MPI_SUCCESS; // function never returns
}
static inline 
int MPI_Allreduce( // Combines values from all processes and distributes the result back to all processes 
	const void* sendbuf,   // starting address of send buffer (choice)
	void* recvbuf,         // starting address of receive buffer (choice)
	int count,             // number of elements in send buffer (non-negative)
	MPI_Datatype datatype, // data type of elements of send buffer (handle)
	MPI_Op op,             // operation (handle)
	MPI_Comm comm          // communicator (handle)
) // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node117.htm#Node117
{
	assert(&comm != &MPI_COMM_NULL);
	assert(sendbuf == MPI_IN_PLACE);
	return MPI_SUCCESS;
}
void fatal_error_(MPI_Comm * comm, int * errcode, ...){
	switch(*errcode){
		case MPI_ERR_COMM : puts("[] *** MPI_ERR_COMM: invalid communicator\n[] *** MPI_ERRORS_ARE_FATAL (will now abort)"); MPI_Abort(*comm, *errcode);
	}
}
static inline
int MPI_Barrier( // Blocks until all processes in the communicator have reached this routine. 
	MPI_Comm comm // [in] communicator (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node100.htm
	return MPI_SUCCESS;
}
static inline 
int MPI_Bcast( // Broadcasts a message from the process with rank "root" to all other processes of the communicator 
	void* buffer,          // starting address of buffer (choice)
	int count,             // number of entries in buffer (non-negative integer)
	MPI_Datatype datatype, // data type of buffer (handle)
	int root,              // rank of broadcast root (integer)
	MPI_Comm comm          // communicator (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node101.htm
	if(comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	assert(0); // TODO: is there something to implement?
	return MPI_SUCCESS; 
}
int MPI_Comm_call_errhandler( // Call the error handler installed on a communicator 
	MPI_Comm comm, // [in] communicator with error handler (handle) 
	int errorcode  // [in] error code (integer) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node223.htm
//	if(comm == MPI_COMM_NULL){
//		return MPI_ERR_COMM;
//	}
	comm->errhandler_->func_(&comm, &errorcode);
	return MPI_SUCCESS;
}
int MPI_Comm_call_errhandler_(
	MPI_Comm comm, 
	int errorcode
){
	MPI_Comm_call_errhandler(comm, errorcode);
	return errorcode;
}
int MPI_Cart_create( // Makes a new communicator to which topology information has been attached 
	MPI_Comm comm_old,  // [in] input communicator (handle) 
	int ndims,          // [in] number of dimensions of cartesian grid (integer) 
	int *dims,          // [in] integer array of size ndims specifying the number of processes in each dimension 
	int *periods,       // [in] logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension 
	int reorder,        // [in] ranking may be reordered (true) or not (false) (logical) 
	MPI_Comm *comm_cart // [out] communicator with new cartesian topology (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node192.htm
	if(dims == NULL) return MPI_ERR_DIMS;
	for(int i = 0; i != ndims; ++i) if(dims[i] <= 0) return MPI_ERR_DIMS;
	*comm_cart = comm_old; // TODO: allocate unique handle
	return MPI_SUCCESS;
}
int MPI_Cart_sub( // Partitions a communicator into subgroups which form lower-dimensional cartesian subgrids 
	MPI_Comm comm,    // [in] communicator with cartesian structure (handle) 
	int *remain_dims, // [in] the ith entry of remain_dims specifies whether the ith dimension is kept in the subgrid (true) or is dropped (false) (logical vector) 
	MPI_Comm *newcomm // [out] communicator containing the subgrid that includes the calling process (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node198.htm
	if(&comm == &MPI_COMM_NULL) return MPI_ERR_COMM;
	*newcomm = comm; // TODO: allocate unique handle
	return MPI_SUCCESS;
}
static inline
int MPI_Comm_compare( // Compares two communicators
  MPI_Comm comm1, // [in] comm1 (handle) 
  MPI_Comm comm2, // [in] comm2 (handle) 
  int *result // [out] integer which is MPI_IDENT if the contexts and groups are the same, MPI_CONGRUENT if different contexts but identical groups, MPI_SIMILAR if different contexts but similar groups, and MPI_UNEQUAL otherwise 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node155.htm
	if(&comm1 == &MPI_COMM_NULL || &comm2 == &MPI_COMM_NULL) return MPI_ERR_COMM;
	if(result == NULL) return MPI_ERR_ARG;
	*result = MPI_CONGRUENT; //MPI_IDENT;
	return MPI_SUCCESS;
}
static inline
int MPI_Comm_create( // Creates a new communicator 
	MPI_Comm comm,    // [in] communicator (handle) 
	MPI_Group group,  // [in] group, which is a subset of the group of comm (handle) 
	MPI_Comm *newcomm // [out] new communicator (handle) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node156.htm
int MPI_Comm_free( // Marks the communicator object for deallocation 
  MPI_Comm *comm // [in] Communicator to be destroyed (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node157.htm
	if(*comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	free(*comm);
	*comm = MPI_COMM_NULL;
	return MPI_SUCCESS;
} 
// Determines the size of the remote group associated with an inter-communictor 
int MPI_Comm_remote_size(MPI_Comm comm, int *size); 
static inline 
int MPI_Comm_size( // Determines the size of the group associated with a communicator
	MPI_Comm comm, // communicator (handle)
	int *size // number of processes in the group of comm (integer)
){
	if(comm == MPI_COMM_NULL){
		int error = MPI_ERR_COMM;
		comm->errhandler_->func_(&comm, &error);
		return error;
	}
	*size = 1;
	return MPI_SUCCESS;
} // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node155.htm
int MPI_Comm_spawn( // Spawn up to maxprocs instances of a single MPI application 
	const char *command,    // [in] name of program to be spawned (string, at root) 
	char *argv[],           // [in] arguments to command (array of strings, at root) 
	int maxprocs,           // maximum number of processes to start(int at root) 
	MPI_Info info,          // a set of key-value pairs telling the runtime system where and how to start the processes (handle, significant only at root) 
	int root,               // rank of process in which previous arguments are examined (integer) 
	MPI_Comm comm,          // intracommunicator containing group of spawning processes (handle) 
	MPI_Comm *intercomm,    // [out] intercommunicator between original group and the newly spawned group (handle) 
	int array_of_errcodes[] // [out] one code per process (array of integer) 
);
int MPI_Comm_split( // Creates new communicators based on colors and keys 
	MPI_Comm comm,     // [in] communicator (handle)
	int color,         // [in] control of subset assignment (integer)
	int key,           // [in] control of rank assigment (integer)
	MPI_Comm *newcomm  // [out] new communicator (handle)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node156.htm
int MPI_Comm_get_errhandler( // Get the error handler attached to a communicator 
	MPI_Comm comm,              // [in] communicator (handle) 
	MPI_Errhandler *errhandler  // [out] handler currently associated with communicator (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node220.htm
	return MPI_SUCCESS;
}
int MPI_Comm_get_name( // Return the print name from the communicator 
	MPI_Comm comm,   // communicator whose name is to be returned (handle)
	char *comm_name, // the name previously stored on the communicator, or an...
	int *resultlen   // length of returned name (integer)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node179.htm
int MPI_Comm_get_parent( // Return the parent communicator for this process 
  MPI_Comm *parent // [out] the parent communicator (handle) 
);
int MPI_Comm_rank( // MPI_Group_rank Returns the rank of this process in the given group 
	MPI_Comm comm, // [in] group (handle) 
	int* rank      // [out] rank of the calling process in group, or MPI_UNDEFINED if the process is not a member (integer) 
){
	if(comm == MPI_COMM_NULL){
		MPI_Comm_call_errhandler(comm, MPI_ERR_COMM);
		return MPI_ERR_COMM;
	}
	*rank = 0;
	return MPI_SUCCESS;
}
int MPI_Get_count( // Gets the number of "top level" elements 
	MPI_Status *status,    // [in] return status of receive operation (Status) 
	MPI_Datatype datatype, // [out] number of received elements (integer)
	int *count             // [in] datatype of each receive buffer element (han 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node51.htm
int MPI_Get_version( // Return the version number of MPI 
	int* version,   // [out] Version of MPI
	int* subversion // [out] Suversion of MPI
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node209.htm
	*version = 3;
	*subversion = 1;
	return MPI_SUCCESS;
}
int MPI_Comm_set_name( // Sets the print name for a communicator 
	MPI_Comm comm,        // [in] communicator to name (handle) 
	const char *comm_name // [in] Name for communicator (string)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node179.htm
inline int MPI_Comm_create( // Creates a new communicator 
	MPI_Comm comm,    // [in]  communicator (handle)
	MPI_Group group,  // [in]  group, which is a subset of the group of comm (handle)
	MPI_Comm *newcomm // [out] new communicator (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node156.htm
	if(&comm == &MPI_COMM_NULL) return MPI_ERR_COMM;
	if(group == MPI_GROUP_NULL) return MPI_ERR_GROUP;
	newcomm = (MPI_Comm*)malloc(sizeof(MPI_Comm));
	return MPI_SUCCESS;
}
static inline
int MPI_Comm_create_errhandler( // Create a communicator error handler 
	MPI_Comm_errhandler_function *comm_errhandler_fn, // [in] user defined error handling procedure (function) 
	MPI_Errhandler *errhandler                        // [out] MPI error handler (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node218.htm
	*errhandler = (struct MPI_Errhandler_impl_*)malloc(sizeof(struct MPI_Errhandler_impl_));
	(*errhandler)->func_ = comm_errhandler_fn;
	return MPI_SUCCESS;
}
int MPI_Comm_create_group( // must be called by all processes in group, which is a subgroup of the group of comm
	MPI_Comm comm,    // intracommunicator (handle)
	MPI_Group group,  // group, which is a subset of the group of comm (handle)
	int tag,          // tag (integer)
	MPI_Comm *newcomm // new communicator (handle)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node156.htm
static inline 
int MPI_Comm_disconnect( // MPI_Comm_disconnect
	MPI_Comm *comm // [in] communicator (handle)
){
	if(*comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	free(*comm);
	*comm = MPI_COMM_NULL;
	return MPI_SUCCESS;
}
static inline 
int MPI_Comm_dup( // Duplicates an existing communicator with all its cached information 
	MPI_Comm comm,    // communicator (handle)
	MPI_Comm *newcomm // copy of comm (handle)
)
{
	if(&comm == &MPI_COMM_NULL) return MPI_ERR_COMM;
	*newcomm = (struct MPI_Comm_impl_*)malloc(sizeof(struct MPI_Comm_impl_));
	(**newcomm).errhandler_ = comm->errhandler_;
	return MPI_SUCCESS;
}
int MPI_Comm_set_errhandler( // Set the error handler for a communicator 
	MPI_Comm comm,             // [in] communicator (handle) 
	MPI_Errhandler errhandler  // [in] new error handler for communicator (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node218.htm
	assert(comm != MPI_COMM_NULL);
	if(comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	comm->errhandler_ = errhandler; //->func_ = errhandler->func_;
	return MPI_SUCCESS;
}
static inline 
int MPI_Gather( // Gathers together values from a group of processes 
	const void *sendbuf, // [in] starting address of send buffer (choice) 
	int sendcnt, // [in] number of elements in send buffer (integer)
	MPI_Datatype sendtype, // [in] data type of send buffer elements (handle)
	void *recvbuf, // [out] address of receive buffer (choice, significant only at root) 
	int recvcnt, // [in] number of elements for any single receive (integer, significant only at root) 
	MPI_Datatype recvtype, // [in] data type of recv buffer elements (significant only at root) (handle) 
	int root, // [in] rank of receiving process (integer) 
	MPI_Comm comm // [in] communicator (handle) 
);
static inline 
int MPI_Get_processor_name( //  the name of the processor on which it was called at the moment of the call. 
	char *name, // A unique specifier for the actual (as opposed to virtual) node.
	int *resultlen // Length (in printable characters) of the result returned in name
){
	if(gethostname(name, MPI_MAX_PROCESSOR_NAME) > 0) 
		return MPI_ERR_UNKNOWN;
	*resultlen = strlen(name);
	return MPI_SUCCESS;
}  // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node210.htm#Node215
int MPI_Group_free( // Frees a group 
	MPI_Group *group // group (handle)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node153.htm
int MPI_Group_incl( // Produces a group by reordering an existing group and taking only listed members 
 	MPI_Group group,  // [in] group (handle) 
 	int n,              // [in] number of elements in array ranks (and size of newgroup ) (integer)
 	const int ranks[], 
 	MPI_Group *newgroup // [in] ranks of processes in group to appear in newgroup (array of integers) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node152.htm
int MPI_Group_excl( // Produces a group by reordering an existing group and taking only unlisted members 
  MPI_Group group,    //[in] group (handle) 
  int n,              // [in] number of elements in array ranks (integer) 
  const int ranks[],         // [in] array of integer ranks in group not to appear in newgroup 
  MPI_Group *newgroup // [out] new group derived from above, preserving the order defined by group (handle) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node152.htm
int MPI_Group_range_excl( // Produces a group by excluding ranges of processes from an existing group 
	MPI_Group group,    // [in] group (handle) 
	int n,              // [in] number of elements in array ranks (integer) 
	int ranges[][3],    // [in] a one-dimensional array of integer triplets of the form (first rank, last rank, stride), indicating the ranks in group of processes to be excluded from the output group newgroup . 
	MPI_Group *newgroup // [out] new group derived from above, preserving the order in group (handle) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node152.htm
int MPI_Group_range_incl( // Creates a new group from ranges of ranks in an existing group 
	MPI_Group group,    // [in] group (handle) 
	int n,              // [in] number of triplets in array ranges (integer) 
	int ranges[][3],    // [in] a one-dimensional array of integer triplets, of the form (first rank, last rank, stride) indicating ranks in group or processes to be included in newgroup. 
	MPI_Group *newgroup // [out] new group derived from above, in the order defined by ranges (handle) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node152.htm
int MPI_Group_rank( // Determines the rank of the calling process in the communicator 
	MPI_Group group, // group (handle)
	int *rank        // rank of the calling process in group, or MPI_UNDEFINED if the process is not a member (integer)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node151.htm
int MPI_Group_size( // Returns the size of a group 
	MPI_Group group, // [in] group (handle)
	int *size        // [out] number of processes in the group (integer) 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node151.htm
int MPI_Group_compare( // Compares two groups 
	MPI_Group group1, // [in] group1 (handle) 
	MPI_Group group2, // [in] group2 (handle) 
	int *result       // [out] integer which is MPI_IDENT if the order and members of the two groups are the same, MPI_SIMILAR if only the members are the same, and MPI_UNEQUAL otherwise 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node151.htm#Node151
int MPI_Group_intersection( // Produces a group as the intersection of two existing groups 
	MPI_Group group1,   // [in] first group (handle) 
	MPI_Group group2,   // [in] second group (handle) 
	MPI_Group *newgroup // [out] intersection group (handle) 
);
int MPI_Group_translate_ranks( // Translates the ranks of processes in one group to those in another group 
	MPI_Group group1,   // [in] group1 (handle) 
	int n,              // [in] number of ranks in ranks1 and ranks2 arrays (integer) 
	const int ranks1[], // [in] array of zero or more valid ranks in group1 
	MPI_Group group2,   // [in] group2 (handle) 
	int ranks2[]        // [out] array of corresponding ranks in group2, MPI_UNDEFINED when no correspondence exists. 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node151.htm
int MPI_Group_union( // Produces a group by combining two groups 
	MPI_Group group1,   // first group (handle)
	MPI_Group group2,   // second group (handle)
	MPI_Group *newgroup // union group (handle)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node152.htm
int MPI_Error_class( // Converts an error code into an error class 
	int errorcode,  // Error code returned by an MPI routine
	int *errorclass // Error class associated with errorcode
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node222.htm#Node222
int MPI_Error_string( // Return a string for a given error code 
	int errorcode, // [in] Error code returned by an MPI routine or an MPI error class 
	char *string,  // [out] Text that corresponds to the errorcode 
	int *resultlen // [out] Length of string 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node221.htm#Node221
int MPI_Errhandler_free( // Frees an MPI-style errorhandler 
	MPI_Errhandler *errhandler // [in-out] MPI error handler (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node221.htm
	free(*errhandler);
	errhandler = MPI_ERRHANDLER_NULL;
	return MPI_SUCCESS;
}
static inline 
int MPI_Finalize( // Terminates MPI execution environment 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node28.htm
	free(MPI_COMM_WORLD->errhandler_); // TODO revise
	free(MPI_COMM_WORLD);
	free(MPI_COMM_SELF->errhandler_); // TODO revise
	free(MPI_COMM_SELF);

	MPI_Errhandler_free(&MPI_ERRORS_RETURN);
	MPI_Errhandler_free(&MPI_ERRORS_ARE_FATAL);

	free(MPI_INT);
	free(MPI_FLOAT_INT);

	return MPI_SUCCESS;
} 
int MPI_Finalized( // Indicates whether MPI_Finalize has been called 
	int *flag // [out] true if MPI was finalized (logical)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node227.htm
static inline 
int MPI_Igather( // Gathers together values from a group of processes 
	const void *sendbuf,   // [in] starting address of send buffer (choice) 
	int sendcnt,           // [in] number of elements in send buffer (integer)
	MPI_Datatype sendtype, // [in] data type of send buffer elements (handle)
	void *recvbuf,         // [out] address of receive buffer (choice, significant only at root) 
	int recvcnt,           // [in] number of elements for any single receive (integer, significant only at root) 
	MPI_Datatype recvtype, // [in] data type of recv buffer elements (significant only at root) (handle) 
	int root,              // [in] rank of receiving process (integer) 
	MPI_Comm comm,         // [in] communicator (handle) 
	MPI_Request *request   // [out] communicatio request (handle)
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node130.htm
	if(comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	assert(0); // TODO implementation
	return MPI_SUCCESS;
}
static inline 
int MPI_Init( // Initialize the MPI execution environment 
	int *argc,   // [in] Pointer to the number of arguments 
	char ***argv // [in] Pointer to the argument vector 
){
	MPI_Comm_create_errhandler(&fatal_error_, &MPI_ERRORS_ARE_FATAL);
	MPI_Comm_create_errhandler(&no_op_error_, &MPI_ERRORS_RETURN);
//	MPI_ERRORS_ARE_FATAL = (MPI_Comm)malloc(sizeof(struct MPI_Comm_impl_));//{&fatal_error_};
//	static struct MPI_Errhandler_impl_ MPI_ERRORS_RETURN    = {&no_op_error_};

//	MPI_COMM_NULL = (MPI_Comm)malloc(sizeof(struct MPI_Comm_impl_));
//	MPI_COMM_NULL->errhandler_ = (struct MPI_Errhandler_impl_*)malloc(sizeof(struct MPI_Errhandler_impl_)); 
//	MPI_Comm_set_errhandler(MPI_COMM_NULL, MPI_ERRORS_ARE_FATAL);

	MPI_COMM_WORLD = (MPI_Comm)malloc(sizeof(struct MPI_Comm_impl_));
	MPI_COMM_WORLD->errhandler_ = (struct MPI_Errhandler_impl_*)malloc(sizeof(struct MPI_Errhandler_impl_)); 
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
//	MPI_COMM_WORLD->errhandler_ = MPI_ERRORS_ARE_FATAL;

	MPI_COMM_SELF = (MPI_Comm)malloc(sizeof(struct MPI_Comm_impl_));
	MPI_COMM_SELF->errhandler_ = (struct MPI_Errhandler_impl_*)malloc(sizeof(struct MPI_Errhandler_impl_)); 
	MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_ARE_FATAL);
//	MPI_COMM_SELF->errhandler_ = MPI_ERRORS_ARE_FATAL;

	MPI_BYTE = 0;//(MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
//	MPI_BYTE->lb = 0;
//	MPI_BYTE->extent = sizeof(unsigned char);

	MPI_CHAR = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_CHAR->count = 1;
	MPI_CHAR->blocklens = (int*)malloc(1*sizeof(int));
	MPI_CHAR->indices = (MPI_Aint*)malloc(1*sizeof(MPI_Aint));
	MPI_CHAR->old_types = (MPI_Datatype*)malloc(1*sizeof(MPI_Datatype));;
	MPI_CHAR->blocklens[0] = 1;
	MPI_CHAR->indices[0] = 0;
	MPI_CHAR->old_types[0] = MPI_BYTE;

	MPI_INT = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_INT->count = 1;
	MPI_INT->blocklens = (int*)malloc(1*sizeof(int));
	MPI_INT->indices = (MPI_Aint*)malloc(1*sizeof(MPI_Aint));
	MPI_INT->old_types = (MPI_Datatype*)malloc(1*sizeof(MPI_Datatype));;
	MPI_INT->blocklens[0] = sizeof(int);
	MPI_INT->indices[0] = 0;
	MPI_INT->old_types[0] = MPI_BYTE;

	MPI_DOUBLE = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_DOUBLE->count = 1;
	MPI_DOUBLE->blocklens = (int*)malloc(1*sizeof(int));
	MPI_DOUBLE->indices = (MPI_Aint*)malloc(1*sizeof(MPI_Aint));
	MPI_DOUBLE->old_types = (MPI_Datatype*)malloc(1*sizeof(MPI_Datatype));;

	MPI_FLOAT_INT = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_FLOAT_INT->count = 1;
	MPI_FLOAT_INT->blocklens = (int*)malloc(1*sizeof(int));
	MPI_FLOAT_INT->indices = (MPI_Aint*)malloc(1*sizeof(MPI_Aint));
	MPI_FLOAT_INT->old_types = (MPI_Datatype*)malloc(1*sizeof(MPI_Datatype));;
	MPI_FLOAT_INT->blocklens[0] = sizeof(float) + sizeof(int);
	MPI_FLOAT_INT->indices[0] = 0;
	MPI_FLOAT_INT->old_types[0] = MPI_BYTE;


/*
	MPI_INT = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_INT->lb = 0;
	MPI_INT->extent = sizeof(int);

	MPI_DOUBLE = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_DOUBLE->lb = 0;
	MPI_DOUBLE->extent = sizeof(double);

	MPI_BYTE = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_BYTE->lb = 0;
	MPI_BYTE->extent = sizeof(unsigned char);

	MPI_FLOAT_INT = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	MPI_FLOAT_INT->lb = 0;
	MPI_FLOAT_INT->extent = sizeof(float) + sizeof(int);
*/

	return MPI_SUCCESS;
} // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node28.htm
int MPI_Init_thread( // Initialize the MPI execution environment 
	int *argc,    // [in] Pointer to the number of arguments 
	char ***argv, // [in] Pointer to the argument vector 
	int required, // [in] Level of desired thread support 
	int *provided // [out] Level of provided thread support 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm
int MPI_Initialized( // Indicates whether MPI_Init has been called. 
  int *flag // [out] Flag is true if MPI_Init or MPI_Init_thread has been called and false otherwise. 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node225.htm
int MPI_Intercomm_create( // Creates an intercommuncator from two intracommunicators 
	MPI_Comm local_comm,     // [in] Local (intra)communicator 
	int local_leader,        // [in] Rank in local_comm of leader (often 0) 
	MPI_Comm peer_comm,      // [in] Communicator used to communicate between a designated process in the other communicator. Significant only at the process in local_comm with rank local_leader. 
	int remote_leader,       // [in] Rank in peer_comm of remote leader (often 0) 
	int tag,                 // [in] Message tag to use in constructing intercommunicator; if multiple MPI_Intercomm_creates are being made, they should use different tags (more precisely, ensure that the local and remote leaders are using different tags for each MPI_intercomm_create). 
	MPI_Comm *newintercomm   // [out] Created intercommunicator 
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node168.htm
int MPI_Iprobe( // Nonblocking test for a message
	int source,        // [in] source rank, or MPI_ANY_SOURCE (integer) 
	int tag,           // [in] tag value or MPI_ANY_TAG (integer) 
	MPI_Comm comm,     // [in] communicator (handle) 
	int *flag,         // [out] True if a message with the specified source, tag...
	MPI_Status *status // [out] status object (Status) 
);
int MPI_Irecv( // Start a nonblocking receive
	void* buf,             // [out] initial address of receive buffer (choice) 
	int count,             // [in] number of elements in receive buffer (non-negative integer)
	MPI_Datatype datatype, // [in] datatype of each receive buffer element (handle)
	int source,            // [in] rank of source or MPI_ANY_SOURCE (integer)
	int tag,               // [in] message tag or MPI_ANY_TAG (integer)
	MPI_Comm comm,         // [in] communicator (handle)
	MPI_Request *request   // [out] communication request (handle)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node63.htm#Node63
int MPI_Is_thread_main( //  This function can be called by a thread to determine if it is the main thread (the thread that called MPI_INIT or MPI_INIT_THREAD).
	int *flag // true if calling thread is main thread, false otherwise (logical)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm
// Returns the upper bound on the amount of space needed to pack a message
int MPI_Pack_size(
	int incount, 
	MPI_Datatype datatype, 
	MPI_Comm comm, 
	int *size
);
int MPI_Mprobe(
	int source,           // rank of source or MPI_ANY_SOURCE (integer)
	int tag,              // message tag or MPI_ANY_TAG (integer)
	MPI_Comm comm,        // communicator (handle)
	MPI_Message *message, // returned message (handle)
	MPI_Status *status    // status object (Status)
);
int MPI_Mrecv(
	void* buf,             // initial address of receive buffer (choice)
	int count,             // number of elements in receive buffer (non-negati)
	MPI_Datatype datatype, // datatype of each receive buffer element (handle)
	MPI_Message *message,  // message (handle)
	MPI_Status *status     // status object (Status)
); 
// Blocking test for a message
int MPI_Probe( // like MPI_IMPROBE except that it is a blocking call that returns only after a matching message has been found.
	int source,        // [in] source rank, or MPI_ANY_SOURCE (integer) 
	int tag,           // [in] tag value or MPI_ANY_TAG (integer) 
	MPI_Comm comm,     // [in] communicator (handle) 
	MPI_Status *status // [out] status object (Status) 
);
int MPI_Query_thread( //  The following function can be used to query the current level of thread support.
	int *provided // provided level of thread support (integer)
); // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm#Node303
static inline 
int MPI_Send( // Performs a blocking send 
	const void *buf,       // [in] initial address of send buffer (choice) 
	int count,             // [in] number of elements in send buffer (nonnegat...) 
	MPI_Datatype datatype, // [in] datatype of each send buffer element (handle) 
	int dest,              // [in] rank of destination (integer) 
	int tag,               // [in] message tag (integer) 
	MPI_Comm comm          // [in] communicator (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node47.htm#Node47
	int rank = -1; MPI_Comm_rank(comm, &rank);
	assert(rank != dest);
	int size = -1; MPI_Comm_size(comm, &size);
	assert(dest < size); // Invalid rank has value 1 but must be nonnegative and less than 1 
	return MPI_SUCCESS;
}
int MPI_Sendrecv( // Sends and receives a message 
	const void *sendbuf,   // [in] initial address of send buffer (choice) 
	int sendcount,         // [in] number of elements in send buffer (integer) 
	MPI_Datatype sendtype, // [in] type of elements in send buffer (handle) 
	int dest,              // [in] rank of destination (integer) 
	int sendtag,           // [in] send tag (integer) 
	void *recvbuf,         // [out] initial address of receive buffer (choice) 
	int recvcount,         // [in] number of elements in receive buffer (integer) 
	MPI_Datatype recvtype, // [in] type of elements in receive buffer (handle) 
	int source,            // [in] rank of source (integer) 
	int recvtag,           // [in] receive tag (integer) 
	MPI_Comm comm,         // [in] communicator (handle) 
	MPI_Status *status     // [out] status object (Status). This refers to the receive operation. 
);
static inline
int MPI_Recv( // Blocking receive for a message
	void *buf,             // [out] initial address of receive buffer (choice) 
	int count,             // [in] maximum number of elements in receive buffer...
	MPI_Datatype datatype, // [in] datatype of each receive buffer element...
	int source,            // [in] rank of source (integer) 
	int tag,               // [in] message tag (integer) 
	MPI_Comm comm,         // [in] communicator (handle) 
	MPI_Status *status     // [out] status object (Status) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node50.htm#Node50
	int rank = -1; MPI_Comm_rank(comm, &rank);
	assert(rank != source);
	assert(0);
	return MPI_SUCCESS;
}
static inline int MPI_Type_commit( // Commits the datatype 
	MPI_Datatype *datatype // [in] datatype (handle) 
){
	return MPI_SUCCESS;
}
static inline 
int MPI_Type_contiguous( // Creates a contiguous datatype 
	int count, // [in] replication count (nonnegative integer) 
	MPI_Datatype oldtype, // [in] old datatype (handle) 
	MPI_Datatype *newtype // [out] new datatype (handle) 
);
static inline 
int MPI_Type_dup( // MPI_Type_dup
	MPI_Datatype oldtype, // [in] datatype (handle) 
	MPI_Datatype *newtype // [out] copy of type (handle) 
){
	int c;
	*newtype = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	(*newtype)->count = oldtype->count;
	(*newtype)->blocklens = (int*)malloc(oldtype->count*sizeof(int));
	(*newtype)->indices = (MPI_Aint*)malloc(oldtype->count*sizeof(MPI_Aint));
	(*newtype)->old_types = (MPI_Datatype*)malloc(oldtype->count*sizeof(MPI_Datatype));
	for(c = 0 ; c != oldtype->count; ++c){
		(*newtype)->blocklens[c] = oldtype->blocklens[c];
		(*newtype)->indices[c] = oldtype->indices[c];
		MPI_Type_dup(oldtype->old_types[c], (*newtype)->old_types + c);
	}
	return MPI_SUCCESS;
}
static inline
int MPI_Type_extent( // Returns the extent of a datatype, deprecated->MPI_Type_get_extent
	MPI_Datatype datatype, // [in] datatype (handle) 
	MPI_Aint *extent // [out] datatype extent (integer) 
){
	assert(0);
	return MPI_SUCCESS;
}
static inline int MPI_Type_free( // Frees the datatype 
	MPI_Datatype *datatype // [in] datatype that is freed (handle) 
){
	int c;
//	*newtype = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
//	(*newtype)->count = oldtype->count;
//	; = (MPI_Datatype*)malloc(oldtype->count*sizeof(MPI_Datatype));
	for(c = 0 ; c != (*datatype)->count; ++c){
		MPI_Type_free((*datatype)->old_types + c);
	}
	free((*datatype)->old_types);
	free((*datatype)->indices);
	free((*datatype)->blocklens);
	free((*datatype));
	return MPI_SUCCESS;
}
static inline
int MPI_Type_get_extent( // Get the lower bound and extent for a Datatype 
	MPI_Datatype datatype, // [in] datatype to get information on (handle) 
	MPI_Aint *lb, // [out] lower bound of datatype (integer) 
	MPI_Aint *extent // [out] extent of datatype (integer) 
){
	if(!lb || !extent) return MPI_ERR_ARG;
	if(datatype == MPI_BYTE){
		*lb = 0;
		*extent = 1;
		return MPI_SUCCESS;
	}
	if(datatype->count == 1){
		MPI_Aint oldlb = 0;
		MPI_Aint oldextent = 0;
		int s = MPI_Type_get_extent(datatype->old_types[0], &oldlb, &oldextent);
		assert(s == MPI_SUCCESS);
		*lb = 0;
		*extent = datatype->blocklens[0]*oldextent;
	}else assert(0);
	return MPI_SUCCESS;
}
static inline
int MPI_Type_size( // Return the number of bytes occupied by entries in the datatype 
	MPI_Datatype datatype, // [in] datatype (handle) 
	int *size // [out] datatype size (integer) 
){
//	assert(size);
	if(datatype == MPI_BYTE){
		*size = 1;
		return MPI_SUCCESS;
	}
	if(datatype == MPI_INT){
		*size = sizeof(int);
		return MPI_SUCCESS;
	}
	if(datatype->count == 1){
		int oldsize = -1;
		MPI_Type_size(datatype->old_types[0], &oldsize);
		*size = datatype->blocklens[0]*oldsize;
	}else assert(0);
//	*size = datatype->extent;
/*	switch(datatype){
		MPI_INT      : *size = sizeof(int); break;
		MPI_FLOAT    : *size = sizeof(float); break;
		MPI_DOUBLE   : *size = sizeof(double); break;
		MPI_FLOAT_INT: *size = sizeof(float) + sizeof(int); break;
		default: assert(0);
	}*/
	return MPI_SUCCESS;
}
static inline
int MPI_Type_ub( // Returns the upper bound of a datatype 
	MPI_Datatype datatype, //  [in] datatype (handle) 
	MPI_Aint *displacement
){
	if(!displacement) return MPI_ERR_ARG;
	if(!datatype) return MPI_ERR_TYPE;
	MPI_Aint lb = -1;
	MPI_Aint extent = -1;
	MPI_Type_get_extent(datatype, &lb, &extent);
	*displacement = lb + extent;
	return MPI_SUCCESS;
}
static inline
int MPI_Type_vector( // Creates a vector (strided) datatype 
	int count, // [in] number of blocks (nonnegative integer) 
	int blocklength, // [in] number of elements in each block (nonnegative integer) 
	int stride, //  [in] number of elements between start of each block (integer) 
	MPI_Datatype old_type, // [in] old datatype (handle) 
	MPI_Datatype *newtype_p // [out] new datatype (handle) 
);
int MPI_Wait( // Waits for an MPI request to complete 
	MPI_Request *request, // [in] request (handle) 
	MPI_Status *status    // [out] status object (Status). May be MPI_STATUS_IGNORE. 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node64.htm
	assert(0); // TODO implement
	return MPI_SUCCESS;
}
static inline
double MPI_Wtime( // Returns an elapsed time on the calling processor 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node37.htm
    struct timespec tw;
    clock_gettime(CLOCK_MONOTONIC, &tw);
	return 1.0*tw.tv_sec + 1e-9*tw.tv_nsec;;
}
static inline
double MPI_Wtick( // Returns the resolution of MPI_Wtime 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node37.htm
	return 1e-9;
}

inline int MPI_Status_set_cancelled(MPI_Status *status, int flag)
{
    return MPI_SUCCESS;
}

inline int MPI_Test_cancelled(const MPI_Status *status, int *flag)
{
    if (flag) *flag = -1; // -1 true
    return MPI_SUCCESS;
}

typedef void (MPI_User_function)(void *a, void *b, int *len, MPI_Datatype *);
inline int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op)
{
    return MPI_SUCCESS;
}

inline int MPI_Op_free(MPI_Op *op)
{
    if (op) *op = MPI_NO_OP;
    return MPI_SUCCESS;
}


//const int MPI_MAX_PROCESSOR_NAME = 128;
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
static inline 
int MPI_Gather( // Gathers together values from a group of processes 
	const void *sendbuf, // [in] starting address of send buffer (choice) 
	int sendcnt, // [in] number of elements in send buffer (integer)
	MPI_Datatype sendtype, // [in] data type of send buffer elements (handle)
	void *recvbuf, // [out] address of receive buffer (choice, significant only at root) 
	int recvcnt, // [in] number of elements for any single receive (integer, significant only at root) 
	MPI_Datatype recvtype, // [in] data type of recv buffer elements (significant only at root) (handle) 
	int root, // [in] rank of receiving process (integer) 
	MPI_Comm comm // [in] communicator (handle) 
){ // http://mpi-forum.org/docs/mpi-3.1/mpi31-report/node103.htm
	if(comm == MPI_COMM_NULL) return MPI_ERR_COMM;
	assert(root == 0);
	int sendsize = -1;
	int recvsize = -1;
#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif
	MPI_Type_size(sendtype, &sendsize);
	MPI_Type_size(recvtype, &recvsize);
	memcpy((char*)recvbuf, (const char*)sendbuf, max(sendcnt*sendsize, recvcnt*recvsize));
	return MPI_SUCCESS;
}

static inline 
int MPI_Type_contiguous( // Creates a contiguous datatype 
	int count, // [in] replication count (nonnegative integer) 
	MPI_Datatype oldtype, // [in] old datatype (handle) 
	MPI_Datatype *newtype // [out] new datatype (handle) 
){
//	MPI_Aint oldlb;
//	MPI_Aint oldextent;
//	int s = MPI_Type_get_extent(oldtype, &oldlb, &oldextent);
	*newtype = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	if(oldtype->count == 1){
		(*newtype)->count = 1;
		(*newtype)->blocklens = (int*)malloc(sizeof(int));
		(*newtype)->blocklens[0] = count;
		(*newtype)->indices = (MPI_Aint*)malloc(sizeof(MPI_Aint));
		(*newtype)->indices[0] = 0;
		(*newtype)->old_types = (MPI_Datatype*)malloc(count*sizeof(MPI_Datatype));
		(*newtype)->old_types[0] = MPI_BYTE;
	}else assert(0);
//	if(count < 0) return MPI_ERR_INTERN;
//	if(!newtype) return MPI_ERR_INTERN;
//	(*newtype)->lb = 0;
//	if(s != MPI_SUCCESS) return MPI_ERR_TYPE;
//	(*newtype)->extent = count*oldextent;
	return MPI_SUCCESS;
}
static inline
int MPI_Type_struct( // Creates a struct datatype 
	int count, // [in] number of blocks (integer) -- also number of entries in arrays array_of_types , array_of_displacements and array_of_blocklengths 
	int blocklens[], // [in] number of elements in each block (array) 
	MPI_Aint indices[], // [in] byte displacement of each block (array) 
	MPI_Datatype old_types[], // [in] type of elements in each block (array of handles to datatype objects)
	MPI_Datatype *newtype // [out] new datatype (handle) 
){
	*newtype = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	(*newtype)->count = count;
	(*newtype)->blocklens = (int*)malloc(count*sizeof(int));
	(*newtype)->indices = (MPI_Aint*)malloc(count*sizeof(MPI_Aint));
	(*newtype)->old_types = (MPI_Datatype*)malloc(count*sizeof(MPI_Datatype));
	memcpy((*newtype)->blocklens, blocklens, count*sizeof(int));
	memcpy((*newtype)->indices, indices, count*sizeof(MPI_Aint));
	
	return MPI_SUCCESS;
}
static inline
int MPI_Type_vector( // Creates a vector (strided) datatype 
	int count, // [in] number of blocks (nonnegative integer) 
	int blocklength, // [in] number of elements in each block (nonnegative integer) 
	int stride, //  [in] number of elements between start of each block (integer) 
	MPI_Datatype old_type, // [in] old datatype (handle) 
	MPI_Datatype *newtype_p // [out] new datatype (handle) 
){
	assert(0);
	if(count < 0) return MPI_ERR_INTERN;
	*newtype_p = (MPI_Datatype)malloc(sizeof(struct MPI_Datatype_impl_));
	if(!newtype_p) return MPI_ERR_INTERN;
//	(*newtype)->lb = 0;
	MPI_Aint oldlb;
	MPI_Aint oldextent;
	int s = MPI_Type_get_extent(old_type, &oldlb, &oldextent);
//	(*newtype)
	if(s != MPI_SUCCESS) return MPI_ERR_TYPE;
//	(*newtype)
//	(*newtype)->extent = count*oldextent;
	return MPI_SUCCESS;
}



inline int MPI_Get_address(const void *location, MPI_Aint *address)
{
    return MPI_SUCCESS;
}


// Memory handling 

//#include <stdlib.h>

int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr)
{
    //if (baseptr) *(void **)baseptr = malloc(size);
    if (baseptr) *(void **)baseptr = 0;
    return MPI_SUCCESS;
}

inline int MPI_Free_mem(void *base)
{
    //free(base);
    return MPI_SUCCESS;
}

//struct MPI_Request {
//   bool operator!=(const MPI_Request &o) { return false; }
//};
//static struct MPI_Request MPI_REQUEST_NULL;

int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status);
//{
//    if (flag) *flag = -1; // true
//    return MPI_SUCCESS;
//}

inline int MPI_Cancel(MPI_Request *request)
{
    return MPI_SUCCESS;
}

inline int MPI_Request_free(MPI_Request *request)
{
    return MPI_SUCCESS;
}

inline int MPI_Start(MPI_Request *request)
{
    return MPI_SUCCESS;
}

inline int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])
{
    return MPI_SUCCESS;
}

inline int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount, int array_of_indices[], MPI_Status array_of_statuses[])
{
    if (outcount) *outcount = 0;
    return MPI_SUCCESS;
}

//static MPI_Status* MPI_STATUSES_IGNORE = NULL;

enum {
    MPI_GRAPH,
    MPI_CART,
    MPI_DIST_GRAPH
};

enum {
    MPI_PROC_NULL,
    MPI_ANY_SOURCE,
    MPI_ANY_TAG,
    MPI_UNDEFINED,
    MPI_BSEND_OVERHEAD,
    MPI_KEYVAL_INVALID,
    MPI_LOCK_EXCLUSIVE,
    MPI_LOCK_SHARED,
    MPI_ROOT
}; 

enum {
    MPI_MODE_APPEND,
    MPI_MODE_CREATE,
    MPI_MODE_DELETE_ON_CLOSE,
    MPI_MODE_EXCL,
    MPI_MODE_NOCHECK,
    MPI_MODE_NOPRECEDE,
    MPI_MODE_NOPUT,
    MPI_MODE_NOSTORE,
    MPI_MODE_NOSUCCEED,
    MPI_MODE_RDONLY,
    MPI_MODE_RDWR,
    MPI_MODE_SEQUENTIAL,
    MPI_MODE_UNIQUE_OPEN,
    MPI_MODE_WRONLY
    
};


inline int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{
    return MPI_SUCCESS;
}

int MPI_Comm_test_inter(MPI_Comm comm, int *flag);
//{
//    if (flag) *flag = -1; // -1 = true
//    return MPI_SUCCESS;
//}

inline int MPI_Topo_test(MPI_Comm comm, int *status)
{
    if (status) *status = MPI_UNDEFINED;
    return MPI_SUCCESS;
}

// const on last arg should not be there
inline int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm, const MPI_Comm *newcomm)
{
    return MPI_SUCCESS;
}

// const on last arg should not be there
inline int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm, const MPI_Comm *newcomm)
{
    return MPI_SUCCESS;
}

inline int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    return MPI_SUCCESS;
}



#endif
