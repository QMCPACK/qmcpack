<!--- &2>/dev/null
    xdg-open $0; exit
--->
[comment]: # (Comment)

# [Boost].MPI3
*Alfredo A. Correa*  
<alfredo.correa@gmail.com>

[Boost].MPI3 is a C++ library wrapper for standard MPI3.

[Boost].MPI3 is not an official Boost library.
However Boost.MPI3 is designed following the principles of Boost and the STL.

[Boost].MPI3 is not a derivative of Boost.MPI and it is unrelated to the, now deprecated, official MPI-C++ interface.
It adds features which were missing in Boost.MPI (which only covers MPI-1), with an iterator-based interface and MPI-3 features (RMA and Shared memory).
[Boost].MPI3 is written from scratch in C++14.

[Boost].MPI3 depends and has been compiled against Boost +1.53 and one of the MPI library implementations, OpenMPI +1.9, MPICH +3.2.1 or MVAPICH, using the following compilers gcc +5.4.1, clang +6.0, PGI 18.04. 
The current version of the library (wrapper) is `0.71`, (programmatically accesible from `./version.hpp`).

## Introduction

MPI is a large library for run-time parallelism where several paradigms coexist.
It was is originally designed as standardized and portable message-passing system to work on a wide variety of parallel computing architectures.

The last standard, MPI-3, uses a combination of techniques to achieve parallelism, Message Passing (MP), (Remote Memory Access (RMA) and Shared Memory (SM).
We try here to give a uniform interface and abstractions for these features by means of wrapper function calls and concepts brought familiar to C++ and the STL.

## Motivation: The problem with the standard interface 

A typical C-call for MP looks like this,

```c++
int status_send = MPI_Send(&numbers, 10, MPI_INT, 1, 0, MPI_COMM_WORLD);
assert(status_send == MPI_SUCCESS);
... // concurrently with 
int status_recv = MPI_Recv(&numbers, 10, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
assert(status_recv == MPI_SUCCESS);
```

In principle this call can be made from a C++ program. 
However there are obvious drawbacks from using this standard interface.

Here we enumerate some of problems,

* Function calls have many arguments (e.g. 6 or 7 arguments in average)
* Many mandatory arguments are redundant or could easily have a default natural value ( e.g. message tags are not always necessary).
* Use of raw pointers and sizes, (e.g. `&number` and `1`)
* Argument is type-erased by `void*`.
* Only primitive types (e.g. `MPI_INT`) can be passed.
* Consistency between pointer types and data-types is responsibility of the user.
* Only contiguous memory blocks can be used with this interface.
* Error codes are stored and had to be checked after each function call.
* Use of handles (such as `MPI_COMM_WORLD`), handles do not have a well defined semantics.

A call of this type would be an improvement:

```c++
world.send(numbers.begin(), numbers.end(), 1);
... // concurrently with 
world.receive(numbers.begin(), numbers.end(), 0); 
```

For other examples, see here: [http://mpitutorial.com/tutorials/mpi-send-and-receive/](http://mpitutorial.com/tutorials/mpi-send-and-receive/)

MPI used to ship with a C++-style interfaces.
It turns out that this interface was a very minimal change over the C version, and for good reasons it was dropped.

The Boost.MPI3 library was designed to use simultaneously (interleaved) with the standard C interface of MPI. 
In this way, changes to existing code can be made incrementally.
Mixing the standard C interface with the Boost.MPI3 is not complicated but requires more knowledge of the library internals than the one provided in this document.

## Installation

The library is "header-only"; no separate compilation is necessary.
Most functions are inline or template functions.
In order to compile it requires an MPI distribution (e.g. OpenMPI or MPICH2) and the corresponding compiler-wrapper (`mpic++` or `mpicxx`).
Currently the library requieres C++14 (usually activated with the compiler option `-std=c++14`) and Boost. In particular it depends on Boost.Serialization and may require linking to this library if values passed are not basic types (`-lboost_serialization`). A typical compilation/run command looks like this:

```bash
$ mpic++ -std=c++14 -O3 mpi3/test/communicator_send.cpp -o communicator_send.x -lboost_serialization
$ mpirun -n 8 ./communicator_send.x
```

In a system such as Red Hat, the dependencies can by installed by

```bash
$ dnf install gcc-c++ boost-devel openmpi-devel mpich-devel
```

The library is tested frequently against `openmpi` and `mpich`, and less frequently with `mvapich2`.

## Testing

The library has a basic `ctest` based testing system.

```c++
cd mpi3/test
mkdir build; cd build
cmake .. && make && ctest
```

## Initialization

Like MPI, Boost.MPI3 requires some global library initialization.
The library includes `mpi3/main.hpp` which wraps around this initialization steps and *simulates* a main function. 
In this way, a parallel program looks very much like normal programs, except that the main function has a third argument with the default global communicator passed in.

```c++
#include "mpi3/version.hpp"
#include "mpi3/main.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3; 
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){
	if(world.rank() == 0) cout << mpi3::version() << '\n';
	return 0;
}
```

Here `world` is a communicator object that is a wrapper over MPI communicator handle.

Changing the `main` program to this syntax in existing code can be too intrusive. 
For this reason a more traditional initialization is also possible.
The alternative initialization is done by instantiating the `mpi3::environment` object (from with the global communicator `.world()` is extracted).

```c++
#include "mpi3/environment.hpp"
int main(int argc, char** argv){
	mpi3::environment env(argc, argv);
	auto world = env.world(); // communicator is extracted from the environment 
    // ... code here
	return 0;
}
```

## Communicators

In the last example, `world` is a global communicator (not necessarely the same as `MPI_COMM_WORLD`, but a copy of it).
There is no global communicator variable `world` that can be accessed directly in a nested function.
The idea behind this is to avoid using the global communicators in nested functions of the program unless they are explicitly passed in the function call.
Communicators are usually passed by reference to nested functions.
Even in traditional MPI it is a mistake to assume that the `MPI_COMM_WORLD` is the only available communicator.

`mpi3::communicator` represent communicators with value-semantics.
This means that `mpi3::communicator` can be copied or passed by reference.
A communicator and their copies are different entities that compare equal.
Communicators can be empty, in a state that is analogous to `MPI_COMM_NULL` but with proper value semantics.

Like in MPI communicators can be duplicated (copied into a new instance) or split.
They can be also compared. 

```c++
mpi3::communicator world2 = world;
assert( world2 == world );
mpi3::communicator hemisphere = world/2;
mpi3::communicator interleaved = world%2;
```

This program for example splits the global communicator in two sub-communicators one of size 2 (including process 0 and 1) and one with size 6 (including 2, 3, ... 7);

```c++
#include "mpi3/main.hpp"
#include "mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){
    assert(world.size() == 8); // this program can only be run in 8 processes
    mpi3::communicator comm = (world <= 1);
    assert(!comm || (comm && comm.size() == 2));
    return 0;
}
```

Communicators give also index access to individual `mpi3::processes` ranging from `0` to `comm.size()`. 
For example, `world[0]` referrers to process 0 or the global communicator.
An `mpi3::process` is simply a rank inside a communicator.
This concept doesn't exist explicit in the standard C interface, but it simplifies the syntax for message passing.

Splitting communicators can be done more traditionally via the `communicator::split` member function. 

Communicators are used to pass messages and to create memory windows.
A special type of communicator is a shared-communicator `mpi3::shared_communicator`.

## Message Passing

This section describes the features related to the message passing (MP) functions in the MPI library.
In C-MPI information is passed via pointers to memory.
This is expected in a C-based interface and it is also very efficient.
In Boost.MPI, information is passed exclusively by value semantics. 
Although there are optimizations that amortize the cost, we decided to generalize the pointer interface and leave the value-based message passing for a higher-level syntax. 

Here we replicate the design of STL to process information, that is, aggregated data is passed mainly via iterators. (Pointer is a type of iterator).

For example in STL data is copied between ranges in this way.
```c++
std::copy(origin.begin(), origin.end(), destination.begin());
```

The caller of function copy doesn't need to worry about he type of the `origin` and `destination` containers, it can mix pointers and iterators and the function doesn't need more redundant information than the information passed. 
The programmer is responsible for managing the memory and making sure that design is such that the algorithm can access the data referred by the passed iterators.

Contiguous iterators (to built-in types) are particularity efficient because they can be mapped to pointers at compile time. This in turn is translated into a MPI primitive function call.
The interface for other type of iterators or contiguous iterators to non-build-in type are simulated, mainly via buffers and serialization.
The idea behind this is that generic message passing function calls can be made to work with arbitrary data types.

The main interface for message passing in Boost.MPI3 are member functions of the communicator.
For example `communicator::send`, `::receive` and `::barrier`. 
The functions `::rank` and `::size` allows each process to determine their unique identity inside the communicator.

```c++
int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
    assert(world.size() == 2);
	if(world.rank() == 0){
	   std::vector<double> v = {1.,2.,3.};
	   world.send(v.begin(), v.end(), 1); // send to rank 1
	}else if(world.rank() == 1){
	   std::vector<double> v(3);
	   world.receive(v.begin(), v.end(), 0); // receive from rank 1
	   assert( v == std::vector{1.,2.,3.} );
	}
	world.barrier(); // synchronize execution here
	return 0;
}
```

Other important functions are `::gather`, `::broadcast` and `::accumulate`. 
This syntax has a more or less obvious (but simplified) mapping to the standard C-MPI interface.
In Boost.MPI3 however all, these functions have reasonable defaults that make the function call shorted and less prone to errors and with the C-MPI interface.

For more examples, look into `./mpi3/tests/`, `./mpi3/examples/` and `./mpi3/exercises/`.

The interface described above is iterator based and is a direct generalization of the C-interface which works with pointers.
If the iterators are contiguous and the associated value types are primitive MPI types, the function is directly mapped to the C-MPI call.

Alternatively, value-based interface can be used.
We will show the terse syntax, using the process objects.

```c++
int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
    assert(world.size() == 2);
	if(world.rank() == 0){
	   double v = 5.;
	   world[1] << v;
	}else if(world.rank() == 1){
	   double v = -1.;
	   world[0] >> v;
	   assert(v == 5.);
	}
	return 0;
}
```

## Remote Memory Access

Remote Memory (RM) is handled by `mpi3::window` objects. 
`mpi3::window`s are created by `mpi3::communicator` via a collective (member) functions.
Since `mpi3::window`s represent memory, it cannot be copied (but can be moved). 

```c++
mpi3::window w = world.make_window(begin, end);
```

Just like in the MPI interface, local access and remote access is synchronized by a `window::fence` call.
Read and write remote access is performed via put and get functions.

```c++
w.fence();
w.put(begin, end, rank);
w.fence();
```

This is minimal example using `put` and `get` functions.

```c++
#include "mpi3/main.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3; using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<double> darr(world.rank()?0:100);
	mpi3::window<double> w = world.make_window(darr.data(), darr.size());
	w.fence();
	if(world.rank() == 0){
		std::vector<double> a = {5., 6.};
		w.put(a.begin(), a.end(), 0);
	}
	world.barrier();
	w.fence();
	std::vector<double> b(2);
	w.get(b.begin(), b.end(), 0);
	w.fence();
	assert( b[0] == 5.);
	world.barrier();

	return 0;
}
```

In this example, memory from process 0 is shared across the communicator, and accessible through a common window.
Process 0 writes (`window::put`s) values in the memory (this can be done locally or remotely). 
Later all processes read from this memory. 
`put` and `get` functions take at least 3 arguments (and at most 4).
The first two is a range of iterators, while the third is the destination/source process rank (called "target_rank"). 

Relevant examples and test are located in For more examples, look into `./mpi3/tests/`, `./mpi3/examples/` and `./mpi3/exercises/`.

`mpi3::window`s may carry type information (as `mpi3::window<double>`) or not (`mpi3::window<>`)

## Shared Memory

Shared memory (SM) uses the underlying capability of the operating system to share memory from process within the same node. 
Historically shared memory has an interface similar to that of remove access.
Only communicators that comprise a single node can be used to create a share memory window.
A special type of communicator can be created by splitting a given communicator.

`mpi3::shared_communicator node = world.split_shared();`

If the job is launched in single node, `node` will be equal (congruent) to `world`.
Otherwise the global communicator will be split into a number of (shared) communicators equal to the number of nodes.

`mpi3::shared_communicator`s can create `mpi3::shared_window`s. 
These are special type of memory windows.

```c++
#include "mpi3/main.hpp"

namespace mpi3 = boost::mpi3; using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::shared_communicator node = world.split_shared();
	mpi3::shared_window<int> win = node.make_shared_window<int>(node.rank()==0?1:0);

	assert(win.base() != nullptr and win.size<int>() == 1);

	win.lock_all();
	if(node.rank()==0) *win.base<int>(0) = 42;
	for (int j=1; j != node.size(); ++j){
		if(node.rank()==0) node.send_n((int*)nullptr, 0, j);//, 666);
	    else if(node.rank()==j) node.receive_n((int*)nullptr, 0, 0);//, 666);
	}
	win.sync();

	int l = *win.base<int>(0);
	win.unlock_all();

	int minmax[2] = {-l,l};
	node.all_reduce_n(&minmax[0], 2, mpi3::max<>{});
	assert( -minmax[0] == minmax[1] );
	cout << "proc " << node.rank() << " " << l << std::endl;

	return 0;
}
```

For more examples, look into `./mpi3/tests/`, `./mpi3/examples/` and `./mpi3/exercises/`.

# Beyond MP, RMA and SHM

MPI provides a very low level abstraction to inter-process communication.
Higher level of abstractions can be constructed on top of MPI and by using the wrapper the works is simplified considerably.

## Mutex

Mutexes can be implemented fairly simply on top of RMA.
Mutexes are used similarly than in threaded code, 
it prevents certain blocks of code to be executed by more than one process (rank) at a time.

```c++
#include "mpi3/main.hpp"
#include "mpi3/mutex.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3; using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::mutex m(world);
	{
		m.lock();
		cout << "locked from " << world.rank() << '\n';
		cout << "never interleaved " << world.rank() << '\n';
		cout << "forever blocked " << world.rank() << '\n';
		cout << std::endl;
		m.unlock();
	}
	return 0;
}
```

(Recursive mutexes are not implemented yet)

Mutexes themselves can be used to implement atomic operations on data.

# Ongoing work

We are implementing memory allocators for remote memory, atomic classes and asynchronous remote function calls.
Higher abstractions and use patterns will be implemented, specially those that fit into the patterns of the STL algorithms and containers.

# Conclusion

The goal is to provide a type-safe, efficient, generic interface for MPI.
We achieve this by leveraging template code and classes that C++ provides.
Typical low-level use patterns become extremely simple, and that exposes higher-level patterns. 
