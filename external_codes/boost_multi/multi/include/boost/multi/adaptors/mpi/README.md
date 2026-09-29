<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# Multi MPI Adaptor

_© Alfredo A. Correa, 2025_

(documentation in progress)

MPI libraries are a standard for parallel computing in C and C++.
The Multi MPI Adaptor provides ways to interface arrays with MPI library implementations.
This feature helps use arrays in MPI-based programs by streamlining the communication of array contents between different MPI processes.

The functions in the adaptor do not access the array data, which ensures that it is compatible with GPU-aware implementations of MPI.

## Contents
[[_TOC_]]

## Interfaces

The message-passing interface of MPI generally works with messages to communicate data.
An MPI message is described in 3 parts: a buffer, a count, and a datatype.

For example, the function to send data from one process to another (destination) via a communicator is:

```cpp
int MPI_Send(const void* buffer, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator);
```

MPI library implementations provide dozens of functions like this to send, receive, and process data between MPI processes.
This adaptor doesn't try to replace these functions;
instead, it works by providing a function that generates messages, explicitly calculating the datatypes to describe the arrays.

The MPI user-defined datatypes are generated from the Multi arrays, including their type information and layout information (strides).

The basic usage consists of creating a `multi::mpi::message` object from the array elements.
A message is created by passing a reference to the elements of an array, `multi::mpi::message(my_array.elements())`.
(Elements are not copied in the process.)

The message can then be later decomposed into a buffer, a count, and a datatype for use in the MPI functions.

In this example, which runs in 2 processes, creates an array that is communicated from process 0 to process 1:

```cpp
// compile with `mpic++ -std=c++17 example.cpp -o example.x`
// run with     `mpirun -n 2 example.x`
#include <boost/multi/mpi/adaptors.hpp>

int main() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto const A = multi::array<int, 2>({
		{1, 2, 3},
		{4, 5, 6}
	});

	if(world_rank == 0) {
		auto const& A_msg = multi::mpi::message(A.elements());
		MPI_Send(A_msg.buffer(), A_msg.count(), A_msg.datatype(), 1, 0, MPI_COMM_WORLD);
	} else if(world_rank == 1) {
		multi::array<int, 2> B({2, 3});

		auto&& B_msg = multi::mpi::message(B.elements());
		MPI_Recv(B_msg.buffer(), B_msg.count(), B_msg.datatype(), 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		assert(B == A);
	}
}
```

The examples in this documentation are two-dimensional for simplicity and illustration purposes.
The usage is general to an arbitrary number of dimensions.

## Advanced usage

### Subarrays

A subarray can also be communicated, for example, a small 2x2 block of the original array.
Replacing `message(AA.elements())` with `message(A({0, 2}, {0, 2}).elements())` and `message(B.elements())` with `message(B({0, 2}, {0, 2}).elements())` will result in a communication of a subset of elements.

```cpp
	...
	if(world_rank == 0) {
		auto const& msg = multi::mpi::message(A({0, 2}, {0, 2}).elements());
		MPI_Send(msg.buffer(), msg.count(), msg.datatype(), 1, 0, MPI_COMM_WORLD);
	} else if(world_rank == 1) {
		multi::array<int, 2> B({2, 3});

		auto&& msg = multi::mpi::message(B({0, 2}, {0, 2}).elements());
		MPI_Recv(msg.buffer(), msg.count(), msg.datatype(), 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		assert(B({0, 2}, {0, 2}) == A({0, 2}, {0, 2})); // only the 2x2 block is communicated
	}
```

### Rearrangement of elements

It is essential to understand that, due to the way MPI works, the array's value is not what is communicated but only its fundamental elements (in a canonical order).
We emphasize this detail by passing the `.elements()` range for message construction, not the array per se.

A consequence of this is that the user has to ensure consistency in the shape of the receiving end, as in the previous example.
Communicating a 2x3 array and receiving a 2x2 array will be an error because they have different numbers of elements.

Similarly, a 2x3 array can be communicated into a 3x2 array, although the elements will be rearranged, which is typically not desired.

Still, a reasonable use of rearrangement of elements could involve transposition of the array during communication.
The key to rearranging the elements is that the layouts can be different in different processes.

```cpp
int main() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto const A = multi::array<int, 2>({
		{1, 2, 3},
		{4, 5, 6}
	});

	if(world_rank == 0) {
		auto const& A_msg = multi::mpi::message(A.elements());
		MPI_Send(A_msg.buffer(), A_msg.count(), A_msg.datatype(), 1, 0, MPI_COMM_WORLD);
	} else if(world_rank == 1) {
		multi::array<int, 2> B({3, 2});

		auto&& BT_msg = multi::mpi::message(B.tranposed().elements());
		MPI_Recv(BT_msg.buffer(), BT_msg.count(), BT_msg.datatype(), 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		assert(B == multi::array<int, 2>{
			{1, 4},
			{2, 5},
			{3, 6}
		});
	}
}
```

In this example, the layouts of `A` and the transpose of `B` are different since the elements are arranged differently in memory.
However, they will be equal because the logical element arrangement is the same.
The result is that `B` is the transposed version of `A`.

### Reduction and other operations

Since the message represents a set of elements, reductions and other MPI computations can be used directly element-wise.

```cpp
int main() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto const A = multi::array<int, 2>({
		{1, 2, 3},
		{4, 5, 6}
	});

	auto const& A_msg = multi::mpi::message(A.elements());
	MPI_Reduce(MPI_SUM, A_msg.buffer(), A_msg.count(), A_msg.datatype(), MPI_COMM_WORLD);
}
```

### Iteration and skeletons

If subelements of an array need to be communicated repeatedly, it is wasteful to produce a new message each time.
The key is that all subarrays of a larger array have the same layout.

Suppose we want to communicate the rows of an array in random order.

```cpp
int main() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	multi::array<int, 2> A({
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12}
	});

	if(rank == 0) {
		std::vector<int> perm(A.size()); std::iota(perm.begin(), perm.end(), 0);
		for(auto&& row_index: perm) {
			auto const& msg = multi::mpi::message(A[row_index].elements());
			MPI_Send(msg.buffer(), msg.count(), msg.datatype(), 1, 0, MPI_COMM_WORLD);
		}
	} else if(rank == 1) {
		multi::array<int, 2> B({3, 4});
		for(auto&& row: B) {
			auto&& msg = multi::mpi::message(row.elements());
			MPI_Recv(msg.buffer(), msg.count(), msg.datatype(), 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}
```

The MPI datatype computation is unfortunately repeated for each message inside both loops.
The alternative is to use MPI skeletons, which only contain datatype information.
The loops can be replaced with this code:

```cpp
		...
		auto sk = multi::mpi::skeleton<int>(A.front().layout());
		for(auto&& row_index: perm) {
			MPI_Send(A[row_index].base(), sk.count(), sk.datatype(), 1, 0, MPI_COMM_WORLD);
		}
		...
```

## GPU-aware MPI

Messages and skeletons can be generated for arrays on the GPU.
If the MPI implementation is GPU-aware (e.g., Spectrum MPI), it can communicate the GPU array elements with the same interface.

In the examples above, the arrays can be replaced with `multi::array<int, 2, thrust::cuda::allocator<int> >`, and the communication will automatically use GPU hardware.

## Serialization

The advantage of using datatypes and messages is that data doesn't need to be copied explicitly into a local buffer.
However, in certain cases it is possible that communication using datatypes is slower than compacting the array data into a contiguous buffer.

A completely different alternative to the use of message and datatypes is to use serialization.

```cpp
#include <sstream>

int main() {
	auto const A = multi::array<int, 2>({
		{1, 2, 3},
		{4, 5, 6}
	});

	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0) {
		std::ostringstream oss;
		boost::archive:::binary_oarchive oa(oss);
		oa << A();
		
		MPI_Send(A_msg.str().data(), A_msg.str().data(), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
	} else if(rank == 1) {
		auto B = multi::array<int, 2>({2, 3});

		MPI_Status status; MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
		int count; MPI_Get_count(&status, MPI_CHAR, &count);

		std::string buffer(count);
		MPI_Recv(buffer.data(), bufer.size(), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		std::istringstream iss(std::move(buffer));
		boost::archive::binary_iarchive ia(iss);
		ia >> B();

		assert(B == A);
	}
}
```
