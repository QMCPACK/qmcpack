C++ wrapper on top of HDF5 library C interfaces

Users only need `hdf_archive` class to open/close and read/write files.

### I/O Modes of `hdf_archive`

`hdf_archive` supports three I/O modes determined at construction:

1. **Serial mode** (`hdf_archive()`):
   - Standard serial POSIX I/O without MPI-IO property lists.
   - For single-rank runs or when only a single designated rank (e.g., rank 0) accesses the file.
   - Supports both reading and writing.

2. **Parallel independent mode** (`hdf_archive(comm, false)` or `hdf_archive(comm)`):
   - All ranks in `comm` independently open and read the file concurrently.
   - **Predetermined read-only**: Mutating calls (`create()`, opening `H5F_ACC_RDWR`, creating groups in `push()`, `write()`, `append()`, `flush()`) will throw `std::runtime_error`.
   - When parallel HDF5 (`ENABLE_PHDF5`) is available and `comm.size() > 1`, uses the HDF5 MPI-IO driver (`H5Pset_fapl_mpio`) with independent dataset transfer.
   - When parallel HDF5 is not available, gracefully decays to serial POSIX read-only access across all ranks.

3. **Parallel collective mode** (`hdf_archive(comm, true)`):
   - All ranks in `comm` participate collectively in parallel I/O.
   - Supports both reading and writing.
   - Requires parallel HDF5 (`ENABLE_PHDF5` and `comm.size() > 1`). Configures collective metadata operations and collective dataset transfers (`H5FD_MPIO_COLLECTIVE`).
   - Throws `std::runtime_error` upon construction if parallel HDF5 is not available.

`hdf_wrapper_functions.h` wraps raw C functions of HDF5.

`hdf_datatype.h` handles the mapping between C type and HDF5 native type.

`hdf_dataspace.h` handles HDF5 multidimentional dataspace.

`hdf_dataproxy` is a template class to support any kind of datatype written to HDF5 file as a single dataset.
Its specialization are 
STL containers, including vector, bitset and string, in `hdf_stl.h`;
OhmmsPETE containers, including Vector, Matrix and Array, in `hdf_pete.h`;
Afredo's multi container for multidimentional arrays, in `hdf_multi.h`.
Users are recommended to include the corresponding header if a non-STL data container is used.

`hdf_hyperslab.h` supports hyperslab selection in filespace. In production.


`hdf_hyperslab` reads from  and writes into data containers which requires `../type_traits/container_traits.h`
to support features like resizing containers.
`container_traits` has specializations for STL vector in `container_traits.h`,
 OhmmsPETE containers in `container_traits_ohmms.h`
 and Afredo's multidimentional arrays `container_traits_multi.h`
When using `hdf_hyperslab`, users are required to include the corresponding header if a non-STL data container is used.

Although users need to include a few headers to operate a feature with full functionality, it reduces header file entanglement and saves compilation time.

A bit more about multidimensional data. Take a datatype in memory `Matrix<TinyVector<std::complex<double>, 3>>` as an example.
The dataset on the file has a rank of 2 (Matrix) + 1 (TinyVector) + 1 (std::complex) + 0 (double) = 4
The Matrix container holds `TinyVector<std::complex<double>, 3>` value type.
`getShape<TinyVector<std::complex<double>, 3>>` gives the dimensions of a matrix.
