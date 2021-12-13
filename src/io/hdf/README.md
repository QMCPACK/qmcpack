C++ wrapper on top of HDF5 library C interfaces

Users only need `hdf_archive` class to open/close and read/write files.

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

`hdf_double_hyperslab.h` supports hyperslab selection in both filespace and memory space. Not completed yet due to limited demand.

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
