#ifndef KERNEL_BUFFER_HELPER_HPP
#define KERNEL_BUFFER_HELPER_HPP
template <typename T>
inline __device__ T* shared_memory_proxy()
{
    extern __shared__ unsigned char memory[];
    return reinterpret_cast<T*>(memory);
}
#endif
