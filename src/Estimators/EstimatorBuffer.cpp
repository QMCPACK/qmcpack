#include "EstimatorBuffer.hpp"

namespace qmcplusplus
{

// template<typename T>
// EstimatorBuffer<T>::EstimatorBuffer(std::size_t size) : data_(size) {}

template<typename T>
EstimatorBuffer<T>::EstimatorBuffer(std::size_t size) : data_(size) {}

template<typename T>
void EstimatorBuffer<T>::reserve(std::size_t size) { data_.reserve(size); }

}
