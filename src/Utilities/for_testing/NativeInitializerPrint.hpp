
#include <iostream>

namespace qmcplusplus
{

/** This wrapper is to allow us to leave the user facing operator<< for classes alone
 */
template<typename OBJECT>
class NativePrint
{
public:
  using type_t = OBJECT;
  NativePrint(const OBJECT& obj) : obj_(obj) {}
  const OBJECT& get_obj() const { return obj_; }
private:
  OBJECT obj_;
};

template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const NativePrint<TinyVector<T, D>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for(int i = 0 ; i < D - 1; ++i)
    out << std::setw(12) << std::setprecision(10) << vec[i] << ", ";
  out << std::setw(12) << std::setprecision(10) << vec[D - 1] << " }";
  return out;
}

template<class T>
std::ostream& operator<<(std::ostream& out, const NativePrint<std::vector<T>>& np_vec)
{
  out << "{ ";
  auto vec = np_vec.get_obj();
  for(T& t : vec)
    out << std::setprecision(10) << t << ", ";
  out << " }";
  return out;
}

} // namespace qmcplusplus
