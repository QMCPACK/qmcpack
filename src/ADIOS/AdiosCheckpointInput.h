#include <adios_read.h>
#include <string>
#include <vector>
#include <Message/Communicate.h>

namespace ADIOS
{


class AdiosCheckpointInput
{

  ADIOS_SELECTION* sel;
  ADIOS_FILE* adios_file_handle;
  std::vector<std::string> read_queue;
  std::vector<void*> buffers;
  std::vector<int> sizes;


public:
  const enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;


  AdiosCheckpointInput(std::string&);
  ~AdiosCheckpointInput();

  template <typename T>
  void getVector(const std::string& var_name, std::vector<T>& buffer);

  template <typename T>
  T getScalar(const std::string& var_name);

  template <>
  std::string getScalar(const std::string& var_name);

  void clear();

  template <typename T>
  std::vector<T> retrieveVector(const std::string& var_name);

  void performReads();

  void queueRead(const std::string& var_name);

};

}
