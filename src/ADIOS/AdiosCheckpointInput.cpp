//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <ADIOS/AdiosCheckpointInput.h>
#include <Configuration.h>
#include <adios_read.h>

namespace ADIOS
{

AdiosCheckpointInput::AdiosCheckpointInput(std::string& checkpoint_file)
{
  adios_read_init(method, OHMMS::Controller->getMPI(), "verbose=3;abort_on_error");
  adios_file_handle = adios_read_open_file(checkpoint_file.c_str(), method,
                      comm);
  if (adios_file_handle == NULL)
  {
    qmcplusplus::app_error() << adios_errmsg() << std::endl;
  }
  //Select for local values
  sel = adios_selection_writeblock (OHMMS::Controller->rank());
}

AdiosCheckpointInput::~AdiosCheckpointInput()
{
  adios_read_close(adios_file_handle);
  adios_read_finalize_method(method);
  for (int i = 0; i < buffers.size(); i++)
    free(buffers[i]);
}

void AdiosCheckpointInput::queueRead(const std::string& var_name)
{
  read_queue.push_back(var_name);
}

void AdiosCheckpointInput::performReads()
{
  for (int i = 0; i < read_queue.size(); i++)
  {
    const char* var_name = read_queue[i].c_str();
    //grab the information from the metadata about this variable
    ADIOS_VARINFO* adios_inq = adios_inq_var(adios_file_handle, var_name);
    //Calculate the total size of the vector based on the metadata
    uint64_t total_size = adios_type_size(adios_inq->type, NULL);
    for (int i = 0; i < adios_inq->ndim; i++)
      total_size *= adios_inq->dims[i];
    //allocate and then store for when retrieve
    void* buff = malloc(total_size);
    buffers.push_back(buff);
    //schedule the read
    adios_schedule_read(adios_file_handle, sel, var_name, 0, 1, buff);
    //Don't need the information about the variable anymore
    adios_free_varinfo(adios_inq);
  }
  //perform the read
  adios_perform_reads(adios_file_handle, 1);
}

template <typename T>
std::vector<T> AdiosCheckpointInput::retrieveVector(const std::string& var_name)
{
  /*
    Retrieve one of the buffers that we have read from the file system.
  */
  void* buff_to_return;
  int buff_size;
  bool buff_found = false;
  for (int i = 0; i < read_queue.size(); i++)
  {
    if (read_queue[i] == var_name)
    {
      buff_found = true;
      buff_to_return = buffers[i];
      buff_size = sizes[i];
      //delete this value from the queues
      buffers.erase(buffers.begin() + i);
      sizes.erase(sizes.begin() + i);
      read_queue.erase(read_queue.begin() + i);
      break;
    }
  }
  if (!buff_found)
  {
    qmcplusplus::app_error() << "Checkpoint file does not contain: " << var_name << std::endl;
    APP_ABORT("qmcapp");
    return ;
  }
  T* raw_buff = static_cast<T*>(buff_to_return);
  std::vector<T> vec_buff(raw_buff, raw_buff + buff_size);
  //free some memory
  free(raw_buff);
  return vec_buff;
}

void AdiosCheckpointInput::clear()
{
  read_queue.clear();
  sizes.clear();
  for (int i = 0; i < buffers.size(); i++)
  {
    free(buffers[i]);
  }
  sizes.clear();
}

template <>
std::string AdiosCheckpointInput::getScalar(const std::string& var_name)
{
  //Scalars are all stored in the metadata so we can read them without disk access
  ADIOS_VARINFO* adios_inq = adios_inq_var(adios_file_handle, var_name.c_str());
  int string_size = adios_type_size(adios_inq->type, adios_inq->value);
  std::string string_val = static_cast<char*>(*adios_inq->value);
  adios_free_varinfo(adios_inq);
  return string_val;
}


template <typename T>
T AdiosCheckpointInput::getScalar(const std::string& var_name)
{
  //Scalars are all stored in the metadata so we can read them without disk access
  ADIOS_VARINFO* adios_inq = adios_inq_var(adios_file_handle, var_name.c_str());
  if (adios_type_size(adios_inq->type, adios_inq->value) != sizeof(T))
  {
    qmcplusplus::app_error() << "Data type does not match data type found in file: " << std::endl;
    return ;
  }
  T scalar_value = static_cast<T>(*adios_inq->value);
  adios_free_varinfo(adios_inq);
  return scalar_value;
}



template <typename T>
void AdiosCheckpointInput::getVector(const std::string& var_name, std::vector<T>& buffer)
{
  ADIOS_VARINFO* adios_inq = adios_inq_var(adios_file_handle, var_name.c_str());
  //Check to make sure the T is the same size as the data type on the disk
  if (adios_type_size(adios_inq->type, NULL) != sizeof(T))
  {
    qmcplusplus::app_error() << "Data type does not match data type found in file: " << std::endl;
    APP_ABORT("qmcapp");
    return ;
  }
  uint64_t total_size = 1;
  for (int i = 0; i < adios_inq->ndim; i++)
    total_size *= adios_inq->dims[i];
  //make sure we have enough space to copy all the data from the file
  buffer.reserve(total_size);
  //schedule the read
  adios_schedule_read(adios_file_handle, sel, var_name.c_str(), 0, 1, &(buffer[0]));
  //perform the read
  adios_perform_reads(adios_file_handle, 1);
  //Don't need the information about the variable anymore
  adios_free_varinfo(adios_inq);
}
}
