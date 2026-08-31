//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/DeepQMC/DeepQMCBridge.h"

extern "C"
{
#include <Python.h>
}

#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifndef DEEPQMC_INFER_BRIDGE_SOURCE_DIR
#define DEEPQMC_INFER_BRIDGE_SOURCE_DIR ""
#endif

namespace qmcplusplus
{
namespace
{
class PyObjectHandle
{
public:
  explicit PyObjectHandle(PyObject* obj = nullptr) : obj_(obj) {}
  ~PyObjectHandle() { Py_XDECREF(obj_); }

  PyObjectHandle(const PyObjectHandle&)            = delete;
  PyObjectHandle& operator=(const PyObjectHandle&) = delete;

  PyObjectHandle(PyObjectHandle&& other) noexcept : obj_(other.obj_) { other.obj_ = nullptr; }
  PyObjectHandle& operator=(PyObjectHandle&& other) noexcept
  {
    if (this != &other)
    {
      Py_XDECREF(obj_);
      obj_       = other.obj_;
      other.obj_ = nullptr;
    }
    return *this;
  }

  PyObject* get() const { return obj_; }
  PyObject* release()
  {
    PyObject* tmp = obj_;
    obj_          = nullptr;
    return tmp;
  }
  void reset()
  {
    Py_XDECREF(obj_);
    obj_ = nullptr;
  }

private:
  PyObject* obj_;
};

std::string formatPythonError()
{
  if (!PyErr_Occurred())
    return "unknown Python error";

  PyObject* type_raw  = nullptr;
  PyObject* value_raw = nullptr;
  PyObject* trace_raw = nullptr;
  PyErr_Fetch(&type_raw, &value_raw, &trace_raw);
  PyErr_NormalizeException(&type_raw, &value_raw, &trace_raw);
  PyObjectHandle type(type_raw);
  PyObjectHandle value(value_raw);
  PyObjectHandle trace(trace_raw);

  PyObjectHandle value_str(PyObject_Str(value.get() ? value.get() : Py_None));
  const char* value_cstr = value_str.get() ? PyUnicode_AsUTF8(value_str.get()) : nullptr;
  if (value_cstr && value_cstr[0] != '\0')
    return value_cstr;
  PyObjectHandle type_str(PyObject_Str(type.get() ? type.get() : Py_None));
  const char* type_cstr = type_str.get() ? PyUnicode_AsUTF8(type_str.get()) : nullptr;
  return type_cstr ? type_cstr : "unknown Python error";
}

[[noreturn]] void throwPythonError(const std::string& context)
{ throw std::runtime_error(context + ": " + formatPythonError()); }

void ensurePythonInitialized()
{
  if (!Py_IsInitialized())
    Py_Initialize();
}

void prependSysPath(const std::string& path)
{
  if (path.empty())
    return;

  PyObject* borrowed_sys_path = PySys_GetObject("path");
  if (!borrowed_sys_path || !PyList_Check(borrowed_sys_path))
    throw std::runtime_error("Python sys.path is not available as a list");

  PyObjectHandle py_path(PyUnicode_FromString(path.c_str()));
  if (!py_path.get())
    throwPythonError("creating Python path string failed");
  if (PyList_Insert(borrowed_sys_path, 0, py_path.get()) != 0)
    throwPythonError("inserting DeepQMC bridge path into sys.path failed");
}

PyObjectHandle makeFloatList(const std::vector<DeepQMCBridge::RealType>& values)
{
  PyObjectHandle list(PyList_New(static_cast<Py_ssize_t>(values.size())));
  if (!list.get())
    throwPythonError("allocating Python list failed");

  for (std::size_t i = 0; i < values.size(); ++i)
  {
    PyObject* item = PyFloat_FromDouble(values[i]);
    if (!item)
      throwPythonError("creating Python float failed");
    // PyList_SET_ITEM steals the new reference.
    PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(i), item);
  }
  return list;
}

std::vector<DeepQMCBridge::RealType> extractFloatList(PyObject* obj, const char* name)
{
  if (!PyList_Check(obj))
  {
    std::ostringstream msg;
    msg << "DeepQMC Python bridge returned non-list for " << name;
    throw std::runtime_error(msg.str());
  }

  const Py_ssize_t size = PyList_Size(obj);
  std::vector<DeepQMCBridge::RealType> values;
  values.reserve(static_cast<std::size_t>(size));
  for (Py_ssize_t i = 0; i < size; ++i)
  {
    PyObject* item     = PyList_GetItem(obj, i); // borrowed
    const double value = PyFloat_AsDouble(item);
    if (PyErr_Occurred())
      throwPythonError(std::string("converting DeepQMC result list '") + name + "' to double failed");
    values.push_back(static_cast<DeepQMCBridge::RealType>(value));
  }
  return values;
}

class PythonDeepQMCBridge : public DeepQMCBridge
{
public:
  PythonDeepQMCBridge(std::string model_path, std::string python_module_path)
      : model_path_(std::move(model_path)), python_module_path_(std::move(python_module_path))
  {
    if (model_path_.empty())
      throw std::runtime_error("DeepQMC wavefunction requires a non-empty model path");

    ensurePythonInitialized();
    PyGILState_STATE gil_state = PyGILState_Ensure();
    try
    {
      prependSysPath(DEEPQMC_INFER_BRIDGE_SOURCE_DIR);
      prependSysPath(python_module_path_);

      PyObject* modules = PyImport_GetModuleDict();
      if (modules && PyDict_DelItemString(modules, "deepqmc_infer_bridge") != 0)
        PyErr_Clear();

      python_module_ = PyObjectHandle(PyImport_ImportModule("deepqmc_infer_bridge"));
      if (!python_module_.get())
        throwPythonError("importing deepqmc_infer_bridge failed");

      PyObjectHandle python_class(PyObject_GetAttrString(python_module_.get(), "DeepQMCInferBridge"));
      if (!python_class.get() || !PyCallable_Check(python_class.get()))
        throwPythonError("loading DeepQMCInferBridge class failed");

      python_instance_ = PyObjectHandle(PyObject_CallFunction(python_class.get(), "s", model_path_.c_str()));
      if (!python_instance_.get())
        throwPythonError("constructing DeepQMCInferBridge failed");
    }
    catch (...)
    {
      python_instance_.reset();
      python_module_.reset();
      PyGILState_Release(gil_state);
      throw;
    }
    PyGILState_Release(gil_state);
  }

  BatchResult evaluateLogBatch(const std::vector<RealType>& ion_coords,
                               const std::vector<RealType>& electron_coords,
                               int mol_idx,
                               int batch_size,
                               int n_elec) const override
  {
    if (batch_size < 0 || n_elec < 0)
      throw std::runtime_error("DeepQMC batch_size and n_elec must be non-negative");
    const std::size_t expected_electron_size =
        static_cast<std::size_t>(batch_size) * static_cast<std::size_t>(n_elec) * OHMMS_DIM;
    if (electron_coords.size() != expected_electron_size)
      throw std::runtime_error("DeepQMC electron coordinate size does not match batch_size * n_elec * 3");

    BatchResult batch_result;
    std::lock_guard<std::mutex> lock(call_mutex_);
    PyGILState_STATE gil_state = PyGILState_Ensure();
    try
    {
      PyObjectHandle py_ions(makeFloatList(ion_coords));
      PyObjectHandle py_electrons(makeFloatList(electron_coords));
      PyObjectHandle result(PyObject_CallMethod(python_instance_.get(), "compute_log_gl", "OOiii", py_ions.get(),
                                                py_electrons.get(), mol_idx, batch_size, n_elec));
      if (!result.get())
        throwPythonError("DeepQMCInferBridge.compute_log_gl failed");
      if (!PyTuple_Check(result.get()) || PyTuple_Size(result.get()) != 3)
        throw std::runtime_error("DeepQMCInferBridge.compute_log_gl must return a 3-tuple");

      batch_result.log_values      = extractFloatList(PyTuple_GetItem(result.get(), 0), "log_values");
      batch_result.grad_log_values = extractFloatList(PyTuple_GetItem(result.get(), 1), "grad_log_values");
      batch_result.lap_log_values  = extractFloatList(PyTuple_GetItem(result.get(), 2), "lap_log_values");
    }
    catch (...)
    {
      PyGILState_Release(gil_state);
      throw;
    }
    PyGILState_Release(gil_state);
    return batch_result;
  }

  ~PythonDeepQMCBridge()
  {
    std::lock_guard<std::mutex> lock(call_mutex_);
    if (Py_IsInitialized())
    {
      PyGILState_STATE gil_state = PyGILState_Ensure();
      python_instance_.reset();
      python_module_.reset();
      PyGILState_Release(gil_state);
    }
    else
    {
      python_instance_.release();
      python_module_.release();
    }
  }

private:
  std::string model_path_;
  std::string python_module_path_;
  PyObjectHandle python_module_;
  PyObjectHandle python_instance_;
  /** Serialize calls into this Python bridge instance.
   *
   * DeepQMCWF clones share one bridge, so independent OpenMP crowds may reach
   * the same Python object concurrently.  The GIL protects Python C API access,
   * but Python extension/runtime code used underneath DeepQMC may release the
   * GIL.  This instance-local mutex avoids reentrant calls into the Python
   * bridge without introducing any thread creation or execution policy; it only
   * protects this external runtime object.
   */
  mutable std::mutex call_mutex_;
};

class UnavailableDeepQMCBridge : public DeepQMCBridge
{
public:
  explicit UnavailableDeepQMCBridge(std::string reason) : reason_(std::move(reason)) {}

  BatchResult evaluateLogBatch(const std::vector<RealType>&, const std::vector<RealType>&, int, int, int) const override
  { throw std::runtime_error("DeepQMC inference bridge is not available: " + reason_); }

private:
  std::string reason_;
};
} // namespace

std::unique_ptr<const DeepQMCBridge> makePythonDeepQMCBridge(std::string model_path, std::string python_module_path)
{ return std::make_unique<PythonDeepQMCBridge>(std::move(model_path), std::move(python_module_path)); }

std::unique_ptr<const DeepQMCBridge> makeUnavailableDeepQMCBridge(std::string reason)
{ return std::make_unique<UnavailableDeepQMCBridge>(std::move(reason)); }

} // namespace qmcplusplus
