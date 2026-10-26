##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""Calculate memory used by the Nexus host process and its children"""

#===========================================#
# The following is adapted from recipe 8.2  #
# in the Python Cookbook by Jean Brouwers   #
#===========================================#

import os


def _VmB(virt_mem_key: str, pid: int | None = None) -> float:
    """Read ``/proc/<pid>/status`` and get the value on the line with ``virt_mem_key``."""
    if pid is None:
        pid = os.getpid()

    proc_status = f"/proc/{pid}/status"
    try:
        with open(proc_status) as t:
            v = t.read()
    except Exception:  # noqa: BLE001
        return 0.0  # non-Linux?

    # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    try:
        i = v.index(virt_mem_key)
        v = v[i:].split(maxsplit=3)
    except ValueError:
        return 0.0

    if len(v) < 3:
        return 0.0  # invalid format?

    value = float(v[1])
    units = v[2].lower()
    if units == "kb":
        value = value * 1024.0
    elif units == "mb":
        value = value * 1024.0 * 1024.0

    return value


def get_children(pid: int) -> list | list[int]:
    """Read ``/proc/<pid>/task/<pid>/children`` and return a list of children"""
    proc_children = f"/proc/{pid}/task/{pid}/children"
    try:
        with open(proc_children,'r') as t:
            v = t.read()
    except Exception:  # noqa: BLE001
        return []

    children = [int(c) for c in v.split()]
    return children


def all_children(pid: int | None = None) -> list | list[int]:
    """Recursive function to get all children of a process."""
    if not pid:
        pid = os.getpid()

    child_list = get_children(pid)
    all_list = child_list[:]
    for child_pid in child_list:
        child2 = all_children(child_pid)
        all_list.extend(child2)
    return all_list


def memory(since: float = 0.0, *, children: bool = False) -> tuple[float, float]:
    """Return memory usage in bytes."""
    self_mem = _VmB("VmSize:") - since
    child_mem = 0.0
    if children:
        for child_pid in all_children():
            child_mem += _VmB("VmSize:", pid=child_pid)
    return self_mem, child_mem


def resident(since: float = 0.0, *, children: bool = False) -> tuple[float, float]:
    """Return resident memory usage in bytes."""
    self_mem = _VmB("VmRSS:") - since
    child_mem = 0.0
    if children:
        for child_pid in all_children():
            child_mem += _VmB("VmRSS:", pid=child_pid)
    return self_mem, child_mem


def stacksize(since: float = 0.0, *, children: bool = False) -> tuple[float, float]:
    """Return stack size in bytes."""
    self_mem = _VmB("VmStk:") - since
    child_mem = 0.0
    if children:
        for child_pid in all_children():
            child_mem += _VmB("VmStk:", pid=child_pid)
    return self_mem, child_mem
