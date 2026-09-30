#!/usr/bin/env python3

"""Helper to launch a CTest GPU test on the GPU assigned by CTest's resource scheduler.

The first argument names the backend-specific device visibility variable, such
as ROCR_VISIBLE_DEVICES or CUDA_VISIBLE_DEVICES. The launcher reads CTest's GPU
resource assignment, maps its logical device ID through any visibility mask
provided by a job scheduler, restricts the child process to that device, and
then replaces itself with the requested test command.
"""

import os
import re
import sys


def get_gpu_id():
    """Return the GPU ID assigned to resource group zero by CTest."""
    assignment = os.environ.get("CTEST_RESOURCE_GROUP_0_GPUS")
    if assignment is None:
        sys.stderr.write(
            "GPU test has no CTest resource assignment. "
            "Configure with QMC_CTEST_NUM_GPUS and run the test through ctest.\n"
        )
        sys.exit(1)

    match = re.search(r"(?:^|;)id:([^,;]+),slots:[0-9]+(?:;|$)", assignment)
    if match is None:
        sys.stderr.write(f"Invalid CTest GPU resource assignment: {assignment!r}\n")
        sys.exit(1)
    return match.group(1)


def map_visible_gpu(gpu_id, visibility_variable):
    """Map a CTest logical ID through a visibility mask from a job launcher."""
    gpu_index = int(gpu_id)
    visible_devices = os.environ.get(visibility_variable)
    if visible_devices:
        device_ids = [device.strip() for device in visible_devices.split(",")]
        if gpu_index >= len(device_ids):
            sys.stderr.write(
                f"CTest assigned GPU {gpu_id}, but {visibility_variable} only exposes "
                f"{len(device_ids)} device(s).\n"
            )
            sys.exit(1)
        return device_ids[gpu_index]
    return gpu_id


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("GPU test resource launcher requires a visibility variable and command to run.\n")
        return 1

    visibility_variable = sys.argv[1]
    gpu_id = map_visible_gpu(get_gpu_id(), visibility_variable)
    os.environ[visibility_variable] = gpu_id

    os.execvpe(sys.argv[2], sys.argv[2:], os.environ)


if __name__ == "__main__":
    sys.exit(main())
