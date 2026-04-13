#! /usr/bin/env python3

try:
    import pytest
    pytest_version = tuple(map(int, pytest.__version__.split(".")))
    if pytest_version >= (6, 2, 4):
        print("True")
    else:
        print("False")
except ImportError:
    print("False")
