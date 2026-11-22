import subprocess
import sys
from os import PathLike
from pathlib import Path

bin_dir = Path(__file__).parent / "bin"


def run(script_name: PathLike):
    script_path = bin_dir / script_name
    result = subprocess.run([script_path] + sys.argv[1:])
    sys.exit(result.returncode)


def eshdf():
    run("eshdf")


def nxs_redo():
    run("nxs-redo")


def nxs_sim():
    run("nxs-sim")


def nxs_test():
    run("nxs-test")


def qdens():
    run("qdens")


def qdens_radial():
    run("qdens-radial")


def qmca():
    run("qmca")


def qmc_fit():
    run("qmc-fit")
