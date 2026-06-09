import subprocess
import sys
from os import PathLike
from pathlib import Path

bin_dir = Path(__file__).parent/"bin"

def run(script_name: PathLike):
    script_path = bin_dir/script_name
    result = subprocess.run([script_path] + sys.argv[1:])
    sys.exit(result.returncode)
#end def run

def eshdf():
    run("eshdf")
#end def eshdf

def nxs_redo():
    run("nxs-redo")
#end def nxs_redo

def nxs_sim():
    run("nxs-sim")
#end def nxs_sim

def nxs_test():
    run("nxs-test")
#end def nxs_test

def qdens():
    run("qdens")
#end def qdens

def qdens_radial():
    run("qdens-radial")
#end def qdens_radial

def qmca():
    run("qmca")
#end def qmca

def qmc_fit():
    run("qmc-fit")
#end def qmc_fit
