#!/bin/bash

echo --- START environment setup `date`

command -v spack >/dev/null 2>&1 || { echo "Error: spack not found on PATH." >&2; exit 1; }

# serial   : single install
# 8up      : 8 installs
# par48    : install -j 48 
# makefile : make -j

parallelmode=par48
#parallelmode=makefile

install_environment () {
case "$parallelmode" in
    serial )
	echo --- Serial install
	spack install
	;;
    par48 )
	echo --- spack install -j 48
	spack install -j 48
	;;
    8up )
	echo --- Running 8 installs simultaneously
        spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4 &
	sleep 1; spack install -j 4
	;;
    makefile )
	echo --- Install via parallel make
	spack concretize
	spack env depfile >Makefile
	make -j 48 SPACK_COLOR=always --output-sync=recurse
	;;
    * )
	echo Unknown parallelmode
	exit 1
	;;
esac
}

here=`pwd`

if [ -e `dirname "$0"`/ornl_versions.sh ]; then
    echo --- Contents of ornl_versions.sh
    cat `dirname "$0"`/ornl_versions.sh
    echo --- End of contents of ornl_versions.sh
    source `dirname "$0"`/ornl_versions.sh
else
    echo Did not find version numbers script ornl_versions.sh
    exit 1
fi

plat=`lscpu|grep Vendor|sed 's/.*ID:[ ]*//'`
case "$plat" in
    GenuineIntel )
	ourplatform=Intel
	;;
    AuthenticAMD )
	ourplatform=AMD
	;;
    * )
	# ARM support should be trivial, but is not yet done
	echo Unknown platform
	exit 1
	;;
esac

echo --- Installing for $ourplatform architecture
ourhostname=`hostname|sed 's/\..*//g'`
echo --- Host is $ourhostname
echo --- Using spack `which spack`
echo --- SPACK_USER_CONFIG_PATH=$SPACK_USER_CONFIG_PATH

for havempi in mpi nompi
do	       

theenv=envgccnew${havempi}
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:true"
spack env activate $theenv

spack add gcc@${gcc_vnew}
spack add hwloc
spack add git
spack add ninja
spack add cmake@${cmake_vnew}^gcc@${gcc_vnew}
spack add libxml2
spack add boost@${boost_vnew}^gcc@${gcc_vnew}
spack add util-linux-uuid^gcc@${gcc_vnew}
spack add python^gcc@${gcc_vnew}
if [ "$havempi" == "mpi" ]; then
    spack add openmpi@${ompi_vnew}^gcc@${gcc_vnew}
fi

if [ "$havempi" == "mpi" ]; then
    spack add hdf5@${hdf5_vnew} +fortran +hl +mpi ^gcc@${gcc_vnew}
else
    spack add hdf5@${hdf5_vnew} +fortran +hl ~mpi ^gcc@${gcc_vnew}
fi

spack add fftw@${fftw_vnew} -mpi ^gcc@${gcc_vnew} #Avoid MPI for simplicity
spack add openblas threads=openmp ^gcc@${gcc_vnew} 
#PK: amdlibflame@5.2 will fail with gcc 15.2.0 as of 2026-03-24.f Link error
#if [ "$ourplatform" == "AMD" ]; then
#spack add amdblis; spack add amdlibflame; #spack add amd-aocl
#fi
if [ "$ourplatform" == "Intel" ]; then
spack add intel-oneapi-mkl
fi

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
if [ "$havempi" == "mpi" ]; then
    spack add py-mpi4py
fi

spack add py-numpy@${numpy_vnew}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}

if [ "$havempi" == "mpi" ]; then
# Complete install only for MPI environments:
spack add quantum-espresso +mpi +qmcpack
export CMAKE_BUILD_PARALLEL_LEVEL=8 # For PySCF
#spack add py-pyscf  # Does not work with gccnew=15.2
#spack add dftd4
#spack add rmgdft@develop # Use develop version to avoid vendored SCALPACK compilation bugs (20250324)
#spack add rmgdft

#Luxury options for actual science use:
spack add py-requests # for pseudo helper
spack add py-ase      # full Atomic Simulation Environment
spack add libffi
spack add graphviz +pangocairo # NEXUS requires optional PNG support in dot
spack add py-pydot    # NEXUS optional

#spack add py-spglib   # NEXUS optional  Forces numpy<2 currently + scikit-build issue
#spack add py-seekpath # NEXUS optional
#spack add py-pycifrw  # NEXUS optional
#NOT IN SPACK spack add py-cif2cell # NEXUS optional
fi

install_environment
unset CMAKE_BUILD_PARALLEL_LEVEL
spack env deactivate

done

for havempi in mpi nompi
do	       

theenv=envgccold${havempi}
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:true"
spack env activate $theenv

spack add gcc@${gcc_vold}
spack add hwloc
spack add git
spack add ninja
spack add cmake@${cmake_vold}^gcc@${gcc_vold}
spack add libxml2
spack add boost@${boost_vold}^gcc@${gcc_vold}
spack add util-linux-uuid^gcc@${gcc_vold}
spack add python^gcc@${gcc_vold}
if [ "$havempi" == "mpi" ]; then
    spack add openmpi@${ompi_vnew}^gcc@${gcc_vold}
fi

if [ "$havempi" == "mpi" ]; then
    spack add hdf5@${hdf5_vold} +fortran +hl +mpi ^gcc@${gcc_vold}
else
    spack add hdf5@${hdf5_vold} +fortran +hl ~mpi ^gcc@${gcc_vold}
fi

spack add fftw@${fftw_vold} -mpi ^gcc@${gcc_vold} #Avoid MPI for simplicity
spack add openblas threads=openmp ^gcc@${gcc_vold}
if [ "$ourplatform" == "AMD" ]; then
spack add amdblis; spack add amdlibflame; #spack add amd-aocl
fi
if [ "$ourplatform" == "Intel" ]; then
spack add intel-oneapi-mkl
fi

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
if [ "$havempi" == "mpi" ]; then
    spack add py-mpi4py
fi

spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}

if [ "$havempi" == "mpi" ]; then
# Complete install only for MPI environments:
spack add quantum-espresso +mpi +qmcpack
export CMAKE_BUILD_PARALLEL_LEVEL=8 # For PySCF
spack add py-pyscf
#spack add dftd4
#spack add rmgdft@develop # Use develop version to avoid vendored SCALPACK compilation bugs (20250324)
#spack add rmgdft

#Luxury options for actual science use:
spack add py-requests # for pseudo helper
spack add py-ase      # full Atomic Simulation Environment
spack add libffi
spack add graphviz +pangocairo # NEXUS requires optional PNG support in dot
spack add py-pydot    # NEXUS optional

#spack add py-spglib   # NEXUS optional 
#spack add py-seekpath # NEXUS optional
#spack add py-pycifrw  # NEXUS optional
#NOT IN SPACK spack add py-cif2cell # NEXUS optional
fi

install_environment
unset CMAKE_BUILD_PARALLEL_LEVEL
spack env deactivate
done

for havempi in mpi nompi
do	    

theenv=envclangnew${havempi}
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:true"
spack env activate $theenv

spack add llvm@${llvm_vnew}
#spack add gcc@${gcc_vold}
spack add hwloc
spack add git
spack add ninja
#spack add cmake@${cmake_vnew}^gcc@${gcc_vold}
spack add cmake@${cmake_vnew}
spack add libxml2
#spack add boost@${boost_vnew}^gcc@${gcc_vold}
spack add boost@${boost_vnew}
#spack add util-linux-uuid^gcc@${gcc_vold}
spack add util-linux-uuid
#spack add python^gcc@${gcc_vold}
spack add python
if [ "$havempi" == "mpi" ]; then
#spack add openmpi@${ompi_vnew}^gcc@${gcc_vold}
spack add openmpi@${ompi_vnew}
fi

if [ "$havempi" == "mpi" ]; then
#spack add hdf5@${hdf5_vnew} +fortran +hl +mpi ^gcc@${gcc_vold}
spack add hdf5@${hdf5_vnew} +fortran +hl +mpi 
else
#spack add hdf5@${hdf5_vnew} +fortran +hl ~mpi ^gcc@${gcc_vold}
spack add hdf5@${hdf5_vnew} +fortran +hl ~mpi 
fi
#spack add fftw@${fftw_vnew} -mpi ^gcc@${gcc_vold} #Avoid MPI for simplicity
spack add fftw@${fftw_vnew} -mpi  #Avoid MPI for simplicity
#spack add openblas threads=openmp ^gcc@${gcc_vold}
spack add openblas threads=openmp 
#PK: amdlibflame@5.2 will fail as of 2026-03-24
#if [ "$ourplatform" == "AMD" ]; then
#spack add amdblis; spack add amdlibflame; #spack add amd-aocl
#fi
if [ "$ourplatform" == "Intel" ]; then
spack add intel-oneapi-mkl
fi
 
spack add py-lxml
spack add py-matplotlib
spack add py-pandas
if [ "$havempi" == "mpi" ]; then
	spack add py-mpi4py
fi
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}
install_environment
spack env deactivate
done

# Build LLVM offload with preferred GCC since CUDA may not support new GCC
# Build with new CMake
# TO DO: Match chosen cuda with version installed on system
for havempi in mpi nompi
do	      
theenv=envclangoffload${havempi}
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:true"
spack env activate $theenv

spack add gcc@${gcc_vllvmoffload}
spack add cuda@${cuda_voffload} +allow-unsupported-compilers
spack add llvm@${llvm_voffload} openmp=project +libomptarget +libomptarget_debug ^cuda@${cuda_voffload}

spack add hwloc # Try hwloc to solve duplicate 2026-03-24
spack add git
spack add ninja
spack add cmake@${cmake_vnew}^gcc@${gcc_vllvmoffload}
spack add libxml2
spack add boost@${boost_vold}^gcc@${gcc_vllvmoffload}
spack add util-linux-uuid^gcc@${gcc_vllvmoffload}
spack add python^gcc@${gcc_vllvmoffload}
spack add openmpi@${ompi_vnew}^gcc@${gcc_vllvmoffload}
if [ "$havempi" == "mpi" ]; then
spack add hdf5@${hdf5_vold} +fortran +hl +mpi ^gcc@${gcc_vllvmoffload}
else
spack add hdf5@${hdf5_vold} +fortran +hl ~mpi ^gcc@${gcc_vllvmoffload}
fi
spack add fftw@${fftw_vold} -mpi ^gcc@${gcc_vllvmoffload} #Avoid MPI for simplicity
spack add openblas threads=openmp ^gcc@${gcc_vllvmoffload}
#PK: amdlibflame@5.2 will fail 2026-03-24. Link error
#if [ "$ourplatform" == "AMD" ]; then
#spack add amdblis; spack add amdlibflame; #spack add amd-aocl
#fi
if [ "$ourplatform" == "Intel" ]; then
spack add intel-oneapi-mkl
fi

spack add py-lxml
spack add py-matplotlib
spack add py-pandas

if [ "$havempi" == "mpi" ]; then
spack add py-mpi4py
fi

spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}
install_environment
spack env deactivate
done

if [ "$ourplatform" == "AMD" ]; then

for havempi in mpi nompi
do	      

theenv=envamdclang${havempi}
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:true"
spack env activate $theenv

#Use older likely offload compatible version of GCC
spack add gcc@${gcc_vllvmoffload}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2
spack add boost@${boost_vold}^gcc@${gcc_vllvmoffload}
spack add util-linux-uuid^gcc@${gcc_vllvmoffload}
spack add python^gcc@${gcc_vllvmoffload}

if [ "$havempi" == "mpi" ]; then
spack add openmpi@${ompi_vnew}^gcc@${gcc_vllvmoffload}
fi
#spack add hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl +mpi
if [ "$havempi" == "mpi" ]; then
spack add hdf5@${hdf5_vold} +fortran +hl +mpi
else
spack add hdf5@${hdf5_vold} +fortran +hl ~mpi
fi
spack add fftw@${fftw_vold} -mpi ^gcc@${gcc_vllvmoffload} #Avoid MPI for simplicity
spack add openblas threads=openmp ^gcc@${gcc_vllvmoffload}
#MAR26spack add amdblis^gcc@${gcc_vold}; spack add amdlibflame^gcc@${gcc_vold}; #spack add amd-aocl

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
if [ "$havempi" == "mpi" ]; then
spack add py-mpi4py
fi
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}
if [ "$havempi" == "mpi" ]; then
spack add quantum-espresso +mpi +qmcpack
fi
install_environment
spack env deactivate
done
fi

#MAR26#if [ "$ourplatform" == "Intel" ]; then
#MAR26#theenv=envinteloneapinompi
#MAR26#echo --- Setting up $theenv `date`
#MAR26#spack env create $theenv
#MAR26##spack -e $theenv config add "concretizer:unify:when_possible"
#MAR26#spack -e $theenv config add "concretizer:unify:true"
#MAR26#spack env activate $theenv
#MAR26#
#MAR26#spack add gcc@${gcc_vintel}
#MAR26#spack add git
#MAR26#spack add ninja
#MAR26#spack add cmake@${cmake_vnew}
#MAR26#spack add libxml2%gcc@${gcc_vintel}
#MAR26#spack add boost@${boost_vnew}%gcc@${gcc_vintel}
#MAR26#spack add util-linux-uuid%gcc@${gcc_vintel}
#MAR26#spack add python%gcc@${gcc_vintel}
#MAR26#spack add hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
#MAR26#spack add fftw@${fftw_vnew}%gcc@${gcc_vintel} -mpi #Avoid MPI for simplicity
#MAR26#
#MAR26#spack add py-lxml
#MAR26#spack add py-matplotlib
#MAR26#spack add py-pandas
#MAR26#spack add py-numpy@${numpy_vold}
#MAR26#spack add py-scipy
#MAR26#spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
#MAR26#install_environment
#MAR26#spack env deactivate
#MAR26#
#MAR26#theenv=envinteloneapimpi
#MAR26#echo --- Setting up $theenv `date`
#MAR26#spack env create $theenv
#MAR26##spack -e $theenv config add "concretizer:unify:when_possible"
#MAR26#spack -e $theenv config add "concretizer:unify:true"
#MAR26#spack env activate $theenv
#MAR26#
#MAR26#spack add gcc@${gcc_vintel}
#MAR26#spack add git
#MAR26#spack add ninja
#MAR26#spack add cmake@${cmake_vnew}
#MAR26#spack add libxml2%gcc@${gcc_vintel}
#MAR26#spack add boost@${boost_vnew}%gcc@${gcc_vintel}
#MAR26#spack add util-linux-uuid%gcc@${gcc_vintel}
#MAR26#spack add python%gcc@${gcc_vintel}
#MAR26#spack add hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
#MAR26#spack add fftw@${fftw_vnew}%gcc@${gcc_vintel} -mpi #Avoid MPI for simplicity
#MAR26#
#MAR26#spack add py-lxml
#MAR26#spack add py-matplotlib
#MAR26#spack add py-pandas
#MAR26#spack add py-numpy@${numpy_vold}
#MAR26#spack add py-scipy
#MAR26#spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
#MAR26#install_environment
#MAR26#spack env deactivate
#MAR26#fi
#MAR26
#MAR26# CAUTION: Removing build deps reveals which spack packages do not have correct runtime deps specified and may result in breakage
#MAR26#echo --- Removing build deps
#MAR26#for f in `spack env list`
#MAR26#do
#MAR26#    spack env activate $f
#MAR26#    spack gc --yes-to-all
#MAR26#    echo --- Software for environment $f
#MAR26#    spack env status
#MAR26#    spack find
#MAR26#    spack env deactivate
#MAR26#done

echo --- Making loads files
for f in `spack env list`
do
    spack env activate $f
    spack module tcl refresh -y
    spack env loads
    spack env deactivate
done

echo --- FINISH environment setup `date`
