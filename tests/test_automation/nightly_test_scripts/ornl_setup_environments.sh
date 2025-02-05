#!/bin/bash

echo --- START environment setup `date`

# serial   : single install
# 8up      : 8 installs
# par48    : install -j 48 
# makefile : make -j

parallelmode=par48

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
	make -j
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

theenv=envgccnewmpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vnew}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vnew}
spack add boost@${boost_vnew}%gcc@${gcc_vnew}
spack add util-linux-uuid%gcc@${gcc_vnew}
spack add python%gcc@${gcc_vnew}
spack add openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack add hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
spack add fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vnew} threads=openmp
#spack add blis%gcc@${gcc_vnew} threads=openmp
#spack add libflame%gcc@${gcc_vnew} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-mpi4py
spack add py-numpy@${numpy_vnew}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
spack add quantum-espresso@7.4 +mpi +qmcpack
export CMAKE_BUILD_PARALLEL_LEVEL=8 # For PySCF
spack add py-pyscf
#spack add rmgdft #Fails to compile with GCC14 due to bug in vendored SCALAPACK 

#Luxury options for actual science use:
spack add py-requests # for pseudo helper
spack add py-ase      # full Atomic Simulation Environment
#spack add graphviz +libgd # NEXUS requires optional PNG support in dot
spack add libffi
spack add graphviz +pangocairo # NEXUS requires optional PNG support in dot
spack add py-pydot    # NEXUS optional
spack add py-spglib   # NEXUS optional 
spack add py-seekpath # NEXUS optional
spack add py-pycifrw  # NEXUS optional
#NOT IN SPACK spack add py-cif2cell # NEXUS optional

install_environment
spack env deactivate

theenv=envgccnewnompi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vnew}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vnew}
spack add boost@${boost_vnew}%gcc@${gcc_vnew}
spack add util-linux-uuid%gcc@${gcc_vnew}
spack add python%gcc@${gcc_vnew}
#spack add openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack add hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl ~mpi
spack add fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vnew} threads=openmp
#spack add blis%gcc@${gcc_vnew} threads=openmp
#spack add libflame%gcc@${gcc_vnew} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
#spack add py-mpi4py
spack add py-numpy@${numpy_vnew}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl ~mpi

install_environment
spack env deactivate

theenv=envgccoldnompi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vold}
spack add cmake@${cmake_vold}%gcc@${gcc_vold}
spack add libxml2@${libxml2_v}%gcc@${gcc_vold}
spack add boost@${boost_vold}%gcc@${gcc_vold}
spack add util-linux-uuid%gcc@${gcc_vold}
spack add python%gcc@${gcc_vold}
#spack add openmpi@${ompi_vnew}%gcc@${gcc_vold}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl ~mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vold} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vold} threads=openmp
#spack add blis%gcc@${gcc_vold} threads=openmp
#spack add libflame%gcc@${gcc_vold} threads=openmp
spack add git
spack add ninja

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
#spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vold} +fortran +hl ~mpi
install_environment
spack env deactivate

theenv=envgccoldmpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vold}
spack add cmake@${cmake_vold}%gcc@${gcc_vold}
spack add libxml2@${libxml2_v}%gcc@${gcc_vold}
spack add boost@${boost_vold}%gcc@${gcc_vold}
spack add util-linux-uuid%gcc@${gcc_vold}
spack add python%gcc@${gcc_vold}
spack add openmpi@${ompi_vnew}%gcc@${gcc_vold}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl +mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vold} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vold} threads=openmp
#spack add blis%gcc@${gcc_vold} threads=openmp
#spack add libflame%gcc@${gcc_vold} threads=openmp
spack add git
spack add ninja

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vold} +fortran +hl +mpi
#spack add quantum-espresso@7.4 +mpi +qmcpack
#export CMAKE_BUILD_PARALLEL_LEVEL=8 # For PySCF
#spack add py-pyscf
spack add rmgdft

#Luxury options for actual science use:
spack add py-requests # for pseudo helper
spack add py-ase      # full Atomic Simulation Environment
#spack add graphviz +libgd # NEXUS requires optional PNG support in dot
spack add libffi
spack add graphviz +pangocairo # NEXUS requires optional PNG support in dot
spack add py-pydot    # NEXUS optional
spack add py-spglib   # NEXUS optional 
spack add py-seekpath # NEXUS optional
spack add py-pycifrw  # NEXUS optional
#NOT IN SPACK spack add py-cif2cell # NEXUS optional

install_environment
spack env deactivate

theenv=envclangnewmpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vnew}
spack add llvm@${llvm_vnew}%gcc@${gcc_vnew}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vnew}
spack add boost@${boost_vnew}%gcc@${gcc_vnew}
spack add util-linux-uuid%gcc@${gcc_vnew}
spack add python%gcc@${gcc_vnew}
spack add openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack add hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
spack add fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vnew} threads=openmp
#spack add blis%gcc@${gcc_vnew} threads=openmp
#spack add libflame%gcc@${gcc_vnew} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
install_environment
spack env deactivate


# Build LLVM offload with preferred GCC since CUDA may not support new GCC
# Build with new CMake
# TO DO: Match chosen cuda with version installed on system
theenv=envclangoffloadmpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vllvmoffload}
spack add cuda@${cuda_voffload} +allow-unsupported-compilers
spack add llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload} targets=all
#spack add llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload} targets=all cuda_arch=70

spack add hwloc
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vllvmoffload}
spack add boost@${boost_vold}%gcc@${gcc_vllvmoffload}
spack add util-linux-uuid%gcc@${gcc_vllvmoffload}
spack add python%gcc@${gcc_vllvmoffload}
spack add openmpi@${ompi_vnew}%gcc@${gcc_vllvmoffload}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl +mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vllvmoffload} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vllvmoffload} threads=openmp
#spack add blis%gcc@${gcc_vllvmoffload} threads=openmp
#spack add libflame%gcc@${gcc_vllvmoffload} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl +mpi
install_environment
spack env deactivate

theenv=envclangoffloadnompi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vllvmoffload}
spack add cuda@${cuda_voffload} +allow-unsupported-compilers
spack add llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload} targets=all

spack add hwloc
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vllvmoffload}
spack add boost@${boost_vold}%gcc@${gcc_vllvmoffload}
spack add util-linux-uuid%gcc@${gcc_vllvmoffload}
spack add python%gcc@${gcc_vllvmoffload}
#spack add openmpi@${ompi_vnew}%gcc@${gcc_vllvmoffload}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl ~mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vllvmoffload} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vllvmoffload} threads=openmp
#spack add blis%gcc@${gcc_vllvmoffload} threads=openmp
#spack add libflame%gcc@${gcc_vllvmoffload} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
#spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl ~mpi
install_environment
spack env deactivate


if [ "$ourplatform" == "AMD" ]; then
theenv=envamdclangmpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

#Use older likely offload compatible version of GCC
spack add gcc@${gcc_vllvmoffload}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vllvmoffload}
spack add boost@${boost_vold}%gcc@${gcc_vllvmoffload}
spack add util-linux-uuid%gcc@${gcc_vllvmoffload}
spack add python%gcc@${gcc_vllvmoffload}
spack add openmpi@${ompi_vnew}%gcc@${gcc_vllvmoffload}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl +mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vllvmoffload} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vllvmoffload} threads=openmp
#spack add blis%gcc@${gcc_vllvmoffload} threads=openmp
#spack add libflame%gcc@${gcc_vllvmoffload} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl +mpi
spack add quantum-espresso@7.4 +mpi +qmcpack

#spack add gcc@${gcc_vnew}
#spack add git
#spack add ninja
#spack add cmake@${cmake_vnew}
#spack add libxml2@${libxml2_v}%gcc@${gcc_vnew}
#spack add boost@${boost_vnew}%gcc@${gcc_vnew}
#spack add util-linux-uuid%gcc@${gcc_vnew}
#spack add python%gcc@${gcc_vnew}
#spack add openmpi@${ompi_vnew}%gcc@${gcc_vnew}
#spack add hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
#spack add fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
#spack add openblas%gcc@${gcc_vnew} threads=openmp
##spack add blis%gcc@${gcc_vnew} threads=openmp
##spack add libflame%gcc@${gcc_vnew} threads=openmp
#
#spack add py-lxml
#spack add py-matplotlib
#spack add py-pandas
#spack add py-mpi4py
#spack add py-scipy
#spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl +mpi
#spack add quantum-espresso@7.4 +mpi +qmcpack

#spack add rmgdft
install_environment
spack env deactivate

theenv=envamdclangnompi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vllvmoffload}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vllvmoffload}
spack add boost@${boost_vold}%gcc@${gcc_vllvmoffload}
spack add util-linux-uuid%gcc@${gcc_vllvmoffload}
spack add python%gcc@${gcc_vllvmoffload}
#spack add openmpi@${ompi_vnew}%gcc@${gcc_vllvmoffload}
spack add hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl ~mpi
spack add fftw@${fftw_vold}%gcc@${gcc_vllvmoffload} -mpi #Avoid MPI for simplicity
spack add openblas%gcc@${gcc_vllvmoffload} threads=openmp
#spack add blis%gcc@${gcc_vllvmoffload} threads=openmp
#spack add libflame%gcc@${gcc_vllvmoffload} threads=openmp

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
#spack add py-mpi4py
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vold}%gcc@${gcc_vllvmoffload} +fortran +hl ~mpi
install_environment
spack env deactivate
fi


if [ "$ourplatform" == "Intel" ]; then
theenv=envinteloneapinompi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vintel}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vintel}
spack add boost@${boost_vnew}%gcc@${gcc_vintel}
spack add util-linux-uuid%gcc@${gcc_vintel}
spack add python%gcc@${gcc_vintel}
spack add hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
spack add fftw@${fftw_vnew}%gcc@${gcc_vintel} -mpi #Avoid MPI for simplicity

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
install_environment
spack env deactivate

theenv=envinteloneapimpi
echo --- Setting up $theenv `date`
spack env create $theenv
spack -e $theenv config add "concretizer:unify:when_possible"
spack env activate $theenv

spack add gcc@${gcc_vintel}
spack add git
spack add ninja
spack add cmake@${cmake_vnew}
spack add libxml2@${libxml2_v}%gcc@${gcc_vintel}
spack add boost@${boost_vnew}%gcc@${gcc_vintel}
spack add util-linux-uuid%gcc@${gcc_vintel}
spack add python%gcc@${gcc_vintel}
spack add hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
spack add fftw@${fftw_vnew}%gcc@${gcc_vintel} -mpi #Avoid MPI for simplicity

spack add py-lxml
spack add py-matplotlib
spack add py-pandas
spack add py-numpy@${numpy_vold}
spack add py-scipy
spack add py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vintel} +fortran +hl ~mpi
install_environment
spack env deactivate
fi

# CAUTION: Removing build deps reveals which spack packages do not have correct runtime deps specified and may result in breakage
#echo --- Removing build deps
#for f in `spack env list`
#do
#    spack env activate $f
#    spack gc --yes-to-all
#    echo --- Software for environment $f
#    spack env status
#    spack find
#    spack env deactivate
#done

echo --- Making loads files
for f in `spack env list`
do
    spack env activate $f
    spack module tcl refresh -y
    spack env loads
    spack env deactivate
done

echo --- FINISH environment setup `date`
