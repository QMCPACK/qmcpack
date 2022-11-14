#!/bin/bash

#
# Installs compilers & libraries for QMCPACK testing via SPACK 
#
#          *** DESTROYS EXISTING INSTALLATION ***
#
#

# If you are behind a firewall:
# git config --global url."https://".insteadOf git://
# or in .gitconfig
#[url "https://"]
#	insteadOf = git://

echo --- START initial setup `date`

# Bug avoidance 20200902
#if [ -e /usr/lib/aomp/bin/clang ]; then
#    echo AOMP Clang install detected. This will break llvm install.
#    echo Suggest temporarily: sudo chmod og-rx /usr/lib/aomp
#    exit 1
#fi

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

if [ -e $HOME/apps/spack ]; then
    rm -r -f $HOME/apps/spack
fi
if [ -e $HOME/.spack ]; then
    rm -r -f $HOME/.spack
fi
mkdir $HOME/.spack

# Setup build multiplicity and preferred directories for spack
# Choose the fastest filesytem. Don't abuse shared nodes.
case "$ourhostname" in
    nitrogen )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 128

EOF
# Use flash /scratch for builds
	rm -r -f /scratch/$USER/spack_build_stage
	mkdir /scratch/$USER/spack_build_stage

	#Use system installed SSL. See https://spack.readthedocs.io/en/latest/getting_started.html#openssl
	cat >>$HOME/.spack/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1k
          prefix: /usr
          buildable: False
EOF
#	cat >>$HOME/.spack/packages.yaml<<EOF
#packages:
#EOF
# Experiment to See if clang builds and other problems clear 20200827
#	cat >>$HOME/.spack/packages.yaml<<EOF
#    all:
#        target: [x86_64]
#EOF
;;
    sulfur )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 96

EOF
#  concretizer: clingo
# Use flash /scratch for builds
	rm -r -f /scratch/$USER/spack_build_stage
	mkdir /scratch/$USER/spack_build_stage
	#Use system installed SSL. See https://spack.readthedocs.io/en/latest/getting_started.html#openssl
	#Use system RHEL gcc8 compiler for cmake to solve bootstrap problems
	cat >>$HOME/.spack/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1k
          prefix: /usr
          buildable: False
EOF
	;;
    *)
	echo "*** WARNING: Unknown host in initial ourhostname case statement. No custom onfiguration"
	;;
esac

cat >$HOME/.spack/modules.yaml<<EOF
modules:
  prefix_inspections::
    bin:
    - PATH
    man:
    - MANPATH
    share/man:
    - MANPATH
    share/aclocal:
    - ACLOCAL_PATH
    lib/pkgconfig:
    - PKG_CONFIG_PATH
    lib64/pkgconfig:
    - PKG_CONFIG_PATH
    share/pkgconfig:
    - PKG_CONFIG_PATH
    '':
    - CMAKE_PREFIX_PATH
EOF

cat >$HOME/.spack/spack.yaml<<EOF
spack:
  concretization:
    unify:  true
EOF

if [ ! -e $HOME/apps ]; then
mkdir $HOME/apps
fi

cd $HOME/apps
#DEBUG git clone --depth 1 https://github.com/spack/spack.git 
git clone https://github.com/spack/spack.git 

cd $HOME/apps/spack
# For reproducibility, use a specific version of Spack
# Prefer to use tagged releases https://github.com/spack/spack/releases
#git checkout 675210bd8bd1c5d32ad1cc83d898fb43b569ed74
#commit 675210bd8bd1c5d32ad1cc83d898fb43b569ed74 (HEAD -> develop, origin/develop, origin/HEAD)
#Author: Erik Schnetter <schnetter@gmail.com>
#Date:   Mon Jan 10 15:33:14 2022 -0500
#
#automake: New version 1.16.5 (#28299)

echo --- Git version and last log entry
git log -1

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin

## Consider using a GCC toolset on Red Hat systems to use
## recent compilers with better architecture support.
## e.g. dnf install gcc-toolset-11
#if [ -e /opt/rh/gcc-toolset-11/root/bin/gcc ]; then
#    echo --- Added gcc-toolset-11 to path for RHEL provided GCC11 compilers
#    export PATH=/opt/rh/gcc-toolset-11/root/bin/:$PATH
#else
#if [ -e /opt/rh/gcc-toolset-10/root/bin/gcc ]; then
#    echo --- Added gcc-toolset-10 to path for RHEL provided GCC10 compilers
#    export PATH=/opt/rh/gcc-toolset-10/root/bin/:$PATH
#fi
#fi

export DISPLAY="" 
export SPACK_ROOT=$HOME/apps/spack
export PATH=$SPACK_ROOT/bin:$PATH
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Changing RMGDFT boost dependency
sed -i 's/^ .* depends.*boost@1.61.*//g' $HOME/apps/spack/var/spack/repos/builtin/packages/rmgdft/package.py
echo --- Spack list
spack find
echo --- Spack compilers
spack compilers
echo --- Spack compiler add
spack compiler find
echo --- Spack compilers
spack compilers
echo --- Modules list
module list
echo --- End listings

echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
echo --- load gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
module list
spack compiler find
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vold}
spack install gcc@${gcc_vold}
echo --- load gcc@${gcc_vold}
spack load gcc@${gcc_vold}
module list
spack compiler find
spack unload gcc@${gcc_vold}
#echo --- gcc@master
#spack install gcc@master
#echo --- load gcc@master
#spack load gcc@master
#module list
#spack compiler find
#spack unload gcc@master
#echo --- gcc@${gcc_vcuda}
#spack install gcc@${gcc_vcuda}
#spack load gcc@${gcc_vcuda}
#spack compiler find
#spack unload gcc@${gcc_vcuda}
if [ "$ourplatform" == "Intel" ]; then
echo --- gcc@${gcc_vintel}
spack install gcc@${gcc_vintel}
spack load gcc@${gcc_vintel}
spack compiler find
spack unload gcc@${gcc_vintel}
fi
echo --- gcc@${gcc_vnvhpc}
spack install gcc@${gcc_vnvhpc}
spack load gcc@${gcc_vnvhpc}
spack compiler find
spack unload gcc@${gcc_vnvhpc}
echo --- llvm@${llvm_vnew}
spack install llvm@${llvm_vnew}
spack load llvm@${llvm_vnew}
spack compiler find
spack unload llvm@${llvm_vnew}
#echo --- llvm@main
#spack install llvm@main +cuda cuda_arch=70
#spack load llvm@main
#spack compiler find
#spack unload llvm@main
echo --- Cleanup
spack gc --yes-to-all
echo --- Spack compilers
spack compilers
echo --- Modules list
module list
echo --- End listings
echo --- FINISH initial setup `date`
bash $HOME/.cron_jobs/ornl_setup_environments.sh

