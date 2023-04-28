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
# Bug avoidance 20230404
if [ -e /opt/rocm-*/bin/rocminfo ]; then
    echo Spack LLVM16 installation will fail with ROCm present
    echo Suggest temporarily: sudo chmod og-rx /opt/rocm-*
    exit 1
fi
# Bug avoidance 20200902
#if [ -e /usr/lib/aomp/bin/clang ]; then
#    echo AOMP Clang install detected. This will break LLVM install.
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
    nitrogen2 )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 96
  connect_timeout: 120

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
;;
    nitrogen )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 128
  connect_timeout: 120

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
;;
    sulfur )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 96
  connect_timeout: 120

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

#cat >$HOME/.spack/spack.yaml<<EOF
#spack:
#  concretization:
#    unify:  true
#EOF
cat >$HOME/.spack/spack.yaml<<EOF
spack:
  concretization:
    unify:  when_possible
EOF

if [ ! -e $HOME/apps ]; then
mkdir $HOME/apps
fi

cd $HOME/apps

git clone https://github.com/spack/spack.git

cd $HOME/apps/spack

# For reproducibility, use a specific version of Spack
# Prefer to use tagged releases https://github.com/spack/spack/releases
git checkout 6edfc070926f7934eda5882de20ac3fb193e310a
#commit 6edfc070926f7934eda5882de20ac3fb193e310a (HEAD -> develop, origin/develop, origin/HEAD)
#Author: Alec Scott <hi@alecbcs.com>
#Date:   Mon Apr 24 12:04:31 2023 -0700
#
#    megadock: add v4.1.1 (#37154)

echo --- Git version and last log entry
git log -1

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin

# Consider using a GCC toolset on Red Hat systems to use
# recent compilers with better architecture support.
# e.g. dnf install gcc-toolset-11
#if [ -e /opt/rh/gcc-toolset-11/root/bin/gcc ]; then
#    echo --- Added gcc-toolset-11 to path for RHEL provided GCC11 compilers
#    export PATH=/opt/rh/gcc-toolset-11/root/bin/:$PATH
#fi
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
echo --- Bootstrap
spack bootstrap now

#echo --- Changing RMGDFT boost dependency
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

echo --- gcc@${gcc_vnew} `date`
spack install gcc@${gcc_vnew}
echo --- load gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
module list
spack compiler find
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vold}  `date`
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
echo --- gcc@${gcc_vintel}  `date`
spack install gcc@${gcc_vintel}
spack load gcc@${gcc_vintel}
spack compiler find
spack unload gcc@${gcc_vintel}
fi
echo --- gcc@${gcc_vnvhpc}  `date`
spack install gcc@${gcc_vnvhpc}
spack load gcc@${gcc_vnvhpc}
spack compiler find
spack unload gcc@${gcc_vnvhpc}
echo --- llvm@${llvm_vnew}  `date`
spack install llvm@${llvm_vnew}
spack load llvm@${llvm_vnew}
spack compiler find
spack unload llvm@${llvm_vnew}
#echo --- llvm@main
#spack install llvm@main +cuda cuda_arch=70
#spack load llvm@main
#spack compiler find
#spack unload llvm@main
#echo --- Cleanup
#spack gc --yes-to-all
echo --- gcc@${gcc_vllvmoffload} for offload  `date`
spack install gcc@${gcc_vllvmoffload}
spack load gcc@${gcc_vllvmoffload}
spack compiler find
spack unload gcc@${gcc_vllvmoffload}

echo --- llvm@${llvm_voffload} for offload  `date`
spack install gcc@${gcc_vllvmoffload}
spack install cuda@${cuda_voffload} +allow-unsupported-compilers
spack install llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload} ~libcxx +compiler-rt ~lldb ~gold ~omp_as_runtime targets=all
spack load llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload}
spack compiler find
spack unload llvm@${llvm_voffload}%gcc@${gcc_vllvmoffload}

echo --- Spack compilers  `date`
spack compilers
echo --- Modules list
module list
echo --- End listings
echo --- FINISH initial setup `date`
bash $HOME/.cron_jobs/ornl_setup_environments.sh

echo --- REMEMBER REMEMBER
echo If ROCm installed, sudo chmod og+rx /opt/rocm-*
