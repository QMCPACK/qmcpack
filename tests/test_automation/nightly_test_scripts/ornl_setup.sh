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

export SPACK_ROOT=$HOME/apps/spack
export SPACK_USER_CONFIG_PATH=$HOME/apps/spack_user_config  # Avoid using $HOME/.spack
if [ ! -e $SPACK_USER_CONFIG_PATH ]; then
    mkdir $SPACK_USER_CONFIG_PATH
else
    rm -r -f $SPACK_USER_CONFIG_PATH/*
fi

if [ -e $HOME/.spack ]; then
    rm -r -f $HOME/.spack 
fi


# Pin the version of the package repository for reproducibility
# Docs: https://spack.readthedocs.io/en/latest/repositories.html#updating-and-pinning
# The package repo will be at something like ~/.spack/package_repos/fncqgg4/

cat >>$SPACK_USER_CONFIG_PATH/repos.yaml<<EOF
repos:
  builtin:
    commit: ba81a382694253e885f2f3012992ed8d7df2fe08

#commit ba81a382694253e885f2f3012992ed8d7df2fe08 (HEAD -> develop, origin/develop)
#Author: Alec Scott <scott112@llnl.gov>
#Date:   Wed Apr 8 16:29:47 2026 -0700
#
#    helm: new package (#4179)
#
#    Signed-off-by: Alec Scott <alec@llnl.gov>
EOF


# Setup build multiplicity and preferred directories for spack
# Choose the fastest filesytem. Don't abuse shared nodes.
case "$ourhostname" in
    nitrogen2 )
	cat >$SPACK_USER_CONFIG_PATH/config.yaml<<EOF
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
	cat >>$SPACK_USER_CONFIG_PATH/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1k
          prefix: /usr
EOF
;;
    nitrogen )
	cat >$SPACK_USER_CONFIG_PATH/config.yaml<<EOF
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
	cat >>$SPACK_USER_CONFIG_PATH/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1k
          prefix: /usr
EOF
;;
    sulfur )
	cat >$SPACK_USER_CONFIG_PATH/config.yaml<<EOF
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
	cat >>$SPACK_USER_CONFIG_PATH/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1k
          prefix: /usr
EOF
	;;
    *)
	echo "*** WARNING: Unknown host in initial ourhostname case statement. No custom onfiguration"
	;;
esac

cat >$SPACK_USER_CONFIG_PATH/modules.yaml<<EOF
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

cat >$SPACK_USER_CONFIG_PATH/spack.yaml<<EOF
spack:
  concretization:
     unify:  true
EOF

if [ ! -e $HOME/apps ]; then
mkdir $HOME/apps
fi

cd $HOME/apps

if [ -e spack ]; then
    cd spack
    echo -- Resetting existing spack git repo
    git checkout -f
    git clean -fd
    git checkout develop
    git pull
    cd ..
else
    git clone https://github.com/spack/spack.git
fi

if [ ! -e spack/CHANGELOG.md ]; then
    echo "--- FAILED TO FIND spack/CHANGELOG.md . BAD CLONE or I/O PROBLEMS. ABORTING"
    exit 1	 
fi

cd $SPACK_ROOT
# For reproducibility, use a specific version of Spack
# Prefer to use tagged releases https://github.com/spack/spack/releases
git checkout 45b9069e4997f4b453b2770886b1e8ba790980f4
#commit 45b9069e4997f4b453b2770886b1e8ba790980f4 (HEAD -> develop, origin/develop, origin/HEAD)
#Author: Harmen Stoppels <me@harmenstoppels.nl>
#Date:   Wed Apr 8 14:38:32 2026 +0200
#
#    new_installer.py: sub_process < /dev/null (#52221)

# Limit overly strong rmg boost dependency to allow concretizer:unify:true
#sed -ibak 's/boost@1.61.0:1.82.0/boost@1.61.0:1.82.0", when="@:6.1.2/g' var/spack/repos/spack_repo/builtin/packages/rmgdft/package.py

echo --- Git version and last log entry
git log -1

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin

# Consider using a GCC toolset on Red Hat systems to use
# recent compilers with better architecture support.
# e.g. dnf install gcc-toolset-14
#if [ -e /opt/rh/gcc-toolset-14/enable ]; then
#    echo --- Using gcc-toolset-14 for newer compilers
#    source /opt/rh/gcc-toolset-14/enable 
#fi

export DISPLAY="" 
# SPACK_ROOT & SPACK_USER_CONFIG_PATH already set
export PATH=$SPACK_ROOT/bin:$PATH
. $SPACK_ROOT/share/spack/setup-env.sh
echo --- Bootstrap
spack bootstrap now
spack gpg init

echo --- Setup/use mirror in /scratch/$USER/spack_mirror
if [ ! -e /scratch/$USER/spack_mirror ];then
    mkdir /scratch/$USER/spack_mirror
    #spack mirror create -d /scratch/$USER/spack_mirror #  Will error due to no packages
else
    echo --- Mirror already exists. Carefully consider if refresh needed after major spack updates.
    du -ksh /scratch/$USER/spack_mirror
fi
spack mirror add localmirror /scratch/$USER/spack_mirror
spack buildcache update-index /scratch/$USER/spack_mirror
spack mirror set --autopush localmirror  # enable automatic push for an existing mirror
spack mirror set --unsigned localmirror  # disable signing and verification

echo --- Spack list
spack find
echo --- Spack compilers
spack compilers
echo --- Spack compiler find
spack compiler find
echo --- Spack compilers
spack compilers
echo --- Modules list
module list
echo --- End listings

#DEBUGecho --- Screen llvm compilations
#DEBUGecho --- Trying without cuda_arch
#DEBUGecho --- llvm@22 for offload  `date`
#DEBUGspack spec llvm@22 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack install llvm@22 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack load llvm@22 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack unload llvm
#DEBUGecho --- llvm@21 for offload  `date`
#DEBUGspack spec llvm@21 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack install llvm@21 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack load llvm@21 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack unload llvm
#DEBUGecho --- llvm@20 for offload  `date`
#DEBUGspack spec llvm@20 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack install llvm@20 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack load llvm@20 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack unload llvm
#DEBUGecho --- llvm@19 for offload  `date`
#DEBUGspack spec llvm@19 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack install llvm@19 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack load llvm@19 +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
#DEBUGspack unload llvm

# 2026-03-28: Install LLVM first. Will use system gcc and avoid problems/spack bugs with newer gcc, binutils etc.
echo --- llvm@${llvm_voffload} for offload  `date`

spack spec llvm@${llvm_voffload} openmp=project +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
spack install llvm@${llvm_voffload} openmp=project +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
spack load llvm@${llvm_voffload} openmp=project +libomptarget +libomptarget_debug ^cuda@${cuda_voffload} +allow-unsupported-compilers
spack unload llvm

echo --- llvm@${llvm_vnew}  `date`
if [ "${llvm_vnew}" == "${llvm_voffload}" ]; then
    echo Skipping: Offload and New LLVM versions are identical
else
    spack install llvm@${llvm_vnew} openmp=project
    spack load llvm@${llvm_vnew} openmp=project
    spack compiler find
    spack unload llvm@${llvm_vnew}
fi

echo --- gcc@${gcc_vold}  `date`
if [ "${gcc_vllvmoffload}" == "${gcc_vold}" ]; then
    echo Skipping: Already available via offload version
else
    spack install gcc@${gcc_vold}
    spack load gcc@${gcc_vold}
    module list
    spack compiler find
    spack unload gcc@${gcc_vold}
fi

echo --- gcc@${gcc_vnew} `date`
if [ "${gcc_vllvmoffload}" == "${gcc_vnew}" ] || [ "${gcc_vold}" == "${gcc_vnew}" ] ; then
    echo Skipping: Already available via offload or old versions
else
    spack install gcc@${gcc_vnew}
    spack load gcc@${gcc_vnew}
    module list
    spack compiler find
    spack unload gcc@${gcc_vnew}
fi

echo --- Spack compilers  `date`
spack compilers
echo --- Modules list
module list
echo --- End listings
echo --- FINISH initial setup `date`
bash $HOME/.cron_jobs/ornl_setup_environments.sh
