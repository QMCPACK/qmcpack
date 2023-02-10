# Build stage with Spack pre-installed and ready to be used
FROM spack/centos-stream as builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)

RUN which tar
RUN tar -cf tarfilename.tar -T /dev/null

RUN spack install gcc@11.3.0
RUN spack load gcc@11.3.0
RUN spack compiler add $(spack location -i gcc@11.3.0)
RUN spack compilers

RUN yum install texinfo -y

RUN mkdir /opt/spack-environment \
&& export gcc_vnew=11.3.0\
&& export cmake_vnew=3.24.2 \
&& export libxml2_v=2.9.13 \
&& export boost_vnew=1.79.0 \
&& export hdf5_vnew=1.12.2 \
&& export fftw_vnew=3.3.10 \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - gcc@${gcc_vnew}%gcc@${gcc_vnew}" \
&&   echo "  - texinfo" \
&&   echo "  - binutils" \
&&   echo "  - git" \
&&   echo "  - ninja" \
&&   echo "  - cmake@${cmake_vnew}" \
&&   echo "  - libxml2@${libxml2_v}%gcc@${gcc_vnew}" \
&&   echo "  - boost@${boost_vnew}%gcc@${gcc_vnew}" \
&&   echo "  - util-linux-uuid%gcc@${gcc_vnew}" \
&&   echo "  - python%gcc@${gcc_vnew}" \
&&   echo "  - hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl ~mpi" \
&&   echo "  - fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi" \
&&   echo "  - openblas%gcc@${gcc_vnew} threads=openmp" \
&&   echo "  - py-lxml" \
&&   echo "  - py-matplotlib" \
&&   echo "  - py-pandas" \
&&   echo "  - py-scipy" \
&&   echo "  - py-h5py ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl ~mpi" \
&&   echo "  concretizer:" \
&&   echo "    unify: true" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "    build_stage:" \
&&   echo "      - $HOME" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml


# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
FROM quay.io/centos/centos:stream

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

# building from source is missing ld (the linker)
RUN yum install -y gcc-toolset-11-binutils
RUN yum install -y rsync
RUN export PATH=/opt/rh/gcc-toolset-11/root/bin/:$PATH
RUN yum install -y libgcc.i686 glibc-devel.i686 libstdc++-devel.i686
RUN ranlib /opt/view/lib/gcc/x86_64-pc-linux-gnu/11.3.0/libgcc.a

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c", "$*", "--" ]
CMD [ "/bin/bash" ]

