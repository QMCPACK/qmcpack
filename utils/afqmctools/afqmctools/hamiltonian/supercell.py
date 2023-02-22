import time
import os
import sys
import math
import h5py
import numpy
from functools import reduce
from pyscf import lib
from pyscf.pbc import tools, df
from mpi4py import MPI
from afqmctools.utils.parallel import fair_share, bisect
from afqmctools.utils.pyscf_utils import (
        load_from_pyscf_chk
        )

from afqmctools.utils.io import (
        format_fixed_width_floats,
        format_fixed_width_strings
        )


def write_hamil_supercell(comm, scf_data, hamil_file, chol_cut,
                          verbose=True, cas=None, ortho_ao=False,
                          maxvecs=20, exxdiv='ewald', nelec=None):

    nproc = comm.size
    rank = comm.rank

    tstart = time.process_time()()

    # Unpack pyscf data.
    # 1. core (1-body) Hamiltonian.
    hcore = scf_data['hcore']
    # 2. Rotation matrix to orthogonalised basis.
    X = scf_data['X']
    # 3. Pyscf cell object.
    cell = scf_data['cell']
    # 4. kpoints
    kpts = scf_data['kpts']
    # 5. MO occupancies.
    Xocc = scf_data['Xocc']
    # 6. Number of MOs per kpoint.
    nmo_pk = scf_data['nmo_pk']
    nkpts = len(kpts)
    nao = cell.nao_nr()

    # Update for GDF
    mydf = df.FFTDF(cell, kpts)
    nmo_max = numpy.max(nmo_pk)
    assert(nmo_max is not None)

    kconserv = tools.get_kconserv(cell, kpts)
    ik2n, nmo_tot = setup_basis_map(Xocc, nmo_max, nkpts, nmo_pk, ortho_ao)
    if comm.rank == 0:
        h5file = h5py.File(hamil_file, "w")
        h5grp = h5file.create_group("Hamiltonian")
    else:
        h5file = None

    if rank == 0:
        h1size = write_one_body(hcore, X, nkpts, nmo_pk, nmo_max, ik2n, h5grp)
    comm.barrier()
    numv = 0

    if rank == 0 and verbose:
        print(" # Time to reach cholesky: {:.2e} s".format(time.process_time()()-tstart))
        sys.stdout.flush()
    tstart = time.process_time()()

    # Setup parallel partition of work.
    part = Partition(comm, maxvecs, nmo_tot, nmo_max, nkpts)
    if comm.rank == 0 and verbose:
        print(" # Each kpoint is distributed accross {} mpi tasks."
              .format(part.nproc_pk))
        sys.stdout.flush()
    maxvecs = part.maxvecs
    # Set up mapping for shifted FFT grid.
    gmap, Qi, ngs = generate_grid_shifts(cell)
    if rank == 0 and verbose:
        app_mem = (nkpts*nkpts*nmo_max*nmo_max*2*ngs*16) / 1024.0**3
        print(" # Approx. total memory required: {:.2e} GB.".format(app_mem))
        sys.stdout.flush()

    if rank == 0 and verbose:
        print(" # Generating orbital products.")
        sys.stdout.flush()

    # Left and right pair densities. Only left contains Coulomb kernel.
    t0 = time.process_time()()
    Xaoik, Xaolj = gen_orbital_products(cell, mydf, X, nmo_pk, ngs,
                                        part, kpts, nmo_max)
    t1 = time.process_time()()
    if part.rank == 0 and verbose:
        print(" # Time to generate orbital products: {:.2e} s".format(t1-t0))
        sys.stdout.flush()
    # Finally perform Cholesky decomposition.
    solver = Cholesky(part, kconserv, gtol_chol=chol_cut)
    cholvecs = solver.run(comm, Xaoik, Xaolj, part, kpts,
                          nmo_pk, nmo_max, Qi, gmap)
    num_chol = cholvecs.shape[-1]
    comm.barrier()
    # Maximum number of integrals to include in single block in h5 file.
    max_block_size = int(2e6)
    # Write cholesky vectors to file.
    v2cnts = write_cholesky(h5file, comm, cholvecs, part, ik2n, nmo_pk,
                            chol_cut, max_block_size)
    if rank == 0:
        # Write remaining descriptor arrays to file.
        write_info(h5grp, cell, v2cnts, h1size, nmo_tot,
                   num_chol, kpts, ortho_ao, Xocc, nmo_pk,
                   ik2n, exxdiv, nelec=nelec)
        h5file.close()

def write_rhoG_supercell(comm, scf_data, hdf_file, Gcut,
                        verbose=True,
                        ortho_ao=False, write_real=False, phdf=False):

    nproc = comm.size
    rank = comm.rank

    tstart = time.process_time()()

    # Unpack pyscf data.
    # 1. Rotation matrix to orthogonalised basis.
    X = scf_data['X']
    # 2. Pyscf cell object.
    cell = scf_data['cell']
    # 3. kpoints
    kpts = scf_data['kpts']
    # 4. Number of MOs per kpoint.
    nmo_pk = scf_data['nmo_pk']
    nkpts = len(kpts)
    nao = cell.nao_nr()

    quit()

    # Update for GDF
    mydf = df.FFTDF(cell, kpts)
    nmo_max = numpy.max(nmo_pk)
    assert(nmo_max is not None)

    kconserv = tools.get_kconserv(cell, kpts)
    ik2n, nmo_tot = setup_basis_map(Xocc, nmo_max, nkpts, nmo_pk, ortho_ao)
    if comm.rank == 0:
        h5file = h5py.File(hamil_file, "w")
        h5grp = h5file.create_group("rhoG")
    else:
        h5file = None

    numv = 0

    # Setup parallel partition of work.
    part = Partition(comm, 0, nmo_tot, nmo_max, nkpts)
    if comm.rank == 0 and verbose:
        print(" # Each kpoint is distributed accross {} mpi tasks."
              .format(part.nproc_pk))
        sys.stdout.flush()
    # Set up mapping for shifted FFT grid.
    if rank == 0 and verbose:
        app_mem = (nmo_max*nmo_max*ngs*16) / 1024.0**3
        print(" # Approx. local memory required: {:.2e} GB.".format(app_mem))
        sys.stdout.flush()

    if rank == 0 and verbose:
        print(" # Generating orbital products.")
        sys.stdout.flush()

    for k in range(part.nkk):
        k1 = part.n2k1[k]
        k2 = part.n2k2[k]
        i0 = part.ij0 // nmo_pk[k2]
        iN = part.ijN // nmo_pk[k2]
        if part.ijN % nmo_pk[k2] != 0:
            iN += 1
        if iN > nmo_pk[k1]:
            iN = nmo_pk[k1]
        pij = part.ij0 % nmo_pk[k2]
        n_ = min(part.ijN, nmo_pk[k1]*nmo_pk[k2]) - part.ij0
        X_t = X[k1][:,i0:iN].copy()
        Xaoik[k,:,0:n_] = mydf.get_mo_pairs_G((X_t,X[k2]),
                                            (kpts[k1],kpts[k2]),
                                            kpts[k2]-kpts[k1],
                                            compact=False)[:,pij:pij+n_]
        X_t = None
        Xaolj[k,:,:] = Xaoik[k,:,:]
        coulG = tools.get_coulG(cell, kpts[k2]-kpts[k1], mesh=mydf.mesh)
        Xaoik[k,:,:] *= (coulG*cell.vol/ngs**2).reshape(-1,1)

def setup_basis_map(Xocc, nmo_max, nkpts, nmo_pk, ortho_ao):
    # setup basic mapping, use eigenvalues later
    ik2n = -1*numpy.ones((nmo_max,nkpts), dtype=numpy.int32)
    # ik2n[:,:]=-1
    cnt = 0
    # if ortho_ao:
    for ki in range(nkpts):
        for i in range(nmo_pk[ki]):
            ik2n[i,ki] = cnt
            cnt += 1
    nmo_tot = cnt
    return ik2n, nmo_tot


def write_one_body(hcore, X, nkpts, nmo_pk, nmo_max, ik2n, h5grp, gtol=1e-8):
    h1e_X = numpy.zeros((nkpts,nmo_max*(nmo_max+1)//2),dtype=numpy.complex128)
    for ki in range(nkpts):
        h1 = numpy.dot(X[ki][:,0:nmo_pk[ki]].T.conj(),
                       numpy.dot(hcore[ki], X[ki][:,0:nmo_pk[ki]]))
        ij = 0
        for i in range(nmo_pk[ki]):
            for j in range(0, i+1):
                h1e_X[ki,ij] = h1[i,j]
                ij += 1
        h1 = None
    h1size = write_h1(h5grp, h1e_X, nmo_pk, ik2n, gtol)
    return h1size


def write_cholesky(h5file, comm, cholvecs, part, ik2n, nmo_pk, gtol, max_ints):
    nproc = comm.size
    rank = comm.rank
    numv = cholvecs.shape[-1]
    if comm.rank == 0:
        h5grp_v2 = h5file.create_group("Hamiltonian/Factorized")
        # upper bound
        maxngrp = nproc * (((part.nkk*part.nij*numv)-1) // max_ints + 1)
        v2cnts = numpy.zeros(maxngrp, dtype=numpy.int32)
        h5nblk = 0
        vcnt, ngrp = write_cholMat_singleWriter(h5grp_v2, h5nblk,
                                                cholvecs, nmo_pk, part.n2k1,
                                                part.n2k2, part.ij0, ik2n,
                                                gtol, max_ints)
        if ngrp > 0:
            v2cnts[h5nblk:h5nblk+ngrp] = vcnt[0:ngrp]
            h5nblk += ngrp

        for np in range(1,nproc):
            cv = comm.recv(source=np, tag=np)
            n2k1_ = comm.recv(source=np, tag=nproc+np)
            n2k2_ = comm.recv(source=np, tag=2*nproc+np)
            ij0_ = comm.recv(source=np, tag=3*nproc+np)
            vcnt, ngrp = write_cholMat_singleWriter(h5grp_v2, h5nblk,
                                                    cv, nmo_pk, n2k1_,
                                                    n2k2_, ij0_, ik2n,
                                                    gtol, max_ints)
            if ngrp > 0:
                v2cnts[h5nblk:h5nblk+ngrp] = vcnt[0:ngrp]
                h5nblk += ngrp
            cv = None
            comm.barrier()   # to avoid a rain of messages on head node
        v2cnts = v2cnts[0:h5nblk]
        dummy = h5grp_v2.create_dataset("block_sizes", data=v2cnts[:],
                                        chunks=True, maxshape=(None,),
                                        compression="gzip" )
    else:
        for np in range(1,nproc):
            if rank == np:
                comm.send(cholvecs, dest=0, tag=rank)
                comm.send(part.n2k1, dest=0, tag=nproc+rank)
                comm.send(part.n2k2, dest=0, tag=2*nproc+rank)
                comm.send(part.ij0, dest=0, tag=3*nproc+rank)
            comm.barrier()
            v2cnts = None
    return v2cnts

def write_info(h5grp, cell, v2cnts, h1size, nmo_tot, numv,
               kpts, ortho_ao, Xocc, nmo_pk, ik2n, exxdiv, nelec=None):
    nkpts = len(kpts)
    if nelec is not None:
        nup, ndown = nelec
    else:
        nup = nkpts*(cell.nelectron+cell.spin)//2
        ndown = nkpts*cell.nelectron-nup
    nelectron = nup + ndown

    h2size = numpy.sum(v2cnts)
    dims = numpy.array([h1size, h2size, v2cnts.size,
                        nmo_tot, nup, ndown, 0, numv])
    h5grp.create_dataset("dims", data=dims)
    h5grp.create_dataset("ComplexIntegrals",
                         data=numpy.array([1]),
                         dtype=numpy.int32)

    if ortho_ao:
        occ = numpy.arange(0,nup+ndown)
        occ[nup:nup+ndown] += (nmo_tot-nup)
    else:
        # write occupation numbers of ground state
        # gamma for now
        occ = numpy.zeros(nup+ndown)
        cnt=0
        for ki in range(nkpts):
            for i in range(nmo_pk[ki]):
                if Xocc[ki][i] > 0.9:
                    occ[cnt] = ik2n[i,ki]
                    cnt += 1
        if cnt != nup:
            print(" # WARNING: Number of up electrons does not match "
                  "occupation string. ")
        cnt = 0
        for ki in range(nkpts):
            for i in range(nmo_pk[ki]):
                if Xocc[ki][i] > 1.9:
                    occ[nup+cnt] = ik2n[i,ki]
                    cnt += 1
        if cnt != ndown:
            print(" # WARNING: Number of down electrons does not match "
                  "occupation string. ")
    h5grp.create_dataset("occups", data=occ)

    # zero electron energies, including madelung term
    e0 = nkpts * cell.energy_nuc()
    if exxdiv=='ewald':
        madelung = tools.pbc.madelung(cell, kpts)
        e0 += madelung*nelectron * -.5
        print(" # Adding ewald correction to the energy: "
              "{:13.8e}".format(-0.5*madelung*nelectron))
    h5grp.create_dataset("Energies", data=numpy.array([e0, 0]))


def generate_grid_shifts(cell):
    mesh = cell.mesh
    ngs = numpy.prod(mesh)
    # define mapping to shift fft by Q=kk-ki+kl-kj
    g1 = numpy.arange(ngs, dtype=numpy.int32)
    g1 = g1.reshape(mesh,order='C')
    gmap = numpy.zeros((27,ngs), numpy.int32)
    # all 27 possible shifts
    Qi = numpy.zeros((27,3),dtype=numpy.float64)

    kvecs = cell.reciprocal_vectors()
    ii = 0
    for nx in range(-1,2):
        for ny in range(-1,2):
            for nz in range(-1,2):
                Qi[ii,:] = numpy.dot(numpy.array([nx,ny,nz]),kvecs)
                g0 = numpy.zeros(mesh, dtype=numpy.int32)
                for ix in range(mesh[0]):
                    for iy in range(mesh[1]):
                        for iz in range(mesh[2]):
                            g0[ix,iy,iz] = g1[(ix+nx)%mesh[0],
                                              (iy+ny)%mesh[1],
                                              (iz+nz)%mesh[2]]
                g0 = g0.reshape(-1, order='C')
                gmap[ii,:] = g0[:]
                ii += 1
    return gmap, Qi, ngs


class Partition(object):

    def __init__(self, comm, maxvecs, nmo_tot, nmo_max, nkpts, kp_sym=False):
        # to do:
        # keep usage of n2k1,n2k2 but assign entire rows of the (k1,k2) matrix
        # also add extra partitioning of i*j matrix, along nproc_pk cores
        # algorithm is the same
        self.maxvecs = maxvecs * nmo_tot
        self.rank = comm.rank
        self.size = comm.size

        if comm.size <= nkpts:
            self.kkbounds = numpy.zeros(comm.size+1, dtype=numpy.int32)
            if kp_sym:
                work = nkpts
            else:
                work = nkpts*nkpts
            for i in range(comm.size):
                self.kkbounds[i], self.kkbounds[i+1] = fair_share(work,
                                                                  comm.size,
                                                                  i)
            self.kk0, self.kkN = fair_share(work, comm.size, comm.rank)
            self.nproc_pk = 1
        else:
            if (comm.size > nkpts) and (comm.size%nkpts != 0):
                print(" # If nproc > nkpts, nproc must evenly divide number of"
                      "  k-points.")
                sys.exit()
            self.nproc_pk = comm.size // nkpts
            if kp_sym:
                self.kk0 = comm.rank // self.nproc_pk
                self.kkN = self.kk0 + 1
            else:
                self.mykpt = comm.rank // self.nproc_pk
                self.kk0 = self.mykpt * nkpts
                self.kkN = self.kk0 + nkpts

        self.ij0, self.ijN = fair_share(nmo_max*nmo_max, self.nproc_pk, comm.rank%self.nproc_pk)
        self.nij = self.ijN - self.ij0
        self.nkk = self.kkN - self.kk0
        if not kp_sym:
            self.n2k1 = numpy.zeros(self.nkk, dtype=numpy.int32)
            self.n2k2 = numpy.zeros(self.nkk, dtype=numpy.int32)
            cnt = 0
            for k in range(self.kk0, self.kkN):
                self.n2k1[cnt] = k // nkpts
                self.n2k2[cnt] = k % nkpts
                cnt += 1

class PartitionOld(object):

    def __init__(self, comm, maxvecs, nmo_tot, nmo_max, nkpts):
        # to do:
        # keep usage of n2k1,n2k2 but assign entire rows of the (k1,k2) matrix
        # also add extra partitioning of i*j matrix, along nproc_pk cores
        # algorithm is the same
        self.maxvecs = maxvecs * nmo_tot
        self.rank = comm.rank
        self.size = comm.size

        if comm.size <= nkpts:
            self.kkbounds = numpy.zeros(comm.size+1, dtype=numpy.int32)
            for i in range(comm.size):
                self.kkbounds[i], self.kkbounds[i+1] = fair_share(nkpts*nkpts,
                                                                  comm.size,
                                                                  i)
            self.kk0, self.kkN = fair_share(nkpts*nkpts, comm.size, comm.rank)
            self.nproc_pk = 1
        else:
            if (comm.size > nkpts) and (comm.size%nkpts != 0):
                print(" # If nproc > nkpts, nproc must evenly divide number of"
                      "  k-points.")
                sys.exit()
            self.nproc_pk = comm.size // nkpts
            if comm.rank == 0:
                print(" # Each kpoint is distributed accross {} mpi tasks."
                      .format(self.nproc_pk))
                sys.stdout.flush()
            self.mykpt = comm.rank // self.nproc_pk
            self.kk0 = self.mykpt * nkpts
            self.kkN = self.kk0 + nkpts

        self.ij0, self.ijN = fair_share(nmo_max*nmo_max,
                                        self.nproc_pk,
                                        comm.rank%self.nproc_pk)
        self.nij = self.ijN - self.ij0
        self.nkk = self.kkN - self.kk0
        self.n2k1 = numpy.zeros(self.nkk, dtype=numpy.int32)
        self.n2k2 = numpy.zeros(self.nkk, dtype=numpy.int32)
        cnt = 0
        for k in range(self.kk0, self.kkN):
            self.n2k1[cnt] = k // nkpts
            self.n2k2[cnt] = k % nkpts
            cnt += 1

def gen_orbital_products(cell, mydf, X, nmo_pk, ngs, part, kpts, nmo_max):
    # left (includes v(G)) and right pair densities
    nkpts = len(kpts)
    try:
        Xaoik = numpy.zeros((part.nkk,ngs,part.nij),
                            dtype=numpy.complex128)
        Xaolj = numpy.zeros((part.nkk,ngs,part.nij),
                            dtype=numpy.complex128)
    except:
        mem = part.nkk*ngs*part.nij*16*2 / 1024.0**3
        print(" # Problems allocating memory. Trying to allocate: {}"
              " GB.".format(mem))
        sys.exit()

    for k in range(part.nkk):
        k1 = part.n2k1[k]
        k2 = part.n2k2[k]
        i0 = part.ij0 // nmo_pk[k2]
        iN = part.ijN // nmo_pk[k2]
        if part.ijN % nmo_pk[k2] != 0:
            iN += 1
        if iN > nmo_pk[k1]:
            iN = nmo_pk[k1]
        pij = part.ij0 % nmo_pk[k2]
        n_ = min(part.ijN, nmo_pk[k1]*nmo_pk[k2]) - part.ij0
        X_t = X[k1][:,i0:iN].copy()
        Xaoik[k,:,0:n_] = mydf.get_mo_pairs_G((X_t,X[k2]),
                                            (kpts[k1],kpts[k2]),
                                            kpts[k2]-kpts[k1],
                                            compact=False)[:,pij:pij+n_]
        X_t = None
        Xaolj[k,:,:] = Xaoik[k,:,:]
        coulG = tools.get_coulG(cell, kpts[k2]-kpts[k1], mesh=mydf.mesh)
        Xaoik[k,:,:] *= (coulG*cell.vol/ngs**2).reshape(-1,1)
    t1 = time.process_time()()
    return Xaoik, Xaolj


class Cholesky(object):

    def __init__(self, part, kconserv, gtol_chol=1e-5, verbose=True):
        try:
            self.maxres_buff = numpy.zeros(5*part.size, dtype=numpy.float64)
        except MemoryError:
            print(" # Problems allocating memory for auxiliary structures for "
                  "Cholesky solver.")
            sys.exit()
        self.kconserv = kconserv
        self.gtol_chol = gtol_chol
        self.verbose = verbose

    def find_k3k4(self, nproc):
        vmax = 0
        for i in range(nproc):
            if self.maxres_buff[i*5+4] > vmax:
                vmax = self.maxres_buff[i*5+4]
                k3 = int(self.maxres_buff[i*5])
                k4 = int(self.maxres_buff[i*5+1])
                i3 = int(self.maxres_buff[i*5+2])
                i4 = int(self.maxres_buff[i*5+3])
        return k3, k4, i3, i4, vmax

    def generate_diagonal(self, Xaoik, Xaolj, part, nmo_pk):
        maxv = 0
        try:
            residual = numpy.zeros((part.nkk,part.nij), dtype=numpy.float64)
        except MemoryError:
            print(" # Problems allocating memory for residuals array.")
            sys.exit()
        for k in range(part.nkk):
            k1 = part.n2k1[k]
            k2 = part.n2k2[k]
            for ij in range(part.nij):
                if (ij+part.ij0) >= nmo_pk[k1]*nmo_pk[k2]:
                    break
                intg = numpy.dot(Xaoik[k,:,ij], Xaolj[k,:,ij].conj())
                if (intg.real < 0) | (abs(intg.imag) > 1e-9):
                    i = (ij+part.ij0) // part.nmo_pk[k2]
                    j = (ij+part.ij0) % part.nmo_pk[k2]
                    print(" # ERROR: Negative or complex diagonal term: "
                          "{} {} {} {} {:13.8e}".format(k1,i,k2,j,intg))
                residual[k,ij] = intg.real
                if abs(intg) > maxv:
                    maxv = abs(intg)
                    k1max = k1
                    k2max = k2
                    i1max = (ij+part.ij0) // nmo_pk[k2]
                    i2max = (ij+part.ij0) % nmo_pk[k2]

        return residual, k1max, k2max, i1max, i2max, maxv

    def run(self, comm, Xaoik, Xaolj, part, kpts, nmo_pk, nmo_max, Qi, gmap):
        # Setup residual matrix.
        ngs = Xaoik.shape[1]
        nkpts = len(kpts)
        t0 = time.process_time()()
        residual, k1max, k2max, i1max, i2max, maxv = (
                self.generate_diagonal(Xaoik, Xaolj, part, nmo_pk)
                )
        t1 = time.process_time()()
        if part.rank == 0 and self.verbose:
            print(" # Time to generate diagonal (initial residual):"
                  " {:.2e}".format(t1-t0))
            sys.stdout.flush()
        comm.Allgather(numpy.array([k1max,k2max,i1max,i2max,maxv],
                       dtype=numpy.float64),self.maxres_buff)
        vmax = 0
        k3, k4, i3, i4, vmax = self.find_k3k4(comm.size)
        try:
            done = numpy.zeros((nkpts,nkpts,nmo_max,nmo_max), numpy.int32)
            maxresidual = numpy.zeros(part.maxvecs, dtype=numpy.float64)
            cholvecs = numpy.zeros((part.nkk,part.nij,part.maxvecs),
                                        dtype=numpy.complex128)
            Xkl = numpy.zeros(ngs, dtype=numpy.complex128)
            Xkl0 = numpy.zeros(ngs+part.maxvecs, dtype=numpy.complex128)
            Vbuff = numpy.zeros(part.maxvecs, dtype=numpy.complex128)
        except MemoryError:
            print(" # Problems allocating memory for auxiliary structures for "
                  "Cholesky solver.")
        done[k3,k4,i3,i4] = 1

        tstart = time.process_time()()
        if comm.rank == 0:
            sys.stdout.flush()

        more = True   # for parallel
        numv = 0
        if comm.rank == 0:
            header = ["iteration", "max_residual", "total_time",
                      "time_k3k4", "time_comp_cholv", "time_buff"]
            if self.verbose:
                print(format_fixed_width_strings(header))
        while more:
            t0 = time.process_time()()
            # stop condition
            if numv >= part.maxvecs:
                print(" Too many vectors needed to converge. Increase maximum "
                      "number of vectors.")
                break
            kkmax = k3*nkpts + k4
            if comm.size <= nkpts:
                ipr = bisect(part.kkbounds[1:comm.size+1], kkmax)
                if (kkmax >= part.kk0) & (kkmax < part.kkN):
                    assert(comm.rank == ipr)
            else:
                i34 = i3*nmo_pk[k4] + i4
                for i in range(part.nproc_pk):
                    ij0_, ijN_ = fair_share(nmo_max*nmo_max, part.nproc_pk,i)
                    if i34 < ijN_:
                        ipr = k3*part.nproc_pk + i
                        break

            # bcast Xaolj/CV[k3,k4,i34,0:numv]
            if comm.rank == ipr:
                assert(((i3*nmo_pk[k4]+i4) >= part.ij0) &
                       ((i3*nmo_pk[k4]+i4) < part.ijN))
                Xkl0[0:ngs] = Xaolj[kkmax-part.kk0,0:ngs,i3*nmo_pk[k4]+i4-part.ij0]
                Xkl0[ngs:ngs+numv] = (
                        cholvecs[kkmax-part.kk0,i3*nmo_pk[k4]+i4-part.ij0,0:numv]
                        )
                Vbuff[0:numv] = (
                        cholvecs[kkmax-part.kk0,i3*nmo_pk[k4]+i4-part.ij0,0:numv]
                        )
                comm.Bcast(Xkl0[0:ngs+numv], root=ipr)
            else:
                comm.Bcast(Xkl0[0:ngs+numv], root=ipr)
                Vbuff[0:numv] = Xkl0[ngs:ngs+numv]

            # add new Cholesky vector
            # 1. evaluate new column (ik|imax,kmax)

            t1 = time.process_time()()
            tadd = 0.0

            for k in range(part.nkk):
                k1 = part.n2k1[k]
                k2 = part.n2k2[k]
                if k3 == self.kconserv[k1,k2,k4]:
                    q1 = kpts[k2] - kpts[k1] + kpts[k3] - kpts[k4]
                    if numpy.sum(abs(q1)) > 1e-9:
                        t_ = time.process_time()()
                        ip = -1
                        for ii in range(27):
                            if numpy.sum(numpy.linalg.norm(q1-Qi[ii,:])) < 1e-12:
                                ip = ii
                                break
                        if ip < 0:
                            print(" # Could not find Q: {} {} ".format(q1,Qi))
                            sys.exit()
                        for ix in range(ngs):
                            Xkl[ix] = Xkl0[gmap[ip,ix]]
                        tadd += time.process_time()() - t_
                    else:
                        Xkl[0:ngs] = Xkl0[0:ngs]
                    n_ = min(nmo_pk[k1]*nmo_pk[k2], part.ijN) - part.ij0
                    cholvecs[k,0:n_,numv] = numpy.dot(Xaoik[k,:,0:n_].T,
                                                      Xkl.conj())

            t2 = time.process_time()()

            # 2. subtract projection along previous components
            cholvecs[:,:,numv] -= numpy.dot(cholvecs[:,:,0:numv],
                                            Vbuff[0:numv].conj())
            cholvecs[:,:,numv] /= math.sqrt(vmax)

            # update residual
            residual -= (cholvecs[:,:,numv]*cholvecs[:,:,numv].conj()).real

            maxv = 0
            for k in range(part.nkk):
                k1 = part.n2k1[k]
                k2 = part.n2k2[k]
                for ij in range(part.nij):
                    if (ij+part.ij0) >= nmo_pk[k1]*nmo_pk[k2]:
                        break
                    if abs(residual[k,ij]) > maxv :
                        maxv = abs(residual[k,ij])
                        k1max = k1
                        k2max = k2
                        i1max = (ij+part.ij0) // nmo_pk[k2]
                        i2max = (ij+part.ij0) % nmo_pk[k2]

            t3 = time.process_time()()

            # assemble full CV on head node
            comm.Allgather(numpy.array([k1max,k2max,i1max,i2max,maxv],
                                       dtype=numpy.float64), self.maxres_buff)
            k3, k4, i3, i4, vmax = self.find_k3k4(comm.size)

            t4 = time.process_time()()

            # only root keeps track of residual and I/O
            if part.rank == 0:

                maxresidual[numv] = abs(vmax)

                # print and evaluate stop condition
                output = [vmax, t4-t0, t3-t2, t2-t1, t1-t0]
                if self.verbose:
                    print("{:17d} ".format(numv)+format_fixed_width_floats(output))
                # print("{:8d}  {:13.8e}".format(numv, vmax))
                tstart = time.process_time()()

                if numv%100 == 0:
                    sys.stdout.flush()

            numv += 1
            if vmax < self.gtol_chol:
                more = False

            if (done[k3,k4,i3,i4] > 0) and more:
                print("Error in Cholesky decomposition. "
                      "done[imax,kmax]>0.",k3,k4,i3,i4,vmax)
                sys.exit()
            done[k3,k4,i3,i4]=1

            t6 = time.process_time()()
        return cholvecs[:,:,:numv]

def write_h1(h5grp, intgs, npk, ik2n, gtol=1e-6):

    nkpts = intgs.shape[0]
    assert(npk.shape[0] == nkpts)
    assert(ik2n.shape[1] == nkpts)
    cnt=0
    for ki in range(nkpts):
        ij = 0
        for i in range(npk[ki]):
            for j in range(0, i+1):
                if abs(intgs[ki,ij])>gtol:
                    cnt+=1
                ij+=1

    H1_indx = h5grp.create_dataset("H1_indx", (2*cnt,), dtype='int32')
    if intgs.dtype == numpy.dtype(complex):
        H1 = h5grp.create_dataset("H1", (cnt,2), dtype=numpy.float64)
        cnt=0
        for ki in range(nkpts):
            ij = 0
            for i in range(npk[ki]):
                for j in range(0, i+1):
                    if abs(intgs[ki,ij])>gtol:
                        H1_indx[2*cnt] = ik2n[i,ki]
                        H1_indx[2*cnt+1] = ik2n[j,ki]
                        H1[cnt,0] = intgs[ki,ij].real
                        H1[cnt,1] = intgs[ki,ij].imag
                        cnt+=1
                    ij+=1
    else:
        H1 = h5grp.create_dataset("H1", (cnt,), dtype=numpy.float64)
        cnt=0
        for ki in range(nkpts):
            ij = 0
            for i in range(npk[ki]):
                for j in range(0, i+1):
                    if abs(intgs[ki,ij])>gtol:
                        H1_indx[2*cnt] = ik2n[i,ki]
                        H1[cnt] = intgs[ki,ij]
                        cnt+=1
                    ij+=1
    return cnt

# input: CVMat[nkk,nij,nvec]
def write_cholMat_singleWriter(h5grp, h5nblk, CV, nopk, n2k1, n2k2, ij0, ik2n,
                               gtol=1e-6, MaxIntgs=2000000):

    assert(len(ik2n.shape)==2)
    assert(len(n2k1.shape)==1)
    assert(len(n2k2.shape)==1)
    assert(len(CV.shape)==3)

    nkpts = len(nopk)
    nmo_max = ik2n.shape[0]
    nkk = CV.shape[0]
    nij = CV.shape[1]
    nvec = CV.shape[2]
    nmo_tot=numpy.sum(nopk)

    gamma = (CV.dtype == numpy.dtype(float))

    assert(n2k1.shape[0]==nkk)
    assert(n2k2.shape[0]==nkk)
    assert(ik2n.shape[1]==nkpts)
    assert(ij0+nij <= nmo_max*nmo_max)

    cnt=0
    for k in range(nkk):
        k1 = n2k1[k]
        k2 = n2k2[k]
        for ij in range(nij):
            if (ij+ij0) >= nopk[k1]*nopk[k2]:
                break
            for nv in range(nvec):
                if abs(CV[k,ij,nv]/nkpts) > gtol:
                    cnt+=1

    vcnt = numpy.array([0])
    if cnt==0:
        return vcnt,0

    nblk = cnt//MaxIntgs
    if cnt%MaxIntgs > 0:
        nblk+=1
    vcnt = numpy.zeros(nblk,dtype="int32")

    if gamma:
        V2 = numpy.zeros((MaxIntgs,),dtype=numpy.float64)
    else:
        V2 = numpy.zeros((MaxIntgs,2),dtype=numpy.float64)
    V2_indx = numpy.zeros(2*MaxIntgs,dtype="int32")

    fct = 1.0/math.sqrt(nkpts*1.0)
    cnt=0
    ngrp = 0
    if gamma:
        for k in range(nkk):
            k1 = n2k1[k]
            k2 = n2k2[k]
            for ij in range(nij):
                if (ij+ij0) >= nopk[k1]*nopk[k2]:
                    break
                for nv in range(nvec):
                    if abs(CV[k,ij,nv]/nkpts) > gtol:
                        i = (ij+ij0) // nopk[k2]
                        j = (ij+ij0) % nopk[k2]
                        V2_indx[2*cnt] = ik2n[i,k1]*nmo_tot+ik2n[j,k2]
                        V2_indx[2*cnt+1] = nv
                        V2[cnt] = CV[k,ij,nv]*fct
                        cnt+=1
                        if cnt==MaxIntgs:
                            ix = str(h5nblk+ngrp)
                            h5dset = h5grp.create_dataset("vals_"+ix,
                                                          data=V2[0:cnt],
                                                          chunks=True,
                                                          maxshape=(None,),
                                                          compression="gzip")
                            h5dset_indx = (
                                    h5grp.create_dataset("index_"+ix,
                                                         data=V2_indx[0:2*cnt],
                                                         chunks=True,
                                                         maxshape=(None,),
                                                         compression="gzip")
                                    )
                            vcnt[ngrp]=cnt
                            cnt=0
                            ngrp+=1
    else:
        for k in range(nkk):
            k1 = n2k1[k]
            k2 = n2k2[k]
            for ij in range(nij):
                if (ij+ij0) >= nopk[k1]*nopk[k2]:
                    break
                for nv in range(nvec):
                    if abs(CV[k,ij,nv]/nkpts) > gtol:
                        i = (ij+ij0) // nopk[k2]
                        j = (ij+ij0) % nopk[k2]
                        V2_indx[2*cnt] = ik2n[i,k1]*nmo_tot+ik2n[j,k2]
                        V2_indx[2*cnt+1] = nv
                        V2[cnt,0] = CV[k,ij,nv].real*fct
                        V2[cnt,1] = CV[k,ij,nv].imag*fct
                        cnt+=1
                        if cnt==MaxIntgs:
                            ix = str(h5nblk+ngrp)
                            h5dset = h5grp.create_dataset("vals_"+ix,
                                                          data=V2[0:cnt,0:2],
                                                          chunks=True,
                                                          maxshape=(None,2),
                                                          compression="gzip")
                            h5dset_indx = (
                                h5grp.create_dataset("index_"+ix,
                                                     data=V2_indx[0:2*cnt],
                                                     chunks=True,
                                                     maxshape=(None,),
                                                     compression="gzip")
                                )
                            vcnt[ngrp]=cnt
                            cnt=0
                            ngrp+=1


    if cnt>0:
        ix = str(h5nblk+ngrp)
        if gamma:
            h5dset = h5grp.create_dataset("vals_"+ix,
                                          data=V2[0:cnt],
                                          chunks=True,
                                          maxshape=(None,),
                                          compression="gzip")
        else:
            h5dset = h5grp.create_dataset("vals_"+ix,
                                          data=V2[0:cnt,0:2],
                                          chunks=True,
                                          maxshape=(None,2),
                                          compression="gzip")
        h5dset_indx = h5grp.create_dataset("index_"+ix,
                                           data=V2_indx[0:2*cnt],
                                           chunks=True,
                                           maxshape=(None,),
                                           compression="gzip")
        vcnt[ngrp]=cnt
        ngrp+=1

    assert(ngrp==nblk)
    return vcnt,ngrp
