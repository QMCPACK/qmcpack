import h5py
import math
import numpy
import sys
import time
import os
from pyscf import lib
from pyscf.pbc import tools, df
from afqmctools.hamiltonian.supercell import generate_grid_shifts, Partition
from afqmctools.utils.parallel import bisect, fair_share
from afqmctools.utils.io import (
        format_fixed_width_floats,
        format_fixed_width_strings,
        to_qmcpack_complex
        )
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk

def alloc_helper(shape, dtype=numpy.float64, name='array', verbose=False):
    """Numpy array allocator helper.

    Parameters
    ----------
    shape : tuple
        Array shape.
    dtype : numpy data type
        Array dtype. Default: numpy.float64
    name : string
        Array name. Default: array.
    verbose : bool
        Print array information. Default: False.

    Returns
    -------
    array : :class:`numpy.ndarray`
        Array.
    """
    mem = numpy.prod(shape)*numpy.dtype(dtype).itemsize / (1024.0**3)
    if verbose:
        print(" # Allocating array: {:s} with size {:.2e} GB".format(name, mem))
        print(" # Shape {}.".format(shape))
    try:
        return numpy.zeros(shape, dtype=dtype)
    except MemoryError:
        print(" # Problems allocating array: {:s} with size {:.2e}".format(name, mem))
        sys.exit()



def write_hamil_kpoints(comm, scf_data, hamil_file, chol_cut,
                        verbose=True, cas=None, max_vecs=20,
                        ortho_ao=False, exxdiv='ewald', nelec=None,
                        phdf=False):
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
    nao = cell.nao_nr()
    qk_to_k2, kminus = construct_qk_maps(cell, kpts)


    h5file = FileHandler(comm,hamil_file,"w",phdf)
    h5file.grp = h5file.h5f.create_group("Hamiltonian")
    h5file.grp_v2 = h5file.h5f.create_group("Hamiltonian/KPFactorized")
    if h5file.error:
        sys.exit()
    write_basic(comm, cell, kpts, hcore, h5file, X, nmo_pk,
                qk_to_k2, kminus, verbose=verbose, nelec=nelec)

    if comm.rank == 0 and verbose:
        print(" # Time to reach Cholesky: {:13.8e} s.".format(time.process_time()()-tstart))
        sys.stdout.flush()
    tstart = time.process_time()()

    solver = KPCholesky(comm, cell, kpts, max_vecs, nmo_pk,
                        qk_to_k2, kminus, gtol_chol=chol_cut,
                        verbose=verbose)
    solver.run(comm, X, h5file)
    if comm.rank == 0 and verbose:
        print(" # Time to perform Cholesky: {:13.8e} s.".format(time.process_time()()-tstart))
        sys.stdout.flush()

    comm.barrier()
    if not phdf and comm.rank==0:
        comm.barrier()
        nkpts = len(kpts)
        for rk in range(1,comm.size):
            h5f2 = FileHandler(comm,"rank"+str(rk)+"_"+hamil_file,"r",False)
            for Q in range(nkpts):
                if Q > kminus[Q]:
                    continue
                Ldim = h5f2.h5f["/Hamiltonian/KPFactorized/Ldim"+str(Q)][:]
                nkk = Ldim[0]
                nij = Ldim[1]
                kk0 = Ldim[2]
                ij0 = Ldim[3]
                ijN = Ldim[4]
                numv = Ldim[5]
                LQ2 = h5f2.h5f["/Hamiltonian/KPFactorized/L"+str(Q)][:]
                LQ2 = numpy.reshape(LQ2,(nkk,nij*numv,2))
                h5file.grp_v2["L"+str(Q)][kk0:kk0+nkk,ij0*numv:ijN*numv,:] = LQ2[:,:,:]
            h5f2.close()
            os.remove("rank"+str(rk)+"_"+hamil_file)
        h5file.close()
    else:
        h5file.close()
        comm.barrier()

    comm.barrier()

def write_rhoG_kpoints(comm, scf_data, hdf_file, Gcut,
                        verbose=True, 
                        ortho_ao=False, phdf=False):
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
    nao = cell.nao_nr()
    qk_to_k2, kminus = construct_qk_maps(cell, kpts)


    h5file = FileHandler(comm,hdf_file,"w",phdf)
    h5file.grp = h5file.h5f.create_group("rhoG")
    if h5file.error:
        sys.exit()
    write_basic(comm, cell, kpts, hcore, h5file, X, nmo_pk,
                qk_to_k2, kminus, verbose=verbose, nelec=nelec)

    solver = KPCholesky(comm, cell, kpts, 0, nmo_pk,
                        qk_to_k2, kminus, gtol_chol=0,
                        verbose=verbose)
#    solver.run(comm, X, h5file)

    comm.barrier()
    if not phdf and comm.rank==0:
        comm.barrier()
        nkpts = len(kpts)
        for rk in range(1,comm.size):
            h5f2 = FileHandler(comm,"rank"+str(rk)+"_"+hdf_file,"r",False)
            for Q in range(nkpts):
                if Q > kminus[Q]:
                    continue
                Ldim = h5f2.h5f["/rhoG/dim"+str(Q)][:]
                nkk = Ldim[0]
                nij = Ldim[1]
                kk0 = Ldim[2]
                ij0 = Ldim[3]
                ijN = Ldim[4]
                numG = Ldim[5]
                LQ2 = h5f2.h5f["/rhoG/rho"+str(Q)][:]
#                LQ2 = numpy.reshape(LQ2,(nkk,nij*numv,2))
#                h5file.grp["rho"+str(Q)][kk0:kk0+nkk,ij0*numv:ijN*numv,:] = LQ2[:,:,:]
            h5f2.close()
            os.remove("rank"+str(rk)+"_"+hdf_file)
        h5file.close()
    else:
        h5file.close()
        comm.barrier()

    comm.barrier()


class FileHandler:
    def __init__(self, comm, filename,ftype="w",phdf=False,
                g1="Hamiltonian",g2="KPFactorized"):
        self.phdf=phdf
        if phdf:
            try:
                self.h5f = h5py.File(filename, ftype,
                                   driver='mpio', comm=comm)
                self.h5f.atomic = False
                self.error = 0
            except:
                if comm.rank == 0:
                    print("Parallel hdf5 required.")
                self.error= 1
        else:
            try:
                if comm.rank == 0:
                    self.h5f = h5py.File(filename, ftype)
                else:
                    self.h5f = h5py.File("rank"+str(comm.rank)+"_"+filename, ftype)
                self.error = 0
            except:
                if comm.rank == 0:
                    print("Error creating hdf5 file.")
                self.error= 1
#        if ftype=="w":
#            self.grp = self.h5f.create_group(g1)
#            self.grp_v2 = self.h5f.create_group(g1+"/"+g2)

    def close(self):
        self.h5f.close()

def write_basic(comm, cell, kpts, hcore, h5file, X, nmo_pk, qk_to_k2, kminus,
                exxdiv='ewald', verbose=False, nelec=None):
    """Write basic system information including one-body Hamiltonian to file.
    """
    nkpts = len(kpts)
    nmo_tot = numpy.sum(nmo_pk)
    if nelec is not None:
        nup, ndown = nelec
    else:
        nup = nkpts * (cell.nelectron+cell.spin) // 2
        ndown = nkpts*cell.nelectron - nup
    nelectron = nup + ndown

    # zero electron energies, including madelung term
    if comm.rank == 0:
        e0 = nkpts * cell.energy_nuc()
        if exxdiv == 'ewald':
            madelung = tools.pbc.madelung(cell, kpts)
            emad = -0.5*nelectron*madelung
            e0 += emad
            if verbose:
                print(" # Adding ewald correction to the energy:"
                      " {}".format(emad))
        sys.stdout.flush()

    comm.barrier()

    dims_ = h5file.grp.create_dataset("dims", (8,), dtype=numpy.int32)
    h5file.grp.create_dataset("ComplexIntegrals",
                              data=numpy.array([1]),
                              dtype=numpy.int32)
    et_ = h5file.grp.create_dataset("Energies", (2,), dtype=numpy.float64)
    kp_ = h5file.grp.create_dataset("KPoints", (nkpts,3), dtype=numpy.float64)
    nmo_pk_ = h5file.grp.create_dataset("NMOPerKP", (nkpts,), dtype=numpy.int32)
    q_ = h5file.grp.create_dataset("QKTok2", (nkpts,nkpts), dtype=numpy.int32)
    kminus_ = h5file.grp.create_dataset("MinusK", (nkpts,), dtype=numpy.int32)
    if comm.rank == 0:
        dims_[:] = numpy.array([0, 0, nkpts, nmo_tot, nup, ndown, 0, 0])
        et_[:] = numpy.array([e0, 0])
        kp_[:,:] = kpts[:,:]
        nmo_pk_[:] = nmo_pk[:]
        kminus_[:] = kminus[:]
        q_[:,:] = qk_to_k2[:,:]

    comm.barrier()

    for ki in range(nkpts):
        if comm.rank == 0:
            h1 = numpy.dot(X[ki][:,0:nmo_pk[ki]].T.conj(),
                           numpy.dot(hcore[ki], X[ki][:,0:nmo_pk[ki]]))
            write_kp_h1(h5file.grp, ki, nmo_pk[ki], h1)
        else:
            write_kp_h1(h5file.grp, ki, nmo_pk[ki], None)

    comm.barrier()

class KPCholesky(object):

    def __init__(self, comm, cell, kpts, maxvecs, nmo_pk, QKToK2, kminus,
                 gtol_chol=1e-5, verbose=True):
        self.QKToK2 = QKToK2
        self.kminus = kminus
        nmo_max = numpy.max(nmo_pk)
        nmo_tot = numpy.sum(nmo_pk)
        if comm.rank == 0 and verbose:
            print(" # Total number of orbitals: {}".format(nmo_tot))
            sys.stdout.flush()
        assert(nmo_max is not None)
        nkpts = len(kpts)
        self.nkpts = nkpts
        self.kpts = kpts
        self.cell = cell
        self.part = Partition(comm, maxvecs, nmo_tot, nmo_max, nkpts,
                              kp_sym=True)
        maxvecs = maxvecs*nmo_max
        self.maxvecs = maxvecs
        self.nmo_max = nmo_max
        self.gmap, self.Qi, self.ngs = generate_grid_shifts(cell)
        if comm.rank == 0 and verbose:
            print(" # Each kpoint is distributed accross {} mpi tasks."
                  .format(self.part.nproc_pk))
            sys.stdout.flush()
        try:
            self.maxres_buff = numpy.zeros(5*comm.size, dtype=numpy.float64)
        except MemoryError:
            print(" # Problems allocating memory for auxiliary structures for "
                  "Cholesky solver.")
            sys.exit()
        self.df = df.FFTDF(cell,kpts)
        tstart = time.process_time()()

        self.nmo_pk = nmo_pk
        self.gtol_chol = gtol_chol
        self.verbose = verbose

    def find_k3k4(self, nproc):
        vmax = 0
        for i in range(nproc):
            if self.maxres_buff[i*5+4] > vmax :
                vmax = self.maxres_buff[i*5+4]
                k3 = int(self.maxres_buff[i*5])
                k4 = int(self.maxres_buff[i*5+1])
                i3 = int(self.maxres_buff[i*5+2])
                i4 = int(self.maxres_buff[i*5+3])
        return k3, k4, i3, i4, vmax

    def generate_orbital_products(self, Q, X, Xaoik, Xaolj):
        kpts = self.kpts
        cell = self.cell
        ngs = self.ngs
        for K in range(self.part.nkk):
            k1 = K + self.part.kk0
            k2 = self.QKToK2[Q][k1]
            if self.part.ij0 > self.nmo_pk[k1]*self.nmo_pk[k2]:
                continue
            i0 = self.part.ij0 // self.nmo_pk[k2]
            iN = self.part.ijN // self.nmo_pk[k2]
            if self.part.ijN % self.nmo_pk[k2] != 0:
                iN += 1
            if iN > self.nmo_pk[k1]:
                iN = self.nmo_pk[k1]
            pij = self.part.ij0 % self.nmo_pk[k2]
            n_ = min(self.part.ijN,self.nmo_pk[k1]*self.nmo_pk[k2]) - self.part.ij0
            X_t = X[k1][:,i0:iN].copy()
            Xaoik[K,:,0:n_] = self.df.get_mo_pairs_G((X_t,X[k2].copy()),
                                                  (kpts[k1],kpts[k2]),
                                                  (kpts[k2]-kpts[k1]),
                                                   compact=False)[:,pij:pij+n_]
            Xaolj[K,:,:] = Xaoik[K,:,:]
            coulG = tools.get_coulG(cell, kpts[k2]-kpts[k1], mesh=self.df.mesh)
            Xaoik[K,:,:] *= (coulG*cell.vol/ngs**2).reshape(-1,1)

    def generate_diagonal(self, Q, Xaoik, Xaolj):
        maxv = 0
        residual = numpy.zeros((self.part.nkk,self.part.nij),dtype=numpy.float64)
        k1max = -1
        k2max = -1
        i1max = -1
        i2max = -1
        for k in range(self.part.nkk):
            k1 = k + self.part.kk0
            k2 = self.QKToK2[Q][k1]
            for ij in range(self.part.nij):
                if (ij+self.part.ij0) >= self.nmo_pk[k1]*self.nmo_pk[k2]:
                    break
                intg = numpy.dot(Xaoik[k,:,ij],
                                 Xaolj[k,:,ij].conj())
                if (intg.real < 0) | (abs(intg.imag) > 1e-9):
                    i = (ij+self.part.ij0) // self.nmo_pk[k2]
                    j = (ij+self.part.ij0) % self.nmo_pk[k2]
                    print("ERROR!!! Negative or complex diagonal term: "
                          "{} {} {} {} {:13.8e}".format(k1,i,k2,j,intg))
                residual[k,ij] = intg.real
                if abs(intg) > maxv:
                    maxv = abs(intg)
                    k1max = k1
                    k2max = k2
                    i1max = (ij+self.part.ij0) // self.nmo_pk[k2]
                    i2max = (ij+self.part.ij0) % self.nmo_pk[k2]

        return residual, k1max, k2max, i1max, i2max, maxv

    def run(self, comm, X, h5file):
        # Unpack for convenience.
        ngs = self.ngs
        nmo_max = self.nmo_max
        nkpts = self.nkpts
        part = self.part
        nmo_pk = self.nmo_pk
        QKToK2 = self.QKToK2
        kpts = self.kpts
        # Setup residual matrix.
        mem = 2*16*nkpts*ngs*nmo_max**2 / 1024**3
        if self.verbose and comm.rank == 0:
            print(" # Approx total memory required for orbital products: "
                  "{:.2e} GB.".format(mem))
        mem_verbose = comm.rank == 0 and self.verbose > 1
        if mem_verbose:
            print(" # Approx memory per MPI task for auxiliary structures.")
        Xaoik = alloc_helper((part.nkk,ngs,part.nij), dtype=numpy.complex128,
                             name='Xaoik', verbose=mem_verbose)
        Xaolj = alloc_helper((part.nkk,ngs,part.nij), dtype=numpy.complex128,
                             name='Xaolj', verbose=mem_verbose)
        done = alloc_helper((nkpts,nmo_max,nmo_max), numpy.int32,
                            name='done', verbose=mem_verbose)
        maxresidual = alloc_helper((part.maxvecs,), dtype=numpy.float64,
                                   name='maxresidual', verbose=mem_verbose)
        cholvecs = alloc_helper((part.nkk,part.nij,part.maxvecs), dtype=numpy.complex128,
                                  name='cholvecs', verbose=mem_verbose)
        Xkl = alloc_helper((ngs,), dtype=numpy.complex128,
                           name='xkl', verbose=mem_verbose)
        Xkl0 = alloc_helper((ngs+part.maxvecs,), dtype=numpy.complex128,
                            name='xkl0', verbose=mem_verbose)
        Vbuff = alloc_helper((part.maxvecs,), dtype=numpy.complex128,
                             name='Vbuff', verbose=mem_verbose)
        num_cholvecs = alloc_helper((nkpts,), dtype=numpy.int32,
                                    name='num_cholvecs', verbose=mem_verbose)
        header = ["iteration", "max_residual", "total_time",
                  "time_k3k4", "time_comp_cholv", "time_buff"]

        for Q in range(nkpts):
            t0 = time.process_time()()
            if Q > self.kminus[Q]:
                continue
            if comm.rank == 0 and self.verbose:
                print(" # Calculating factorization for momentum: {}".format(Q))
                print(" # Generating orbital products")
                sys.stdout.flush()

            t1 = time.process_time()()

            maxresidual[:] = 0
            done[:,:,:] = 0
            cholvecs[:,:,:] = 0
            start = time.time()
            self.generate_orbital_products(Q, X, Xaoik, Xaolj)
            if self.verbose and comm.rank == 0:
                print(" # Time to generate orbital products: "
                      "{:13.8e}".format(time.time()-start))
            residual, k1max, k2max, i1max, i2max, maxv = (
                    self.generate_diagonal(Q, Xaoik, Xaolj)
                    )

            buff = numpy.array([k1max,k2max,i1max,i2max,maxv],
                               dtype=numpy.float64)
            comm.Allgather(buff, self.maxres_buff)
            k3, k4, i3, i4, vmax = self.find_k3k4(comm.size)
            done[k3,i3,i4] = 1

            tstart = time.process_time()()
            if comm.rank == 0:
                sys.stdout.flush()

            cnt = 0
            vmaxold = vmax
            more = True   # for parallel
            numv = 0
            while more:

                t0 = time.process_time()()
                # stop condition
                if comm.rank == 0 and self.verbose:
                    if numv == 0:
                        print(format_fixed_width_strings(header))
                if numv >= self.maxvecs:
                    print(" Too many vectors needed to converge. "
                          "Increase maximum number of vectors.")
                    break

                if comm.size <= nkpts:
                    ipr = bisect(part.kkbounds[1:comm.size+1],k3)
                    if (k3 >= part.kk0) and (k3 < part.kkN):
                        assert(comm.rank==ipr)
                else:
                    i34 = i3*nmo_pk[k4] + i4
                    for i in range(part.nproc_pk):
                        ij0_, ijN_ = fair_share(nmo_max*nmo_max,
                                                part.nproc_pk,
                                                i)
                        if i34 < ijN_:
                            ipr = k3*part.nproc_pk + i
                            break

                if comm.rank == ipr:
                    assert((i3*nmo_pk[k4]+i4 >= part.ij0))
                    assert((i3*nmo_pk[k4]+i4 < part.ijN))
                    Xkl0[0:ngs] = Xaolj[k3-part.kk0,
                                        0:ngs,
                                        i3*nmo_pk[k4]+i4-part.ij0]
                    Xkl0[ngs:ngs+numv] = cholvecs[k3-part.kk0,
                                                  i3*nmo_pk[k4]+i4-part.ij0,
                                                  0:numv]
                    Vbuff[0:numv] = cholvecs[k3-part.kk0,
                                             i3*nmo_pk[k4]+i4-part.ij0,
                                             0:numv]
                    comm.Bcast(Xkl0[0:ngs+numv], root=ipr)
                else:
                    comm.Bcast(Xkl0[0:ngs+numv], root=ipr)
                    Vbuff[0:numv] = Xkl0[ngs:ngs+numv]

                # add new Cholesky vector
                # 1. evaluate new column (ik|imax,kmax)

                t1 = time.process_time()()
                tadd = 0.0
                for k in range(part.nkk):
                    k1 = k + part.kk0
                    k2 = QKToK2[Q][k1]
                    if part.ij0 > nmo_pk[k1]*nmo_pk[k2]:
                        continue
                    if numpy.sum(abs(kpts[k2]-kpts[k1]+kpts[k3]-kpts[k4])) > 1e-9:
                        t_ = time.process_time()()
                        q1 = kpts[k2]-kpts[k1]+kpts[k3]-kpts[k4]
                        ip = -1
                        for ii in range(27):
                            if numpy.sum(numpy.linalg.norm(q1-self.Qi[ii,:])) < 1e-12:
                                ip = ii
                                break
                        if ip < 0:
                            print("Could not find Q: {} {}".format(q1,self.Qi))
                            sys.exit()
                        for ix in range(ngs):
                            Xkl[ix] = Xkl0[self.gmap[ip,ix]]
                        tadd += time.process_time()() - t_
                    else:
                        Xkl[0:ngs] = Xkl0[0:ngs]
                    n_ = min(nmo_pk[k1]*nmo_pk[k2], part.ijN) - part.ij0
                    cholvecs[k,0:n_,numv] = numpy.dot(Xaoik[k,:,0:n_].T,
                                                      Xkl.conj())
                t2 = time.process_time()()

                # 2. substract projection along previous components
                cholvecs[:,:,numv] -= numpy.dot(cholvecs[:,:,0:numv],
                                                Vbuff[0:numv].conj())
                cholvecs[:,:,numv] /= math.sqrt(vmax)

                # update residual
                residual -= (cholvecs[:,:,numv]*cholvecs[:,:,numv].conj()).real

                maxv = 0
                k1max = -1
                k2max = -1
                i1max = -1
                i2max = -1
                for k in range(part.nkk):
                    k1 = k + part.kk0
                    k2 = self.QKToK2[Q][k1]
                    # if ij0 > nmo_pk[k1]*nmo_pk[k2]:
                        # continue
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
                buff = numpy.array([k1max,k2max,i1max,i2max,maxv],
                                   dtype=numpy.float64)
                comm.Allgather(buff, self.maxres_buff)
                k3, k4, i3, i4, vmax = self.find_k3k4(comm.size)

                t4 = time.process_time()()

                # only root keeps track of residual and I/O
                if comm.rank == 0:

                    maxresidual[numv] = abs(vmax)
                    # print and evaluate stop condition
                    output = [vmax, t4-t0, t3-t2, t2-t1, t1-t0]
                    if self.verbose:
                        print("{:17d} ".format(numv)+format_fixed_width_floats(output))
                    tstart = time.process_time()()

                    if numv%100 == 0:
                        sys.stdout.flush()

                numv += 1
                if vmax < self.gtol_chol:
                    more = False

                if done[k3,i3,i4] > 0 and more:
                    # TODO: What is this?
                    print("Error in Cholesky decomposition. "
                          "done[imax,kmax]>0.",k3,k4,i3,i4,vmax)
                    sys.exit()
                done[k3,i3,i4] = 1

                t6 = time.process_time()()

            comm.barrier()
            num_cholvecs[Q] = numv

            if h5file.phdf or comm.rank==0:
                LQ = h5file.grp_v2.create_dataset("L"+str(Q),
                                              (nkpts,nmo_max*nmo_max*numv,2),
                                              dtype=numpy.float64)
                # cholvecs[nkk,nij,maxvecs]
                for kk in range(part.nkk):
                    T_ = to_qmcpack_complex(cholvecs[kk,:,0:numv].copy())
                    T_ = numpy.reshape(T_,(-1,2))/math.sqrt(nkpts*1.0)
                    LQ[kk+part.kk0,part.ij0*numv:part.ijN*numv,:] = T_
                    T_ = None
            else:
                Ldim = h5file.grp_v2.create_dataset("Ldim"+str(Q),
                                data=numpy.array([part.nkk,part.nij,part.kk0,
                                                part.ij0,part.ijN,numv],dtype=numpy.int32))
                LQ = h5file.grp_v2.create_dataset("L"+str(Q),
                                              (part.nkk,part.nij*numv,2),
                                              dtype=numpy.float64)
                # cholvecs[nkk,nij,maxvecs]
                for kk in range(part.nkk):
                    T_ = to_qmcpack_complex(cholvecs[kk,:,0:numv].copy())
                    T_ = numpy.reshape(T_,(-1,2))/math.sqrt(nkpts*1.0)
                    LQ[kk,:,:] = T_
                    T_ = None
            comm.barrier()

        h5file.grp.create_dataset("NCholPerKP", data=num_cholvecs)
        comm.barrier()


def construct_qk_maps(cell, kpts):
    # For a given value of Q (any vector in the 1BZ), there are nkpts kpoint
    # pairs that map to Q such that: (k1 - k2) + G  = Q , where G belongs to the
    # reciprocal lattice cell.  The mapping is taken: K=k1, and Q=(k1 - k2) + G
    nkpts = len(kpts)
    Qpts = kpts - kpts[0]

    QKToK2 = numpy.zeros((nkpts,nkpts),dtype=numpy.int32) - 1
    kminus = numpy.zeros((nkpts,),dtype=numpy.int32) - 1

    shifts = []
    kvecs = cell.reciprocal_vectors()
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                K = numpy.dot(numpy.array([i,j,k]), kvecs)
                shifts.append(K)

    for iq, Q in enumerate(Qpts):
        for ia, ka in enumerate(kpts):
            for ic, kc in enumerate(kpts):
                q = ka - kc
                if numpy.abs(numpy.dot(q-Q,q-Q)) < 1e-10:
                    if QKToK2[iq,ia] >= 0:
                        print("Error: More than 1 solution in "
                              "construct_QKmaps.")
                        sys.exit()
                    QKToK2[iq,ia]=ic
                else:
                    for G in shifts:
                        q2 = q + G
                        if numpy.abs(numpy.dot(q2-Q,q2-Q)) < 1e-10:
                            if QKToK2[iq,ia] >= 0:
                                print("Error: More than 1 solution in"
                                      " construct_QKmaps.")
                                sys.exit()
                            QKToK2[iq,ia]=ic
                            break
            if QKToK2[iq,ia] < 0:
                print("Error: Could not construct QK mapping.")
                sys.exit()

    for iq, Q in enumerate(Qpts):
        for iqm, Qm in enumerate(Qpts):
            q = Q + Qm
            for G in shifts:
                if numpy.abs(numpy.dot(q-G,q-G)) < 1e-10:
                    if kminus[iq] >= 0:
                        print("Error: More than 1 solution to Q + Qm = G. ")
                        sys.exit()
                    kminus[iq] = iqm
                    break
        if kminus[iq] < 0:
            print("Error: Could not solve Q + Qm = G.")
            sys.exit()

    return QKToK2, kminus

def write_kp_h1(h5grp, ki, nmo, h1):
    if h1 is not None:
        assert(h1.shape[0]==nmo)
        assert(h1.shape[1]==nmo)
    H1 = h5grp.create_dataset("H1_kp"+str(ki), (nmo,nmo,2), dtype=numpy.float64)
    if h1 is not None:
        for i in range(h1.shape[0]):
            for j in range(h1.shape[1]):
                H1[i,j,0] = h1[i,j].real
                H1[i,j,1] = h1[i,j].imag
