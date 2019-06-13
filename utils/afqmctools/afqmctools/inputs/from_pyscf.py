import sys
try:
    from pyscf.lib.chkfile import load_mol
    from pyscf.pbc.lib.chkfile import load_cell
except ImportError:
    print("Cannot find pyscf.")
    sys.exit()
from afqmctools.hamiltonian.supercell import write_hamil_supercell
from afqmctools.hamiltonian.kpoint import write_hamil_kpoints
from afqmctools.hamiltonian.mol import write_hamil_mol
from afqmctools.utils.qmcpack_utils import write_xml_input
from afqmctools.utils.pyscf_utils import (
        load_from_pyscf_chk,
        load_from_pyscf_chk_mol
        )
from afqmctools.wavefunction.mol import write_wfn_mol
from afqmctools.wavefunction.pbc import write_wfn_pbc


def write_qmcpack(comm, chkfile, hamil_file, threshold,
                  ortho_ao=False, gdf=False, kpoint=False, verbose=False,
                  cas=None, qmc_input=None, wfn_file=None,
                  write_hamil=True, ndet_max=None):
    """Dispatching routine dependent on options.
    """
    try:
        obj = load_cell(chkfile)
        pbc = True
    except NameError:
        pbc = False

    if pbc:
        if comm.rank == 0 and verbose:
            print(" # Generating Hamiltonian and wavefunction from pyscf cell"
                  " object.")
        scf_data = load_from_pyscf_chk(chkfile, orthoAO=ortho_ao)
        if comm.rank == 0:
            nelec = write_wfn_pbc(scf_data, ortho_ao, wfn_file, verbose=verbose,
                                  ndet_max=ndet_max)
        if write_hamil:
            if kpoint:
                write_hamil_kpoints(comm, scf_data, hamil_file, threshold,
                                    verbose=verbose, cas=cas,
                                    ortho_ao=ortho_ao)
            else:
                write_hamil_supercell(comm, scf_data, hamil_file, threshold,
                                      verbose=verbose, cas=cas,
                                      ortho_ao=ortho_ao)
    else:
        if comm.rank == 0 and verbose:
            print(" # Generating Hamiltonian and wavefunction from pysc mol"
                  " object.")
        if comm.size > 1:
            if comm.rank == 0:
                print(" # Error molecular integral generation must be done "
                      "in serial.")
            sys.exit()
        scf_data = load_from_pyscf_chk_mol(chkfile)
        if write_hamil:
            write_hamil_mol(scf_data, hamil_file, threshold,
                            verbose=verbose, cas=cas, ortho_ao=ortho_ao)
        write_wfn_mol(scf_data, ortho_ao, wfn_file)

    if comm.rank == 0 and qmc_input is not None:
        write_xml_input(qmc_input, hamil_file, wfn_file=wfn_file)
