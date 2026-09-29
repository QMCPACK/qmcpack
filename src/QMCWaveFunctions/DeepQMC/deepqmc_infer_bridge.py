######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2026 QMCPACK developers.
##
## File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
##
## File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#######################################################################################

"""Python side of QMCPACK's experimental DeepQMC inference bridge.

The C++ WaveFunctionComponent calls DeepQMCInferBridge.compute_log_gl with a
batch of walker electron coordinates.  The returned quantities are in the
QMCPACK WaveFunctionComponent convention: log(psi), grad log(psi), and per
particle laplacian log(psi).

Import order is intentional.  The standalone miniapp found that importing
DeepQMC before JAX avoids initialization problems.
"""

import deepqmc  # noqa: F401  keep before jax
import jax
import jax.numpy as jnp
import jax_dataclasses as jdc
import os
import pickle
import sys
import types
from pathlib import Path
from hydra import compose, initialize_config_dir
from hydra.utils import instantiate
# PySCF imports h5py at module import time.  In embedded QMCPACK, binary h5py
# wheels can collide with the HDF5 library already linked into QMCPACK.  The He
# inference prototype does not need PySCF checkpoint I/O, so avoid importing the
# real h5py module just to satisfy PySCF's import-time dependency.
_h5py_stub = types.ModuleType('h5py')
_h5py_stub.version = types.SimpleNamespace(version='3.0.0')
_h5py_stub.File = type('File', (), {})
sys.modules.setdefault('h5py', _h5py_stub)

from deepqmc.app import instantiate_ansatz
from deepqmc.molecule import Molecule
from deepqmc.hamil import MolecularHamiltonian
from deepqmc.parallel import replicate_on_devices, scatter_electrons_to_devices
from deepqmc.types import PhysicalConfiguration


def _find_hydra_config(model_path):
    """Find a resolved DeepQMC Hydra config associated with a checkpoint.

    DeepQMC checkpoints store parameters/optimizer/sampler state, but the Haiku
    module tree is defined by the training Hydra config. Prefer an explicit
    environment override, then common DeepQMC workdir layouts:

      workdir/chkpt-*.pt
      workdir/training/chkpt-*.pt
      workdir/evaluation/chkpt-*.pt
    """
    env_path = os.environ.get('DEEPQMC_HYDRA_CONFIG')
    if env_path:
        path = Path(env_path).expanduser()
        if not path.exists():
            raise FileNotFoundError(f'DEEPQMC_HYDRA_CONFIG does not exist: {path}')
        return path

    checkpoint = Path(model_path).expanduser()
    search_dirs = [checkpoint.parent, *checkpoint.parents]
    for directory in search_dirs:
        candidate = directory / '.hydra' / 'config.yaml'
        if candidate.exists():
            return candidate
    return None


def _instantiate_problem_from_hydra_config(config_path):
    """Instantiate MolecularHamiltonian and ansatz from a saved Hydra config."""
    from omegaconf import OmegaConf

    cfg = OmegaConf.load(config_path)
    hamiltonian = instantiate(cfg.hamil, _recursive_=True, _convert_='all')
    ansatz_cfg = instantiate(cfg.ansatz, _recursive_=True, _convert_='all')
    return hamiltonian, instantiate_ansatz(hamiltonian, ansatz_cfg)


def _instantiate_default_he_problem():
    """Fallback prototype default used by older local He checkpoints."""
    mol = Molecule(coords=[[0.0, 0.0, 0.0]], charges=[2], charge=0, spin=0, unit='bohr')
    hamiltonian = MolecularHamiltonian(mol=mol)

    deepqmc_dir = os.path.dirname(deepqmc.__file__)
    config_dir = os.path.join(deepqmc_dir, 'conf/ansatz')
    with initialize_config_dir(version_base=None, config_dir=config_dir):
        cfg = compose(config_name='default')

    ansatz_cfg = instantiate(cfg, _recursive_=True, _convert_='all')
    return hamiltonian, instantiate_ansatz(hamiltonian, ansatz_cfg)


def _load_checkpoint_without_h5py(path):
    """Load a DeepQMC CheckpointStore checkpoint without importing deepqmc.log.

    deepqmc.log imports h5py for training logs.  Embedded QMCPACK already links
    HDF5, which can conflict with binary h5py wheels, but inference only needs
    the pickled TrainState params.  This mirrors the relevant part of
    CheckpointStore.load / deserialize_train_state while avoiding h5py import.
    """
    with Path(path).open('rb') as f:
        step, train_state = pickle.load(f)

    if train_state.sampler['elec'].get('r', None) is not None:
        if train_state.sampler['elec']['r'].ndim == 6:
            return step, train_state
    if train_state.sampler['elec'].get('tau', None) is not None:
        if train_state.sampler['elec']['tau'].ndim == 3:
            train_state.sampler['elec']['tau'] = train_state.sampler['elec']['tau'].mean(axis=-1)
    train_state.sampler['elec']['tau'] = jnp.repeat(
        train_state.sampler['elec']['tau'][..., None], jax.device_count(), axis=-1
    )
    params, opt = replicate_on_devices((train_state.params, train_state.opt))
    sampler = train_state.sampler
    sampler['elec'] = scatter_electrons_to_devices(sampler['elec'])
    sampler['elec']['tau'] = jnp.squeeze(sampler['elec']['tau'], axis=-1)
    sampler['nuc'], sampler['update_nuc_counter'] = replicate_on_devices(
        (sampler['nuc'], sampler['update_nuc_counter'])
    )
    return step, type(train_state)(sampler, params, opt)


class DeepQMCInferBridge:
    def __init__(self, model_path):
        hydra_config = _find_hydra_config(model_path)
        if hydra_config is not None:
            self.H, self.ansatz = _instantiate_problem_from_hydra_config(hydra_config)
        else:
            self.H, self.ansatz = _instantiate_default_he_problem()

        step, train_state = _load_checkpoint_without_h5py(Path(model_path))
        params = train_state.params

        def drop_first_two_dims(x):
            if hasattr(x, 'ndim') and x.ndim >= 2:
                return x[0, 0]
            return x

        self.params = jax.tree_util.tree_map(drop_first_two_dims, params)

        ansatz = self.ansatz
        params = self.params

        def single_log_grad_lap(r, R, mol_idx):
            phys_conf = PhysicalConfiguration(R=R, r=r, mol_idx=jnp.array([mol_idx]))

            def log_psi_fn(r_flat):
                return ansatz.apply(params, jdc.replace(phys_conf, r=r_flat.reshape(-1, 3))).log

            r_flat = r.flatten()
            log_val = log_psi_fn(r_flat)
            grad_flat = jax.grad(log_psi_fn)(r_flat)
            hess = jax.hessian(log_psi_fn)(r_flat)
            diag = jnp.diag(hess).reshape(-1, 3)
            lap_per_electron = jnp.sum(diag, axis=1)
            return log_val, grad_flat.reshape(-1, 3), lap_per_electron

        @jax.jit
        def _compute_log_gl_batch(r_batch, R, mol_idx):
            vals, grads, laps = jax.vmap(lambda r: single_log_grad_lap(r, R, mol_idx))(r_batch)
            return vals, grads, laps

        self._compute_log_gl_batch = _compute_log_gl_batch

    def compute_log_gl(self, nuclear_coords, electron_coords, mol_idx, batch_size, n_elec):
        R = jnp.array(nuclear_coords).reshape(-1, 3)
        r_batch = jnp.array(electron_coords).reshape(batch_size, n_elec, 3)
        vals, grads, laps = self._compute_log_gl_batch(r_batch, R, mol_idx)
        jax.block_until_ready((vals, grads, laps))
        return vals.tolist(), grads.reshape(-1).tolist(), laps.reshape(-1).tolist()
