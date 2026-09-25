##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  error_handler.py                                                  #
#    PWscf parse + recovery                                          #
#                                                                    #
#    generate_pwscf(..., error_handling=True) attaches this chain.   #
#    Spawn can later fan the same handlers out and pick first        #
#    success.                                                        #
#====================================================================#


import os
import re

from ...pwscf_input import PwscfInput
from ..chain_error_handler import ChainErrorHandler, as_simulation, drop_bookkeeping, scan_file


# Tags mainly obtained from AiiDA PwBaseWorkChain
PWSCF_ERROR_MARKERS = (
    ('Maximum CPU time exceeded', 'walltime'),
    ('Program stopped by user request', 'user_stop'),
    ('convergence NOT achieved after', 'electronic_convergence'),
    ('history already reset at previous step: stopping', 'bfgs_history'),
    ('charge is wrong', 'wrong_charge'),
    ('not orthogonal operation', 'symmetry_not_orthogonal'),
    ('Cholesky failed in invchol', 'cholesky'),
    ('problems computing cholesky', 'cholesky'),
    ('too many bands are not converged', 'unconverged_bands'),
    ('S matrix not positive definite', 's_matrix_not_positive_definite'),
    ('zhegvd failed', 'zhegvd_failed'),
    ('[Q, R] = qr(X, 0) failed', 'qr_failed'),
    ('eigenvectors failed to converge', 'unconverged_eigenvectors'),
    ('dexx is negative', 'negative_dexx'),
    ('Error in routine broyden', 'broyden_failure'),
    (
        'Not enough space allocated for radial FFT: try restarting with a larger cell_factor',
        'not_enough_fft_space',
        ),
    ('some nodes have no k-points', 'npools_too_high'),
    ('probably because G_par is NOT a reciprocal lattice vector', 'gpar_error'),
    )

PWSCF_ABORT_MARKERS = (
    'slurmstepd: error',
    'CANCELLED AT',
    'Job step aborted',
    'srun: error',
    'srun: Job step aborted',
    'DUE TO TIME LIMIT',
    'oom-kill',
    'Out Of Memory',
    'Caught signal',
    'Segmentation fault',
    'Fatal Error',
    'mpirun: command not found',
    'mpiexec: command not found',
    )

UNRECOVERABLE_TAGS = frozenset((
    'wrong_charge',
    'symmetry_not_orthogonal',
    'gpar_error',
    'negative_dexx',
    'npools_too_high',
    'command_not_found',
    ))

DIAGO_TAGS = frozenset((
    'cholesky',
    'unconverged_bands',
    's_matrix_not_positive_definite',
    'zhegvd_failed',
    'qr_failed',
    'unconverged_eigenvectors',
    'broyden_failure',
    ))

QE_MIXING_DEFAULT = 0.7
QE_CELL_FACTOR_DEFAULT = 2.0
DIAGO_ORDER = ('david', 'cg')

NBND_RE = re.compile(br'number of Kohn-Sham states\s*=\s*(\d+)')


def parse_pwscf_text(text):
    """Return a unique list of error tags found in a stdout string."""
    tags = []
    seen = set()
    if not text:
        return tags
    #end if
    for marker, tag in PWSCF_ERROR_MARKERS:
        if marker in text and tag not in seen:
            tags.append(tag)
            seen.add(tag)
        #end if
    #end for
    return tags
#end def parse_pwscf_text


class PwscfErrorHandler(ChainErrorHandler):
    """PWscf recovery for ``generate_pwscf(..., error_handling=True)``.

    Chain: one handler patch per attempt.  Pass this object (or an
    int for ``max_runs``) as ``error_handling``.
    """

    unrecoverable = UNRECOVERABLE_TAGS

    def __init__(
        self,
        start=None,
        max_runs=3,
        handlers=None,
        mixing_factor=0.7,
        mixing_min=0.05,
        mixing_default=QE_MIXING_DEFAULT,
        diago_order=DIAGO_ORDER,
        nbnd_factor=0.05,
        nbnd_minimum=4,
        ):
        ChainErrorHandler.__init__(self, start=start, max_runs=max_runs, handlers=handlers)
        self.mixing_factor  = float(mixing_factor)
        self.mixing_min     = float(mixing_min)
        self.mixing_default = float(mixing_default)
        self.diago_order    = tuple(diago_order)
        self.nbnd_factor    = float(nbnd_factor)
        self.nbnd_minimum   = int(nbnd_minimum)
        self.tried_diago    = set()
    #end def __init__

    def default_handlers(self):
        return (
            self.handle_walltime,
            self.handle_electronic_convergence,
            self.handle_diagonalization,
            self.handle_nbnd,
            self.handle_fft_space,
            self.handle_retry,
            )
    #end def default_handlers

    def reset(self):
        ChainErrorHandler.reset(self)
        self.tried_diago = set()
    #end def reset

    def parse(self, sim):
        '''Overlapping with parsing in PWscf.check_sim_status. Can be combined in the future.
        In short term, check_sim_status can reverse read last N lines from stdout to reduce load.
        It currently parses for errors from stdout and stderr, plus number of bands. 
        Number of bands is useful for error handling where nbnd is increased for convergence.'''
        sim = as_simulation(sim)
        locdir  = getattr(sim, 'locdir', None)
        outfile = getattr(sim, 'outfile', None)
        errfile = getattr(sim, 'errfile', None)
        out_path = None
        err_path = None
        if locdir and outfile:
            out_path = os.path.join(locdir, outfile)
        #end if
        if locdir and errfile:
            err_path = os.path.join(locdir, errfile)
        #end if

        markers = [m for m, _ in PWSCF_ERROR_MARKERS]
        found, nbnd_match = scan_file(out_path, markers, extra_re=NBND_RE)
        tags = []
        seen = set()
        marker_to_tag = {m: t for m, t in PWSCF_ERROR_MARKERS}
        for marker in found:
            tag = marker_to_tag.get(marker)
            if tag and tag not in seen:
                tags.append(tag)
                seen.add(tag)
            #end if
        #end for
        abort_found, _ = scan_file(err_path, PWSCF_ABORT_MARKERS)
        for marker in abort_found:
            tag = 'command_not_found' if 'command not found' in marker else 'abort'
            if tag not in seen:
                tags.append(tag)
                seen.add(tag)
            #end if
        #end for

        job_done = False
        if out_path and os.path.isfile(out_path):
            done_found, _ = scan_file(out_path, ('JOB DONE',))
            job_done = bool(done_found)
        #end if

        nbnd = None
        if nbnd_match is not None:
            nbnd = int(nbnd_match.group(1))
        #end if
        return dict(errors=tags, job_done=job_done, nbnd=nbnd)
    #end def parse

    def params_from_input(self, inp):
        params = {}
        electrons = getattr(inp, 'electrons', None)
        if electrons is not None:
            if 'mixing_beta' in electrons:
                params['mixing_beta'] = electrons.mixing_beta
            #end if
            if 'diagonalization' in electrons:
                params['diagonalization'] = electrons.diagonalization
            #end if
            if 'electron_maxstep' in electrons:
                params['electron_maxstep'] = electrons.electron_maxstep
            #end if
        #end if
        control = getattr(inp, 'control', None)
        if control is not None and 'restart_mode' in control:
            params['restart_mode'] = control.restart_mode
        #end if
        system = getattr(inp, 'system', None)
        if system is not None and 'nbnd' in system:
            params['nbnd'] = system.nbnd
        #end if
        cell = getattr(inp, 'cell', None)
        if cell is not None and 'cell_factor' in cell:
            params['cell_factor'] = cell.cell_factor
        #end if
        return params
    #end def params_from_input

    def apply_patch(self, inp, patch):
        '''Apply the patch to the input file.
        Ideally these are tightly organized with error markers'''
        patch = drop_bookkeeping(patch)
        if not patch:
            return
        #end if
        if 'mixing_beta' in patch or 'diagonalization' in patch or 'electron_maxstep' in patch:
            electrons = getattr(inp, 'electrons', None)
            if electrons is None:
                electrons = PwscfInput.section_types['electrons']()
                inp.electrons = electrons
            #end if
            if 'mixing_beta' in patch:
                electrons.mixing_beta = patch['mixing_beta']
            #end if
            if 'diagonalization' in patch:
                electrons.diagonalization = patch['diagonalization']
            #end if
            if 'electron_maxstep' in patch:
                electrons.electron_maxstep = patch['electron_maxstep']
            #end if
        #end if
        if 'restart_mode' in patch:
            inp.control.restart_mode = patch['restart_mode']
        #end if
        if 'nbnd' in patch:
            inp.system.nbnd = patch['nbnd']
        #end if
        if 'cell_factor' in patch:
            cell = getattr(inp, 'cell', None)
            if cell is None:
                cell = PwscfInput.section_types['cell']()
                inp.cell = cell
            #end if
            cell.cell_factor = patch['cell_factor']
        #end if
    #end def apply_patch

    def handle_walltime(self, params, products):
        errors = products.get('errors') or []
        if 'walltime' not in errors and 'user_stop' not in errors:
            return None
        #end if
        return dict(restart_mode='restart')
    #end def handle_walltime

    def handle_electronic_convergence(self, params, products):
        errors = products.get('errors') or []
        if 'electronic_convergence' not in errors:
            return None
        #end if
        beta = float(params.get('mixing_beta', self.mixing_default))
        nxt  = max(self.mixing_min, beta * self.mixing_factor)
        patch = dict(restart_mode='from_scratch')
        if nxt < beta:
            patch['mixing_beta'] = nxt
        #end if
        maxstep = params.get('electron_maxstep')
        if maxstep is not None and int(maxstep) < 80:
            patch['electron_maxstep'] = 80
        #end if
        if 'mixing_beta' not in patch and 'electron_maxstep' not in patch:
            return None
        #end if
        return patch
    #end def handle_electronic_convergence

    def handle_diagonalization(self, params, products):
        errors = products.get('errors') or []
        if not any(tag in DIAGO_TAGS for tag in errors):
            return None
        #end if
        current = params.get('diagonalization', DIAGO_ORDER[0])
        self.tried_diago.add(current)
        for algo in self.diago_order:
            if algo not in self.tried_diago:
                self.tried_diago.add(algo)
                return dict(diagonalization=algo, restart_mode='from_scratch')
            #end if
        #end for
        return None
    #end def handle_diagonalization

    def handle_nbnd(self, params, products):
        errors = products.get('errors') or []
        if 'unconverged_bands' not in errors:
            return None
        #end if
        nbnd = params.get('nbnd')
        if nbnd is None:
            nbnd = products.get('nbnd')
        #end if
        if nbnd is None:
            return None
        #end if
        nbnd = int(nbnd)
        bump = max(int(nbnd * self.nbnd_factor), self.nbnd_minimum)
        return dict(nbnd=nbnd + bump, restart_mode='from_scratch')
    #end def handle_nbnd

    def handle_fft_space(self, params, products):
        errors = products.get('errors') or []
        if 'not_enough_fft_space' not in errors:
            return None
        #end if
        factor = float(params.get('cell_factor', QE_CELL_FACTOR_DEFAULT))
        return dict(cell_factor=2.0 * factor, restart_mode='from_scratch')
    #end def handle_fft_space

    def handle_retry(self, params, products):
        errors = products.get('errors') or []
        if not errors or 'abort' in errors:
            return dict()
        #end if
        return None
    #end def handle_retry
#end class PwscfErrorHandler
