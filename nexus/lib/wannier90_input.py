from simulation import SimulationInput

class Wannier90Input(SimulationInput):
    """
    Class for generating and processing Wannier90 input files.
    Based on Wannier90 3.1.0 manual.
    """
    
    available_versions = [(3, 1, 0)
                          # Example of how to add a new version
                          #, (4, 0, 0)
                          ]
    # Required keywords
    required_v310 = ['num_wann', 'num_bands']

    # All known keywords with their types
    integers_v310 = [
        'num_wann', 'num_bands', 'num_iter', 'num_print_cycles',
        'num_dump_cycles', 'num_cg_steps', 'conv_window', 'conv_noise_num',
        'num_excluded_bands', 'mp_grid', 'kmesh_tol', 'num_spins',
        'spinor_mode', 'nfermi', 'timing_level', 'num_guide_cycles',
        'num_no_guide_iter', 'search_shells', 'conv_window', 'bands_num_points',
        'fermi_surface_num_points', 'ws_search_size', 'tran_num_bb', 'tran_num_ll',
        'tran_num_rr', 'tran_num_cc', 'tran_num_lc', 'tran_num_cr',
        'tran_num_cell_ll', 'tran_num_cell_rr', 'tran_num_bandc'
    ]

    floats_v310 = [
        'conv_tol', 'conv_noise_amp', 'conv_noise_amp_global',
        'conv_window_threshold', 'dis_win_min', 'dis_win_max',
        'dis_froz_min', 'dis_froz_max', 'fermi_energy',
        'kmesh_spacing', 'zone_fraction', 'translation_centre_frac',
        'wannier_plot_radius', 'wannier_plot_scale', 'wannier_plot_spinor_mode',
        'hr_cutoff', 'dist_cutoff', 'precond_threshold', 'omega_invariant',
        'scale_factor', 'energy_offset', 'trial_step', 'fixed_step',
        'dis_mix_ratio', 'dis_conv_tol', 'symmetrize_eps', 'slwf_lambda',
        'tran_win_min', 'tran_win_max', 'tran_energy_step',
        'tran_group_threshold', 'ws_distance_tol'
    ]

    bools_v310 = [
        'guiding_centres', 'use_bloch_phases', 'write_xyz', 'gamma_only',
        'write_vdw', 'write_hr_diag', 'write_rmn', 'write_bvec',
        'write_bstate', 'write_r_pos', 'write_proj_site', 'shell_list',
        'write_force', 'write_tb', 'write_u_matrices', 'search_shells',
        'skip_B1_tests', 'postproc_setup', 'auto_projections', 'precond',
        'spinors', 'site_symmetry', 'slwf_constrain', 'translate_home_cell',
        'use_ws_distance', 'transport', 'tran_write_ht', 'tran_read_ht',
        'tran_use_same_lead', 'wannier_plot', 'bands_plot', 'fermi_surface_plot',
        'write_hr', 'write_xyz', 'write_rmn', 'write_tb'
    ]

    strings_v310 = [
        'wannier_mode', 'restart', 'restart_mode', 'iprint',
        'length_unit', 'wvfn_formatted', 'spin', 'spinors',
        'devel_flag', 'optimisation', 'bloch_outer_window',
        'translation_mode', 'translation_vector_frac',
        'wannier_plot_format', 'wannier_plot_mode',
        'wannier_plot_supercell', 'wannier_plot_list',
        'hr_plot_format', 'dist_plot_format', 'transport_mode',
        'one_dim_axis', 'dist_cutoff_mode', 'bands_plot_format',
        'bands_plot_mode', 'fermi_surface_plot_format',
        'wannier_plot_spinor_mode'
    ]

    # Special multi-value parameters
    arrays_v310 = [
        'unit_cell_cart',      # Unit cell vectors in Cartesian coordinates
        'atoms_cart',          # Atomic positions in Cartesian coordinates
        'atoms_frac',          # Atomic positions in fractional coordinates
        'mp_grid',             # Dimensions of Monkhorst-Pack grid
        'kpoints',             # List of k-points
        'projections',         # Projection functions for initial guess
        'exclude_bands',       # Bands to exclude
        'select_projections',  # Projections to select
        'kpoint_path',         # K-point path for band structure
        'bands_plot_project',  # WFs to project bands onto
        'slwf_centres',        # Centers for selective localization
        'dis_spheres',         # Spheres for k-dependent disentanglement
        'nnkpts'               # Explicit nearest-neighbor k-points
    ]

    
    blocks_list_v310 = [
        'unit_cell_cart',
        'atoms_cart',
        'atoms_frac',
        'mp_grid',
        'kpoints',
        'projections',
    ]

    versioned_inputs = {(3, 1, 0): {
        'required': required_v310,
        'integers': integers_v310,
        'floats': floats_v310,
        'bools': bools_v310,
        'strings': strings_v310,
        'blocks_list': blocks_list_v310
    }, 
    # Example of how to add a new version
    # (4, 0, 0): {
    #     'required': required_v400.update(required_v310),
    #     'integers': integers_v400.update(integers_v310),
    #     'floats': floats_v400.update(floats_v310),
    #     'bools': bools_v400.update(bools_v310),
    #     'strings': strings_v400.update(strings_v310),
    #     'blocks': blocks_list_v400.update(blocks_list_v310)
    # }
    }


    def __init__(self, version:str='3.1.0'):
        """Initialize an empty Wannier90 input."""
        self.params = {}
        self.blocks = {}
        self.version = tuple([int(i) for i in version.split('.')])
        self.check_version()

    def check_version(self):
        """Check the version of the Wannier90 input."""
        closest_version = None
        for version in self.available_versions:
            if self.version >= version:
                closest_version = version

        if closest_version is None:
            raise ValueError(f"Unsupported version: {self.version}")
        else:
            self.version = closest_version
            for key, value in self.versioned_inputs[closest_version].items():
                self[key] = value

    def read(self, filename):
        """Read a Wannier90 input file."""
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        current_block = None
        for line in lines:
            line = line.strip()
            lower_line = line.lower()
            if not line or line.startswith('!') or line.startswith('#'): 
                continue
                
            if lower_line.startswith('begin'):
                current_block = line.split()[1].lower()
                self.blocks[current_block] = []
                continue
            elif lower_line.startswith('end'):
                current_block = None
                continue
                
            if current_block:
                self.blocks[current_block].append(line)
            else:
                # Handle parameter lines
                # Parameters can be in an interchangeable format:
                #   key = value
                #   key : value
                if '=' in line:
                    key, value = [x.strip() for x in line.split('=', 1)]
                    self.params[key.lower()] = value
                if ':' in line:
                    key, value = [x.strip() for x in line.split(':', 1)]
                    self.params[key.lower()] = value

        for block in self.blocks:
            self.blocks[block] = '\n'.join(self.blocks[block])

    def write_text(self):
        """Generate Wannier90 input file contents as a string."""
        text = []
        
        # Write parameters first
        for key, value in self.params.items():
            text.append(f"{key} = {value}")
        
        # Add a blank line before blocks if there were parameters
        if self.params:
            text.append('')
            
        # Write blocks with proper begin/end markers and spacing
        for block_name, content in self.blocks.items():
            text.append(f"begin {block_name}")
            text.append(content)
            text.append(f"end {block_name}")
            text.append('')  # Add blank line after each block
            
        # Remove trailing newline
        if text and not text[-1]:
            text.pop()
        return '\n'.join(text)

    def write(self, filepath):
        """Write Wannier90 input to a file."""
        with open(filepath, 'w') as f:
            f.write(self.write_text())

    def validate(self):
        """Validate the input parameters."""
        # Check required parameters
        for param in self.required:
            if param not in self.params:
                raise ValueError(f"Required parameter '{param}' is missing")
        
        # Validate parameter types
        for param, value in self.params.items():
            if param in self.integers:
                try:
                    int(value)
                except ValueError:
                    raise ValueError(f"Parameter '{param}' must be an integer")
            elif param in self.floats:
                try:
                    float(value)
                except ValueError:
                    raise ValueError(f"Parameter '{param}' must be a float")
            elif param in self.bools:
                if value.lower() not in ['true', 'false', 't', 'f']:
                    raise ValueError(f"Parameter '{param}' must be a boolean")

    def set_param(self, key, value):
        """Set a parameter value."""
        self.params[key.lower()] = value

    def get_param(self, key):
        """Get a parameter value."""
        return self.params.get(key.lower())

    def add_block(self, name, content):
        """Add a block of data."""
        self.blocks[name.lower()] = content

    def get_block(self, name):
        """Get a block of data."""
        return self.blocks.get(name.lower())

def format_unit_cell_cart(axes, units='bohr'):
    """Format unit cell vectors block.
    
    Parameters
    ----------
    axes : array_like
        3x3 array of unit cell vectors
    units : str, optional
        Units for the cell vectors ('Bohr' or 'A')
    """
    allowed_units = ['bohr', 'a']
    if units.lower() not in allowed_units:
        raise ValueError(f"Invalid units: {units}. Must be one of: {', '.join(allowed_units)}")
    else:
        units = units.capitalize()
    #end if 
    block = [units]
    for v in axes:
        block.append('{: 8.6f}  {: 8.6f}  {: 8.6f}'.format(*v))
    return '\n'.join(block)


def format_atoms_cart(atoms, positions, units='bohr'):
    """Format atomic positions in Cartesian coordinates.
    
    Parameters
    ----------
    atoms : list
        List of atomic symbols
    positions : array_like
        Array of atomic positions
    units : str, optional
        Units for positions ('Bohr' or 'A')
    """
    allowed_units = ['bohr', 'a']
    if units.lower() not in allowed_units:
        raise ValueError(f"Invalid units: {units}. Must be one of: {', '.join(allowed_units)}")
    else:
        units = units.capitalize()
    #end if 
    block = [units]
    for atom, pos in zip(atoms, positions):
        block.append('{:2s}  {: 8.6f}  {: 8.6f}  {: 8.6f}'.format(atom, *pos))
    return '\n'.join(block)


def format_atoms_frac(atoms, positions):
    """Format atomic positions in fractional coordinates.
    
    Parameters
    ----------
    atoms : list
        List of atomic symbols
    positions : array_like
        Array of fractional atomic positions
    """
    block = []
    for atom, pos in zip(atoms, positions):
        block.append('{:2s}  {: 8.6f}  {: 8.6f}  {: 8.6f}'.format(atom, *pos))
    return '\n'.join(block)


def format_kpoints(kpoints):
    """Format explicit k-points list.
    
    Parameters
    ----------
    kpoints : array_like
        Array of k-points
    """
    block = []
    for k in kpoints:
        block.append('{: 8.6f}  {: 8.6f}  {: 8.6f}'.format(*k))
    return '\n'.join(block)


def format_kpoint_path(path):
    """Format k-point path for band structure.
    
    Parameters
    ----------
    path : list or dict
        K-point path specification, either:
        - list of tuples: [(label1, k1), (label2, k2), ...]
        - dict from seekpath: {'path': [...], 'point_coords': {...}}
    """
    block = []
    if isinstance(path, dict):
        # Handle seekpath format
        kpath = path['path']
        coords = path['point_coords']
        for i in range(len(kpath)):
            start, end = kpath[i][0], kpath[i][1]
            line = '{:s}  {: 8.6f}  {: 8.6f}  {: 8.6f}    {:s}  {: 8.6f}  {: 8.6f}  {: 8.6f}'.format(
                start, *coords[start], end, *coords[end])
            block.append(line)
    else:
        # Handle direct format
        for i in range(len(path)):
            start, end = path[i][0], path[i][1]
            label1, k1 = start
            label2, k2 = end
            line = '{:s}  {: 8.6f}  {: 8.6f}  {: 8.6f}    {:s}  {: 8.6f}  {: 8.6f}  {: 8.6f}'.format(
                label1, *k1, label2, *k2)
            block.append(line)
    return '\n'.join(block)


def format_projections(projections):
    """Format orbital projections.
    
    Parameters
    ----------
    projections : list
        List of projection specifications, e.g. ['Ga:s,p', 'As:p']
    """
    if isinstance(projections, str):
        projections = [projections]
    return '\n'.join(projections)


def generate_wannier90_input(
        prefix          = 'wannier',
        system         = None,
        num_wann       = None,
        num_bands      = None,
        unit_cell_cart = None,
        atoms_cart     = None,
        atoms_frac     = None,
        mp_grid        = None,
        kpoints        = None,
        kpoint_path    = None,
        projections    = None,
        units          = None,
        **kwargs
        ):
    """
    Generate a Wannier90 input file (.win)
    
    Parameters
    ----------
    prefix : str, optional
        Prefix for the Wannier90 calculation
    system : PhysicalSystem object, optional
        System containing structure information
    num_wann : int
        Number of Wannier functions
    num_bands : int
        Number of bands to include in the wannierization
    unit_cell_cart : array_like, optional
        Unit cell vectors in Cartesian coordinates
    atoms_cart : tuple or list, optional
        Atomic positions in Cartesian coordinates, format: [(atom1, pos1), ...]
    atoms_frac : tuple or list, optional
        Atomic positions in fractional coordinates, format: [(atom1, pos1), ...]
    mp_grid : tuple of int, optional
        Monkhorst-Pack grid dimensions
    kpoints : array_like, optional
        Explicit k-points list
    kpoint_path : list or dict, optional
        K-point path for band structure, can be:
        - list of tuples: [(label1, k1), (label2, k2), ...]
        - dict from seekpath: {'path': [...], 'point_coords': {...}}
    projections : list, optional
        Orbital projections, e.g. ['Ga:s,p', 'As:p']
    units : str, optional
        Units for positions and cell vectors ('bohr' or 'a')
    **kwargs : dict
        Additional Wannier90 parameters
    """
    
    win = Wannier90Input()
    
    # Set required parameters
    if num_wann is None or num_bands is None:
        win.error('num_wann and num_bands are required parameters')
    win.set_param('num_wann', num_wann)
    win.set_param('num_bands', num_bands)
    
    # Handle structure information
    if system is not None:
        # Get structure from PhysicalSystem
        if units is not None:
            system_units = system.change_units(units).copy()
        else:
            system_units = system
            units = system.generation_info.structure.units

        axes = system_units.structure.axes
        pos = system_units.structure.pos
        elem = system_units.structure.elem
        win.add_block('unit_cell_cart', format_unit_cell_cart(axes, units))
        win.add_block('atoms_cart', format_atoms_cart(elem, pos, units))
    else:
        available_units = ['Bohr', 'A']
        available_units = [u.lower() for u in available_units]
        assert units.lower() in available_units, f'units must be one of: {", ".join(available_units)}'
        # Use provided structure information
        if unit_cell_cart is not None:
            win.add_block('unit_cell_cart', format_unit_cell_cart(unit_cell_cart, units))
        if atoms_cart is not None:
            atoms, positions = zip(*atoms_cart)
            win.add_block('atoms_cart', format_atoms_cart(atoms, positions, units))
        elif atoms_frac is not None:
            atoms, positions = zip(*atoms_frac)
            
            win.add_block('atoms_frac', format_atoms_frac(atoms, positions))

    # Handle k-points
    if mp_grid is not None:
        win.set_param('mp_grid', '{} {} {}'.format(*mp_grid))
        if kpoints is not None:
            win.add_block('kpoints', format_kpoints(kpoints))
    
    # Handle k-point path
    if kpoint_path is not None:
        win.add_block('kpoint_path', format_kpoint_path(kpoint_path))

    # Handle projections
    if projections is not None:
        win.add_block('projections', format_projections(projections))

    # Handle all other parameters
    for key, value in kwargs.items():
        if isinstance(value, (list, tuple)):
            # Handle array-type parameters
            if all(isinstance(x, str) for x in value):
                win.set_param(key, ' '.join(value))
            else:
                win.set_param(key, ' '.join(map(str, value)))
        else:
            win.set_param(key, value)

    return win


if __name__ == "__main__":
    # Using direct structure input
    print('Using direct structure input')

    win = generate_wannier90_input(
        prefix = 'silicon',
        unit_cell_cart = [[5.43, 0.0, 0.0],
                        [0.0, 5.43, 0.0],
                        [0.0, 0.0, 5.43]],
        atoms_cart = [('Si', (0.0, 0.0, 0.0)),
                    ('Si', (2.715, 2.715, 2.715))],
        num_wann = 8,
        num_bands = 12,
        mp_grid = (4,4,4),
        projections = ['Si:sp3'],
        # kpoints = [[0.0, 0.0, 0.0],
        #           [0.5, 0.5, 0.5]],
        # TODO: Check later on with seekpath
        # Seekpath format
        kpoint_path = {'path': [['G', 'X'], 
                                ['X', 'K']],
                       'point_coords': {'G': [0.0, 0.0, 0.0],
                                        'X': [0.5, 0.5, 0.5],
                                        'K': [0.375, 0.375, 0.75]}},
        band_num_points = 100,
        units = 'A'
    )

    # Get input file contents as string
    win_text = win.write_text()
    print(win_text)    

    print('Using PhysicalSystem input')
    from physical_system import generate_physical_system    
    from structure import Structure
    system = generate_physical_system(
        structure = Structure(
            axes = [[5.43, 0.0, 0.0],
                    [0.0, 5.43, 0.0],
                    [0.0, 0.0, 5.43]],
            pos = [(0.0, 0.0, 0.0), (2.715, 2.715, 2.715)],
            elem = ['Si', 'Si'], 
            units = 'A'),
        net_charge = 0,
        net_spin = 0,
        
    )

    win = generate_wannier90_input(
        system = system,
        num_wann = 8,
        num_bands = 12,
        mp_grid = (4,4,4),
        projections = ['Si:sp3'],
        # Direct kpoint path
        kpoint_path = [[('G', [0.0, 0.0, 0.0]), ('X', [0.5, 0.5, 0.5])],
                       [('X', [0.5, 0.5, 0.5]), ('K', [0.375, 0.375, 0.75])]],
        band_num_points = 100)  

    win_text = win.write_text()
    print(win_text)  