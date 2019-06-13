import h5py
import xml.etree.ElementTree as et
from afqmctools.wavefunction.mol import read_qmcpack_wfn

def write_xml_input(qmc_in, hamil_file, wfn_file, series=0, options=None):
    """Generate template input file from hamiltonian and wavefunction.

    Parameters
    ----------
    qmc_in : string
        Input xml file name to write to.
    hamil_file : string
        HDF5 file containing Hamiltonian.
    wfn_file : string
        HDF5 file containing Wavefunction.
    series : int
        Simulation series number. Optional. Default 0.
    options : dict
        Dictionary of input options to add to basic template. Default None.
        Example:
            options = {
                'execute': {
                    'nWalkers': 10,
                    'blocks': 1000,
                    'steps': 4
                    'Estimator': {
                        back_propagation': {
                            'ortho': 1,
                            'nskip': 8,
                            'block_size': 100,
                            'nsteps': 1000
                        }
                    }
                },
                'Wavefunction': {
                    'nodes': 4
                },
                'Propagator': {
                    'nodes': 4,
                    'vbias_bound': 40
                }
            }
    """
    walker_types = ['NONE', 'CLOSED', 'COLLINEAR', 'NONCOLLINEAR']
    if wfn_file is None:
        wfn_file = hamil_file
    with h5py.File(wfn_file, 'r') as  fh5:
        try:
            dims = fh5['Wavefunction/PHMSD/dims'][:]
            wfn_type = 'PHMSD'
        except KeyError:
            dims = fh5['Wavefunction/NOMSD/dims'][:]
            wfn_type = 'NOMSD'
        nalpha = dims[1]
        nbeta = dims[2]
        nmo = dims[0]
        walker_type = walker_types[dims[3]]

    base = '''<simulation method="afqmc">
    <project id="qmc" series="{:d}"/>
    <AFQMCInfo name="info0">
        <parameter name="NMO">{:d}</parameter>
        <parameter name="NAEA">{:d}</parameter>
        <parameter name="NAEB">{:d}</parameter>
    </AFQMCInfo>
    <Hamiltonian name="ham0" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">{:s}</parameter>
    </Hamiltonian>
    <Wavefunction name="wfn0" type="{:s}" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">{:s}</parameter>
      <parameter name="cutoff">1e-8</parameter>
    </Wavefunction>
    <WalkerSet name="wset0" type="shared">
      <parameter name="walker_type">{:s}</parameter>
    </WalkerSet>
    <Propagator name="prop0" info="info0">
      <parameter name="hybrid">yes</parameter>
    </Propagator>
    <execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">
   </execute>
</simulation>'''.format(series, nmo, nalpha, nbeta,
                        hamil_file, wfn_type, wfn_file,
                        walker_type)
    basic = {
        'execute': {
            'ncores': 1,
            'timestep': 0.005,
            'blocks': 10000,
            'steps': 10,
            'nWalkers': 10
            }
        }
    if options is not None:
        for g, d in options.items():
            if basic.get(g) is not None:
                for k, v in d.items():
                    basic[g][k] = v
            else:
                basic[g] = d
    tree = et.ElementTree(et.fromstring(base))
    root = tree.getroot()
    for g, d in basic.items():
        for k, v in d.items():
            if k == 'Estimator':
                for e, o in v.items():
                    add_estimator(root, e, o)
            else:
                add_param(root, g, k, v)
    indent(root)
    tree.write(qmc_in)

def add_param(root, block, name, val):
    node = root.find(block)
    param = et.Element('parameter', name=name)
    param.text = str(val)
    node.append(param)

def add_estimator(root, name, vals):
    node = root.find('execute')
    param = et.Element('Estimator', name=name)
    node.append(param)
    for k, v in vals.items():
        add_param(root, 'execute/Estimator', k, v)

def indent(elem, ilevel=0):
    # Stackoverflow.
    i = "\n" + ilevel*"  "
    if len(elem) > 0:
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, ilevel+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
