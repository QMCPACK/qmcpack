
import testing
from testing import value_eq,object_eq


associated_files = dict()

def get_files():
    return testing.collect_unit_test_file_paths('xmlreader',associated_files)
#end def get_files



def test_files():
    filenames = [
        'vmc.in.xml',
        ]
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files



def test_import():
    import xmlreader
    from xmlreader import XMLreader,readxml
#end def test_import



def test_read():
    from generic import obj
    from xmlreader import readxml,XMLelement

    ref = obj(
        simulation = obj(
            project = obj(
                id              = 'vmc',
                series          = '0',
                application = obj(
                    class_          = 'serial',
                    name            = 'qmcpack',
                    role            = 'molecu',
                    version         = '1.0',
                    )
                ),
            qmc = obj(
                method          = 'vmc',
                move            = 'pbyp',
                parameter1 = obj(
                    name            = 'walkers',
                    text            = '1',
                    ),
                parameter2 = obj(
                    name            = 'warmupSteps',
                    text            = '20',
                    ),
                parameter3 = obj(
                    name            = 'blocks',
                    text            = '200',
                    ),
                parameter4 = obj(
                    name            = 'steps',
                    text            = '10',
                    ),
                parameter5 = obj(
                    name            = 'subSteps',
                    text            = '2',
                    ),
                parameter6 = obj(
                    name            = 'timestep',
                    text            = '0.4',
                    ),
                ),
            qmcsystem = obj(
                hamiltonian = obj(
                    name            = 'h0',
                    target          = 'e',
                    type            = 'generic',
                    estimator = obj(
                        name            = 'KEcorr',
                        psi             = 'psi0',
                        source          = 'e',
                        type            = 'chiesa',
                        ),
                    pairpot1 = obj(
                        name            = 'ElecElec',
                        source          = 'e',
                        target          = 'e',
                        type            = 'coulomb',
                        ),
                    pairpot2 = obj(
                        name            = 'IonIon',
                        source          = 'ion0',
                        target          = 'ion0',
                        type            = 'coulomb',
                        ),
                    pairpot3 = obj(
                        format          = 'xml',
                        name            = 'PseudoPot',
                        source          = 'ion0',
                        type            = 'pseudo',
                        wavefunction    = 'psi0',
                        pseudo = obj(
                            elementType     = 'C',
                            href            = 'C.BFD.xml',
                            ),
                        ),
                    pairpot4 = obj(
                        ecut            = '60.0',
                        name            = 'MPC',
                        physical        = 'no',
                        source          = 'e',
                        target          = 'e',
                        type            = 'MPC',
                        ),
                    ),
                particleset1 = obj(
                    name            = 'e',
                    random          = 'yes',
                    group1 = obj(
                        mass            = '1.0',
                        name            = 'u',
                        size            = '32',
                        parameter1 = obj(
                            name            = 'charge',
                            text            = '-1',
                            ),
                        parameter2 = obj(
                            name            = 'mass',
                            text            = '1.0',
                            ),
                        ),
                    group2 = obj(
                        mass            = '1.0',
                        name            = 'd',
                        size            = '32',
                        parameter1 = obj(
                            name            = 'charge',
                            text            = '-1',
                            ),
                        parameter2 = obj(
                            name            = 'mass',
                            text            = '1.0',
                            ),
                        ),
                    ),
                particleset2 = obj(
                    name            = 'ion0',
                    group = obj(
                        mass            = '21894.7135906',
                        name            = 'C',
                        size            = '16',
                        attrib = obj(
                            condition       = '0',
                            datatype        = 'posArray',
                            name            = 'position',
                            text            = \
'''6.74632230        6.74632230        0.00000000
                     1.68658057        1.68658057        1.68658058
                     3.37316115        3.37316115        0.00000000
                     5.05974173        5.05974172        1.68658058
                     0.00000000        3.37316115        3.37316115
                     1.68658058        5.05974172        5.05974173
                     3.37316115        6.74632230        3.37316115
                     5.05974173        8.43290287        5.05974173
                     3.37316115        0.00000000        3.37316115
                     5.05974172        1.68658057        5.05974173
                     6.74632230        3.37316115        3.37316115
                     8.43290287        5.05974172        5.05974173
                    10.11948345       10.11948345        6.74632230
                     5.05974173        5.05974172        8.43290287
                     6.74632230        6.74632230        6.74632230
                     8.43290287        8.43290287        8.43290287''',
                            ),
                        parameter1 = obj(
                            name            = 'charge',
                            text            = '4',
                            ),
                        parameter2 = obj(
                            name            = 'valence',
                            text            = '4',
                            ),
                        parameter3 = obj(
                            name            = 'atomicnumber',
                            text            = '6',
                            ),
                        parameter4 = obj(
                            name            = 'mass',
                            text            = '21894.7135906',
                            ),
                        ),
                    ),
                simulationcell = obj(
                    parameter1 = obj(
                        name            = 'lattice',
                        text            = \
'''6.74632230        6.74632230        0.00000000
                  0.00000000        6.74632230        6.74632230
                  6.74632230        0.00000000        6.74632230''',
                        units           = 'bohr',
                        ),
                    parameter2 = obj(
                        name            = 'bconds',
                        text            = 'p p p',
                        ),
                    parameter3 = obj(
                        name            = 'LR_dim_cutoff',
                        text            = '15',
                        ),
                    ),
                wavefunction = obj(
                    name            = 'psi0',
                    target          = 'e',
                    determinantset = obj(
                        slaterdeterminant = obj(
                            determinant1 = obj(
                                group           = 'u',
                                id              = 'updet',
                                size            = '32',
                                sposet          = 'spo_ud',
                                ),
                            determinant2 = obj(
                                group           = 'd',
                                id              = 'downdet',
                                size            = '32',
                                sposet          = 'spo_ud',
                                ),
                            ),
                        ),
                    sposet_builder = obj(
                        href            = 'MISSING.h5',
                        meshfactor      = '1.0',
                        precision       = 'float',
                        source          = 'ion0',
                        tilematrix      = '2 0 0 0 2 0 0 0 2',
                        truncate        = 'no',
                        twistnum        = '0',
                        type            = 'bspline',
                        version         = '0.10',
                        sposet = obj(
                            name            = 'spo_ud',
                            size            = '32',
                            spindataset     = '0',
                            type            = 'bspline',
                            ),
                        ),
                    ),
                ),
            ),
        )

    files = get_files()
    x = readxml(files['vmc.in.xml'])
    assert(isinstance(x,XMLelement))
    x.remove_hidden()
    o = x.to_obj()

    assert(object_eq(o,ref))
#end def test_read
