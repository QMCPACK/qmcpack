
import versions
import testing
from testing import divert_nexus_log,restore_nexus_log
from testing import value_eq,object_eq,check_object_eq

associated_files = dict()

def get_files():
    return testing.collect_unit_test_file_paths('qmcpack_input',associated_files)
#end def get_files


def format_value(v):
    import numpy as np
    s = ''
    if isinstance(v,np.ndarray):
        pad = 12*' '
        s = 'np.array([\n'
        if len(v.shape)==1:
            s += pad
            for vv in v:
                s += format_value(vv)+','
            #end for
            s = s[:-1]
        else:
            for vv in v:
                s += pad + format_value(list(vv))+',\n'
            #end for
            s = s[:-2]
        #end if
        s += '])'
    elif isinstance(v,(str,np.string_)):
        s = "'"+v+"'"
    else:
        s = str(v)
    #end if
    return s
#end def format_value


def make_serial_reference(qi):
    s = qi.serial()
    ref = '    ref = {\n'
    for k in sorted(s.keys()):
        v = s[k]
        ref +="        '{}' : {},\n".format(k,format_value(v))
    #end for
    ref += '        }\n'
    return ref
#end def make_serial_reference


serial_references = dict()


def generate_serial_references():
    import numpy as np
    from qmcpack_input import onerdm

    # reference for read/write
    ref = {
        '_metadata/lattice/units' : 'bohr',
        '_metadata/position/condition' : '0',
        '_metadata/position/datatype' : 'posArray',
        'simulation/calculations/0/blocks' : 70,
        'simulation/calculations/0/checkpoint' : -1,
        'simulation/calculations/0/method' : 'vmc',
        'simulation/calculations/0/move' : 'pbyp',
        'simulation/calculations/0/samplesperthread' : 2,
        'simulation/calculations/0/steps' : 5,
        'simulation/calculations/0/substeps' : 2,
        'simulation/calculations/0/timestep' : 0.3,
        'simulation/calculations/0/walkers' : 1,
        'simulation/calculations/0/warmupsteps' : 20,
        'simulation/calculations/1/blocks' : 80,
        'simulation/calculations/1/checkpoint' : -1,
        'simulation/calculations/1/method' : 'dmc',
        'simulation/calculations/1/move' : 'pbyp',
        'simulation/calculations/1/nonlocalmoves' : 'yes',
        'simulation/calculations/1/steps' : 5,
        'simulation/calculations/1/timestep' : 0.02,
        'simulation/calculations/1/warmupsteps' : 2,
        'simulation/calculations/2/blocks' : 600,
        'simulation/calculations/2/checkpoint' : -1,
        'simulation/calculations/2/method' : 'dmc',
        'simulation/calculations/2/move' : 'pbyp',
        'simulation/calculations/2/nonlocalmoves' : 'yes',
        'simulation/calculations/2/steps' : 5,
        'simulation/calculations/2/timestep' : 0.005,
        'simulation/calculations/2/warmupsteps' : 10,
        'simulation/project/application/class_' : 'serial',
        'simulation/project/application/name' : 'qmcapp',
        'simulation/project/application/role' : 'molecu',
        'simulation/project/application/version' : 1.0,
        'simulation/project/id' : 'qmc',
        'simulation/project/series' : 0,
        'simulation/qmcsystem/hamiltonians/h0/estimators/KEcorr/name' : 'KEcorr',
        'simulation/qmcsystem/hamiltonians/h0/estimators/KEcorr/psi' : 'psi0',
        'simulation/qmcsystem/hamiltonians/h0/estimators/KEcorr/source' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/estimators/KEcorr/type' : 'chiesa',
        'simulation/qmcsystem/hamiltonians/h0/estimators/SpinDensity/grid' : np.array([
            72,44,44]),
        'simulation/qmcsystem/hamiltonians/h0/estimators/SpinDensity/name' : 'SpinDensity',
        'simulation/qmcsystem/hamiltonians/h0/estimators/SpinDensity/type' : 'spindensity',
        'simulation/qmcsystem/hamiltonians/h0/name' : 'h0',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/ElecElec/name' : 'ElecElec',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/ElecElec/source' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/ElecElec/target' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/ElecElec/type' : 'coulomb',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/IonIon/name' : 'IonIon',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/IonIon/source' : 'ion0',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/IonIon/target' : 'ion0',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/IonIon/type' : 'coulomb',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/ecut' : 60.0,
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/name' : 'MPC',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/physical' : False,
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/source' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/target' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/MPC/type' : 'MPC',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/format' : 'xml',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/name' : 'PseudoPot',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/pseudos/O/elementtype' : 'O',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/pseudos/O/href' : 'O.opt.xml',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/pseudos/V/elementtype' : 'V',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/pseudos/V/href' : 'V.opt.xml',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/source' : 'ion0',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/type' : 'pseudo',
        'simulation/qmcsystem/hamiltonians/h0/pairpots/PseudoPot/wavefunction' : 'psi0',
        'simulation/qmcsystem/hamiltonians/h0/target' : 'e',
        'simulation/qmcsystem/hamiltonians/h0/type' : 'generic',
        'simulation/qmcsystem/particlesets/e/groups/d/charge' : -1,
        'simulation/qmcsystem/particlesets/e/groups/d/mass' : 1.0,
        'simulation/qmcsystem/particlesets/e/groups/d/name' : 'd',
        'simulation/qmcsystem/particlesets/e/groups/d/size' : 200,
        'simulation/qmcsystem/particlesets/e/groups/u/charge' : -1,
        'simulation/qmcsystem/particlesets/e/groups/u/mass' : 1.0,
        'simulation/qmcsystem/particlesets/e/groups/u/name' : 'u',
        'simulation/qmcsystem/particlesets/e/groups/u/size' : 200,
        'simulation/qmcsystem/particlesets/e/name' : 'e',
        'simulation/qmcsystem/particlesets/e/random' : True,
        'simulation/qmcsystem/particlesets/ion0/groups/O/atomicnumber' : 8,
        'simulation/qmcsystem/particlesets/ion0/groups/O/charge' : 6,
        'simulation/qmcsystem/particlesets/ion0/groups/O/mass' : 29164.3928678,
        'simulation/qmcsystem/particlesets/ion0/groups/O/name' : 'O',
        'simulation/qmcsystem/particlesets/ion0/groups/O/position' : np.array([
            [0.00978311, 1.81708472, 1.78656736],
            [10.85992161, 15.33331378, -1.78656738],
            [2.75326234, 11.04571415, -2.495713],
            [8.1164424, 6.10468435, 2.495713],
            [2.71381355, 6.02493499, 2.55909075],
            [8.15589117, 11.12546351, -2.55909075],
            [5.45729278, 15.41306313, -1.72318961],
            [5.41241194, 1.73733537, 1.72318961],
            [10.87948783, 1.81708472, 1.78656736],
            [21.72962633, 15.33331378, -1.78656738],
            [13.62296706, 11.04571415, -2.495713],
            [18.98614712, 6.10468435, 2.495713],
            [13.58351827, 6.02493499, 2.55909075],
            [19.02559589, 11.12546351, -2.55909075],
            [16.3269975, 15.41306313, -1.72318961],
            [16.28211666, 1.73733537, 1.72318961],
            [0.00978311, 10.39228397, 1.78656736],
            [10.85992161, 6.75811453, -1.78656738],
            [-2.7336961, 11.04571415, 6.06884775],
            [13.60340084, 6.10468435, -6.06884775],
            [8.20077199, 6.02493499, -6.00547],
            [2.66893273, 11.12546351, 6.00547],
            [5.45729278, 6.83786388, -1.72318961],
            [5.41241194, 10.31253462, 1.72318961],
            [10.87948783, 10.39228397, 1.78656736],
            [21.72962633, 6.75811453, -1.78656738],
            [8.13600862, 11.04571415, 6.06884775],
            [24.47310556, 6.10468435, -6.06884775],
            [19.07047671, 6.02493499, -6.00547],
            [13.53863745, 11.12546351, 6.00547],
            [16.3269975, 6.83786388, -1.72318961],
            [16.28211666, 10.31253462, 1.72318961]]),
        'simulation/qmcsystem/particlesets/ion0/groups/O/size' : 32,
        'simulation/qmcsystem/particlesets/ion0/groups/O/valence' : 6,
        'simulation/qmcsystem/particlesets/ion0/groups/V/atomicnumber' : 23,
        'simulation/qmcsystem/particlesets/ion0/groups/V/charge' : 13,
        'simulation/qmcsystem/particlesets/ion0/groups/V/mass' : 92861.5851912,
        'simulation/qmcsystem/particlesets/ion0/groups/V/name' : 'V',
        'simulation/qmcsystem/particlesets/ion0/groups/V/position' : np.array([
            [2.45778327, 8.39460555, 0.22661828],
            [8.41192147, 8.75579295, -0.22661828],
            [5.2012625, 13.04339257, -4.0556621],
            [5.66844224, 4.10700593, 4.0556621],
            [13.32748799, 8.39460555, 0.22661828],
            [19.28162619, 8.75579295, -0.22661828],
            [16.07096722, 13.04339257, -4.0556621],
            [16.53814696, 4.10700593, 4.0556621],
            [7.94474171, 8.39460555, -8.33794247],
            [2.92496303, 8.75579295, 8.33794247],
            [5.2012625, 4.46819332, -4.0556621],
            [5.66844224, 12.68220518, 4.0556621],
            [18.81444643, 8.39460555, -8.33794247],
            [13.79466775, 8.75579295, 8.33794247],
            [16.07096722, 4.46819332, -4.0556621],
            [16.53814696, 12.68220518, 4.0556621]]),
        'simulation/qmcsystem/particlesets/ion0/groups/V/size' : 16,
        'simulation/qmcsystem/particlesets/ion0/groups/V/valence' : 13,
        'simulation/qmcsystem/particlesets/ion0/name' : 'ion0',
        'simulation/qmcsystem/simulationcell/bconds' : np.array([
            'p','p','p']),
        'simulation/qmcsystem/simulationcell/lattice' : np.array([
            [21.73940944, 0.0, 0.0],
            [5.48695844, 8.57519925, -8.56456075],
            [-5.48695844, 8.57519925, 8.56456075]]),
        'simulation/qmcsystem/simulationcell/lr_dim_cutoff' : 15,
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/downdet/group' : 'd',
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/downdet/id' : 'downdet',
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/downdet/size' : 200,
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/downdet/sposet' : 'spo_d',
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/updet/group' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/updet/id' : 'updet',
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/updet/size' : 200,
        'simulation/qmcsystem/wavefunctions/psi0/determinantset/slaterdeterminant/determinants/updet/sposet' : 'spo_u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/coefficients/coeff' : np.array([
            -1.488295706,-1.406709163,-1.232298155,-0.9391459067,-0.5575491618,-0.2186131788,-0.1463697747,-0.09781208605,-0.06418209044,-0.03977101442,-0.02226362717,-0.009458557456,-0.002401473122]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/coefficients/id' : 'eO',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/cusp' : 0.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/elementtype' : 'O',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/rcut' : 6.05,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/O/size' : 13,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/coefficients/coeff' : np.array([
            -2.88368129,-2.686350256,-2.500947608,-2.096756839,-1.444128943,-0.7686333881,-0.5720610092,-0.4061081504,-0.2772741837,-0.1767662649,-0.1010035901,-0.047325819,-0.01700847314]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/coefficients/id' : 'eV',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/cusp' : 0.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/elementtype' : 'V',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/rcut' : 6.05,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/correlations/V/size' : 13,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/function' : 'bspline',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/name' : 'J1',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/print' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/source' : 'ion0',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J1/type' : 'One-Body',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/coeff' : np.array([
            0.529300758,0.3529320289,0.2365993762,0.1604582152,0.1128159005,0.08243318778,0.06023602184,0.04310552718,0.02984314449,0.01958170086,0.01186100803,0.006112206499,0.002625360754]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/id' : 'ud',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/rcut' : 6.05,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/size' : 13,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/speciesa' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/speciesb' : 'd',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/coeff' : np.array([
            0.3569086717,0.2751683418,0.2058897032,0.1520886231,0.111693376,0.08181917929,0.05977972383,0.04283213009,0.02968150709,0.01944788064,0.01196129476,0.006271327336,0.002804432275]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/id' : 'uu',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/rcut' : 6.05,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/size' : 13,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/speciesa' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/speciesb' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/function' : 'bspline',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/name' : 'J2',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/print' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/type' : 'Two-Body',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/coefficients/coeff' : np.array([
            -0.004166620907,0.0003869059334,0.01344638104,-7.5215692e-05,-0.006436299048,0.0008791813519,0.007681280497,-0.006673633544,0.0300621195,0.00157665002,-0.001657156134,-0.01142258435,-0.02006687607,0.005271171591,0.01511417522,0.0008942941789,-0.002018984988,0.01595864928,0.005244762096,0.01545262066,-0.006397246289,-0.0072233246,-0.0008063061353,0.00830708478,0.001242024926,-0.0003962016339]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/coefficients/id' : 'udO',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/coefficients/optimize' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/esize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/especies1' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/especies2' : 'd',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/isize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/ispecies' : 'O',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udO/rcut' : 5.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/coefficients/coeff' : np.array([
            0.000658573315,0.005924655484,0.008096696785,0.002998451182,0.001289481835,8.390092052e-05,0.0174934698,0.004082827829,0.001656608224,-0.01638865932,0.002852247319,-0.01043954065,0.006179637761,-0.000652977982,-0.004542989787,-0.0004825008427,0.03569269894,-0.01539236687,0.007843924995,-0.009660462887,-0.01173827315,0.005074028683,0.001248279616,0.008752252359,-0.003457347502,0.0001174638519]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/coefficients/id' : 'udV',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/coefficients/optimize' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/esize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/especies1' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/especies2' : 'd',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/isize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/ispecies' : 'V',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/udV/rcut' : 5.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/coefficients/coeff' : np.array([
            -0.0006976974299,-0.001602461137,0.002262076236,-0.001250356792,-0.002453974076,0.00100226978,-0.008343708726,0.01062739293,0.01589135522,0.007887562739,-0.0005580320441,-0.01523126657,-0.009565046782,-0.0009005995139,0.01105399926,-0.0002575705031,-0.01652920678,0.00747060564,0.01464528142,0.005133083617,0.006916610617,-0.009683594066,0.001290999707,-0.001322800206,0.003931225142,-0.001163411737]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/coefficients/id' : 'uuO',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/coefficients/optimize' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/esize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/especies1' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/especies2' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/isize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/ispecies' : 'O',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuO/rcut' : 5.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/coefficients/coeff' : np.array([
            0.004388200165,0.001900643263,-0.01549468789,-0.002564479476,0.002118937653,0.0007437421471,-0.0085007067,0.009637603236,-0.01717900977,0.00186285366,-0.006121695671,0.01831402072,0.006890778761,0.003340289515,-0.001491823024,-0.001123033117,-0.008713157223,0.02100098414,-0.03224060809,-0.002479213835,0.001387768485,0.006636471962,0.0004745014561,0.001629700016,-0.001615344115,-0.0001680854702]),
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/coefficients/id' : 'uuV',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/coefficients/optimize' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/coefficients/type' : 'Array',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/esize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/especies1' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/especies2' : 'u',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/isize' : 3,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/ispecies' : 'V',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/correlations/uuV/rcut' : 5.0,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/function' : 'polynomial',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/name' : 'J3',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/print' : True,
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/source' : 'ion0',
        'simulation/qmcsystem/wavefunctions/psi0/jastrows/J3/type' : 'eeI',
        'simulation/qmcsystem/wavefunctions/psi0/name' : 'psi0',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/href' : '../scf/pwscf_output/pwscf.pwscf.h5',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/meshfactor' : 1.0,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/precision' : 'float',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/source' : 'ion0',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_d/name' : 'spo_d',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_d/size' : 200,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_d/spindataset' : 1,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_d/type' : 'bspline',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_u/name' : 'spo_u',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_u/size' : 200,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_u/spindataset' : 0,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/sposets/spo_u/type' : 'bspline',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/tilematrix' : np.array([
            [2, 0, 0],
            [0, 1, -1],
            [0, 1, 1]]),
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/truncate' : False,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/twistnum' : 0,
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/type' : 'bspline',
        'simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/version' : 0.1,
        'simulation/qmcsystem/wavefunctions/psi0/target' : 'e',
        }

    serial_references['VO2_M1_afm.in.xml'] = ref


    # reference for generate/read
    ref = ref.copy()
    ref['simulation/project/application/name'] = 'qmcpack' # name has been updated
    for k in list(ref.keys()):
        if 'jastrow' in k or 'metadata' in k:
            del ref[k]
        #end if
    #end for
    #  generated initial jastrow rather than optimized one
    ref['simulation/project/driver_version'] =  'legacy'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/coeff'] = np.array([0.44140587,0.26944819,0.15547533,0.08413778,0.04227037,0.01951441,0.00820536,0.00312028])
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/id'] = 'ud'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/coefficients/type'] = 'Array'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/rcut'] = 6.651925584744773
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/size'] = 8
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/speciesa'] = 'u'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/ud/speciesb'] = 'd'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/coeff'] = np.array([0.31314348,0.21502161,0.13496763,0.07727679,0.04023251,0.01897712,0.00807963,0.00309418])
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/id'] = 'uu'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/coefficients/type'] = 'Array'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/rcut'] = 6.651925584744773
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/size'] = 8
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/speciesa'] = 'u'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/correlations/uu/speciesb'] = 'u'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/function'] = 'bspline'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/name'] = 'J2'
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/print'] = True
    ref['simulation/qmcsystem/wavefunctions/psi0/jastrows/J2/type'] = 'Two-Body'

    serial_references['VO2_M1_afm.in.xml gen_read'] = ref


    # reference for generate
    #   legacy driver
    ref = ref.copy()
    ref['simulation/qmcsystem/simulationcell/bconds'] = tuple('ppp')
    for k in list(ref.keys()):
        if 'SpinDensity' in k:
            del ref[k]
        #end if
    #end for
    ref['simulation/project/application/version'] = '1.0'
    ref['simulation/qmcsystem/wavefunctions/psi0/sposet_builders/bspline/version'] = '0.10'
    ref['simulation/qmcsystem/hamiltonians/h0/estimators/0/grid'] = (72, 44, 44)
    ref['simulation/qmcsystem/hamiltonians/h0/estimators/0/name'] = 'SpinDensity'
    ref['simulation/qmcsystem/hamiltonians/h0/estimators/0/type'] = 'spindensity'

    serial_references['VO2_M1_afm.in.xml gen'] = ref

    #   batched driver
    ref = ref.copy()
    ref['simulation/project/driver_version'] =  'batched'
    for k in list(ref.keys()):
        if 'estimator' in k or 'calculations' in k or 'MPC' in k:
            del ref[k]
        #end if
    #end for
    ref['simulation/calculations/0/blocks'] = 70
    ref['simulation/calculations/0/method'] = 'vmc'
    ref['simulation/calculations/0/move'] = 'pbyp'
    ref['simulation/calculations/0/steps'] = 5
    ref['simulation/calculations/0/substeps'] = 2
    ref['simulation/calculations/0/timestep'] = 0.3
    ref['simulation/calculations/0/usedrift'] = False
    ref['simulation/calculations/0/walkers_per_rank'] = 64
    ref['simulation/calculations/0/warmupsteps'] = 20
    ref['simulation/calculations/1/blocks'] = 80
    ref['simulation/calculations/1/method'] = 'dmc'
    ref['simulation/calculations/1/move'] = 'pbyp'
    ref['simulation/calculations/1/nonlocalmoves'] = 'yes'
    ref['simulation/calculations/1/steps'] = 5
    ref['simulation/calculations/1/timestep'] = 0.02
    ref['simulation/calculations/1/walkers_per_rank'] = 64
    ref['simulation/calculations/1/warmupsteps'] = 2
    ref['simulation/calculations/2/blocks'] = 600
    ref['simulation/calculations/2/method'] = 'dmc'
    ref['simulation/calculations/2/move'] = 'pbyp'
    ref['simulation/calculations/2/nonlocalmoves'] = 'yes'
    ref['simulation/calculations/2/steps'] = 5
    ref['simulation/calculations/2/timestep'] = 0.005
    ref['simulation/calculations/2/walkers_per_rank'] = 64
    ref['simulation/calculations/2/warmupsteps'] = 10


    serial_references['VO2_M1_afm.in.xml batched gen'] = ref


    # references for afqmc
    ref = {
        'simulation/afqmcinfo/naea' : 5,
        'simulation/afqmcinfo/naeb' : 5,
        'simulation/afqmcinfo/name' : 'info0',
        'simulation/afqmcinfo/nmo' : 9,
        'simulation/execute/blocks' : 1000,
        'simulation/execute/estimator/block_size' : 2,
        'simulation/execute/estimator/name' : 'back_propagation',
        'simulation/execute/estimator/naverages' : 4,
        'simulation/execute/estimator/nsteps' : 200,
        'simulation/execute/estimator/onerdm' : onerdm(),
        'simulation/execute/estimator/ortho' : 1,
        'simulation/execute/ham' : 'ham0',
        'simulation/execute/info' : 'info0',
        'simulation/execute/ncores' : 1,
        'simulation/execute/nwalkers' : 10,
        'simulation/execute/prop' : 'prop0',
        'simulation/execute/steps' : 10,
        'simulation/execute/timestep' : 0.01,
        'simulation/execute/wfn' : 'wfn0',
        'simulation/execute/wset' : 'wset0',
        'simulation/hamiltonian/filename' : 'afqmc.h5',
        'simulation/hamiltonian/filetype' : 'hdf5',
        'simulation/hamiltonian/info' : 'info0',
        'simulation/hamiltonian/name' : 'ham0',
        'simulation/method' : 'afqmc',
        'simulation/project/id' : 'qmc',
        'simulation/project/series' : 0,
        'simulation/propagator/hybrid' : 'yes',
        'simulation/propagator/info' : 'info0',
        'simulation/propagator/name' : 'prop0',
        'simulation/random/seed' : 7,
        'simulation/walkerset/name' : 'wset0',
        'simulation/walkerset/type' : 'shared',
        'simulation/walkerset/walker_type' : 'CLOSED',
        'simulation/wavefunction/cutoff' : 1e-08,
        'simulation/wavefunction/filename' : 'afqmc.h5',
        'simulation/wavefunction/filetype' : 'hdf5',
        'simulation/wavefunction/info' : 'info0',
        'simulation/wavefunction/name' : 'wfn0',
        'simulation/wavefunction/type' : 'NOMSD',
        }


    serial_references['CH4_afqmc.in.xml read'] = ref

    ref = ref.copy()
    serial_references['CH4_afqmc.in.xml write'] = ref

    ref = ref.copy()
    for k in list(ref.keys()):
        if 'estimator' in k:
            del ref[k]
        #end if
    #end for
    ref['simulation/execute/estimators/back_propagation/block_size'] = 2
    ref['simulation/execute/estimators/back_propagation/name'] = 'back_propagation'
    ref['simulation/execute/estimators/back_propagation/naverages'] = 4
    ref['simulation/execute/estimators/back_propagation/nsteps'] = 200
    ref['simulation/execute/estimators/back_propagation/onerdm'] = onerdm()
    ref['simulation/execute/estimators/back_propagation/ortho'] = 1
    ref['simulation/propagator/hybrid'] = True
    serial_references['CH4_afqmc.in.xml compose'] = ref

    ref = ref.copy()
    serial_references['CH4_afqmc.in.xml gen est'] = ref

    ref = ref.copy()
    del ref['simulation/random/seed']
    del ref['simulation/afqmcinfo/naea']
    del ref['simulation/afqmcinfo/naeb']
    del ref['simulation/afqmcinfo/nmo']
    ref['simulation/hamiltonian/filename'] = 'MISSING.h5'
    ref['simulation/wavefunction/filename'] = 'MISSING.h5'
    ref['simulation/execute/blocks'] = 10000
    ref['simulation/execute/timestep'] = 0.005
    for k in list(ref.keys()):
        if 'estimator' in k:
            del ref[k]
        #end if
    #end for
    serial_references['CH4_afqmc.in.xml gen empty'] = ref

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(qi,name):
    from generic import obj
    sr = get_serial_references()[name]
    assert(len(sr)>0)
    sq = qi.serial()
    def remove_metadata(s):
        metadata_keys = []
        for k in s.keys():
            if k.startswith('_metadata'):
                metadata_keys.append(k)
            #end if
        #end for
        for k in metadata_keys:
            del s[k]
        #end for
    #end def remove_metadata
    remove_metadata(sq)
    remove_metadata(sr)
    assert check_object_eq(sq,obj(sr),bypass=True,verbose=True)
#end def check_vs_serial_reference



def test_files():
    filenames = [
        'VO2_M1_afm.in.xml',
        'CH4_afqmc.in.xml',
        ]
    files = get_files()
    assert(set(filenames)==set(files.keys()))
#end def test_files



def test_import():
    from qmcpack_input import QmcpackInput
    from qmcpack_input import simulation,meta,section
#end def test_import



def test_qixml_class_init():
    from generic import obj
    from qmcpack_input import classes

    attr_types = obj(
        tag            = str,
        identifier     = (str,tuple),
        attributes     = list,
        elements       = list,
        text           = str,
        parameters     = list,
        attribs        = list,
        costs          = list,
        h5tags         = list,
        types          = obj,
        write_types    = obj,
        attr_types     = obj,
        precision      = str,
        defaults       = obj,
        collection_id  = str,
        exp_names      = obj,
        params         = list,
        plurals_inv    = obj,
        plurals        = obj,
        expanded_names = obj,
        afqmc_order    = list,
        )
    optional = set(['expanded_names','afqmc_order'])
    assert(len(attr_types)==21)

    def valid_name(s):
        v = True
        v &= isinstance(s,str)
        v &= ' ' not in s and '-' not in s
        v &= s.lower()==s
        return v
    #end def valid_name

    for cls in classes:
        c = obj()
        d = cls.__dict__
        for k,v in d.items():
            if not k.startswith('__'):
                c[k] = v
            #end if
        #end for
        for name,type in attr_types.items():
            if name not in optional:
                assert(name in c)
            #end if
            if name in c:
                v = c[name]
                if v is not None:
                    assert(isinstance(v,type))
                    if isinstance(v,str):
                        assert(valid_name(v))
                    elif isinstance(v,list):
                        for vv in v:
                            assert(valid_name(vv) or vv=='atomic-number')
                        #end for
                    #end if
                #end if
            #end if
        #end for
    #end for

#end def test_qixml_class_init



def test_compose():
    import numpy as np
    from generic import obj
    from qmcpack_input import QmcpackInput
    from qmcpack_input import simulation,meta,section

    qi_comp = QmcpackInput(
        meta(
            lattice  = obj(units='bohr'),
            position = obj(condition='0',datatype='posArray'),
            ),
        simulation(
            project = section(
                id     ='qmc',
                series = 0,
                application = section(
                    name    = 'qmcapp',
                    role    = 'molecu',
                    class_  = 'serial',
                    version = 1.0,
                    ),
                ),
            qmcsystem = section(
                simulationcell = section(
                    lattice = np.array(
                        [[ 21.73940944, 0.00000000,  0.00000000],
                         [  5.48695844, 8.57519925, -8.56456075],
                         [ -5.48695844, 8.57519925,  8.56456075]]),
                    bconds = np.array(tuple('ppp')),
                    lr_dim_cutoff = 15,
                    ),
                particlesets = [
                    section(
                        name = 'e',
                        random = True,
                        groups = [
                            section(
                                name   = 'u',
                                size   = 200,
                                charge = -1,
                                mass   = 1.0,
                                ),
                            section(
                                name   = 'd',
                                size   = 200,
                                charge = -1,
                                mass   = 1.0,
                                ),
                            ],
                        ),
                    section(
                        name = 'ion0',
                        groups = [
                            section(
                                name         = 'V',
                                size         = 16,
                                charge       = 13,
                                valence      = 13,
                                atomicnumber = 23,
                                mass         = 92861.5851912,
                                position     = np.array([
                                    [ 2.45778327,  8.39460555,  0.22661828],
                                    [ 8.41192147,  8.75579295, -0.22661828],
                                    [ 5.20126250, 13.04339257, -4.05566210],
                                    [ 5.66844224,  4.10700593,  4.05566210],
                                    [13.32748799,  8.39460555,  0.22661828],
                                    [19.28162619,  8.75579295, -0.22661828],
                                    [16.07096722, 13.04339257, -4.05566210],
                                    [16.53814696,  4.10700593,  4.05566210],
                                    [ 7.94474171,  8.39460555, -8.33794247],
                                    [ 2.92496303,  8.75579295,  8.33794247],
                                    [ 5.20126250,  4.46819332, -4.05566210],
                                    [ 5.66844224, 12.68220518,  4.05566210],
                                    [18.81444643,  8.39460555, -8.33794247],
                                    [13.79466775,  8.75579295,  8.33794247],
                                    [16.07096722,  4.46819332, -4.05566210],
                                    [16.53814696, 12.68220518,  4.05566210]
                                    ]),
                                ),
                            section(
                                name         = 'O',
                                size         = 32,
                                charge       = 6,
                                valence      = 6,
                                atomicnumber = 8,
                                mass         = 29164.3928678,
                                position     = np.array([
                                    [ 0.00978311,  1.81708472,  1.78656736],
                                    [10.85992161, 15.33331378, -1.78656738],
                                    [ 2.75326234, 11.04571415, -2.49571300],
                                    [ 8.11644240,  6.10468435,  2.49571300],
                                    [ 2.71381355,  6.02493499,  2.55909075],
                                    [ 8.15589117, 11.12546351, -2.55909075],
                                    [ 5.45729278, 15.41306313, -1.72318961],
                                    [ 5.41241194,  1.73733537,  1.72318961],
                                    [10.87948783,  1.81708472,  1.78656736],
                                    [21.72962633, 15.33331378, -1.78656738],
                                    [13.62296706, 11.04571415, -2.49571300],
                                    [18.98614712,  6.10468435,  2.49571300],
                                    [13.58351827,  6.02493499,  2.55909075],
                                    [19.02559589, 11.12546351, -2.55909075],
                                    [16.32699750, 15.41306313, -1.72318961],
                                    [16.28211666,  1.73733537,  1.72318961],
                                    [ 0.00978311, 10.39228397,  1.78656736],
                                    [10.85992161,  6.75811453, -1.78656738],
                                    [-2.73369610, 11.04571415,  6.06884775],
                                    [13.60340084,  6.10468435, -6.06884775],
                                    [ 8.20077199,  6.02493499, -6.00547000],
                                    [ 2.66893273, 11.12546351,  6.00547000],
                                    [ 5.45729278,  6.83786388, -1.72318961],
                                    [ 5.41241194, 10.31253462,  1.72318961],
                                    [10.87948783, 10.39228397,  1.78656736],
                                    [21.72962633,  6.75811453, -1.78656738],
                                    [ 8.13600862, 11.04571415,  6.06884775],
                                    [24.47310556,  6.10468435, -6.06884775],
                                    [19.07047671,  6.02493499, -6.00547000],
                                    [13.53863745, 11.12546351,  6.00547000],
                                    [16.32699750,  6.83786388, -1.72318961],
                                    [16.28211666, 10.31253462,  1.72318961]
                                    ]),
                                ),
                            ],
                        ),
                    ],
                wavefunction = section(
                    name   = 'psi0',
                    target = 'e',
                    sposet_builder = section(
                        type       = 'bspline',
                        href       = '../scf/pwscf_output/pwscf.pwscf.h5',
                        tilematrix = np.array([[2,0,0],
                                               [0,1,-1],
                                               [0,1,1]]),
                        twistnum   = 0,
                        source     = 'ion0',
                        version    = 0.10,
                        meshfactor = 1.0,
                        precision  = 'float',
                        truncate   = False,
                        sposets = [
                            section(
                                type        = 'bspline',
                                name        = 'spo_u',
                                size        = 200,
                                spindataset = 0,
                                ),
                            section(
                                type        = 'bspline',
                                name        = 'spo_d',
                                size        = 200,
                                spindataset = 1,
                                )
                            ],
                        ),
                    determinantset = section(
                        slaterdeterminant = section(
                            determinants = [
                                section(
                                    id     = 'updet',
                                    group  = 'u',
                                    sposet = 'spo_u',
                                    size   = 200,
                                    ),
                                section(
                                    id     = 'downdet',
                                    group  = 'd',
                                    sposet = 'spo_d',
                                    size   = 200,
                                    ),
                                ],
                            ),
                        ),
                    jastrows = [
                        section(
                            type     = 'One-Body',
                            name     = 'J1',
                            function = 'bspline',
                            source   = 'ion0',
                            print    = True,
                            correlations = [
                                section(
                                    elementType = 'O',
                                    size        = 13,
                                    rcut        = 6.05,
                                    cusp        = 0.0,
                                    coefficients = section(
                                        id   = 'eO',
                                        type = 'Array',
                                        coeff = np.array([
                                            -1.488295706, -1.406709163,
                                            -1.232298155, -0.9391459067, 
                                            -0.5575491618, -0.2186131788, 
                                            -0.1463697747, -0.09781208605, 
                                            -0.06418209044, -0.03977101442, 
                                            -0.02226362717, -0.009458557456, 
                                            -0.002401473122])
                                        ),
                                    ),
                                section(
                                    elementType = 'V',
                                    size        = 13,
                                    rcut        = 6.05,
                                    cusp        = 0.0,
                                    coefficients = section(
                                        id   = 'eV',
                                        type = 'Array',
                                        coeff = np.array([
                                            -2.88368129, -2.686350256, 
                                            -2.500947608, -2.096756839, 
                                            -1.444128943, -0.7686333881, 
                                            -0.5720610092, -0.4061081504,
                                            -0.2772741837, -0.1767662649, 
                                            -0.1010035901, -0.047325819, 
                                            -0.01700847314])
                                        ),
                                    ),
                                ],
                            ),
                        section(
                            type     = 'Two-Body',
                            name     = 'J2',
                            function = 'bspline',
                            print    = True,
                            correlations = [
                                section(
                                    speciesA    = 'u',
                                    speciesB    = 'u',
                                    size        = 13,
                                    rcut        = 6.05,
                                    coefficients = section(
                                        id    = 'uu',
                                        type  = 'Array',
                                        coeff = np.array([
                                            0.3569086717, 0.2751683418, 
                                            0.2058897032, 0.1520886231, 
                                            0.111693376, 0.08181917929, 
                                            0.05977972383, 0.04283213009, 
                                            0.02968150709, 0.01944788064, 
                                            0.01196129476, 0.006271327336, 
                                            0.002804432275])
                                        ),
                                    ),
                                section(
                                    speciesA    = 'u',
                                    speciesB    = 'd',
                                    size        = 13,
                                    rcut        = 6.05,
                                    coefficients = section(
                                        id    = 'ud',
                                        type  = 'Array',
                                        coeff = np.array([
                                            0.529300758, 0.3529320289, 
                                            0.2365993762, 0.1604582152, 
                                            0.1128159005, 0.08243318778, 
                                            0.06023602184, 0.04310552718, 
                                            0.02984314449, 0.01958170086, 
                                            0.01186100803, 0.006112206499, 
                                            0.002625360754])
                                        ),
                                    ),
                                ],
                            ),
                        section(
                            type     = 'eeI',
                            name     = 'J3',
                            function = 'polynomial',
                            print    = True,
                            source   = 'ion0',
                            correlations = [
                                section(
                                    ispecies = 'O',
                                    especies1= 'u',
                                    especies2= 'u',
                                    isize    = 3,
                                    esize    = 3,
                                    rcut     = 5.0,
                                    coefficients = section(
                                        id       = 'uuO',
                                        type     = 'Array',
                                        optimize = True,
                                        coeff    = np.array([
                                            -0.0006976974299, -0.001602461137, 
                                            0.002262076236, -0.001250356792, 
                                            -0.002453974076, 0.00100226978, 
                                            -0.008343708726, 0.01062739293, 
                                            0.01589135522, 0.007887562739, 
                                            -0.0005580320441, -0.01523126657, 
                                            -0.009565046782, -0.0009005995139, 
                                            0.01105399926, -0.0002575705031, 
                                            -0.01652920678, 0.00747060564, 
                                            0.01464528142, 0.005133083617, 
                                            0.006916610617, -0.009683594066, 
                                            0.001290999707, -0.001322800206, 
                                            0.003931225142, -0.001163411737])
                                        ),
                                    ),
                                section(
                                    ispecies = 'O',
                                    especies1= 'u',
                                    especies2= 'd',
                                    isize    = 3,
                                    esize    = 3,
                                    rcut     = 5.0,
                                    coefficients = section(
                                        id       = 'udO',
                                        type     = 'Array',
                                        optimize = True,
                                        coeff    = np.array([
                                            -0.004166620907, 0.0003869059334, 
                                            0.01344638104, -7.5215692e-05, 
                                            -0.006436299048, 0.0008791813519, 
                                            0.007681280497, -0.006673633544, 
                                            0.0300621195, 0.00157665002, 
                                            -0.001657156134, -0.01142258435, 
                                            -0.02006687607, 0.005271171591, 
                                            0.01511417522, 0.0008942941789, 
                                            -0.002018984988, 0.01595864928, 
                                            0.005244762096, 0.01545262066, 
                                            -0.006397246289, -0.0072233246, 
                                            -0.0008063061353, 0.00830708478, 
                                            0.001242024926, -0.0003962016339])
                                        ),
                                    ),
                                section(
                                    ispecies = 'V',
                                    especies1= 'u',
                                    especies2= 'u',
                                    isize    = 3,
                                    esize    = 3,
                                    rcut     = 5.0,
                                    coefficients = section(
                                        id       = 'uuV',
                                        type     = 'Array',
                                        optimize = True,
                                        coeff    = np.array([
                                            0.004388200165, 0.001900643263, 
                                            -0.01549468789, -0.002564479476, 
                                            0.002118937653, 0.0007437421471, 
                                            -0.0085007067, 0.009637603236, 
                                            -0.01717900977, 0.00186285366, 
                                            -0.006121695671, 0.01831402072,
                                            0.006890778761, 0.003340289515, 
                                            -0.001491823024, -0.001123033117, 
                                            -0.008713157223, 0.02100098414, 
                                            -0.03224060809, -0.002479213835, 
                                            0.001387768485, 0.006636471962, 
                                            0.0004745014561, 0.001629700016, 
                                            -0.001615344115, -0.0001680854702])
                                        ),
                                    ),
                                section(
                                    ispecies = 'V',
                                    especies1= 'u',
                                    especies2= 'd',
                                    isize    = 3,
                                    esize    = 3,
                                    rcut     = 5.0,
                                    coefficients = section(
                                        id       = 'udV',
                                        type     = 'Array',
                                        optimize = True,
                                        coeff    = np.array([
                                            0.000658573315, 0.005924655484, 
                                            0.008096696785, 0.002998451182, 
                                            0.001289481835, 8.390092052e-05, 
                                            0.0174934698, 0.004082827829, 
                                            0.001656608224, -0.01638865932, 
                                            0.002852247319, -0.01043954065, 
                                            0.006179637761, -0.000652977982, 
                                            -0.004542989787, -0.0004825008427, 
                                            0.03569269894, -0.01539236687, 
                                            0.007843924995, -0.009660462887, 
                                            -0.01173827315, 0.005074028683, 
                                            0.001248279616, 0.008752252359, 
                                            -0.003457347502, 0.0001174638519])
                                        ),
                                    ),
                                ],
                            ),
                        ],
                    ),
                hamiltonian = section(
                    name   = 'h0',
                    type   = 'generic',
                    target = 'e',
                    pairpots = [
                        section( 
                            type   = 'coulomb',
                            name   = 'ElecElec',
                            source = 'e',
                            target = 'e'
                            ),
                        section(
                            type   = 'coulomb',
                            name   = 'IonIon',
                            source = 'ion0',
                            target = 'ion0'
                            ),
                        section(
                            type         = 'pseudo',
                            name         = 'PseudoPot',
                            source       = 'ion0',
                            wavefunction = 'psi0',
                            format       = 'xml',
                            pseudos = [
                                section(
                                    elementType = 'O',
                                    href        = 'O.opt.xml'
                                    ),
                                section(
                                    elementType = 'V',
                                    href        = 'V.opt.xml'
                                    ),
                                ]
                            ),
                        section(
                            type     = 'MPC',
                            name     = 'MPC',
                            source   = 'e',
                            target   = 'e',
                            ecut     = 60.0,
                            physical = False
                            )
                        ],
                    estimators = [
                        section(
                            type = 'spindensity',
                            name = 'SpinDensity',
                            grid = np.array([72,44,44]),
                            ),
                        section(
                            name   = 'KEcorr',
                            type   = 'chiesa',
                            source = 'e',
                            psi    = 'psi0'
                            )
                        ],
                    )
                ),
            calculations = [
                section(
                    method           = 'vmc',
                    move             = 'pbyp',
                    checkpoint       = -1,
                    walkers          = 1,
                    blocks           = 70,
                    steps            = 5,
                    substeps         = 2,
                    timestep         = 0.3,
                    warmupsteps      = 20,
                    samplesperthread = 2,
                    ),
                section(
                    method           = 'dmc',
                    move             = 'pbyp',
                    checkpoint       = -1,
                    blocks           = 80,
                    steps            = 5,
                    timestep         = 0.02,
                    nonlocalmoves    = 'yes',
                    warmupsteps      = 2,
                    ),
                section(
                    method           = 'dmc',
                    move             = 'pbyp',
                    checkpoint       = -1,
                    blocks           = 600,
                    steps            = 5,
                    timestep         = 0.005,
                    nonlocalmoves    = 'yes',
                    warmupsteps      = 10,
                    ),
                ],
            ),
        )
    qi_comp.pluralize()

    check_vs_serial_reference(qi_comp,'VO2_M1_afm.in.xml')


    qi_afqmc = QmcpackInput(
        meta(),
        simulation(
            method = 'afqmc',
            project = section(
                id     = 'qmc',
                series = 0,
                ),
            random = section(
                seed = 7
                ),
            afqmcinfo = section(
                name = 'info0',
                nmo  = 9,
                naea = 5,
                naeb = 5,
                ),
            hamiltonian = section(
                name     = 'ham0',
                info     = 'info0',
                filetype = 'hdf5',
                filename = 'afqmc.h5',
                ),
            wavefunction = section(
                type     = 'NOMSD',
                name     = 'wfn0',
                info     = 'info0',
                filetype = 'hdf5',
                filename = 'afqmc.h5',
                cutoff   = 1e-8,
                ),
            walkerset = section(
                type        = 'shared',
                name        = 'wset0',
                walker_type = 'CLOSED',
                ),
            propagator = section(
                name   = 'prop0',
                info   = 'info0',
                hybrid = True,
                ),
            execute = section(
                info     = 'info0',
                ham      = 'ham0',
                wfn      = 'wfn0',
                wset     = 'wset0',
                prop     = 'prop0',
                blocks   = 1000,
                timestep = 0.01,
                steps    = 10,
                ncores   = 1,
                nwalkers = 10,
                estimators = [
                    section(
                        name       = 'back_propagation',
                        naverages  = 4,
                        block_size = 2,
                        ortho      = 1,
                        nsteps     = 200,
                        onerdm     = section(),
                        )
                    ],
                ),
            )
        )

    check_vs_serial_reference(qi_afqmc,'CH4_afqmc.in.xml compose')

#end def test_compose



def test_generate():
    import numpy as np
    from physical_system import generate_physical_system
    from qmcpack_input import generate_qmcpack_input,spindensity
    from qmcpack_input import back_propagation,onerdm
    
    system = generate_physical_system(
        units    = 'A',
        axes     = '''
                    5.75200000    0.00000000    0.00000000
                    0.00000000    4.53780000    0.00000000
                   -2.90357335    0.00000000    4.53217035
                   ''',
        elem_pos = '''
                   V   1.30060289    4.44223393    0.11992123
                   V   1.54782377    0.09556607    4.41224912
                   V  -0.15118378    2.36446607    2.38600640
                   V   2.99961044    2.17333393    2.14616395
                   O   0.00517700    0.96155982    0.94541073
                   O   2.84324965    3.57624018    3.58675961
                   O  -1.44660967    1.30734018    3.21149591
                   O   4.29503633    3.23045982    1.32067444
                   O   1.43608828    3.18825828    1.35421250
                   O   1.41233837    1.34954172    3.17795785
                   O  -0.01569839    3.61844172    3.62029768
                   O   2.86412504    0.91935828    0.91187267
                   ''',
        tiling = [[ 2,  0,  0],
                  [ 0,  1, -1],
                  [ 0,  1,  1]],
        V = 13,
        O = 6,
        )

    tile_pos = np.array([
        [ 1.30060289e+00,  4.44223393e+00,  1.19921230e-01],
        [ 4.45139712e+00,  4.63336607e+00, -1.19921230e-01],
        [ 2.75238957e+00,  6.90226607e+00, -2.14616395e+00],
        [ 2.99961044e+00,  2.17333393e+00,  2.14616395e+00],
        [ 5.17700000e-03,  9.61559820e-01,  9.45410730e-01],
        [ 5.74682300e+00,  8.11404018e+00, -9.45410740e-01],
        [ 1.45696368e+00,  5.84514018e+00, -1.32067444e+00],
        [ 4.29503633e+00,  3.23045982e+00,  1.32067444e+00],
        [ 1.43608828e+00,  3.18825828e+00,  1.35421250e+00],
        [ 4.31591172e+00,  5.88734172e+00, -1.35421250e+00],
        [ 2.88787496e+00,  8.15624172e+00, -9.11872670e-01],
        [ 2.86412504e+00,  9.19358280e-01,  9.11872670e-01],
        [ 7.05260289e+00,  4.44223393e+00,  1.19921230e-01],
        [ 1.02033971e+01,  4.63336607e+00, -1.19921230e-01],
        [ 8.50438957e+00,  6.90226607e+00, -2.14616395e+00],
        [ 8.75161044e+00,  2.17333393e+00,  2.14616395e+00],
        [ 5.75717700e+00,  9.61559820e-01,  9.45410730e-01],
        [ 1.14988230e+01,  8.11404018e+00, -9.45410740e-01],
        [ 7.20896368e+00,  5.84514018e+00, -1.32067444e+00],
        [ 1.00470363e+01,  3.23045982e+00,  1.32067444e+00],
        [ 7.18808828e+00,  3.18825828e+00,  1.35421250e+00],
        [ 1.00679117e+01,  5.88734172e+00, -1.35421250e+00],
        [ 8.63987496e+00,  8.15624172e+00, -9.11872670e-01],
        [ 8.61612504e+00,  9.19358280e-01,  9.11872670e-01],
        [ 4.20417624e+00,  4.44223393e+00, -4.41224912e+00],
        [ 1.54782377e+00,  4.63336607e+00,  4.41224912e+00],
        [ 2.75238957e+00,  2.36446607e+00, -2.14616395e+00],
        [ 2.99961044e+00,  6.71113393e+00,  2.14616395e+00],
        [ 5.17700000e-03,  5.49935982e+00,  9.45410730e-01],
        [ 5.74682300e+00,  3.57624018e+00, -9.45410740e-01],
        [-1.44660967e+00,  5.84514018e+00,  3.21149591e+00],
        [ 7.19860968e+00,  3.23045982e+00, -3.21149591e+00],
        [ 4.33966163e+00,  3.18825828e+00, -3.17795785e+00],
        [ 1.41233837e+00,  5.88734172e+00,  3.17795785e+00],
        [ 2.88787496e+00,  3.61844172e+00, -9.11872670e-01],
        [ 2.86412504e+00,  5.45715828e+00,  9.11872670e-01],
        [ 9.95617624e+00,  4.44223393e+00, -4.41224912e+00],
        [ 7.29982377e+00,  4.63336607e+00,  4.41224912e+00],
        [ 8.50438957e+00,  2.36446607e+00, -2.14616395e+00],
        [ 8.75161044e+00,  6.71113393e+00,  2.14616395e+00],
        [ 5.75717700e+00,  5.49935982e+00,  9.45410730e-01],
        [ 1.14988230e+01,  3.57624018e+00, -9.45410740e-01],
        [ 4.30539033e+00,  5.84514018e+00,  3.21149591e+00],
        [ 1.29506097e+01,  3.23045982e+00, -3.21149591e+00],
        [ 1.00916616e+01,  3.18825828e+00, -3.17795785e+00],
        [ 7.16433837e+00,  5.88734172e+00,  3.17795785e+00],
        [ 8.63987496e+00,  3.61844172e+00, -9.11872670e-01],
        [ 8.61612504e+00,  5.45715828e+00,  9.11872670e-01]],dtype=float)

    s = system.structure
    if not value_eq(s.pos,tile_pos):
        tile_elem = np.array([
            'V','V','V','V','O','O','O','O','O','O','O','O',
            'V','V','V','V','O','O','O','O','O','O','O','O',
            'V','V','V','V','O','O','O','O','O','O','O','O',
            'V','V','V','V','O','O','O','O','O','O','O','O'],dtype=object)
        s.pos  = tile_pos
        s.elem = tile_elem
    #end if

    #system.structure.order_by_species()

    # legacy drivers
    qi = generate_qmcpack_input(
        input_type      = 'basic',
        system          = system,
        randomsrc       = False,
        pseudos         = ['V.opt.xml','O.opt.xml'],
        spin_polarized  = True,
        twistnum        = 0,
        orbitals_h5     = '../scf/pwscf_output/pwscf.pwscf.h5',
        check_paths     = False,
        estimators      = [spindensity(grid=(72,44,44))],
        qmc = 'dmc',
        # vmc inputs
        vmc_walkers     = 1,
        vmc_warmupsteps = 20,
        vmc_blocks      = 70,
        vmc_steps       = 5,
        vmc_substeps    = 2,
        vmc_timestep    = 0.3,
        vmc_samplesperthread = 2,
        # dmc inputs
        eq_dmc          = True,
        eq_warmupsteps  = 2,
        eq_blocks       = 80,
        eq_steps        =  5,
        eq_timestep     = 0.02,
        # dmc inputs
        warmupsteps     = 10,
        blocks          = 600,
        steps           =  5,
        timestep        = 0.005,
        nonlocalmoves   = 'yes',
        )
    
    qi.pluralize()

    check_vs_serial_reference(qi,'VO2_M1_afm.in.xml gen')


    # batched drivers
    qi = generate_qmcpack_input(
        input_type       = 'basic',
        system           = system,
        randomsrc        = False,
        pseudos          = ['V.opt.xml','O.opt.xml'],
        spin_polarized   = True,
        twistnum         = 0,
        orbitals_h5      = '../scf/pwscf_output/pwscf.pwscf.h5',
        check_paths      = False,
        estimators       = [],
        corrections      = [],
        qmc              = 'dmc',
        driver           = 'batched',
        walkers_per_rank = 64,
        # vmc inputs
        vmc_warmupsteps  = 20,
        vmc_blocks       = 70,
        vmc_steps        = 5,
        vmc_substeps     = 2,
        vmc_timestep     = 0.3,
        # dmc inputs
        eq_dmc           = True,
        eq_warmupsteps   = 2,
        eq_blocks        = 80,
        eq_steps         =  5,
        eq_timestep      = 0.02,
        # dmc inputs     
        warmupsteps      = 10,
        blocks           = 600,
        steps            =  5,
        timestep         = 0.005,
        nonlocalmoves    = 'yes',
        )
    
    qi.pluralize()

    check_vs_serial_reference(qi,'VO2_M1_afm.in.xml batched gen')


    # test afqmc
    qi = generate_qmcpack_input(
        input_type = 'basic_afqmc',
        )
    check_vs_serial_reference(qi,'CH4_afqmc.in.xml gen empty')

    qi = generate_qmcpack_input(
        input_type = 'basic_afqmc',
        seed       = 7,
        nmo        = 9,
        naea       = 5,
        naeb       = 5,
        ham_file   = 'afqmc.h5',
        blocks     = 1000,
        timestep   = 0.01,
        estimators = [
            back_propagation(
                naverages  = 4,
                block_size = 2,
                ortho      = 1,
                nsteps     = 200,
                onerdm     = onerdm(),
                ),
            ],
        )

    check_vs_serial_reference(qi,'CH4_afqmc.in.xml gen est')

#end def test_generate



def test_read():
    from qmcpack_input import QmcpackInput

    files = get_files()

    qi_read = QmcpackInput(files['VO2_M1_afm.in.xml'])
    assert(not qi_read.is_afqmc_input())
    qi_read.pluralize()
    assert(not qi_read.is_afqmc_input())

    # remove extraneous data members for purpose of comparison
    del qi_read._metadata.spo_u
    del qi_read._metadata.spo_d
    spob = qi_read.simulation.qmcsystem.wavefunctions.psi0.sposet_builders
    sposets = spob.bspline.sposets
    del sposets.spo_u.spos
    del sposets.spo_d.spos

    check_vs_serial_reference(qi_read,'VO2_M1_afm.in.xml')


    # test read for afqmc input file
    qi = QmcpackInput(files['CH4_afqmc.in.xml'])
    assert(qi.is_afqmc_input())

    check_vs_serial_reference(qi,'CH4_afqmc.in.xml read')

#end def test_read



def test_write():
    import os
    from qmcpack_input import QmcpackInput

    tpath = testing.setup_unit_test_output_directory('qmcpack_input','test_write')

    files = get_files()

    # test write real space qmc input file
    ref_file   = 'VO2_M1_afm.in.xml'
    write_file = os.path.join(tpath,'write_'+ref_file)

    qi_read = QmcpackInput(files[ref_file])

    text = qi_read.write()
    assert('<simulation>' in text)
    assert('<project' in text)
    assert('<application' in text)
    assert('<qmcsystem>' in text)
    assert('<simulationcell>' in text)
    assert('<particleset' in text)
    assert('<group' in text)
    assert('<attrib' in text)
    assert('<wavefunction' in text)
    assert('<sposet_builder' in text)
    assert('<determinantset>' in text)
    assert('<slaterdeterminant>' in text)
    assert('<determinant' in text)
    assert('<jastrow' in text)
    assert('<hamiltonian' in text)
    assert('<pairpot' in text)
    assert('<pseudo' in text)
    assert('<estimator' in text)
    assert('<qmc' in text)
    assert('<parameter' in text)
    assert('LR_dim_cutoff' in text)
    assert('"posArray"' in text)
    assert('elementType' in text)
    assert('type="spindensity"' in text)

    qi_read.write(write_file)

    qi_write = QmcpackInput(write_file)
    qi_write.pluralize()

    # remove extraneous data members for purpose of comparison
    del qi_write._metadata.spo_u
    del qi_write._metadata.spo_d
    spob = qi_write.simulation.qmcsystem.wavefunctions.psi0.sposet_builders
    sposets = spob.bspline.sposets
    del sposets.spo_u.spos
    del sposets.spo_d.spos

    check_vs_serial_reference(qi_write,ref_file)


    # test write for afqmc input file
    ref_file   = 'CH4_afqmc.in.xml'
    write_file = os.path.join(tpath,'write_'+ref_file)

    qi_read = QmcpackInput(files[ref_file])

    text = qi_read.write()
    assert('<simulation' in text)
    assert('<project' in text)
    assert('<random' in text)
    assert('<AFQMCInfo' in text)
    assert('<Hamiltonian' in text)
    assert('<Wavefunction' in text)
    assert('<WalkerSet' in text)
    assert('<Propagator' in text)
    assert('<execute' in text)
    assert('<Estimator' in text)
    assert('<OneRDM/>' in text)
    assert('<parameter' in text)
    assert('method="afqmc"' in text)
    assert('"NMO"' in text)
    assert('"NAEA"' in text)
    assert('"NAEB"' in text)
    assert('nWalkers' in text)
    assert('name="wfn0"' in text)
    assert('type="NOMSD"' in text)
    assert('afqmc.h5' in text)
    assert('yes' in text)
    assert('name="back_propagation"' in text)

    qi_read.write(write_file)

    qi_write = QmcpackInput(write_file)
    assert(qi_write.is_afqmc_input())

    check_vs_serial_reference(qi_write,'CH4_afqmc.in.xml write')
#end def test_write



def test_get():
    from qmcpack_input import QmcpackInput

    files = get_files()

    qi = QmcpackInput(files['VO2_M1_afm.in.xml'])

    s    = qi.simulation
    p    = s.project
    qs   = s.qmcsystem
    sc   = qs.simulationcell
    pss  = qs.particlesets
    pse  = pss.e
    psi  = pss.ion0
    u    = pse.groups.u
    d    = pse.groups.d
    V    = psi.groups.V
    O    = psi.groups.O
    wf   = qs.wavefunction
    sb   = wf.sposet_builder
    spos = sb.sposets
    spou = spos.spo_u
    spod = spos.spo_d
    ds   = wf.determinantset
    sd   = ds.slaterdeterminant
    du   = sd.determinants.updet
    dd   = sd.determinants.downdet
    js   = wf.jastrows
    J1   = js.J1
    J2   = js.J2
    J3   = js.J3
    h   = qs.hamiltonian
    pps  = h.pairpots
    ee   = pps.ElecElec
    ii   = pps.IonIon
    ecp  = pps.PseudoPot
    ecps = ecp.pseudos
    mpc  = pps.MPC
    ests = h.estimators
    sdens= ests.SpinDensity
    kec  = ests.KEcorr
    vmc  = s.calculations[0]
    dmc1 = s.calculations[1]
    dmc2 = s.calculations[2]

    search_map = dict(
        simulation        = s,
        project           = p,
        qmcsystem         = qs,
        simulationcell    = sc,
        particleset       = pss,
        e                 = pse,
        ion0              = psi,
        u                 = u,
        d                 = d,
        #V                 = V, # can find group or pseudo
        #O                 = O,
        wavefunction      = wf,
        sposet_builder    = sb,
        sposet            = spos,
        spo_u             = spou,
        spo_d             = spod,
        determinantset    = ds,
        slaterdeterminant = sd,
        updet             = du,
        downdet           = dd,
        jastrow           = js,
        J1                = J1,
        J2                = J2,
        J3                = J3,
        hamiltonian       = h,
        pairpot           = pps,
        ElecElec          = ee,
        IonIon            = ii,
        PseudoPot         = ecp,
        pseudo            = ecps,
        MPC               = mpc,
        estimator         = ests,
        SpinDensity       = sdens,
        KEcorr            = kec,
        series            = p.series,
        lattice           = sc.lattice,
        bconds            = sc.bconds,
        lr_dim_cutoff     = sc.lr_dim_cutoff,
        random            = pse.random,
        tilematrix        = sb.tilematrix,
        twistnum          = sb.twistnum,
        meshfactor        = sb.meshfactor,
        precision         = sb.precision,
        truncate          = sb.truncate,
        format            = ecp.format,
        ecut              = mpc.ecut,
        physical          = mpc.physical,
        grid              = sdens.grid,
        psi               = kec.psi,
        )

    missing = 'a b c'.split()


    for k,vref in search_map.items():
        v = qi.get(k)
        assert(v is not None)
        assert(id(v)==id(vref))
    #end for

    for k in missing:
        v = qi.get(k)
        assert(v is None)
    #end for


    qi.pluralize()

    s    = qi.simulation
    p    = s.project
    qs   = s.qmcsystem
    sc   = qs.simulationcell
    pss  = qs.particlesets
    pse  = pss.e
    psi  = pss.ion0
    u    = pse.groups.u
    d    = pse.groups.d
    V    = psi.groups.V
    O    = psi.groups.O
    wfs  = qs.wavefunctions
    wf   = wfs.psi0
    sbs  = wf.sposet_builders
    sb   = sbs.bspline
    spos = sb.sposets
    spou = spos.spo_u
    spod = spos.spo_d
    ds   = wf.determinantset
    sd   = ds.slaterdeterminant
    du   = sd.determinants.updet
    dd   = sd.determinants.downdet
    js   = wf.jastrows
    J1   = js.J1
    J2   = js.J2
    J3   = js.J3
    hs   = qs.hamiltonians
    h    = hs.h0
    pps  = h.pairpots
    ee   = pps.ElecElec
    ii   = pps.IonIon
    ecp  = pps.PseudoPot
    ecps = ecp.pseudos
    mpc  = pps.MPC
    ests = h.estimators
    sdens= ests.SpinDensity
    kec  = ests.KEcorr
    vmc  = s.calculations[0]
    dmc1 = s.calculations[1]
    dmc2 = s.calculations[2]

    search_map = dict(
        simulation        = s,
        project           = p,
        qmcsystem         = qs,
        simulationcell    = sc,
        particleset       = pss,
        e                 = pse,
        ion0              = psi,
        u                 = u,
        d                 = d,
        #V                 = V, # can find group or pseudo
        #O                 = O,
        wavefunction      = wfs,
        sposet_builder    = sbs,
        psi0              = wf,
        sposet            = spos,
        spo_u             = spou,
        spo_d             = spod,
        determinantset    = ds,
        slaterdeterminant = sd,
        updet             = du,
        downdet           = dd,
        jastrow           = js,
        J1                = J1,
        J2                = J2,
        J3                = J3,
        hamiltonian       = hs,
        h0                = h,
        pairpot           = pps,
        ElecElec          = ee,
        IonIon            = ii,
        PseudoPot         = ecp,
        pseudo            = ecps,
        MPC               = mpc,
        estimator         = ests,
        SpinDensity       = sdens,
        KEcorr            = kec,
        series            = p.series,
        lattice           = sc.lattice,
        bconds            = sc.bconds,
        lr_dim_cutoff     = sc.lr_dim_cutoff,
        random            = pse.random,
        tilematrix        = sb.tilematrix,
        twistnum          = sb.twistnum,
        meshfactor        = sb.meshfactor,
        precision         = sb.precision,
        truncate          = sb.truncate,
        format            = ecp.format,
        ecut              = mpc.ecut,
        physical          = mpc.physical,
        grid              = sdens.grid,
        psi               = kec.psi,
        )

    for k,vref in search_map.items():
        v = qi.get(k)
        assert(v is not None)
        assert(id(v)==id(vref))
    #end for

    for k in missing:
        v = qi.get(k)
        assert(v is None)
    #end for


    qi = QmcpackInput(files['CH4_afqmc.in.xml'])

    s  = qi.simulation
    pr = s.project
    r  = s.random
    ai = s.afqmcinfo
    h  = s.hamiltonian
    wf = s.wavefunction
    ws = s.walkerset
    p  = s.propagator
    e  = s.execute
    bp = e.estimator
    dm = bp.onerdm

    search_map = dict(
        simulation   = s,
        project      = pr,
        random       = r,
        afqmcinfo    = ai,
        hamiltonian  = h,
        wavefunction = wf,
        walkerset    = ws,
        propagator   = p,
        execute      = e,
        estimator    = bp,
        onerdm       = dm,
        id           = pr.id,
        series       = pr.series,
        seed         = r.seed,
        cutoff       = wf.cutoff,
        walker_type  = ws.walker_type,
        hybrid       = p.hybrid,
        blocks       = e.blocks,
        timestep     = e.timestep,
        steps        = e.steps,
        ncores       = e.ncores,
        nwalkers     = e.nwalkers,
        naverages    = bp.naverages,
        block_size   = bp.block_size,
        ortho        = bp.ortho,
        nsteps       = bp.nsteps,
        )

    for k,vref in search_map.items():
        v = qi.get(k)
        assert(v is not None)
        assert(id(v)==id(vref))
    #end for

    for k in missing:
        v = qi.get(k)
        assert(v is None)
    #end for

#end def test_get



def test_incorporate_system():
    from physical_system import generate_physical_system
    from qmcpack_input import generate_qmcpack_input

    divert_nexus_log()

    system = generate_physical_system(
        units    = 'A',
        axes     = '''
                    5.75200000    0.00000000    0.00000000
                    0.00000000    4.53780000    0.00000000
                   -2.90357335    0.00000000    4.53217035
                   ''',
        elem_pos = '''
                   V   1.30060289    4.44223393    0.11992123
                   V   1.54782377    0.09556607    4.41224912
                   V  -0.15118378    2.36446607    2.38600640
                   V   2.99961044    2.17333393    2.14616395
                   O   0.00517700    0.96155982    0.94541073
                   O   2.84324965    3.57624018    3.58675961
                   O  -1.44660967    1.30734018    3.21149591
                   O   4.29503633    3.23045982    1.32067444
                   O   1.43608828    3.18825828    1.35421250
                   O   1.41233837    1.34954172    3.17795785
                   O  -0.01569839    3.61844172    3.62029768
                   O   2.86412504    0.91935828    0.91187267
                   ''',
        V = 13,
        O = 6,
        )

    qi = generate_qmcpack_input(
        input_type      = 'basic',
        system          = system,
        pseudos         = ['V.opt.xml','O.opt.xml'],
        spin_polarized  = True,
        twistnum        = 0,
        orbitals_h5     = 'scf.pwscf.h5',
        check_paths     = False,
        qmc             = 'dmc',
        )

    qi_ref = qi.copy()

    shift = 0.1
    s = system.structure
    s.axes += shift
    s.pos  += shift

    qi.incorporate_system(system)

    axes_ref = qi_ref.get('lattice')
    psi_ref  = qi_ref.get('ion0')

    axes     = qi.get('lattice')
    psi      = qi.get('ion0')
    
    assert(value_eq(axes-shift,axes_ref))
    assert(value_eq(psi.groups.V.position-shift,psi_ref.groups.V.position))
    assert(value_eq(psi.groups.O.position-shift,psi_ref.groups.O.position))

    del qi_ref.get('simulationcell').lattice
    del qi_ref.get('qmcsystem').particlesets
    del qi.get('simulationcell').lattice
    del qi.get('qmcsystem').particlesets

    assert(object_eq(qi,qi_ref))

    restore_nexus_log()

#end def test_incorporate_system



def test_generate_kspace_jastrow():
    from qmcpack_input import generate_kspace_jastrow
    kjas = generate_kspace_jastrow(1.0, 2.0, 2, 4)
    expect = '''<jastrow type="kSpace" name="Jk" source="ion0">
   <correlation kc="1.0" type="One-Body" symmetry="isotropic">
      <coefficients id="cG1" type="Array">         
0 0
      </coefficients>
   </correlation>
   <correlation kc="2.0" type="Two-Body" symmetry="isotropic">
      <coefficients id="cG2" type="Array">         
0 0 0 0
      </coefficients>
   </correlation>
</jastrow>
'''
    text = kjas.write()
    assert text == expect
#end def test_generate_kspace_jastrow



def test_excited_state():
    from nexus import generate_physical_system
    from nexus import generate_qmcpack_input

    dia = generate_physical_system(
        units     = 'A',
        axes      = [[ 1.785,  1.785,  0.   ],
                     [ 0.   ,  1.785,  1.785],
                     [ 1.785,  0.   ,  1.785]],
        elem      = ['C','C'],
        pos       = [[ 0.    ,  0.    ,  0.    ],
                     [ 0.8925,  0.8925,  0.8925]],
        tiling    = [3,1,3], 
        kgrid     = (1,1,1), 
        kshift    = (0,0,0), 
        C         = 4
        )
  
  
    # test kp_index, band_index format (format="band")
    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', '0 3 4 4'], #
        pseudos        = ['C.BFD.xml'],
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="36">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 3 4 4
       </occupation>
   </determinant>
   <determinant id="downdet" size="36">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()

    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)


    # test energy_index (format="energy")
    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', '-35 36'], #
        pseudos        = ['C.BFD.xml'],
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="36">
      <occupation mode="excited" spindataset="0" pairs="1" format="energy">             
-35 36
        </occupation>
   </determinant>
   <determinant id="downdet" size="36">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()

    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)

#end def test_excited_state



if versions.seekpath_available:
    def test_symbolic_excited_state():
        from nexus import generate_physical_system
        from nexus import generate_qmcpack_input

        dia = generate_physical_system(
            units     = 'A',
            axes      = [[ 1.785,  1.785,  0.   ],
                         [ 0.   ,  1.785,  1.785],
                         [ 1.785,  0.   ,  1.785]],
            elem      = ['C','C'],
            pos       = [[ 0.    ,  0.    ,  0.    ],
                         [ 0.8925,  0.8925,  0.8925]],
            use_prim  = True,    # Use SeeK-path library to identify prim cell
            tiling    = [2,1,2], 
            kgrid     = (1,1,1), 
            kshift    = (0,0,0), # Assumes we study transitions from Gamma. For non-gamma tilings, use kshift appropriately
            #C         = 4
            )

        qmc_optical = generate_qmcpack_input(
            det_format     = 'old',
            input_type     = 'basic',
            spin_polarized = True,
            system         = dia,
            excitation     = ['up', 'gamma vb x cb'], 
            jastrows       = [],
            qmc            = 'vmc',
            )

        expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 5 3 6
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
        text = qmc_optical.get('slaterdeterminant').write().strip()
        assert(text==expect)


        qmc_optical = generate_qmcpack_input(
            det_format     = 'old',
            input_type     = 'basic',
            spin_polarized = True,
            system         = dia,
            excitation     = ['up', 'gamma vb-1 x cb'], 
            jastrows       = [],
            qmc            = 'vmc',
            )

        expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 4 3 6
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
        text = qmc_optical.get('slaterdeterminant').write().strip()
        assert(text==expect)


        qmc_optical = generate_qmcpack_input(
            det_format     = 'old',
            input_type     = 'basic',
            spin_polarized = True,
            system         = dia,
            excitation     = ['up', 'gamma vb x cb+1'], 
            jastrows       = [],
            qmc            = 'vmc',
            )

        expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 5 3 7
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
        text = qmc_optical.get('slaterdeterminant').write().strip()
        assert(text==expect)

    #end def test_symbolic_excited_state
#end if

