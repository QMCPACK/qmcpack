
from abilities import Callable


class Unit:
    def __init__(self,name,symbol,type,value,shift=0):
        self.name = name
        self.symbol = symbol
        self.type = type
        self.value = value
        self.shift = shift
    #end def __init__
#end class Unit


class UnitConverter:

    unassigned = None


    kb = 1.3806503e-23 #J/K

    count_set = set(['mol'])
    mol = 6.0221415e23

    distance_set = set(['m','A','B','nm','pm','fm','a','b','c'])
    m  = 1.e0
    A  = 1.e-10*m
    B  = .52917720859e-10*m
    nm = 1.e-9*m
    pm = 1.e-12*m
    fm = 1.e-15*m

    time_set = set(['s','ms','ns','ps','fs'])
    s = 1.e0
    ms = 1.e-3*s
    ns = 1.e-9*s
    ps = 1.e-12*s
    fs = 1.e-15*s

    mass_set = set(['kg','me','mp','amu','Da'])
    kg  = 1.e0
    me  = 9.10938291e-31*kg
    mp  = 1.672621777e-27*kg
    amu = 1.660538921e-27*kg
    Da  = amu

    energy_set = set(['J','eV','Ry','Ha','kJ_mol','K','degC','degF'])
    J      = 1.e0
    eV     = 1.60217646e-19*J
    Ry     = 13.6056923*eV
    Ha     = 2*Ry
    kJ_mol = 1000.*J/mol
    K      = J/kb
    degC   = K
    degF   = 9./5.*K

    degC_shift = -273.15
    degF_shift = -459.67 

    charge_set = set(['C','e'])
    C = 1.e0
    e = 1.60217646e-19*C

    pressure_set = set(['Pa','bar','Mbar','GPa','atm'])
    Pa   = 1.e0
    bar  = 1.e5*Pa
    Mbar = 1.e6*bar
    GPa  = 1.e9*Pa
    atm  = 1.01325e5*Pa

    force_set = set(['N','pN'])
    N  = 1.e0
    pN = 1e-12*N

    therm_cond_set = set(['W_mK'])
    W_mK = 1.0

    alatt = unassigned
    blatt = unassigned
    clatt = unassigned

    meter      = Unit('meter'     ,'m' ,'distance',m)
    angstrom   = Unit('Angstrom'  ,'A' ,'distance',A)
    bohr       = Unit('Bohr'      ,'B' ,'distance',B)
    nanometer  = Unit('nanometer' ,'nm','distance',nm)
    picometer  = Unit('picometer' ,'pm','distance',pm)
    femtometer = Unit('femtometer','pm','distance',fm)
    a          = Unit('a'         ,'a' ,'distance',alatt)
    b          = Unit('b'         ,'b' ,'distance',blatt)
    c          = Unit('c'         ,'c' ,'distance',clatt)

    second      = Unit('second'     ,'s' ,'time',s )
    millisecond = Unit('millisecond','ms','time',ms)
    nanosecond  = Unit('nanosecond' ,'ns','time',ns)
    picosecond  = Unit('picosecond' ,'ps','time',ps)
    femtosecond = Unit('femtosecond','fs','time',fs)

    kilogram         = Unit('kilogram'        ,'kg' ,'mass',kg)
    electron_mass    = Unit('electron mass'   ,'me' ,'mass',me)
    proton_mass      = Unit('proton mass'     ,'mp' ,'mass',mp)
    atomic_mass_unit = Unit('atomic mass unit','amu','mass',amu)
    dalton           = Unit('Dalton'          ,'Da' ,'mass',Da)
        
    joule         = Unit('Joule'        ,'J'     ,'energy',J)
    electron_volt = Unit('electron Volt','eV'    ,'energy',eV)
    rydberg       = Unit('Rydberg'      ,'Ry'    ,'energy',Ry)
    hartree       = Unit('Hartree'      ,'Ha'    ,'energy',Ha)
    kJ_mole       = Unit('kJ_mole'      ,'kJ_mol','energy',kJ_mol)
    kelvin        = Unit('Kelvin'       ,'K'     ,'energy',K)
    celcius       = Unit('Celcius'      ,'degC'  ,'energy',degC,degC_shift)
    fahrenheit    = Unit('Fahrenheit'   ,'degF'  ,'energy',degF,degF_shift)

    coulomb         = Unit('Coulomb'        ,'C','charge',C)
    electron_charge = Unit('electron charge','e','charge',e)

    pascal     = Unit('Pascal'    ,'Pa'  ,'pressure',Pa)
    bar        = Unit('bar'       ,'bar' ,'pressure',bar)
    megabar    = Unit('megabar'   ,'Mbar','pressure',Mbar)
    gigapascal = Unit('gigaPascal','Gpa' ,'pressure',GPa)
    atmosphere = Unit('atmosphere','atm' ,'pressure',atm)

    newton     = Unit('Newton'    ,'N' ,'force',N)
    piconewton = Unit('picoNewton','pN','force',pN)

    W_per_mK = Unit('W/(m*K)','W_mK','thermal_cond',W_mK)


    distance_dict = dict([('A',angstrom),\
                          ('B',bohr),\
                          ('a',a),\
                          ('b',b),\
                          ('c',c),\
                          ('m',meter),\
                          ('nm',nanometer),\
                          ('pm',picometer),\
                          ('fm',femtometer),\
                          ])

    mass_dict = dict([    ('kg',kilogram),\
                          ('me',electron_mass),\
                          ('mp',proton_mass),\
                          ('amu',atomic_mass_unit),\
                          ('Da',dalton),\
                          ])
    energy_dict = dict([\
                          ('J',joule),\
                          ('eV',electron_volt),\
                          ('Ry',rydberg),\
                          ('Ha',hartree),\
                          ('kJ_mol',kJ_mole),\
                          ('K',kelvin),\
                          ('degC',celcius),\
                          ('degF',fahrenheit),\
                          ])

    charge_dict = dict([\
                          ('C',coulomb),\
                          ('e',electron_charge),\
                          ])

    pressure_dict = dict([\
                          ('Pa',pascal),\
                          ('bar',bar),\
                          ('Mbar',megabar),\
                          ('GPa',gigapascal),\
                          ('atm',atmosphere),\
                              ])

    force_dict = dict([\
                          ('N',newton),\
                          ('pN',piconewton),\
                              ])

    therm_cond_dict = dict([\
                          ('W_mK',W_per_mK),\
                          ])


    unit_dict = dict([    ('A',angstrom),\
                          ('B',bohr),\
                          ('a',a),\
                          ('b',b),\
                          ('c',c),\
                          ('m',meter),\
                          ('nm',nanometer),\
                          ('pm',picometer),\
                          ('fm',femtometer),\
                          ('kg',kilogram),\
                          ('me',electron_mass),\
                          ('mp',proton_mass),\
                          ('amu',atomic_mass_unit),\
                          ('J',joule),\
                          ('eV',electron_volt),\
                          ('Ry',rydberg),\
                          ('Ha',hartree),\
                          ('kJ_mol',kJ_mole),\
                          ('K',kelvin),\
                          ('degC',celcius),\
                          ('degF',fahrenheit),\
                          ('C',coulomb),\
                          ('e',electron_charge),\
                          ('Pa',pascal),\
                          ('bar',bar),\
                          ('Mbar',megabar),\
                          ('GPa',gigapascal),\
                          ('atm',atmosphere),\
                          ('N',newton),\
                          ('pN',piconewton),\
                          ('W_mK',W_per_mK),\
                          ])


    def __init(self):
        None
    #def __init__

    def convert(value,source_unit,target_unit):
        ui  = UnitConverter.unit_dict[source_unit]
        uo = UnitConverter.unit_dict[target_unit]

        if(ui.type != uo.type):
            print 'ERROR: in UnitConverter.convert()'
            print '   type conversion attempted between'
            print '   '+ui.type+' and '+uo.type
        else:
            value_out = (value*ui.value+ui.shift-uo.shift)/uo.value
        #end if

        return (value_out,target_unit)

    #end def convert
    convert = Callable(convert)


    def convert_scalar_to_all(units,value_orig):
        unit_type = UnitConverter.unit_dict[units].type

        value = dict()
        value['orig'] = value_orig

        for k,u in UnitConverter.unit_dict.iteritems():
            if(u.type == unit_type and u.value!=UnitConverter.unassigned):
                (value[k],utmp) = UnitConverter.convert(value['orig'],units,k)
            #end if
        #end for

        return value
    #end def convert_scalar_to_all
    convert_scalar_to_all = Callable(convert_scalar_to_all)

    
    def convert_array_to_all(units,arr,vectors=None):

        A = dict()
        A['orig'] = arr.copy()

        (n,m) = arr.shape

        if(units=='lattice'):
            if(vectors!=None):
                A['lattice'] = arr
            #end if
            def_unit = 'A'
            arr_use = dot(arr,vectors.A[def_unit])
            ui = UnitConverter.unit_dict[def_unit]
        else:
            arr_use = arr
            ui = UnitConverter.unit_dict[units]
            if(vectors!=None):
                A['lattice'] = dot(arr,linalg.inv(vectors.A[units]))
            #end if
        #end if
            
        At = zeros((n,m));
        for k,u in UnitConverter.unit_dict.iteritems():
            if(u.type == ui.type and u.value!=UnitConverter.unassigned):
                for i in range(n):
                    for j in range(m):
                        At[i,j] = (arr_use[i,j]*ui.value+ui.shift-u.shift)/u.value
                    #end for
                #end for
                A[k] = At.copy()
            #end if
        #end for

        return A
    #end def convert_array_to_all
    convert_array_to_all = Callable(convert_array_to_all)

    def submit_unit(uo):
        UnitConverter.unit_dict[uo.symbol] = uo
    #end def submit_unit
    submit_unit = Callable(submit_unit)
#end class UnitConverter


def convert(value,source_unit,target_unit):
    return UnitConverter.convert(value,source_unit,target_unit)[0]
#end def convert

