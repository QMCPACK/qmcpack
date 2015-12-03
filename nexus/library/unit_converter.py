##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  unit_converter.py                                                 #
#    Support for physical unit conversion of scalars and arrays.     #
#    Several distance, time, mass, energy, charge, pressure, and     #
#    force units are supported.                                      #
#                                                                    #
#  Content summary:                                                  #
#    convert                                                         #
#      User-facing function to convert a scalar or array from one    #
#      unit system to another.                                       #
#                                                                    #
#    UnitConverter                                                   #
#      Class performs unit conversion.  Wrapped by convert function. #
#                                                                    #
#====================================================================#


from numpy import dot,zeros
from numpy.linalg import inv
from generic import obj
from developer import DevBase


class Unit(DevBase):

    unit_dicts = obj(all=obj())

    def __init__(self,type,name,symbol,value,shift=0):
        self.type   = type
        self.name   = name
        self.symbol = symbol
        self.value  = value
        self.shift  = shift

        ud = Unit.unit_dicts
        if type not in ud:
            ud[type] = obj()
        #end if
        td = ud[type]
        td[symbol] = self
        ud.all[symbol] = self
    #end def __init__
#end class Unit


class UnitConverter(DevBase):

    unassigned = None

    kb = 1.3806503e-23 #J/K

    # count
    mol = 6.0221415e23

    # distance
    m  = 1.e0
    A  = 1.e-10*m
    B  = .52917720859e-10*m
    nm = 1.e-9*m
    pm = 1.e-12*m
    fm = 1.e-15*m

    # time
    s = 1.e0
    ms = 1.e-3*s
    ns = 1.e-9*s
    ps = 1.e-12*s
    fs = 1.e-15*s

    # mass
    kg  = 1.e0
    me  = 9.10938291e-31*kg
    mp  = 1.672621777e-27*kg
    amu = 1.660538921e-27*kg
    Da  = amu

    # energy
    J      = 1.e0
    eV     = 1.60217646e-19*J
    Ry     = 13.6056923*eV
    Ha     = 2*Ry
    kJ_mol = 1000.*J/mol
    kcal_mol = 0.04336411531*eV
    K      = J/kb
    degC   = K
    degF   = 9./5.*K

    degC_shift = -273.15
    degF_shift = -459.67 

    # charge
    C = 1.e0
    e = 1.60217646e-19*C

    # pressure
    Pa   = 1.e0
    bar  = 1.e5*Pa
    Mbar = 1.e6*bar
    GPa  = 1.e9*Pa
    atm  = 1.01325e5*Pa

    # force
    N  = 1.e0
    pN = 1e-12*N

    # thermal conductivity
    W_mK = 1.0


    meter            = Unit('distance'  ,'meter'           ,'m'       ,m              )
    angstrom         = Unit('distance'  ,'Angstrom'        ,'A'       ,A              )
    bohr             = Unit('distance'  ,'Bohr'            ,'B'       ,B              )
    nanometer        = Unit('distance'  ,'nanometer'       ,'nm'      ,nm             )
    picometer        = Unit('distance'  ,'picometer'       ,'pm'      ,pm             )
    femtometer       = Unit('distance'  ,'femtometer'      ,'pm'      ,fm             )
                                                                      
    second           = Unit('time'      ,'second'          ,'s'       ,s              )
    millisecond      = Unit('time'      ,'millisecond'     ,'ms'      ,ms             )
    nanosecond       = Unit('time'      ,'nanosecond'      ,'ns'      ,ns             )
    picosecond       = Unit('time'      ,'picosecond'      ,'ps'      ,ps             )
    femtosecond      = Unit('time'      ,'femtosecond'     ,'fs'      ,fs             )
                                                                      
    kilogram         = Unit('mass'      ,'kilogram'        ,'kg'      ,kg             )
    electron_mass    = Unit('mass'      ,'electron mass'   ,'me'      ,me             )
    proton_mass      = Unit('mass'      ,'proton mass'     ,'mp'      ,mp             )
    atomic_mass_unit = Unit('mass'      ,'atomic mass unit','amu'     ,amu            )
    dalton           = Unit('mass'      ,'Dalton'          ,'Da'      ,Da             )
        
    joule            = Unit('energy'    ,'Joule'           ,'J'       ,J              )
    electron_volt    = Unit('energy'    ,'electron Volt'   ,'eV'      ,eV             )
    rydberg          = Unit('energy'    ,'Rydberg'         ,'Ry'      ,Ry             )
    hartree          = Unit('energy'    ,'Hartree'         ,'Ha'      ,Ha             )
    kJ_mole          = Unit('energy'    ,'kJ_mole'         ,'kJ_mol'  ,kJ_mol         )
    kcal_mole        = Unit('energy'    ,'kcal_mole'       ,'kcal_mol',kcal_mol       )
    kelvin           = Unit('energy'    ,'Kelvin'          ,'K'       ,K              )
    celcius          = Unit('energy'    ,'Celcius'         ,'degC'    ,degC,degC_shift)
    fahrenheit       = Unit('energy'    ,'Fahrenheit'      ,'degF'    ,degF,degF_shift)
                                                           
    coulomb          = Unit('charge'    ,'Coulomb'         ,'C'       ,C              )
    electron_charge  = Unit('charge'    ,'electron charge' ,'e'       ,e              )

    pascal           = Unit('pressure'  ,'Pascal'          ,'Pa'      ,Pa             )
    bar              = Unit('pressure'  ,'bar'             ,'bar'     ,bar            )
    megabar          = Unit('pressure'  ,'megabar'         ,'Mbar'    ,Mbar           )
    gigapascal       = Unit('pressure'  ,'gigaPascal'      ,'Gpa'     ,GPa            )
    atmosphere       = Unit('pressure'  ,'atmosphere'      ,'atm'     ,atm            )
                                                                      
    newton           = Unit('force'     ,'Newton'          ,'N'       ,N              )
    piconewton       = Unit('force'     ,'picoNewton'      ,'pN'      ,pN             )
                                                                      
    W_per_mK         = Unit('therm_cond','W/(m*K)'         ,'W_mK'    ,W_mK           )


    distance_dict   = Unit.unit_dicts.distance
    mass_dict       = Unit.unit_dicts.mass
    energy_dict     = Unit.unit_dicts.energy
    charge_dict     = Unit.unit_dicts.charge
    pressure_dict   = Unit.unit_dicts.pressure
    force_dict      = Unit.unit_dicts.force
    therm_cond_dict = Unit.unit_dicts.therm_cond

    unit_dict = Unit.unit_dicts.all


    def __init(self):
        self.error('UnitConverter should not be instantiated')
    #def __init__

    @staticmethod
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

    @staticmethod
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

    
    @staticmethod
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
                A['lattice'] = dot(arr,inv(vectors.A[units]))
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
#end class UnitConverter


def convert(value,source_unit,target_unit):
    return UnitConverter.convert(value,source_unit,target_unit)[0]
#end def convert

