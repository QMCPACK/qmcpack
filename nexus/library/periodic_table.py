##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  periodic_table.py                                                 #
#    Collected periodic table data.                                  #
#                                                                    #
#  Content summary:                                                  #
#    PeriodicTable                                                   #
#      Class representing the periodic table.                        #
#                                                                    #
#    periodic_table or pt                                            #
#      References to a single, complete PeriodicTable instance.      #
#      Used extensively by importing modules for data access.        #
#                                                                    #
#    Element                                                         #
#      Class representing a single element.                          #
#                                                                    #                                        
#====================================================================#


from generic import obj
from developer import DevBase
from unit_converter import UnitConverter


def phys_value_dict(value=None,units=None):    
    vdict = UnitConverter.convert_scalar_to_all(units,value)
    return obj(**vdict)
#end def phys_value_dict




class SimpleElement(DevBase):
    def __init__(self):

        self.atomic_number     = None
        self.name              = None
        self.symbol            = None
        self.group             = None
        self.atomic_weight     = None
        self.atomic_radius     = None 
        self.nuclear_charge    = None 
        self.abundance         = None 
        self.electron_affinity = None 
        self.electronegativity = None 
        self.ionization_energy = None 
        self.ionic_radius      = None 
        self.melting_point     = None 
        self.boiling_point     = None 
        #self. = None
        #self. = None
        #self. = None
        #self. = None

        self.string_rep = None
        self.var_dict = None
    #end def __init__

    def create_var_dict(self):
        self.var_dict = dict()
        self.var_dict['atomic_number'    ] = self.atomic_number    
        self.var_dict['name'             ] = self.name             
        self.var_dict['symbol'           ] = self.symbol           
        self.var_dict['group'            ] = self.group            
        self.var_dict['atomic_weight'    ] = self.atomic_weight    
        self.var_dict['atomic_radius'    ] = self.atomic_radius     
        self.var_dict['nuclear_charge'   ] = self.nuclear_charge    
        self.var_dict['abundance'        ] = self.abundance         
        self.var_dict['electron_affinity'] = self.electron_affinity 
        self.var_dict['electronegativity'] = self.electronegativity 
        self.var_dict['ionization_energy'] = self.ionization_energy 
        self.var_dict['ionic_radius'     ] = self.ionic_radius      
        self.var_dict['melting_point'    ] = self.melting_point     
        self.var_dict['boiling_point'    ] = self.boiling_point     

        self.replace_None()
    #end def create_var_dict

    def replace_None(self):
        none_rep = -1.0
        for k,v in self.var_dict.iteritems():
            if(v==None):
                self.var_dict[k] = none_rep
            #end if
        #end for
        self.atomic_number      = self.var_dict['atomic_number'    ]
        self.name               = self.var_dict['name'             ]
        self.symbol             = self.var_dict['symbol'           ]
        self.group              = self.var_dict['group'            ]
        self.atomic_weight      = self.var_dict['atomic_weight'    ]
        self.atomic_radius      = self.var_dict['atomic_radius'    ]
        self.nuclear_charge     = self.var_dict['nuclear_charge'   ]
        self.abundance          = self.var_dict['abundance'        ]
        self.electron_affinity  = self.var_dict['electron_affinity']
        self.electronegativity  = self.var_dict['electronegativity']
        self.ionization_energy  = self.var_dict['ionization_energy']
        self.ionic_radius       = self.var_dict['ionic_radius'     ]
        self.melting_point      = self.var_dict['melting_point'    ]
        self.boiling_point      = self.var_dict['boiling_point'    ]  
    #end def replace_None

    def create_string_representation(self):
        ind = 4*' '
        
        iformat = '%6i'
        rformat = '%7.5f'

        s = ''
        s += self.symbol+'{\n'
        s += ind + 'atomic_number = '     + str(self.atomic_number)+'\n'
        s += ind + 'name = '              + str(self.name)+'\n'
        s += ind + 'symbol = '            + str(self.symbol)+'\n'
        s += ind + 'group = '             + str(self.group)+'\n'
        s += ind + 'atomic_weight = '     + str(self.atomic_weight)+'\n'
        s += ind + 'atomic_radius = '     + str(self.atomic_radius)+'\n'
        s += ind + 'nuclear_charge = '    + str(self.nuclear_charge)+'\n'
        s += ind + 'abundance = '         + str(self.abundance)+'\n'
        s += ind + 'electron_affinity = ' + str(self.electron_affinity)+'\n'
        s += ind + 'electronegativity = ' + str(self.electronegativity)+'\n'
        s += ind + 'ionization_energy = ' + str(self.ionization_energy)+'\n'
        s += ind + 'ionic_radius = '      + str(self.ionic_radius)+'\n'
        s += ind + 'melting_point = '     + str(self.melting_point)+'\n'
        s += ind + 'boiling_point = '     + str(self.boiling_point)+'\n'
        s += '}\n'

        self.string_rep = s

    #end def create_string_representation
#end class SimpleElement


class Element(SimpleElement):
    def __init__(self,se):
        SimpleElement.__init__(self)
        
        awu = PeriodicTable.atomic_weight_unit      
        aru = PeriodicTable.atomic_radius_unit      
        ncu = PeriodicTable.nuclear_charge_unit     
        eau = PeriodicTable.electron_affinity_unit  
        ieu = PeriodicTable.ionization_energy_units 
        iru = PeriodicTable.ionic_radius_units      
        tcu = PeriodicTable.thermal_cond_units      
        mpu = PeriodicTable.melting_point_units     
        bpu = PeriodicTable.boiling_point_units     

        self.atomic_number = se.atomic_number
        self.name          = se.name
        self.symbol        = se.symbol
        self.group         = PeriodicTable.group_dict[se.group]
        self.abundance     = se.abundance

        self.atomic_weight     = phys_value_dict(se.atomic_weight    , awu)
        self.atomic_radius     = phys_value_dict(se.atomic_radius    , aru)
        self.nuclear_charge    = phys_value_dict(se.nuclear_charge   , ncu)
        self.electron_affinity = phys_value_dict(se.electron_affinity, eau)
        self.ionization_energy = phys_value_dict(se.ionization_energy, ieu)
        self.ionic_radius      = phys_value_dict(se.ionic_radius     , iru)
        self.thermal_cond      = phys_value_dict(se.thermal_cond     , tcu)
        self.melting_point     = phys_value_dict(se.melting_point    , mpu)
        self.boiling_point     = phys_value_dict(se.boiling_point    , bpu)
        
    #end def __init__
#end class Element



class PeriodicTable(DevBase):

    element_set=set([\
       'Ac','Al','Am','Sb','Ar','As','At','Ba','Bk','Be','Bi','B' ,'Br',\
       'Cd','Ca','Cf','C' ,'Ce','Cs','Cl','Cr','Co','Cu','Cm','Dy','Es',\
       'Er','Eu','Fm','F' ,'Fr','Gd','Ga','Ge','Au','Hf','Ha','Hs','He',\
       'Ho','H' ,'In','I' ,'Ir','Fe','Kr','La','Lr','Pb','Li','Lu','Mg',\
       'Mn','Mt','Md','Hg','Mo','Ns','Nd','Ne','Np','Ni','Nb','N' ,'No',\
       'Os','O' ,'Pd','P' ,'Pt','Pu','Po','K' ,'Pr','Pm','Pa','Ra','Rn',\
       'Re','Rh','Rb','Ru','Rf','Sm','Sc','Sg','Se','Si','Ag','Na','Sr',\
       'S' ,'Ta','Tc','Te','Tb','Tl','Th','Tm','Sn','Ti','W' ,'U' ,'V' ,\
       'Xe','Yb','Y' ,'Zn','Zr'])


    element_dict=dict({\
	'Ac':'Actinium',\
	'Al':'Aluminum',\
	'Am':'Americium',\
	'Sb':'Antimony',\
	'Ar':'Argon',\
	'As':'Arsenic',\
	'At':'Astatine',\
	'Ba':'Barium',\
	'Bk':'Berkelium',\
	'Be':'Beryllium',\
	'Bi':'Bismuth',\
	'B':'Boron',\
	'Br':'Bromine',\
	'Cd':'Cadmium',\
	'Ca':'Calcium',\
	'Cf':'Californium',\
	'C' :'Carbon',\
	'Ce':'Cerium',\
	'Cs':'Cesium',\
	'Cl':'Chlorine',\
	'Cr':'Chromium',\
	'Co':'Cobalt',\
	'Cu':'Copper',\
	'Cm':'Curium',\
	'Dy':'Dysprosium',\
	'Es':'Einsteinium',\
	'Er':'Erbium',\
	'Eu':'Europium',\
	'Fm':'Fermium',\
	'F' :'Flourine',\
	'Fr':'Francium',\
	'Gd':'Gadolinium',\
	'Ga':'Gallium',\
	'Ge':'Germanium',\
	'Au':'Gold',\
	'Hf':'Hafnium',\
	'Ha':'Hahnium',\
	'Hs':'Hassium',\
	'He':'Helium',\
	'Ho':'Holmium',\
	'H' :'Hydrogen',\
	'In':'Indium',\
	'I' :'Iodine',\
	'Ir':'Iridium',\
	'Fe':'Iron',\
	'Kr':'Krypton',\
	'La':'Lanthanum',\
	'Lr':'Lawrencium',\
	'Pb':'Lead',\
	'Li':'Lithium',\
	'Lu':'Lutetium',\
	'Mg':'Magnesium',\
	'Mn':'Manganese',\
	'Mt':'Meitnerium',\
	'Md':'Mendelevium',\
	'Hg':'Mercury',\
	'Mo':'Molybdenum',\
	'Ns':'Neilsborium',\
	'Nd':'Neodymium',\
	'Ne':'Neon',\
	'Np':'Neptunium',\
	'Ni':'Nickel',\
	'Nb':'Niobium',\
	'N' :'Nitrogen',\
	'No':'Nobelium',\
	'Os':'Osmium',\
	'O' :'Oxygen',\
	'Pd':'Palladium',\
	'P' :'Phosphorus',\
	'Pt':'Platinum',\
	'Pu':'Plutonium',\
	'Po':'Polonium',\
	'K' :'Potassium',\
	'Pr':'Praseodymium',\
	'Pm':'Promethium',\
	'Pa':'Protactinium',\
	'Ra':'Radium',\
	'Rn':'Radon',\
	'Re':'Rhenium',\
	'Rh':'Rhodium',\
	'Rb':'Rubidium',\
	'Ru':'Ruthenium',\
	'Rf':'Rutherfordium',\
	'Sm':'Samarium',\
	'Sc':'Scandium',\
	'Sg':'Seaborgium',\
	'Se':'Selenium',\
	'Si':'Silicon',\
	'Ag':'Silver',\
	'Na':'Sodium',\
	'Sr':'Strontium',\
	'S' :'Sulfur',\
	'Ta':'Tantalum',\
	'Tc':'Technetium',\
	'Te':'Tellurium',\
	'Tb':'Terbium',\
	'Tl':'Thalium',\
	'Th':'Thorium',\
	'Tm':'Thulium',\
	'Sn':'Tin',\
	'Ti':'Titanium',\
	'W' :'Tungsten',\
	'U' :'Uranium',\
	'V' :'Vanadium',\
	'Xe':'Xenon',\
	'Yb':'Ytterbium',\
	'Y' :'Yttrium',\
	'Zn':'Zinc',\
	'Zr':'Zirconium',\
       })

    group_dict = dict([\
                (0 ,'LanAct'),\
                (1 ,'IA'),\
                (2 ,'IIA'),\
                (3 ,'IIIB'),\
                (4 ,'IVB'),\
                (5 ,'VB'),\
                (6 ,'VIB'),\
                (7 ,'VIIB'),\
                (8 ,'VII'),\
                (9 ,'VII'),\
                (10,'VII'),\
                (11,'IB'),\
                (12,'IIB'),\
                (13,'IIIA'),\
                (14,'IVA'),\
                (15,'VA'),\
                (16,'VIA'),\
                (17,'VIIA'),\
                (18,'0')\
            ])


    atomic_weight_unit      = 'amu'
    atomic_radius_unit      = 'pm'
    nuclear_charge_unit     = 'e'
    electron_affinity_unit  = 'kJ_mol'
    ionization_energy_units = 'eV'
    ionic_radius_units      = 'pm'
    thermal_cond_units      = 'W_mK'
    melting_point_units     = 'degC'
    boiling_point_units     = 'degC'

    def __init__(self):
        self.nelements = None
        self.elements = None

        nelements  = 103
        e = obj()
        for i in range(1,nelements+1):
            e[i] = SimpleElement()
        #end for

        for i in range(1,nelements+1):
            e[i].atomic_number = i
        #end for    


	e[1].symbol='H'
	e[2].symbol='He'
	e[3].symbol='Li'
	e[4].symbol='Be'
	e[5].symbol='B'
	e[6].symbol='C'
	e[7].symbol='N'
	e[8].symbol='O'
	e[9].symbol='F'
	e[10].symbol='Ne'
	e[11].symbol='Na'
	e[12].symbol='Mg'
	e[13].symbol='Al'
	e[14].symbol='Si'
	e[15].symbol='P'
	e[16].symbol='S'
	e[17].symbol='Cl'
	e[18].symbol='Ar'
	e[19].symbol='K'
	e[20].symbol='Ca'
	e[21].symbol='Sc'
	e[22].symbol='Ti'
	e[23].symbol='V'
	e[24].symbol='Cr'
	e[25].symbol='Mn'
	e[26].symbol='Fe'
	e[27].symbol='Co'
	e[28].symbol='Ni'
	e[29].symbol='Cu'
	e[30].symbol='Zn'
	e[31].symbol='Ga'
	e[32].symbol='Ge'
	e[33].symbol='As'
	e[34].symbol='Se'
	e[35].symbol='Br'
	e[36].symbol='Kr'
	e[37].symbol='Rb'
	e[38].symbol='Sr'
	e[39].symbol='Y'
	e[40].symbol='Zr'
	e[41].symbol='Nb'
	e[42].symbol='Mo'
	e[43].symbol='Tc'
	e[44].symbol='Ru'
	e[45].symbol='Rh'
	e[46].symbol='Pd'
	e[47].symbol='Ag'
	e[48].symbol='Cd'
	e[49].symbol='In'
	e[50].symbol='Sn'
	e[51].symbol='Sb'
	e[52].symbol='Te'
	e[53].symbol='I'
	e[54].symbol='Xe'
	e[55].symbol='Cs'
	e[56].symbol='Ba'
	e[57].symbol='La'
	e[58].symbol='Ce'
	e[59].symbol='Pr'
	e[60].symbol='Nd'
	e[61].symbol='Pm'
	e[62].symbol='Sm'
	e[63].symbol='Eu'
	e[64].symbol='Gd'
	e[65].symbol='Tb'
	e[66].symbol='Dy'
	e[67].symbol='Ho'
	e[68].symbol='Er'
	e[69].symbol='Tm'
	e[70].symbol='Yb'
	e[71].symbol='Lu'
	e[72].symbol='Hf'
	e[73].symbol='Ta'
	e[74].symbol='W'
	e[75].symbol='Re'
	e[76].symbol='Os'
	e[77].symbol='Ir'
	e[78].symbol='Pt'
	e[79].symbol='Au'
	e[80].symbol='Hg'
	e[81].symbol='Tl'
	e[82].symbol='Pb'
	e[83].symbol='Bi'
	e[84].symbol='Po'
	e[85].symbol='At'
	e[86].symbol='Rn'
	e[87].symbol='Fr'
	e[88].symbol='Ra'
	e[89].symbol='Ac'
	e[90].symbol='Th'
	e[91].symbol='Pa'
	e[92].symbol='U'
	e[93].symbol='Np'
	e[94].symbol='Pu'
	e[95].symbol='Am'
	e[96].symbol='Cm'
	e[97].symbol='Bk'
	e[98].symbol='Cf'
	e[99].symbol='Es'
	e[100].symbol='Fm'
	e[101].symbol='Md'
	e[102].symbol='No'
	e[103].symbol='Lr'

        for i in range(1,len(e)):
            e[i].name = PeriodicTable.element_dict[e[i].symbol] 
        #end for

	e[1].group = 1
	e[2].group = 18
	e[3].group = 1
	e[4].group = 2
	e[5].group = 13
	e[6].group = 14
	e[7].group = 15
	e[8].group = 16
	e[9].group = 17
	e[10].group = 18
	e[11].group = 1
	e[12].group = 2
	e[13].group = 13
	e[14].group = 14
	e[15].group = 15
	e[16].group = 16
	e[17].group = 17
	e[18].group = 18
	e[19].group = 1
	e[20].group = 2
	e[21].group = 3
	e[22].group = 4
	e[23].group = 5
	e[24].group = 6
	e[25].group = 7
	e[26].group = 8
	e[27].group = 9
	e[28].group = 10
	e[29].group = 11
	e[30].group = 12
	e[31].group = 13
	e[32].group = 14
	e[33].group = 15
	e[34].group = 16
	e[35].group = 17
	e[36].group = 18
	e[37].group = 1
	e[38].group = 2
	e[39].group = 3
	e[40].group = 4
	e[41].group = 5
	e[42].group = 6
	e[43].group = 7
	e[44].group = 8
	e[45].group = 9
	e[46].group = 10
	e[47].group = 11
	e[48].group = 12
	e[49].group = 13
	e[50].group = 14
	e[51].group = 15
	e[52].group = 16
	e[53].group = 17
	e[54].group = 18
	e[55].group = 1
	e[56].group = 2
	e[57].group = 3
	e[58].group = 0
	e[59].group = 0
	e[60].group = 0
	e[61].group = 0
	e[62].group = 0
	e[63].group = 0
	e[64].group = 0
	e[65].group = 0
	e[66].group = 0
	e[67].group = 0
	e[68].group = 0
	e[69].group = 0
	e[70].group = 0
	e[71].group = 0
	e[72].group = 4
	e[73].group = 5
	e[74].group = 6
	e[75].group = 7
	e[76].group = 8
	e[77].group = 9
	e[78].group = 10
	e[79].group = 11
	e[80].group = 12
	e[81].group = 13
	e[82].group = 14
	e[83].group = 15
	e[84].group = 16
	e[85].group = 17
	e[86].group = 18
	e[87].group =1
	e[88].group = 2
	e[89].group = 3
	e[90].group = 0
	e[91].group = 0
	e[92].group = 0
	e[93].group = 0
	e[94].group = 0
	e[95].group = 0
	e[96].group = 0
	e[97].group = 0
	e[98].group = 0
	e[99].group = 0
	e[100].group = 0
	e[101].group = 0
	e[102].group = 0
	e[103].group = 0

	e[1].atomic_weight = 1.00794
	e[2].atomic_weight = 4.002602
	e[3].atomic_weight = 6.941
	e[4].atomic_weight = 9.0122
	e[5].atomic_weight = 10.811
	e[6].atomic_weight = 12.011000
	e[7].atomic_weight = 14.007
	e[8].atomic_weight = 15.999
	e[9].atomic_weight = 18.998
	e[10].atomic_weight = 20.180
	e[11].atomic_weight = 22.990
	e[12].atomic_weight = 24.305
	e[13].atomic_weight = 26.982
	e[14].atomic_weight = 28.086
	e[15].atomic_weight = 30.974
	e[16].atomic_weight = 32.064
	e[17].atomic_weight = 35.453
	e[18].atomic_weight = 39.948
	e[19].atomic_weight = 39.098
	e[20].atomic_weight = 40.08
	e[21].atomic_weight = 44.956
	e[22].atomic_weight = 47.90
	e[23].atomic_weight = 50.942
	e[24].atomic_weight = 51.996
	e[25].atomic_weight = 54.938
	e[26].atomic_weight = 55.845
	e[27].atomic_weight = 58.933
	e[28].atomic_weight = 58.69
	e[29].atomic_weight = 63.546
	e[30].atomic_weight = 65.38
	e[31].atomic_weight = 65.38
	e[32].atomic_weight = 72.61
	e[33].atomic_weight = 74.992
	e[34].atomic_weight = 78.96
	e[35].atomic_weight = 79.904
	e[36].atomic_weight = 83.80
	e[37].atomic_weight = 85.47
	e[38].atomic_weight = 87.956
	e[39].atomic_weight = 88.905
	e[40].atomic_weight = 91.22
	e[41].atomic_weight = 92.906
	e[42].atomic_weight = 95.94
	e[43].atomic_weight = 98.00
	e[44].atomic_weight = 101.07
	e[45].atomic_weight = 102.91
	e[46].atomic_weight = 106.42
	e[47].atomic_weight = 107.87
	e[48].atomic_weight = 112.41
	e[49].atomic_weight = 114.82
	e[50].atomic_weight = 118.69
	e[51].atomic_weight = 121.175
	e[52].atomic_weight = 127.60
	e[53].atomic_weight = 126.90
	e[54].atomic_weight = 131.29
	e[55].atomic_weight = 132.91
	e[56].atomic_weight = 137.33
	e[57].atomic_weight = 138.92
	e[58].atomic_weight = 140.12
	e[59].atomic_weight = 140.91
	e[60].atomic_weight = 144.24
	e[61].atomic_weight = 145.00
	e[62].atomic_weight = 150.36
	e[63].atomic_weight = 151.97
	e[64].atomic_weight = 157.25
	e[65].atomic_weight = 158.924
	e[66].atomic_weight = 162.5
	e[67].atomic_weight = 164.930
	e[68].atomic_weight = 167.26
	e[69].atomic_weight = 169.934
	e[70].atomic_weight = 173.04
	e[71].atomic_weight = 174.97
	e[72].atomic_weight = 178.49
	e[73].atomic_weight = 180.948
	e[74].atomic_weight = 183.85
	e[75].atomic_weight = 186.2
	e[76].atomic_weight = 190.2
	e[77].atomic_weight = 192.2
	e[78].atomic_weight = 195.09
	e[79].atomic_weight = 196.197
	e[80].atomic_weight = 200.59
	e[81].atomic_weight = 204.37
	e[82].atomic_weight = 207.19
	e[83].atomic_weight = 208.980
	e[84].atomic_weight = 209.0
	e[85].atomic_weight = 210.0
	e[86].atomic_weight = 222.0
	e[87].atomic_weight = 223.0
	e[88].atomic_weight = 226.0
	e[89].atomic_weight = 227.028
	e[90].atomic_weight = 204.37
	e[91].atomic_weight = 231.0
	e[92].atomic_weight = 238.03
	e[93].atomic_weight = 237.05
	e[94].atomic_weight = 244.0
	e[95].atomic_weight = 243.0
	e[96].atomic_weight = 245.0
	e[97].atomic_weight = 247.0
	e[98].atomic_weight = 249.0
	e[99].atomic_weight = 254.0
	e[100].atomic_weight = 252.0
	e[101].atomic_weight = 256.0
	e[102].atomic_weight = 254.0
	e[103].atomic_weight = 257

        #atomic radius (in picometers)
	e[1].atomic_radius = 78.000000
	e[2].atomic_radius = 128.000000
	e[3].atomic_radius = 152.000000
	e[4].atomic_radius = 111.300000
	e[5].atomic_radius = 79.500000
	e[6].atomic_radius = 77.200000
	e[7].atomic_radius = 54.900000
	e[8].atomic_radius = 60.400000
	e[9].atomic_radius = 70.900000
	e[10].atomic_radius = 0.000000
	e[11].atomic_radius = 185.800000
	e[12].atomic_radius = 159.900000
	e[13].atomic_radius = 143.200000
	e[14].atomic_radius = 117.600000
	e[15].atomic_radius = 110.500000
	e[16].atomic_radius = 103.500000
	e[17].atomic_radius = 99.400000
	e[18].atomic_radius = 180.000000
	e[19].atomic_radius = 227.200000
	e[20].atomic_radius = 197.400000
	e[21].atomic_radius = 160.600000
	e[22].atomic_radius = 144.800000
	e[23].atomic_radius = 131.100000
	e[24].atomic_radius = 124.900000
	e[25].atomic_radius = 136.700000
	e[26].atomic_radius = 124.100000
	e[27].atomic_radius = 125.300000
	e[28].atomic_radius = 124.600000
	e[29].atomic_radius = 127.800000
	e[30].atomic_radius = 133.500000
	e[31].atomic_radius = 122.100000
	e[32].atomic_radius = 122.500000
	e[33].atomic_radius = 124.500000
	e[34].atomic_radius = 116.000000
	e[35].atomic_radius = 114.500000
	e[36].atomic_radius = 0.000000
	e[37].atomic_radius = 247.500000
	e[38].atomic_radius = 215.100000
	e[39].atomic_radius = 177.600000
	e[40].atomic_radius = 159.000000
	e[41].atomic_radius = 142.900000
	e[42].atomic_radius = 136.300000
	e[43].atomic_radius = 135.200000
	e[44].atomic_radius = 132.500000
	e[45].atomic_radius = 134.500000
	e[46].atomic_radius = 137.600000
	e[47].atomic_radius = 144.500000
	e[48].atomic_radius = 148.900000
	e[49].atomic_radius = 162.600000
	e[50].atomic_radius = 140.500000
	e[51].atomic_radius = 145.000000
	e[52].atomic_radius = 143.200000
	e[53].atomic_radius = 133.100000
	e[54].atomic_radius = 210.000000
	e[55].atomic_radius = 265.500000
	e[56].atomic_radius = 217.400000
	e[57].atomic_radius = 187.000000
	e[58].atomic_radius = 182.500000
	e[59].atomic_radius = 182.000000
	e[60].atomic_radius = 181.400000
	e[61].atomic_radius = 181.000000
	e[62].atomic_radius = 180.200000
	e[63].atomic_radius = 199.500000
	e[64].atomic_radius = 178.700000
	e[65].atomic_radius = 176.300000
	e[66].atomic_radius = 175.200000
	e[67].atomic_radius = 174.300000
	e[68].atomic_radius = 173.400000
	e[69].atomic_radius = 172.400000
	e[70].atomic_radius = 194.000000
	e[71].atomic_radius = 171.800000
	e[72].atomic_radius = 156.400000
	e[73].atomic_radius = 143.000000
	e[74].atomic_radius = 137.000000
	e[75].atomic_radius = 137.100000
	e[76].atomic_radius = 133.800000
	e[77].atomic_radius = 135.700000
	e[78].atomic_radius = 137.300000
	e[79].atomic_radius = 144.200000
	e[80].atomic_radius = 150.300000
	e[81].atomic_radius = 170.000000
	e[82].atomic_radius = 175.000000
	e[83].atomic_radius = 154.500000
	e[84].atomic_radius = 167.300000
	e[85].atomic_radius = 0.000000
	e[86].atomic_radius = 0.000000
	e[87].atomic_radius = 270.000000
	e[88].atomic_radius = 223.000000
	e[89].atomic_radius = 187.800000
	e[90].atomic_radius = 179.800000
	e[91].atomic_radius = 156.100000
	e[92].atomic_radius = 138.500000
	e[93].atomic_radius = 130.000000
	e[94].atomic_radius = 151.300000
	e[95].atomic_radius = 0.000000
	e[96].atomic_radius = 0.000000
	e[97].atomic_radius = 0.000000
	e[98].atomic_radius = 0.000000
	e[99].atomic_radius = 0.000000
	e[100].atomic_radius = 0.000000
	e[101].atomic_radius = 0.000000
	e[102].atomic_radius = 0.000000
	e[103].atomic_radius = 0.000000

        # Nuclear charge (Slater)
        # 0 for those not available
	e[1].nuclear_charge = 1.00
	e[2].nuclear_charge = 1.70
	e[3].nuclear_charge = 1.30
	e[4].nuclear_charge = 1.95
	e[5].nuclear_charge = 2.60
	e[6].nuclear_charge = 3.25
	e[7].nuclear_charge = 3.90
	e[8].nuclear_charge = 4.55
	e[9].nuclear_charge = 5.20
	e[10].nuclear_charge = 5.85
	e[11].nuclear_charge = 2.20
	e[12].nuclear_charge = 2.85
	e[13].nuclear_charge = 3.50
	e[14].nuclear_charge = 4.15
	e[15].nuclear_charge = 4.80
	e[16].nuclear_charge = 5.45
	e[17].nuclear_charge = 6.10
	e[18].nuclear_charge = 6.75
	e[19].nuclear_charge = 2.20
	e[20].nuclear_charge = 2.85
	e[21].nuclear_charge = 3.00
	e[22].nuclear_charge = 3.15
	e[23].nuclear_charge = 3.30
	e[24].nuclear_charge = 3.45
	e[25].nuclear_charge = 3.60
	e[26].nuclear_charge = 3.75
	e[27].nuclear_charge = 3.90
	e[28].nuclear_charge = 4.05
	e[29].nuclear_charge = 4.20
	e[30].nuclear_charge = 4.35
	e[31].nuclear_charge = 5.00
	e[32].nuclear_charge = 5.65
	e[33].nuclear_charge = 6.30
	e[34].nuclear_charge = 6.95
	e[35].nuclear_charge = 7.60
	e[36].nuclear_charge = 8.25
	e[37].nuclear_charge = 2.20
	e[38].nuclear_charge = 2.85
	e[39].nuclear_charge = 3.00
	e[40].nuclear_charge = 3.15
	e[41].nuclear_charge = 3.30
	e[42].nuclear_charge = 3.45
	e[43].nuclear_charge = 3.60
	e[44].nuclear_charge = 3.75
	e[45].nuclear_charge = 3.90
	e[46].nuclear_charge = 4.05
	e[47].nuclear_charge = 4.20
	e[48].nuclear_charge = 4.35
	e[49].nuclear_charge = 5.00
	e[50].nuclear_charge = 5.65
	e[51].nuclear_charge = 6.30
	e[52].nuclear_charge = 6.95
	e[53].nuclear_charge = 7.60
	e[54].nuclear_charge = 8.25
	e[55].nuclear_charge = 2.20
	e[56].nuclear_charge = 2.85
	e[57].nuclear_charge = 2.85
	e[58].nuclear_charge = 2.85
	e[59].nuclear_charge = 2.85
	e[60].nuclear_charge = 2.85
	e[61].nuclear_charge = 2.85
	e[62].nuclear_charge = 2.85
	e[63].nuclear_charge = 2.85
	e[64].nuclear_charge = 2.85
	e[65].nuclear_charge = 2.85
	e[66].nuclear_charge = 2.85
	e[67].nuclear_charge = 2.85
	e[68].nuclear_charge = 2.85
	e[69].nuclear_charge = 2.85
	e[70].nuclear_charge = 2.854
	e[71].nuclear_charge = 3.00
	e[72].nuclear_charge = 3.15
	e[73].nuclear_charge = 3.30
	e[74].nuclear_charge = 4.35
	e[75].nuclear_charge = 3.60
	e[76].nuclear_charge = 3.75
	e[77].nuclear_charge = 3.90
	e[78].nuclear_charge = 4.05
	e[79].nuclear_charge = 4.20
	e[80].nuclear_charge = 4.35
	e[81].nuclear_charge = 5.00
	e[82].nuclear_charge = 5.65
	e[83].nuclear_charge = 6.30
	e[84].nuclear_charge = 6.95
	e[85].nuclear_charge = 7.60
	e[86].nuclear_charge = 8.25
	e[87].nuclear_charge = 2.20
	e[88].nuclear_charge = 1.65
	e[89].nuclear_charge = 1.8
	e[90].nuclear_charge = 1.95
	e[91].nuclear_charge = 1.80
	e[92].nuclear_charge = 1.80
	e[93].nuclear_charge = 1.80
	e[94].nuclear_charge = 1.65
	e[95].nuclear_charge = 4.65
	e[96].nuclear_charge = 1.80
	e[97].nuclear_charge = 1.65
	e[98].nuclear_charge = 1.65
	e[99].nuclear_charge = 1.65
	e[100].nuclear_charge = 1.65
	e[101].nuclear_charge = 1.65
	e[102].nuclear_charge = 1.65
	e[103].nuclear_charge = 1.8
 
	e[1].abundance = 0.880000
	e[2].abundance = 0.000000
	e[3].abundance = 0.006000
	e[4].abundance = 0.000500
	e[5].abundance = 0.001000
	e[6].abundance = 0.090000
	e[7].abundance = 0.030000
	e[8].abundance = 49.400000
	e[9].abundance = 0.030000
	e[10].abundance = 0.000000
	e[11].abundance = 2.640000
	e[12].abundance = 1.940000
	e[13].abundance = 7.570000
	e[14].abundance = 25.800000
	e[15].abundance = 0.090000
	e[16].abundance = 0.050000
	e[17].abundance = 0.190000
	e[18].abundance = 0.000400
	e[19].abundance = 2.400000
	e[20].abundance = 3.390000
	e[21].abundance = 0.000500
	e[22].abundance = 0.410000
	e[23].abundance = 0.010000
	e[24].abundance = 0.020000
	e[25].abundance = 0.090000
	e[26].abundance = 4.700000
	e[27].abundance = 0.004000
	e[28].abundance = 0.010000
	e[29].abundance = 0.010000
	e[30].abundance = 0.010000
	e[31].abundance = 0.001000
	e[32].abundance = 0.000600
	e[33].abundance = 0.000600
	e[34].abundance = 0.000100
	e[35].abundance = 0.000600
	e[36].abundance = 0.000000
	e[37].abundance = 0.030000
	e[38].abundance = 0.010000
	e[39].abundance = 0.003000
	e[40].abundance = 0.020000
	e[41].abundance = 0.002000
	e[42].abundance = 0.001000
	e[43].abundance = 0.000000
	e[44].abundance = 0.000002
	e[45].abundance = 0.000000
	e[46].abundance = 0.000001
	e[47].abundance = 0.000010
	e[48].abundance = 0.000030
	e[49].abundance = 0.000010
	e[50].abundance = 0.001000
	e[51].abundance = 0.000100
	e[52].abundance = 0.000001
	e[53].abundance = 0.000006
	e[54].abundance = 0.000000
	e[55].abundance = 0.000600
	e[56].abundance = 0.030000
	e[57].abundance = 0.002000
	e[58].abundance = 0.004000
	e[59].abundance = 0.000500
	e[60].abundance = 0.002000
	e[61].abundance = 0.000000
	e[62].abundance = 0.000600
	e[63].abundance = 0.000010
	e[64].abundance = 0.000600
	e[65].abundance = 0.000090
	e[66].abundance = 0.000400
	e[67].abundance = 0.000100
	e[68].abundance = 0.000200
	e[69].abundance = 0.000020
	e[70].abundance = 0.000020
	e[71].abundance = 0.000070
	e[72].abundance = 0.000400
	e[73].abundance = 0.000800
	e[74].abundance = 0.006000
	e[75].abundance = 0.000000
	e[76].abundance = 0.000001
	e[77].abundance = 0.000000
	e[78].abundance = 0.000000
	e[79].abundance = 0.000000
	e[80].abundance = 0.000040
	e[81].abundance = 0.000030
	e[82].abundance = 0.002000
	e[83].abundance = 0.000020
	e[84].abundance = 0.000000
	e[85].abundance = 0.000000
	e[86].abundance = 0.000000
	e[87].abundance = 0.000000
	e[88].abundance = 0.000000
	e[89].abundance = 0.000000
	e[90].abundance = 0.001000
	e[91].abundance = 9.0
	e[92].abundance = 0.000300
	e[93].abundance = 0.000000
	e[94].abundance = 0.000000
	e[95].abundance = 0.000000
	e[96].abundance = 0.000000
	e[97].abundance = 0.000000
	e[98].abundance = 0.000000
	e[99].abundance = 0.000000
	e[100].abundance = 0.000000
	e[101].abundance = 0.000000
	e[102].abundance = 0.000000
	e[103].abundance = 0.000000

	# Electron Aff.
	# 0 for those not available
	# Defined as 0 for Elements 2, 25,66 and 72
	e[1].electron_affinity = 72.8
	e[2].electron_affinity = 0.0
	e[3].electron_affinity = 59.6
	e[4].electron_affinity = -18
	e[5].electron_affinity = 26.7
	e[6].electron_affinity = 121.9
	e[7].electron_affinity = -7
	e[8].electron_affinity = 141
	e[9].electron_affinity = 328
	e[10].electron_affinity = -29
	e[11].electron_affinity = 52.9
	e[12].electron_affinity = -21
	e[13].electron_affinity = 44
	e[14].electron_affinity = 133.6
	e[15].electron_affinity = 72
	e[16].electron_affinity = 200.4
	e[17].electron_affinity = 349.0
	e[18].electron_affinity = -35
	e[19].electron_affinity = 48.4
	e[20].electron_affinity = -186
	e[21].electron_affinity = 18.1
	e[22].electron_affinity = 7.6
	e[23].electron_affinity = 50.7
	e[24].electron_affinity = 64.3
	e[25].electron_affinity = 0
	e[26].electron_affinity = 15.7
	e[27].electron_affinity = 63.8
	e[28].electron_affinity = 156
	e[29].electron_affinity = 188.5
	e[30].electron_affinity = 9
	e[31].electron_affinity = 30
	e[32].electron_affinity = 116
	e[33].electron_affinity = 78
	e[34].electron_affinity = 195
	e[35].electron_affinity = 324.7
	e[36].electron_affinity = -39
	e[37].electron_affinity = 46.9
	e[38].electron_affinity = -164
	e[39].electron_affinity = 29.6
	e[40].electron_affinity = 41.1
	e[41].electron_affinity = 86.2
	e[42].electron_affinity = 72.0
	e[43].electron_affinity = 96
	e[44].electron_affinity = 101
	e[45].electron_affinity = 109.7
	e[46].electron_affinity = 53.7
	e[47].electron_affinity = 125.7
	e[48].electron_affinity = -26
	e[49].electron_affinity = 30
	e[50].electron_affinity = 116
	e[51].electron_affinity = 101
	e[52].electron_affinity = 190.2
	e[53].electron_affinity = 295.2
	e[54].electron_affinity = -41
	e[55].electron_affinity = 45.5
	e[56].electron_affinity = -46
	e[57].electron_affinity = 50
	e[58].electron_affinity = 50
	e[59].electron_affinity = 50
	e[60].electron_affinity = 50
	e[61].electron_affinity = 50
	e[62].electron_affinity = 50
	e[63].electron_affinity = 50
	e[64].electron_affinity = 50
	e[65].electron_affinity = 50
	e[66].electron_affinity = 0
	e[67].electron_affinity = 50
	e[68].electron_affinity = 50
	e[69].electron_affinity = 50
	e[70].electron_affinity = 50
	e[71].electron_affinity = 50
	e[72].electron_affinity = 0
	e[73].electron_affinity = 14
	e[74].electron_affinity = 78.6
	e[75].electron_affinity = 14
	e[76].electron_affinity = 106
	e[77].electron_affinity = 151
	e[78].electron_affinity = 205.3
	e[79].electron_affinity = 222.8
	e[80].electron_affinity = -18
	e[81].electron_affinity = 20
	e[82].electron_affinity = 35.1
	e[83].electron_affinity = 91.3
	e[84].electron_affinity = 183
	e[85].electron_affinity = 270
	e[86].electron_affinity = -41
	e[87].electron_affinity = 44
	e[88].electron_affinity = 159
	e[89].electron_affinity = 406
	e[90].electron_affinity = 598.3
	e[91].electron_affinity = 607
	e[92].electron_affinity = 535.6
	e[93].electron_affinity = 0
	e[94].electron_affinity = 0
	e[95].electron_affinity = 0
	e[96].electron_affinity = 0
	e[97].electron_affinity = 0
	e[98].electron_affinity = 0
	e[99].electron_affinity = 50
	e[100].electron_affinity = 0
	e[101].electron_affinity = 0
	e[102].electron_affinity = 0
	e[103].electron_affinity = 0

	# Electronegativity (Pauling)
	# 0 for those not available 
	# Some noble gases defined as zero
	e[1].electronegativity = 2.20
	e[2].electronegativity = 0
	e[3].electronegativity = 0.98
	e[4].electronegativity = 1.57
	e[5].electronegativity = 2.04
	e[6].electronegativity = 2.55
	e[7].electronegativity = 3.04
	e[8].electronegativity = 3.44
	e[9].electronegativity = 3.98
	e[10].electronegativity = 0
	e[11].electronegativity = 0.93
	e[12].electronegativity = 1.31
	e[13].electronegativity = 1.61
	e[14].electronegativity = 1.90
	e[15].electronegativity = 2.19
	e[16].electronegativity = 2.58
	e[17].electronegativity = 3.16
	e[18].electronegativity = 0
	e[19].electronegativity = 0.82
	e[20].electronegativity = 1.00
	e[21].electronegativity = 1.36
	e[22].electronegativity = 1.54
	e[23].electronegativity = 1.63
	e[24].electronegativity = 1.66
	e[25].electronegativity = 1.55
	e[26].electronegativity = 1.83
	e[27].electronegativity = 1.88
	e[28].electronegativity = 1.91
	e[29].electronegativity = 1.90
	e[30].electronegativity = 1.65
	e[31].electronegativity = 1.81
	e[32].electronegativity = 2.01
	e[33].electronegativity = 2.18
	e[34].electronegativity = 2.55
	e[35].electronegativity = 2.96
	e[36].electronegativity = 0
	e[37].electronegativity = 0.82
	e[38].electronegativity = 0.95
	e[39].electronegativity = 1.22
	e[40].electronegativity = 1.33
	e[41].electronegativity = 1.6
	e[42].electronegativity = 2.16
	e[43].electronegativity = 1.9
	e[44].electronegativity = 2.2
	e[45].electronegativity = 2.28
	e[46].electronegativity = 2.20
	e[47].electronegativity = 1.93
	e[48].electronegativity = 1.96
	e[49].electronegativity = 1.78
	e[50].electronegativity = 1.96
	e[51].electronegativity = 2.05
	e[52].electronegativity = 2.1
	e[53].electronegativity = 2.66
	e[54].electronegativity = 2.6
	e[55].electronegativity = 0.79
	e[56].electronegativity = 0.89
	e[57].electronegativity = 1.10
	e[58].electronegativity = 1.12
	e[59].electronegativity = 1.13
	e[60].electronegativity = 1.14
	e[61].electronegativity = 0
	e[62].electronegativity = 1.17
	e[63].electronegativity = 0
	e[64].electronegativity = 1.20
	e[65].electronegativity = 0
	e[66].electronegativity = 1.22
	e[67].electronegativity = 1.23
	e[68].electronegativity = 1.24
	e[69].electronegativity = 1.25
	e[70].electronegativity = 0
	e[71].electronegativity = 1.27
	e[72].electronegativity = 1.3
	e[73].electronegativity = 1.5
	e[74].electronegativity = 2.36
	e[75].electronegativity = 1.9
	e[76].electronegativity = 2.2
	e[77].electronegativity = 2.20
	e[78].electronegativity = 2.28
	e[79].electronegativity = 2.54
	e[80].electronegativity = 2.00
	e[81].electronegativity = 2.04
	e[82].electronegativity = 2.33
	e[83].electronegativity = 2.02
	e[84].electronegativity = 2.0
	e[85].electronegativity = 2.2
	e[86].electronegativity = 0
	e[87].electronegativity = 0.7
	e[88].electronegativity = 0.89
	e[89].electronegativity = 1.1
	e[90].electronegativity = 1.3
	e[91].electronegativity = 1.5
	e[92].electronegativity = 1.38
	e[93].electronegativity = 1.36
	e[94].electronegativity = 1.28
	e[95].electronegativity = 1.3
	e[96].electronegativity = 1.3
	e[97].electronegativity = 1.3
	e[98].electronegativity = 1.3
	e[99].electronegativity = 1.3
	e[100].electronegativity = 1.3
	e[101].electronegativity = 1.3
	e[102].electronegativity = 1.3
	e[103].electronegativity = 1.3

	# ionization energy (in electronvolts].ionization_energy
	e[1].ionization_energy = 13.598
	e[2].ionization_energy = 24.587000
	e[3].ionization_energy = 5.392000
	e[4].ionization_energy = 9.322000
	e[5].ionization_energy = 8.298000
	e[6].ionization_energy = 11.260000
	e[7].ionization_energy = 14.534000
	e[8].ionization_energy = 13.618000
	e[9].ionization_energy = 17.422000
	e[10].ionization_energy = 21.564000
	e[11].ionization_energy = 5.139000
	e[12].ionization_energy = 7.646000
	e[13].ionization_energy = 5.986000
	e[14].ionization_energy = 8.151000
	e[15].ionization_energy = 10.486000
	e[16].ionization_energy = 10.360000
	e[17].ionization_energy = 12.967000
	e[18].ionization_energy = 15.759000
	e[19].ionization_energy = 4.341000
	e[20].ionization_energy = 6.113000
	e[21].ionization_energy = 6.540000
	e[22].ionization_energy = 6.820000
	e[23].ionization_energy = 6.740000
	e[24].ionization_energy = 6.766000
	e[25].ionization_energy = 7.435000
	e[26].ionization_energy = 7.870000
	e[27].ionization_energy = 7.860000
	e[28].ionization_energy = 7.635000
	e[29].ionization_energy = 7.726000
	e[30].ionization_energy = 9.394000
	e[31].ionization_energy = 5.999000
	e[32].ionization_energy = 7.899000
	e[33].ionization_energy = 9.810000
	e[34].ionization_energy = 9.752000
	e[35].ionization_energy = 11.814000
	e[36].ionization_energy = 13.999000
	e[37].ionization_energy = 4.177000
	e[38].ionization_energy = 5.695000
	e[39].ionization_energy = 6.380000
	e[40].ionization_energy = 6.840000
	e[41].ionization_energy = 6.880000
	e[42].ionization_energy = 7.099000
	e[43].ionization_energy = 7.280000
	e[44].ionization_energy = 7.370000
	e[45].ionization_energy = 7.460000
	e[46].ionization_energy = 8.340000
	e[47].ionization_energy = 7.576000
	e[48].ionization_energy = 8.993000
	e[49].ionization_energy = 5.786000
	e[50].ionization_energy = 7.344000
	e[51].ionization_energy = 8.641000
	e[52].ionization_energy = 9.009000
	e[53].ionization_energy = 10.451000
	e[54].ionization_energy = 12.130000
	e[55].ionization_energy = 3.894000
	e[56].ionization_energy = 5.212000
	e[57].ionization_energy = 5.577000
	e[58].ionization_energy = 5.470000
	e[59].ionization_energy = 5.420000
	e[60].ionization_energy = 5.490000
	e[61].ionization_energy = 5.550000
	e[62].ionization_energy = 5.630000
	e[63].ionization_energy = 5.670000
	e[64].ionization_energy = 6.140000
	e[65].ionization_energy = 5.850000
	e[66].ionization_energy = 5.930000
	e[67].ionization_energy = 6.020000
	e[68].ionization_energy = 6.100000
	e[69].ionization_energy = 6.180000
	e[70].ionization_energy = 6.254000
	e[71].ionization_energy = 5.426000
	e[72].ionization_energy = 7.000000
	e[73].ionization_energy = 7.890000
	e[74].ionization_energy = 7.980000
	e[75].ionization_energy = 7.880000
	e[76].ionization_energy = 8.700000
	e[77].ionization_energy = 9.100000
	e[78].ionization_energy = 9.000000
	e[79].ionization_energy = 9.255000
	e[80].ionization_energy = 10.437000
	e[81].ionization_energy = 6.108000
	e[82].ionization_energy = 6.108000
	e[83].ionization_energy = 7.289000
	e[84].ionization_energy = 8.420000
	e[85].ionization_energy = 9.500000
	e[86].ionization_energy = 10.748000
	e[87].ionization_energy = 4.000000
	e[88].ionization_energy = 5.279000
	e[89].ionization_energy = 6.900000
	e[90].ionization_energy = 6.950000
	e[91].ionization_energy = 0.000000
	e[92].ionization_energy = 6.080000
	e[93].ionization_energy = 0.000000
	e[94].ionization_energy = 5.800000
	e[95].ionization_energy = 6.000000
	e[96].ionization_energy = 0.000000
	e[97].ionization_energy = 0.000000
	e[98].ionization_energy = 0.000000
	e[99].ionization_energy = 0.000000
	e[100].ionization_energy = 0.000000
	e[101].ionization_energy = 0.000000
	e[102].ionization_energy = 0.000000
	e[103].ionization_energy = 0.000000


	# Ionic Radius (picometers)
	# Radius for smallest charge where more than one possible
	# Radius for H is for hydride 
	# 0 for those not available or those that don't form ions
	e[1].ionic_radius = 154
	e[2].ionic_radius = 0
	e[3].ionic_radius = 78
	e[4].ionic_radius = 34
	e[5].ionic_radius = 23
	e[6].ionic_radius = 260
	e[7].ionic_radius = 171
	e[8].ionic_radius = 132
	e[9].ionic_radius = 133
	e[10].ionic_radius = 112
	e[11].ionic_radius = 98
	e[12].ionic_radius = 78
	e[13].ionic_radius = 57
	e[14].ionic_radius = 271
	e[15].ionic_radius = 212
	e[16].ionic_radius = 184
	e[17].ionic_radius = 181
	e[18].ionic_radius = 154
	e[19].ionic_radius = 133
	e[20].ionic_radius = 106
	e[21].ionic_radius = 83
	e[22].ionic_radius = 80
	e[23].ionic_radius = 72
	e[24].ionic_radius = 84
	e[25].ionic_radius = 91
	e[26].ionic_radius = 82
	e[27].ionic_radius = 82
	e[28].ionic_radius = 78
	e[29].ionic_radius = 96
	e[30].ionic_radius = 83
	e[31].ionic_radius = 113
	e[32].ionic_radius = 90
	e[33].ionic_radius = 69
	e[34].ionic_radius = 69
	e[35].ionic_radius = 196
	e[36].ionic_radius = 169
	e[37].ionic_radius = 149
	e[38].ionic_radius = 127
	e[39].ionic_radius = 106
	e[40].ionic_radius = 109
	e[41].ionic_radius = 74
	e[42].ionic_radius = 92
	e[43].ionic_radius = 95
	e[44].ionic_radius = 77
	e[45].ionic_radius = 86
	e[46].ionic_radius = 86
	e[47].ionic_radius = 113
	e[48].ionic_radius = 114
	e[49].ionic_radius = 132
	e[50].ionic_radius = 93
	e[51].ionic_radius = 89
	e[52].ionic_radius = 211
	e[53].ionic_radius = 220
	e[54].ionic_radius = 190
	e[55].ionic_radius = 165
	e[56].ionic_radius = 143
	e[57].ionic_radius = 122
	e[58].ionic_radius = 107
	e[59].ionic_radius = 106
	e[60].ionic_radius = 104
	e[61].ionic_radius = 106
	e[62].ionic_radius = 111
	e[63].ionic_radius = 112
	e[64].ionic_radius = 97
	e[65].ionic_radius = 93
	e[66].ionic_radius = 91
	e[67].ionic_radius = 89
	e[68].ionic_radius = 89
	e[69].ionic_radius = 87
	e[70].ionic_radius = 113
	e[71].ionic_radius = 85
	e[72].ionic_radius = 84
	e[73].ionic_radius = 72
	e[74].ionic_radius = 68
	e[75].ionic_radius = 72
	e[76].ionic_radius = 89
	e[77].ionic_radius = 89
	e[78].ionic_radius = 85
	e[79].ionic_radius = 137
	e[80].ionic_radius = 127
	e[81].ionic_radius = 149
	e[82].ionic_radius = 132
	e[83].ionic_radius = 96
	e[84].ionic_radius = 65
	e[85].ionic_radius = 227
	e[86].ionic_radius = 0
	e[87].ionic_radius = 180
	e[88].ionic_radius = 152
	e[89].ionic_radius = 118
	e[90].ionic_radius = 101
	e[91].ionic_radius = 113
	e[92].ionic_radius = 103
	e[93].ionic_radius = 110
	e[94].ionic_radius = 108
	e[95].ionic_radius = 107
	e[96].ionic_radius = 119
	e[97].ionic_radius =118
	e[98].ionic_radius = 117
	e[99].ionic_radius = 116
	e[100].ionic_radius = 115
	e[101].ionic_radius = 114
	e[102].ionic_radius = 113
	e[103].ionic_radius = 112

	# Thermal Conditions (W/mK at 300K)
	# 0 for those not available
	e[1].thermal_cond = 0.1815
	e[2].thermal_cond = 0.152
	e[3].thermal_cond = 84.7
	e[4].thermal_cond = 200
	e[5].thermal_cond = 27
	e[6].thermal_cond = 1960
	e[7].thermal_cond = 0.02598
	e[8].thermal_cond = 0.2674
	e[9].thermal_cond = 0.0279
	e[10].thermal_cond = 0.0493
	e[11].thermal_cond = 141
	e[12].thermal_cond = 156
	e[13].thermal_cond = 273
	e[14].thermal_cond = 148
	e[15].thermal_cond = 0.235
	e[16].thermal_cond = 0.269
	e[17].thermal_cond = 0.0089
	e[18].thermal_cond = 0.0177
	e[19].thermal_cond = 102.4
	e[20].thermal_cond = 200
	e[21].thermal_cond = 15.8
	e[22].thermal_cond = 21.9
	e[23].thermal_cond = 30.7
	e[24].thermal_cond = 93.7
	e[25].thermal_cond = 7.82
	e[26].thermal_cond = 80.2
	e[27].thermal_cond = 100
	e[28].thermal_cond = 90.7
	e[29].thermal_cond = 401
	e[30].thermal_cond = 116
	e[31].thermal_cond = 40.6
	e[32].thermal_cond = 59.9
	e[33].thermal_cond = 50.0
	e[34].thermal_cond = 2.04
	e[35].thermal_cond = 0.122
	e[36].thermal_cond = 0.00949
	e[37].thermal_cond = 58.2
	e[38].thermal_cond = 35.3
	e[39].thermal_cond = 17.2
	e[40].thermal_cond = 22.7
	e[41].thermal_cond = 53.7
	e[42].thermal_cond = 138
	e[43].thermal_cond = 50.6
	e[44].thermal_cond = 117
	e[45].thermal_cond = 150
	e[46].thermal_cond = 71.8
	e[47].thermal_cond = 429
	e[48].thermal_cond = 96.8
	e[49].thermal_cond = 81.6
	e[50].thermal_cond = 66.6
	e[51].thermal_cond = 24.3
	e[52].thermal_cond = 2.35
	e[53].thermal_cond = 0.449
	e[54].thermal_cond = 0.00569
	e[55].thermal_cond = 35.9
	e[56].thermal_cond = 18.4
	e[57].thermal_cond = 13.5
	e[58].thermal_cond = 11.4
	e[59].thermal_cond = 12.5
	e[60].thermal_cond = 16.5
	e[61].thermal_cond = 17.9
	e[62].thermal_cond = 13.3
	e[63].thermal_cond = 13.9
	e[64].thermal_cond = 10.6
	e[65].thermal_cond = 11.1
	e[66].thermal_cond = 10.7
	e[67].thermal_cond = 16.2
	e[68].thermal_cond = 14.3
	e[69].thermal_cond = 16.8
	e[70].thermal_cond = 34.9
	e[71].thermal_cond = 16.4
	e[72].thermal_cond = 23
	e[73].thermal_cond = 57.5
	e[74].thermal_cond = 174
	e[75].thermal_cond = 47.9
	e[76].thermal_cond = 87.6
	e[77].thermal_cond = 147
	e[78].thermal_cond = 71.6
	e[79].thermal_cond = 317
	e[80].thermal_cond = 8.34
	e[81].thermal_cond = 46.1
	e[82].thermal_cond = 35.3
	e[83].thermal_cond = 7.87
	e[84].thermal_cond = 20
	e[85].thermal_cond = 1.7
	e[86].thermal_cond = 0.00364
	e[87].thermal_cond = 15
	e[88].thermal_cond = 18.6
	e[89].thermal_cond = 12
	e[90].thermal_cond = 54.0
	e[91].thermal_cond = 47
	e[92].thermal_cond = 27.6
	e[93].thermal_cond = 6.3
	e[94].thermal_cond = 6.74
	e[95].thermal_cond = 10
	e[96].thermal_cond = 10
	e[97].thermal_cond = 10
	e[98].thermal_cond = 10
	e[99].thermal_cond = 10
	e[100].thermal_cond = 10
	e[101].thermal_cond = 10
	e[102].thermal_cond = 10
	e[103].thermal_cond = 10

	# mpt.m creates e[deg C].melting_point
	e[1].melting_point=-259.14 
	e[2].melting_point=-272.2
	e[3].melting_point=180.54
	e[4].melting_point=1278.000000
	e[5].melting_point=2300.
	e[6].melting_point=3550.000000
	e[7].melting_point=-209.86
	e[8].melting_point=-218.4
	e[9].melting_point=-219.62
	e[10].melting_point=-248.67
	e[11].melting_point=97.81
	e[12].melting_point=648.8
	e[13].melting_point=660.37
	e[14].melting_point=1410.
	e[15].melting_point=44.100000
	e[16].melting_point=112.8
	e[17].melting_point=-100.98
	e[18].melting_point=-189.2
	e[19].melting_point=63.65
	e[20].melting_point=839.000
	e[21].melting_point=1541.
	e[22].melting_point=1660.
	e[23].melting_point=1890. 
	e[24].melting_point=1857. 
	e[25].melting_point=1244.
	e[26].melting_point=1553.
	e[27].melting_point=1495.
	e[28].melting_point=1453.
	e[29].melting_point=1083.4 
	e[30].melting_point=419.58
	e[31].melting_point=29.78
	e[32].melting_point=937.4
	e[33].melting_point=817.00
	e[34].melting_point=217.
	e[35].melting_point=-7.2
	e[36].melting_point=-156.6
	e[37].melting_point=38.89
	e[38].melting_point=769.
	e[39].melting_point=1522
	e[40].melting_point=1852.00
	e[41].melting_point=2468.
	e[42].melting_point=2617.
	e[43].melting_point=2172.
	e[44].melting_point=2310.
	e[45].melting_point=1966 
	e[46].melting_point=1552.
	e[47].melting_point=961.93
	e[48].melting_point=320.9
	e[49].melting_point=156.61
	e[50].melting_point=231.9681
	e[51].melting_point=630.74
	e[52].melting_point=449.5 
	e[53].melting_point=113.5
	e[54].melting_point=-111.9
	e[55].melting_point=28.40 
	e[56].melting_point=725.
	e[57].melting_point=921 
	e[58].melting_point=799 
	e[59].melting_point=931
	e[60].melting_point=1021
	e[61].melting_point=1168
	e[62].melting_point=1077 
	e[63].melting_point=822 
	e[64].melting_point=1313
	e[65].melting_point=1356  
	e[66].melting_point=1356 
	e[67].melting_point=1474
	e[68].melting_point=1529 
	e[69].melting_point=1545   
	e[70].melting_point=819   
	e[71].melting_point=1663  
	e[72].melting_point=2227.0 
	e[73].melting_point=2996
	e[74].melting_point=3410. 
	e[75].melting_point=3180.
	e[76].melting_point=3045. 
	e[77].melting_point=2410.
	e[78].melting_point=1772.
	e[79].melting_point=1064.43
	e[80].melting_point=-38.87
	e[81].melting_point=303.5
	e[82].melting_point=327.502
	e[83].melting_point=271.3
	e[84].melting_point=254.
	e[85].melting_point=302.
	e[86].melting_point=-71.
	e[87].melting_point=27.
	e[88].melting_point=700.
	e[89].melting_point=1050.
	e[90].melting_point=1750.
	e[91].melting_point=1554.000000
	e[92].melting_point=1132.3
	e[93].melting_point=640. 
	e[94].melting_point=641.
	e[95].melting_point=994.
	e[96].melting_point=1340. 
	e[97].melting_point=986.
	e[98].melting_point=900.0000
	
	# bpt.m creates e[deg C].boiling_point
	e[1].boiling_point=-252.87
	e[2].boiling_point=-268.934
	e[3].boiling_point=1347
	e[4].boiling_point=2870.0
	e[5].boiling_point=2550
	e[6].boiling_point=4827.0
	e[7].boiling_point=-195.8
	e[8].boiling_point=-183.962
	e[9].boiling_point=-188.14
	e[10].boiling_point=-246.048
	e[11].boiling_point=882.9
	e[12].boiling_point=1090
	e[13].boiling_point=2467
	e[14].boiling_point=2355
	e[15].boiling_point=280
	e[16].boiling_point=444.674
	e[17].boiling_point=-34.6
	e[18].boiling_point=-185.7
	e[19].boiling_point=774
	e[20].boiling_point=1484
	e[21].boiling_point=2831
	e[22].boiling_point=3287
	e[23].boiling_point=3380
	e[24].boiling_point=2672
	e[25].boiling_point=1962
	e[26].boiling_point=2750
	e[27].boiling_point=2870
	e[28].boiling_point=2732
	e[29].boiling_point=2567
	e[30].boiling_point=907
	e[31].boiling_point=2403
	e[32].boiling_point=2830
	e[33].boiling_point=613.0
	e[34].boiling_point=684.9
	e[35].boiling_point=58.78
	e[36].boiling_point=-152.30
	e[37].boiling_point=688
	e[38].boiling_point=1384
	e[39].boiling_point=3338
	e[40].boiling_point=4377
	e[41].boiling_point=4742
	e[42].boiling_point=4612
	e[43].boiling_point=4877
	e[44].boiling_point=3900
	e[45].boiling_point=3727
	e[46].boiling_point=3140
	e[47].boiling_point=2212
	e[48].boiling_point=765
	e[49].boiling_point=2080
	e[50].boiling_point=2270
	e[51].boiling_point=1750
	e[52].boiling_point=989.8
	e[53].boiling_point=184.35
	e[54].boiling_point=-107.100000
	e[55].boiling_point=678.4
	e[56].boiling_point=1640
	e[57].boiling_point=3457
	e[58].boiling_point=3426
	e[59].boiling_point=3512
	e[60].boiling_point=3068
	e[61].boiling_point=2700
	e[62].boiling_point=1791
	e[63].boiling_point=1597
	e[64].boiling_point=3266
	e[65].boiling_point=3123
	e[66].boiling_point=2562
	e[67].boiling_point=2695
	e[68].boiling_point=2863
	e[69].boiling_point=1947
	e[70].boiling_point=1194
	e[71].boiling_point=3395
	e[72].boiling_point=4602
	e[73].boiling_point=5425
	e[74].boiling_point=5660
	e[75].boiling_point=5627
	e[76].boiling_point=5027
	e[77].boiling_point=4130
	e[78].boiling_point=3827
	e[79].boiling_point=2807
	e[80].boiling_point=356.58
	e[81].boiling_point=1457
	e[82].boiling_point=1740
	e[83].boiling_point=560
	e[84].boiling_point=962
	e[85].boiling_point=337
	e[86].boiling_point=-61.8
	e[87].boiling_point=677
	e[88].boiling_point=1140
	e[86].boiling_point=3200
	e[90].boiling_point=4790
	e[92].boiling_point=3818
	e[93].boiling_point=3902
	e[94].boiling_point=3232
	e[95].boiling_point=2607

        for i in range(1,nelements+1):
            e[i].create_var_dict()
        #end for

        #for i in range(len(e)):
        #    e[i].create_string_representation()
        ##end for


        isotope_masses  =  obj(
            H  = {1:1.00782503207,  2:2.0141017778,  3:3.0160492777},  
            He = {3:3.0160293191,  4:4.00260325415},  
            Li = {6:6.015122795,  7:7.01600455},  
            Be = {9:9.0121822},  
            B  = {10:10.0129370,  11:11.0093054},  
            C  = {12:12.0000000,  13:13.0033548378,  14:14.003241989},  
            N  = {14:14.0030740048,  15:15.0001088982},  
            O  = {16:15.99491461956,  17:16.99913170,  18:17.9991610},  
            F  = {19:18.99840322},  
            Ne = {20:19.9924401754,  21:20.99384668,  22:21.991385114},  
            Na = {23:22.9897692809},  
            Mg = {24:23.985041700,  25:24.98583692,  26:25.982592929},  
            Al = {27:26.98153863},  
            Si = {28:27.9769265325,  29:28.976494700,  30:29.97377017},  
            P  = {31:30.97376163},  
            S  = {32:31.97207100,  33:32.97145876,  34:33.96786690,  36:35.96708076},  
            Cl = {35:34.96885268,  37:36.96590259},  
            Ar = {36:35.967545106,  38:37.9627324,  40:39.9623831225},  
            K  = {39:38.96370668,  40:39.96399848,  41:40.96182576},  
            Ca = {40:39.96259098,  42:41.95861801,  43:42.9587666,  44:43.9554818,  46:45.9536926,  48:47.952534},  
            Sc = {45:44.9559119},  
            Ti = {46:45.9526316,  47:46.9517631,  48:47.9479463,  49:48.9478700,  50:49.9447912},  
            V  = {50:49.9471585,  51:50.9439595},  
            Cr = {50:49.9460442,  52:51.9405075,  53:52.9406494,  54:53.9388804},  
            Mn = {55:54.9380451},  
            Fe = {54:53.9396105,  56:55.9349375,  57:56.9353940,  58:57.9332756},  
            Co = {59:58.9331950},  
            Ni = {58:57.9353429,  60:59.9307864,  61:60.9310560,  62:61.9283451,  64:63.9279660},  
            Cu = {63:62.9295975,  65:64.9277895},  
            Zn = {64:63.9291422,  66:65.9260334,  67:66.9271273,  68:67.9248442,  70:69.9253193},  
            Ga = {69:68.9255736,  71:70.9247013},  
            Ge = {70:69.9242474,  72:71.9220758,  73:72.9234589,  74:73.9211778,  76:75.9214026},  
            As = {75:74.9215965},  
            Se = {74:73.9224764,  76:75.9192136,  77:76.9199140,  78:77.9173091,  80:79.9165213,  82:81.9166994},  
            Br = {79:78.9183371,  81:80.9162906},  
            Kr = {78:77.9203648,  80:79.9163790,  82:81.9134836,  83:82.914136,  84:83.911507,  86:85.91061073},  
            Rb = {85:84.911789738,  87:86.909180527},  
            Sr = {84:83.913425,  86:85.9092602,  87:86.9088771,  88:87.9056121},  
            Y  = {89:88.9058483},  
            Zr = {90:89.9047044,  91:90.9056458,  92:91.9050408,  94:93.9063152,  96:95.9082734},  
            Nb = {93:92.9063781},  
            Mo = {92:91.906811,  94:93.9050883,  95:94.9058421,  96:95.9046795,  97:96.9060215,  98:97.9054082,  100:99.907477},  
            Tc = {97:96.906365,  98:97.907216,  99:98.9062547},  
            Ru = {96:95.907598,  98:97.905287,  99:98.9059393,  100:99.9042195,  101:100.9055821,  102:101.9043493,  104:103.905433},  
            Rh = {103:102.905504},  
            Pd = {102:101.905609,  104:103.904036,  105:104.905085,  106:105.903486,  108:107.903892,  110:109.905153},  
            Ag = {107:106.905097,  109:108.904752},  
            Cd = {106:105.906459,  108:107.904184,  110:109.9030021,  111:110.9041781,  112:111.9027578,  113:112.9044017,  114:113.9033585,  116:115.904756},  
            In = {113:112.904058,  115:114.903878},  
            Sn = {112:111.904818,  114:113.902779,  115:114.903342,  116:115.901741,  117:116.902952,  118:117.901603,  119:118.903308,  120:119.9021947,  122:121.9034390,  124:123.9052739},  
            Sb = {121:120.9038157,  123:122.9042140},  
            Te = {120:119.904020,  122:121.9030439,  123:122.9042700,  124:123.9028179,  125:124.9044307,  126:125.9033117,  128:127.9044631,  130:129.9062244},  
            I  = {127:126.904473},  
            Xe = {124:123.9058930,  126:125.904274,  128:127.9035313,  129:128.9047794,  130:129.9035080,  131:130.9050824,  132:131.9041535,  134:133.9053945,  136:135.907219},  
            Cs = {133:132.905451933},  
            Ba = {130:129.9063208,  132:131.9050613,  134:133.9045084,  135:134.9056886,  136:135.9045759,  137:136.9058274,  138:137.9052472},  
            La = {138:137.907112,  139:138.9063533},  
            Ce = {136:135.907172,  138:137.905991,  140:139.9054387,  142:141.909244},  
            Pr = {141:140.9076528},  
            Nd = {142:141.9077233,  143:142.9098143,  144:143.9100873,  145:144.9125736,  146:145.9131169,  148:147.916893,  150:149.920891},  
            Pm = {145:144.912749,  147:146.9151385},  
            Sm = {144:143.911999,  147:146.9148979,  148:147.9148227,  149:148.9171847,  150:149.9172755,  152:151.9197324,  154:153.9222093},  
            Eu = {151:150.9198502,  153:152.9212303},  
            Gd = {152:151.9197910,  154:153.9208656,  155:154.9226220,  156:155.9221227,  157:156.9239601,  158:157.9241039,  160:159.9270541},  
            Tb = {159:158.9253468},  
            Dy = {156:155.924283,  158:157.924409,  160:159.9251975,  161:160.9269334,  162:161.9267984,  163:162.9287312,  164:163.9291748},  
            Ho = {165:164.9303221},  
            Er = {162:161.928778,  164:163.929200,  166:165.9302931,  167:166.9320482,  168:167.9323702,  170:169.9354643},  
            Tm = {169:168.9342133},  
            Yb = {168:167.933897,  170:169.9347618,  171:170.9363258,  172:171.9363815,  173:172.9382108,  174:173.9388621,  176:175.9425717},  
            Lu = {175:174.9407718,  176:175.9426863},  
            Hf = {174:173.940046,  176:175.9414086,  177:176.9432207,  178:177.9436988,  179:178.9458161,  180:179.9465500},  
            Ta = {180:179.9474648,  181:180.9479958},  
            W  = {180:179.946704,  182:181.9482042,  183:182.9502230,  184:183.9509312,  186:185.9543641},  
            Re = {185:184.9529550,  187:186.9557531},  
            Os = {184:183.9524891,  186:185.9538382,  187:186.9557505,  188:187.9558382,  189:188.9581475,  190:189.9584470,  192:191.9614807},  
            Ir = {191:190.9605940,  193:192.9629264},  
            Pt = {190:189.959932,  192:191.9610380,  194:193.9626803,  195:194.9647911,  196:195.9649515,  198:197.967893},  
            Au = {197:196.9665687},  
            Hg = {196:195.965833,  198:197.9667690,  199:198.9682799,  200:199.9683260,  201:200.9703023,  202:201.9706430,  204:203.9734939},  
            Tl = {203:202.9723442,  205:204.9744275},  
            Pb = {204:203.9730436,  206:205.9744653,  207:206.9758969,  208:207.9766521},  
            Bi = {209:208.9803987},  
            Po = {209:208.9824304,  210:209.9828737},  
            At = {210:209.987148,  211:210.9874963},  
            Rn = {211:210.990601,  220:220.0113940,  222:222.0175777},  
            Fr = {223:223.0197359},  
            Ra = {223:223.0185022,  224:224.0202118,  226:226.0254098,  228:228.0310703},  
            Ac = {227:227.0277521},  
            Th = {230:230.0331338,  232:232.0380553},  
            Pa = {231:231.0358840},  
            U  = {233:233.0396352,  234:234.0409521,  235:235.0439299,  236:236.0455680,  238:238.0507882},  
            Np = {236:236.046570,  237:237.0481734},  
            Pu = {238:238.0495599,  239:239.0521634,  240:240.0538135,  241:241.0568515,  242:242.0587426,  244:244.064204},  
            Am = {241:241.0568291,  243:243.0613811},  
            Cm = {243:243.0613891,  244:244.0627526,  245:245.0654912,  246:246.0672237,  247:247.070354,  248:248.072349},  
            Bk = {247:247.070307,  249:249.0749867},  
            Cf = {249:249.0748535,  250:250.0764061,  251:251.079587,  252:252.081626},  
            Es = {252:252.082980},  
            Fm = {257:257.095105},  
            Md = {258:258.098431,  260:260.10365},  
            No = {259:259.10103},  
            Lr = {262:262.10963},  
            Rf = {265:265.11670},  
            Db = {268:268.12545},  
            Sg = {271:271.13347},  
            Bh = {272:272.13803},  
            Hs = {270:270.13465},  
            Mt = {276:276.15116},  
            Ds = {281:281.16206},  
            Rg = {280:280.16447},  
            Cn = {285:285.17411}
            )


        self.nelements = nelements
        self.simple_elements = e

        self.elements = obj()
        for i in range(1,self.nelements+1):
            elem = self.simple_elements[i]
            element = Element(elem)
            self.elements[elem.symbol] = element
            self[elem.symbol] = element
        #end for

        isotopes = obj()
        for symbol,element in self.elements.iteritems():
            elem_isotopes = obj()
            for mass_number,mass in isotope_masses[symbol].iteritems():
                isotope = element.copy()
                isotope.atomic_weight = phys_value_dict(mass,'amu')
                elem_isotopes[mass_number] = isotope
            #end for
            isotopes[symbol] = elem_isotopes
        #end for
        self.isotopes = isotopes

    #end def __init__

    def show(self):
        for i in range(self.nelements):
            print
            print self.elements[i].string_rep
        #end for
    #end def show
#end class PeriodicTable


pt = PeriodicTable()
periodic_table = pt
ptable = pt




def is_element(name,symbol=False):
    s      = name
    iselem = False
    if isinstance(name,str):
        iselem = name in periodic_table.elements
        if not iselem:
            nlen = len(name)
            if name.find('_')!=-1:
                s,n = name.split('_',1)
            elif nlen>1 and name[1:].isdigit():
                s = name[0:1]
            elif nlen>2 and name[2:].isdigit():
                s = name[0:2]
            #end if
            if len(s)==1:
                s = s.upper()
            elif len(s)==2:
                s = s[0].upper()+s[1].lower()
            #end if
            iselem = s in periodic_table.elements
        #end if
    #end if
    if symbol:
        return iselem,s
    else:
        return iselem
    #end if
#end def is_element










