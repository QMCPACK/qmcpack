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

from __future__ import annotations
from dataclasses import dataclass
from enum import Enum


@dataclass(frozen=True)
class ElementData:
    """Dataclass for storing element data."""
    symbol: str
    atomic_number: int
    atomic_weight: float
    group: int
    isotopes: dict[int, float]

    def __hash__(self):
        return hash((
            self.symbol,
            self.atomic_number,
            self.atomic_weight,
            self.group,
            tuple(self.isotopes.keys()),
            tuple(self.isotopes.values()),
        ))
#end class ElementData


class Elements(ElementData, Enum):
    """Enumeration of all elements in the periodic table.

    Attributes
    ----------
    symbol : str
        In titlecase (H, He, ...)
    atomic_number : int
        Use 0 for a dummy element
        (symbol "Xx", name "Unknown", all properties zero)
    atomic_weight : float
        Average atomic weight in amu [1]_.
    group : int
        Group of the element on the periodic table.
        Lanthanides and Actinides are 0.
    isotopes : dict[int, float]
        A dictionary of the isotopes for the element [2]_.
        This can be accessed as ``Element.Name.isotopes[mass_number]``,
        which yields the relative atomic mass.

    References
    ----------
    .. [1] https://iupac.qmul.ac.uk/AtWt/
    .. [2] https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses

    Examples
    --------
    The fastest way to grab element data is with the following signatures:

    >>> Elements("Hydrogen") is Elements.Hydrogen
    True
    >>> Elements("H") is Elements.Hydrogen
    True
    >>> Elements(1) is Elements.Hydrogen
    True

    If you're unsure of the case of the input you can reliably get
    case-insensitive parsing with the default interface.

    >>> Elements("h") is Elements.Hydrogen
    True
    >>> Elements("hydrogen") is Elements.Hydrogen
    True

    It can also handle leading and trailing whitespace

    >>> Elements("Hydrogen ") is Elements.Hydrogen
    True
    >>> Elements(" H") is Elements.Hydrogen
    True

    If you think the input is up to **one** step from the correct
    signature, (e.g. ``Elements(int(val))``) you can still expect the
    default call to work.

    >>> Elements("1") is Elements.Hydrogen
    True
    >>> Elements(1.0) is Elements.Hydrogen
    True

    Printing an element calls its ``__str__`` method which will return
    just the atomic symbol.

    >>> print(Elements.Hydrogen)
    H

    This also works with f-strings and ``str.format`` calls.

    >>> print(f"{Elements.Hydrogen}")
    H

    If you want to see the data attached to the enum member
    (minus the isotopes), use ``repr``.

    >>> print(repr(Elements.Hydrogen))
    <Elements.Hydrogen: symbol='H', atomic_number=1, atomic_weight=1.008, group=1>

    You can also get this view by just directly typing the element.

    >>> Elements.Hydrogen
    <Elements.Hydrogen: symbol='H', atomic_number=1, atomic_weight=1.008, group=1>
    """

    def __new__(
        cls,
        symbol: str,
        atomic_number: int,
        atomic_weight: float,
        group: int,
        isotopes: dict[int, float]
    ):
        element = ElementData.__new__(cls)
        element._value_ = atomic_number
        return element


    # Override to not print isotopes
    def __repr__(self):
        # This prints on one line, but since it's nearly impossible to
        # grep for this regardless of if it's on one line in the code
        # it is split here for readability.
        return (
            f"<{self.__class__.__name__}.{self.name}: "
            f"symbol='{self.symbol}', "
            f"atomic_number={self.atomic_number}, "
            f"atomic_weight={self.atomic_weight}, "
            f"group={self.group}>"
        )


    def __str__(self) -> str:
        return self.symbol


    @classmethod
    def _missing_(cls, value):
        """Workaround to not having access to ``_add_alias_`` or
        ``_add_value_alias_`` from Python 3.13. This function
        automatically gets called when the traditional lookup fails.
        """
        if isinstance(value, str):
            value = value.strip()
            if value.isalpha(): # `Elements("h")` or `Elements("hydrogen")`
                val_title = value.title()
                if val_title in cls.__members__:
                    return cls.__members__[val_title]

            elif value.isdecimal(): # `Elements("1")`
                try:
                    value = int(value)
                except ValueError: # We don't support `Elements("1.0")`
                    pass

                if isinstance(value, int) and value <= cls.num_elements():
                    return cls._value2member_map_[value]

        raise ValueError(f"Can not determine element for input value: {value}!")


    @staticmethod
    def is_element(
        value: str,
        return_element: bool = False,
    ) -> bool | tuple[bool, Elements]:
        """Robust method that will try to match a wide array of element
        identifier formats, including all that are handled by the parent
        call signature ``Elements(value)``.

        Parameters
        ----------
        value : str
            The string to be checked.
        return_element : bool, default=False
            Return the element that was identified.
            This will return ``value`` if it could not be identified as
            an element.

        Examples
        --------
        Regular elements

        >>> Elements.is_element("H")
        True
        >>> Elements.is_element("He")
        True

        Elements with integer identifiers

        >>> Elements.is_element("H1")
        True
        >>> Elements.is_element("H10")
        True
        >>> Elements.is_element("He1")
        True
        >>> Elements.is_element("He10")
        True

        Elements with integer identifiers and underscores

        >>> Elements.is_element("H_1")
        True
        >>> Elements.is_element("H_10")
        True
        >>> Elements.is_element("He_1")
        True
        >>> Elements.is_element("He_10")
        True

        Elements with integer identifiers and hyphens

        >>> Elements.is_element("H-1")
        True
        >>> Elements.is_element("H-10")
        True
        >>> Elements.is_element("He-1")
        True
        >>> Elements.is_element("He-10")
        True

        You can optionally return the element as well

        >>> Elements.is_element("H_10", return_element=True)
        (True, <Elements.Hydrogen: symbol='H', atomic_number=1, atomic_weight=1.008, group=1>)
        >>> Elements.is_element("He-1", return_element=True)
        (True, <Elements.Helium: symbol='He', atomic_number=2, atomic_weight=4.002602, group=18>)

        If it can't determine the element and you asked it to return the
        element then it will return ``False`` and the supplied value.

        >>> Elements.is_element("Bean", return_element=True)
        (False, 'Bean')
        """

        if isinstance(value, Elements):
            if return_element:
                return True, value
            else:
                return True

        val_len = len(value)
        if "_" in value: # H_1
            value = value.split("_", maxsplit=1)[0]
        elif "-" in value: # H-1
            value = value.split("-", maxsplit=1)[0]
        elif val_len >= 2 and value[1:].isdigit(): # H1 / H10
            value = value[0:1]
        elif val_len >= 3 and value[2:].isdigit(): # He1 / He10
            value = value[0:2]

        value = value.title()
        if value in Elements.__members__:
            if return_element:
                return True, Elements(value)
            else:
                return True
        else:
            if return_element:
                return False, value
            else:
                return False

    @staticmethod
    def num_elements() -> int:
        return len(Elements) - 1 # Dummy atom

    # Name        Symbol  Number   Mass      Group  Isotopes
    Unknown       = "Xx",    0,    0.0,          0, {0:0.0}
    Hydrogen      = "H",     1,    1.0080,       1, {1:1.00782503223,    2:2.01410177812,   3:3.0160492779}
    Helium        = "He",    2,    4.002602,    18, {3:3.0160293201,     4:4.00260325413}
    Lithium       = "Li",    3,    6.94,         1, {6:6.0151228874,     7:7.0160034366}
    Beryllium     = "Be",    4,    9.0121831,    2, {9:9.012183065}
    Boron         = "B",     5,   10.81,        13, {10:10.01293695,    11:11.00930536}
    Carbon        = "C",     6,   12.011,       14, {12:12.0000000,     13:13.00335483507, 14:14.0032419884}
    Nitrogen      = "N",     7,   14.007,       15, {14:14.00307400443, 15:15.00010889888}
    Oxygen        = "O",     8,   15.999,       16, {16:15.99491461957, 17:16.9991317565,  18:17.99915961286}
    Fluorine      = "F",     9,   18.998403162, 17, {19:18.99840316273}
    Neon          = "Ne",   10,   20.1797,      18, {20:19.9924401762,  21:20.993846685,   22:21.991385114}
    Sodium        = "Na",   11,   22.98976928,   1, {23:22.989769282}
    Magnesium     = "Mg",   12,   24.305,        2, {24:23.985041697,   25:24.985836976,   26:25.982592968}
    Aluminum      = "Al",   13,   26.9815384,   13, {27:26.98153853}
    Silicon       = "Si",   14,   28.085,       14, {28:27.97692653465, 29:28.9764946649,  30:29.973770136}
    Phosphorus    = "P",    15,   30.973761998, 15, {31:30.97376199842}
    Sulfur        = "S",    16,   32.06,        16, {32:31.9720711744,  33:32.9714589098,  34:33.967867004,   36:35.96708071}
    Chlorine      = "Cl",   17,   35.45,        17, {35:34.968852682,   37:36.965902602}
    Argon         = "Ar",   18,   39.95,        18, {36:35.967545105,   38:37.96273211,    40:39.9623831237}
    Potassium     = "K",    19,   39.0983,       1, {39:38.9637064864,  40:39.963998166,   41:40.9618252579}
    Calcium       = "Ca",   20,   40.078,        2, {40:39.962590863,   42:41.95861783,    43:42.95876644,    44:43.95548156,     46:45.953689,     48:47.95252276}
    Scandium      = "Sc",   21,   44.955907,     3, {45:44.95590828}
    Titanium      = "Ti",   22,   47.867,        4, {46:45.95262772,    47:46.95175879,    48:47.94794198,    49:48.94786568,     50:49.94478689}
    Vanadium      = "V",    23,   50.9415,       5, {50:49.94715601,    51:50.94395704}
    Chromium      = "Cr",   24,   51.9961,       6, {50:49.94604183,    52:51.94050623,    53:52.94064815,    54:53.93887916}
    Manganese     = "Mn",   25,   54.938043,     7, {55:54.93804391}
    Iron          = "Fe",   26,   55.845,        8, {54:53.93960899,    56:55.93493633,    57:56.93539284,    58:57.93327443}
    Cobalt        = "Co",   27,   58.933194,     9, {59:58.93319429}
    Nickel        = "Ni",   28,   58.6934,      10, {58:57.93534241,    60:59.93078588,    61:60.93105557,    62:61.92834537,     64:63.92796682}
    Copper        = "Cu",   29,   63.546,       11, {63:62.92959772,    65:64.9277897}
    Zinc          = "Zn",   30,   65.38,        12, {64:63.92914201,    66:65.92603381,    67:66.92712775,    68:67.92484455,     70:69.9253192}
    Gallium       = "Ga",   31,   69.723,       13, {69:68.9255735,     71:70.92470258}
    Germanium     = "Ge",   32,   72.630,       14, {70:69.92424875,    72:71.922075826,   73:72.923458956,   74:73.921177761,    76:75.921402726}
    Arsenic       = "As",   33,   74.921595,    15, {75:74.92159457}
    Selenium      = "Se",   34,   78.971,       16, {74:73.922475934,   76:75.919213704,   77:76.919914154,   78:77.91730928,     80:79.9165218,     82:81.9166995}
    Bromine       = "Br",   35,   79.904,       17, {79:78.9183376,     81:80.9162897}
    Krypton       = "Kr",   36,   83.798,       18, {78:77.92036494,    80:79.91637808,    82:81.91348273,    83:82.91412716,     84:83.9114977282,  86:85.9106106269}
    Rubidium      = "Rb",   37,   85.4678,       1, {85:84.9117897379,  87:86.909180531}
    Strontium     = "Sr",   38,   87.62,         2, {84:83.9134191,     86:85.9092606,     87:86.9088775,     88:87.9056125}
    Yttrium       = "Y",    39,   88.905838,     3, {89:88.9058403}
    Zirconium     = "Zr",   40,   91.222,        4, {90:89.9046977,     91:90.9056396,     92:91.9050347,     94:93.9063108,      96:95.9082714}
    Niobium       = "Nb",   41,   92.90637,      5, {93:92.906373}
    Molybdenum    = "Mo",   42,   95.95,         6, {92:91.90680796,    94:93.9050849,     95:94.90583877,    96:95.90467612,     97:96.90601812,    98:97.90540482,  100:99.9074718}
    Technetium    = "Tc",   43,   97.0,          7, {97:96.9063667,     98:97.9072124,     99:98.9062508}
    Ruthenium     = "Ru",   44,  101.07,         8, {96:95.90759025,    98:97.9052868,     99:98.9059341,    100:99.9042143,     101:100.9055769,   102:101.9043441,  104:103.9054275}
    Rhodium       = "Rh",   45,  102.90549,      9, {103:102.905498}
    Palladium     = "Pd",   46,  106.42,        10, {102:101.9056022,  104:103.9040305,   105:104.9050796,   106:105.9034804,    108:107.9038916,   110:109.9051722}
    Silver        = "Ag",   47,  107.8682,      11, {107:106.9050916,  109:108.9047553}
    Cadmium       = "Cd",   48,  112.414,       12, {106:105.9064599,  108:107.9041834,   110:109.90300661,  111:110.90418287,   112:111.90276287,  113:112.90440813, 114:113.90336509,   116:115.90476315}
    Indium        = "In",   49,  114.818,       13, {113:112.90406184, 115:114.903878776}
    Tin           = "Sn",   50,  118.710,       14, {112:111.90482387, 114:113.9027827,   115:114.903344699, 116:115.9017428,    117:116.90295398,  118:117.90160657, 119:118.90331117,   120:119.90220163, 122:121.9034438, 124:123.9052766}
    Antimony      = "Sb",   51,  121.760,       15, {121:120.903812,   123:122.9042132}
    Tellurium     = "Te",   52,  127.60,        16, {120:119.9040593,  122:121.9030435,   123:122.9042698,   124:123.9028171,    125:124.9044299,   126:125.9033109,  128:127.90446128,   130:129.906222748}
    Iodine        = "I",    53,  126.90447,     17, {127:126.9044719}
    Xenon         = "Xe",   54,  131.293,       18, {124:123.905892,   126:125.9042983,   128:127.903531,    129:128.9047808611, 130:129.903509349, 131:130.90508406, 132:131.9041550856, 134:133.90539466, 136:135.907214484}
    Cesium        = "Cs",   55,  132.90545196,   1, {133:132.905451961}
    Barium        = "Ba",   56,  137.327,        2, {130:129.9063207,  132:131.9050611,   134:133.90450818,  135:134.90568838,   136:135.90457573,  137:136.90582714, 138:137.905247}
    Lanthanum     = "La",   57,  138.90547,      3, {138:137.9071149,  139:138.9063563}
    Cerium        = "Ce",   58,  140.116,        0, {136:135.90712921, 138:137.905991,    140:139.9054431,   142:141.9092504}
    Praseodymium  = "Pr",   59,  140.90766,      0, {141:140.9076576}
    Neodymium     = "Nd",   60,  144.242,        0, {142:141.907729,   143:142.90982,     144:143.910093,    145:144.9125793,    146:145.9131226,   148:147.9168993,  150:149.9209022}
    Promethium    = "Pm",   61,  145.0,          0, {145:144.9127559,  147:146.915145}
    Samarium      = "Sm",   62,  150.36,         0, {144:143.9120065,  147:146.9149044,   148:147.9148292,   149:148.9171921,    150:149.9172829,   152:151.9197397,  154:153.9222169}
    Europium      = "Eu",   63,  151.964,        0, {151:150.9198578,  153:152.921238}
    Gadolinium    = "Gd",   64,  157.249,        0, {152:151.9197995,  154:153.9208741,   155:154.9226305,   156:155.9221312,    157:156.9239686,   158:157.9241123,  160:159.9270624}
    Terbium       = "Tb",   65,  158.925354,     0, {159:158.9253547}
    Dysprosium    = "Dy",   66,  162.500,        0, {156:155.9242847,  158:157.9244159,   160:159.9252046,   161:160.9269405,    162:161.9268056,   163:162.9287383,  164:163.9291819}
    Holmium       = "Ho",   67,  164.930329,     0, {165:164.9303288}
    Erbium        = "Er",   68,  167.259,        0, {162:161.9287884,  164:163.9292088,   166:165.9302995,   167:166.9320546,    168:167.9323767,   170:169.9354702}
    Thulium       = "Tm",   69,  168.934219,     0, {169:168.9342179}
    Ytterbium     = "Yb",   70,  173.045,        0, {168:167.9338896,  170:169.9347664,   171:170.9363302,   172:171.9363859,    173:172.9382151,   174:173.9388664,  176:175.9425764}
    Lutetium      = "Lu",   71,  174.96669,      0, {175:174.9407752,  176:175.9426897}
    Hafnium       = "Hf",   72,  178.486,        4, {174:173.9400461,  176:175.9414076,   177:176.9432277,   178:177.9437058,    179:178.9458232,   180:179.946557}
    Tantalum      = "Ta",   73,  180.94788,      5, {180:179.9474648,  181:180.9479958}
    Tungsten      = "W",    74,  183.84,         6, {180:179.9467108,  182:181.94820394,  183:182.95022275,  184:183.95093092,   186:185.9543628}
    Rhenium       = "Re",   75,  186.207,        7, {185:184.9529545,  187:186.9557501}
    Osmium        = "Os",   76,  190.23,         8, {184:183.9524885,  186:185.953835,    187:186.9557474,   188:187.9558352,    189:188.9581442,   190:189.9584437,  192:191.961477}
    Iridium       = "Ir",   77,  192.217,        9, {191:190.9605893,  193:192.9629216}
    Platinum      = "Pt",   78,  195.084,       10, {190:189.9599297,  192:191.9610387,   194:193.9626809,   195:194.9647917,    196:195.96495209,  198:197.9678949}
    Gold          = "Au",   79,  196.966570,    11, {197:196.96656879}
    Mercury       = "Hg",   80,  200.592,       12, {196:195.9658326,  198:197.9667686,   199:198.96828064,  200:199.96832659,   201:200.97030284,  202:201.9706434,  204:203.97349398}
    Thallium      = "Tl",   81,  204.38,        13, {203:202.9723446,  205:204.9744278}
    Lead          = "Pb",   82,  207.2,         14, {204:203.973044,   206:205.9744657,   207:206.9758973,   208:207.9766525}
    Bismuth       = "Bi",   83,  208.98040,     15, {209:208.9803991}
    Polonium      = "Po",   84,  209.0,         16, {209:208.9824308,  210:209.9828741}
    Astatine      = "At",   85,  210.0,         17, {210:209.9871479,  211:210.9874966}
    Radon         = "Rn",   86,  222.0,         18, {211:210.9906011,  220:220.0113941,   222:222.0175782}
    Francium      = "Fr",   87,  223.0,          1, {223:223.019736}
    Radium        = "Ra",   88,  226.0,          2, {223:223.0185023,  224:224.020212,    226:226.0254103,   228:228.0310707}
    Actinium      = "Ac",   89,  227.0,          3, {227:227.0277523}
    Thorium       = "Th",   90,  232.0377,       0, {230:230.0331341,  232:232.0380558}
    Protactinium  = "Pa",   91,  231.03588,      0, {231:231.0358842}
    Uranium       = "U",    92,  238.02891,      0, {233:233.0396355,  234:234.0409523,   235:235.0439301,   236:236.0455682,    238:238.0507884}
    Neptunium     = "Np",   93,  237.0,          0, {236:236.04657,    237:237.0481736}
    Plutonium     = "Pu",   94,  244.0,          0, {238:238.0495601,  239:239.0521636,   240:240.0538138,   241:241.0568517,    242:242.0587428,   244:244.0642053}
    Americium     = "Am",   95,  243.0,          0, {241:241.0568293,  243:243.0613813}
    Curium        = "Cm",   96,  247.0,          0, {243:243.0613893,  244:244.0627528,   245:245.0654915,   246:246.0672238,    247:247.0703541,   248:248.0723499}
    Berkelium     = "Bk",   97,  247.0,          0, {247:247.0703073,  249:249.0749877}
    Californium   = "Cf",   98,  251.0,          0, {249:249.0748539,  250:250.0764062,   251:251.0795886,   252:252.0816272}
    Einsteinium   = "Es",   99,  252.0,          0, {252:252.08298}
    Fermium       = "Fm",  100,  257.0,          0, {257:257.0951061}
    Mendelevium   = "Md",  101,  258.0,          0, {258:258.0984315,  260:260.10365}
    Nobelium      = "No",  102,  259.0,          0, {259:259.10103}
    Lawrencium    = "Lr",  103,  262.0,          0, {262:262.10961}
    Rutherfordium = "Rf",  104,  267.0,          4, {267:267.12179}
    Dubnium       = "Db",  105,  270.0,          5, {268:268.12567}
    Seaborgium    = "Sg",  106,  269.0,          6, {271:271.13393}
    Bohrium       = "Bh",  107,  270.0,          7, {272:272.13826}
    Hassium       = "Hs",  108,  270.0,          8, {270:270.13429}
    Meitnerium    = "Mt",  109,  278.0,          9, {276:276.15159}
    Darmstadtium  = "Ds",  110,  281.0,         10, {281:281.16451}
    Roentgenium   = "Rg",  111,  281.0,         11, {280:280.16514}
    Copernicium   = "Cn",  112,  285.0,         12, {285:285.17712}
    Nihonium      = "Nh",  113,  286.0,         13, {284:284.17873}
    Flerovium     = "Fl",  114,  289.0,         14, {289:289.19042}
    Moscovium     = "Mc",  115,  289.0,         15, {288:288.19274}
    Livermorium   = "Lv",  116,  293.0,         16, {293:293.20449}
    Tennessine    = "Ts",  117,  293.0,         17, {292:292.20746}
    Oganesson     = "Og",  118,  294.0,         18, {294:294.21392}
    # Name        Symbol  Number Mass        Group  Isotopes

    # Add aliases for each element.
    # Once we can use features from Python 3.13, this will be replaced by
    # the Enum sunder `_add_alias_`
    Xx = Unknown
    H  = Hydrogen
    He = Helium
    Li = Lithium
    Be = Beryllium
    B  = Boron
    C  = Carbon
    N  = Nitrogen
    O  = Oxygen
    F  = Fluorine
    Ne = Neon
    Na = Sodium
    Mg = Magnesium
    Al = Aluminum
    Si = Silicon
    P  = Phosphorus
    S  = Sulfur
    Cl = Chlorine
    Ar = Argon
    K  = Potassium
    Ca = Calcium
    Sc = Scandium
    Ti = Titanium
    V  = Vanadium
    Cr = Chromium
    Mn = Manganese
    Fe = Iron
    Co = Cobalt
    Ni = Nickel
    Cu = Copper
    Zn = Zinc
    Ga = Gallium
    Ge = Germanium
    As = Arsenic
    Se = Selenium
    Br = Bromine
    Kr = Krypton
    Rb = Rubidium
    Sr = Strontium
    Y  = Yttrium
    Zr = Zirconium
    Nb = Niobium
    Mo = Molybdenum
    Tc = Technetium
    Ru = Ruthenium
    Rh = Rhodium
    Pd = Palladium
    Ag = Silver
    Cd = Cadmium
    In = Indium
    Sn = Tin
    Sb = Antimony
    Te = Tellurium
    I  = Iodine
    Xe = Xenon
    Cs = Cesium
    Ba = Barium
    La = Lanthanum
    Ce = Cerium
    Pr = Praseodymium
    Nd = Neodymium
    Pm = Promethium
    Sm = Samarium
    Eu = Europium
    Gd = Gadolinium
    Tb = Terbium
    Dy = Dysprosium
    Ho = Holmium
    Er = Erbium
    Tm = Thulium
    Yb = Ytterbium
    Lu = Lutetium
    Hf = Hafnium
    Ta = Tantalum
    W  = Tungsten
    Re = Rhenium
    Os = Osmium
    Ir = Iridium
    Pt = Platinum
    Au = Gold
    Hg = Mercury
    Tl = Thallium
    Pb = Lead
    Bi = Bismuth
    Po = Polonium
    At = Astatine
    Rn = Radon
    Fr = Francium
    Ra = Radium
    Ac = Actinium
    Th = Thorium
    Pa = Protactinium
    U  = Uranium
    Np = Neptunium
    Pu = Plutonium
    Am = Americium
    Cm = Curium
    Bk = Berkelium
    Cf = Californium
    Es = Einsteinium
    Fm = Fermium
    Md = Mendelevium
    No = Nobelium
    Lr = Lawrencium
    Rf = Rutherfordium
    Db = Dubnium
    Sg = Seaborgium
    Bh = Bohrium
    Hs = Hassium
    Mt = Meitnerium
    Ds = Darmstadtium
    Rg = Roentgenium
    Cn = Copernicium
    Nh = Nihonium
    Fl = Flerovium
    Mc = Moscovium
    Lv = Livermorium
    Ts = Tennessine
    Og = Oganesson
#end class Element
