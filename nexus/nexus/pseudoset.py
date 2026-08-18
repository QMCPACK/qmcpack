"""Functionality for gathering, parsing, and validating collections of pseudopotentials."""


from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from os import PathLike
from pathlib import Path
from re import Pattern
from types import MappingProxyType
from typing import ClassVar, Literal

from .developer import DevBase, error, obj, warn
from .periodic_table import Elements
from .physical_system import PhysicalSystem
from .utilities import is_valid_filename


def pp_elem_label(filename,*,guard=False):
    if guard and not is_valid_filename(filename):
        msg = f"Pseudopotential file name {filename} is invalid!"
        raise RuntimeError(msg)

    el = ''
    for c in filename:
        if c=='.' or c=='_' or c=='-':
            break
        #end if
        el+=c
    #end for
    elem_label = el
    is_elem, element = Elements.is_element(el, return_element=True)
    if guard: 
        if not is_elem:
            msg = (
               f"cannot determine element for pseudopotential file: {filename}\n"
                "pseudopotential file names must be prefixed by an atomic symbol or label\n"
                "(e.g. Si, Si1, etc)"
                )
            raise RuntimeError(msg)
        #end if
        return elem_label, element.symbol
    else:
        if isinstance(element, Elements):
            return elem_label, element.symbol, is_elem
        else:
            return elem_label, element, is_elem
    #end if
#end def pp_elem_label


def read_upf_z_valence(file: PathLike) -> int | float:
    """Read Z-valence from a UPF-compliant pseudopotential file."""
    # Bind these to the function so we only compile them once.
    if not (
        hasattr(read_upf_z_valence, "zval_xml_like_pattern")
        and hasattr(read_upf_z_valence, "zval_old_pattern")
        ):
        # Regex:
        # `z[_ ]?valence` -> "z_valence" or "z valence"
        # ` *=? *`        -> "=" or " =" or " = " or " " (whitespace optional)
        # `([\d\.eEdD]+)` -> Capturing group gets any numbers in scientific notation.
        # `\"? *() *\"?`  -> Anything between quotes or not, with optional whitespace around it too.
        # Note: Both of these are similar, but the first is key then value and the second is value then key.
        zval_xml_like_pattern = re.compile(
            pattern = r'z[_ ]?valence *=? *\"? *([\d\.eEdD]+) *\"?',
            flags   = re.IGNORECASE,
            )
        zval_old_pattern = re.compile(
            pattern = r'\"? *([ \d\.eEdD]+) *\"? *=? *z[_ ]?valence',
            flags   = re.IGNORECASE,
            )
        read_upf_z_valence.zval_xml_like_pattern = zval_xml_like_pattern
        read_upf_z_valence.zval_old_pattern = zval_old_pattern
    #end if

    zval = None
    with open(file, "r") as pseudo:
        found_header_start = False
        while not found_header_start:
            line = pseudo.readline()
            if "<PP_HEADER" in line:
                found_header_start = True

        if "/>" in line or "</PP_HEADER>" in line: # One-line header
            zval = re.search(
                pattern = read_upf_z_valence.zval_xml_like_pattern,
                string  = line,
                )
        else:
            # We're at the header, but we don't know where the Z-valence is.
            # Search until we hit a line with a proper end token, or until we hit 200 lines.
            i = 0
            while i < 200:
                i += 1
                line = pseudo.readline().lower()
                if "valence" in line:
                    zval = re.search(
                        pattern = read_upf_z_valence.zval_xml_like_pattern,
                        string  = line,
                        )
                    if zval is None:
                        zval = re.search(
                            pattern = read_upf_z_valence.zval_old_pattern,
                            string  = line,
                            )
                    break
                elif "/>" in line or "</PP_HEADER>" in line:
                    break
                #end if
            #end while
        #end if

    if zval is None:
        error(
           f"Could not find Z valence in file: {file!s}\n"
            "You may need to provide the Z valence manually!"
           )
    else:
        zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        error(
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
    # Round to 8 digits
    if round(zval, 8).is_integer():
        return int(zval)
    else:
        return zval
#end def read_upf_z_valence



def read_qmcpack_xml_z_valence(file: PathLike) -> int | float:
    """Read the Z-valence from a QMCPACK-compatible XML pseudopotential file."""
        # Bind these to the function so we only compile them once.
    if not hasattr(read_qmcpack_xml_z_valence, "zval_pattern"):
        # Regex:
        # `zval  *= *`    -> "zval=" or "zval = " or "zval =" or "zval= "
        # `([\d\.eEdD]+)` -> Capturing group gets any numbers in scientific notation.
        # `\"? *() *\"?`  -> Anything between quotes or not, with optional whitespace around it too.
        read_qmcpack_xml_z_valence.zval_pattern = re.compile(r'zval *= *\" *([\d\.eEdD]+) *\"')

    header_lines = []
    with open(file, "r") as xml:
        header_started = False
        for line in xml:
            if "<header" in line:
                header_started = True

            if header_started:
                header_lines.append(line)

            if "/>" in line or "</header>" in line:
                if line not in header_lines:
                    header_lines.append(line)
                break

    header = " ".join(header_lines)
    zval = re.search(read_qmcpack_xml_z_valence.zval_pattern, header)

    if zval is None:
        error(
           f"Could not find Z valence in file: {file!s}\n"
            "You may need to provide the Z valence manually!"
           )
    else:
        zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        error(
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
    # Round to 8 digits
    if round(zval, 8).is_integer():
        return int(zval)
    else:
        return zval
#end def read_xml_z_valence



def read_potcar_z_valence(file: PathLike) -> int | float:
    """Read the Z-valence from a POTCAR file.

    This function uses the format specifications from the VASP wiki, and
    assumes that the file is a valid POTCAR, so the second line should
    be the Z-valence [1]_.

    References
    ----------
    .. [1] https://vasp.at/wiki/POTCAR#File_format
    """
    if not hasattr(read_potcar_z_valence, "zval_pattern"):
        # Regex:
        # `ZVAL ?= ?`     -> "ZVAL=" or "ZVAL = " or "ZVAL =" or "ZVAL= "
        # `([\d\.eEdD]+)` -> Capturing group gets any numbers in scientific notation.
        read_potcar_z_valence.zval_pattern = re.compile(r"ZVAL ?= ?([\d\.eEdD]+)")
    file = Path(file).resolve()
    with open(file, "r") as potcar:
        potcar.readline() # Skip first line
        z_valence = potcar.readline().strip()
        try:
            zval = float(z_valence)
        except ValueError: # Improperly formatted POTCAR, but try alternative location
            for line in potcar:
                if "ZVAL" in line:
                    zval = re.search(
                        pattern = read_potcar_z_valence.zval_pattern,
                        string  = line
                        )
                    break

            if zval is None:
                error(
                   f"Could not find Z valence in file: {file!s}\n"
                    "You may need to provide the Z valence manually!"
                   )
            else:
                zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        error(
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
    # Round to 8 digits
    if round(zval, 8).is_integer():
        return int(zval)
    else:
        return zval
#end def read_potcar_z_valence





def ppset(label,**codes_pps):
    if label in PseudoSet.labeled_pseudosets:
        msg = f'pseudopotential set label "{label}" is already registered'
        raise ValueError(msg)
    pps_coll = {}
    ref_elements = None
    for code,ppfiles in codes_pps.items():
        missing = set(ppfiles)-PseudoSet.pseudo_files.keys()
        if len(missing)>0:
            msg = f'pseudopotential files "{missing}" are not present in PseudoSet.pseudo_files'
            raise ValueError(msg)
        pps = PseudoSet([PseudoSet.pseudo_files[f] for f in ppfiles])
        code = PseudoSet._check_code_str(code)
        if code not in pps.codes:
            msg = f'pseudopotential files provided for code "{code}" are not compatible with that code'
            raise ValueError(msg)
        elements = set(pps.pseudos)
        if ref_elements is None:
            ref_elements = elements
        elif elements!=ref_elements:
            msg = f'pseudopotential set "{label}" must contain potentials for the same elements for each code'
            raise ValueError(msg)
        pps_coll[code] = pps
    PseudoSet.labeled_pseudosets[label] = pps_coll
#end def ppset



class PseudoSet(DevBase):
    """Object representing a set of pseudopotentials.

    Attributes
    ----------
    pseudos : dict of str: Path
        Dictionary mapping element labels to their pseudos.
    codes : set of one or more of {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"}
        Name of the program(s) that these pseudos are meant for. Due to
        overlap between some programs for file extension, this can
        occasionally contain more than one code.
    Zeff_map : Map of str: int
        A ``dict`` or ``obj`` mapping elements to their effective
        nuclear charges (Z-valences).
    pseudo_dirs : set of Path
        The directories that the pseudopotentials are stored in.
    pseudo_files : dict of str: str (class attribute)
        Dictionary mapping pseudopotential file names to their full paths.
    labeled_pseudosets : dict of str: dict of str: PseudoSet (class attribute)
        Pseudopotential sets registered by label, then compatible code.

    Parameters
    ----------
    pseudos : list of str/Path or map of str/Elements to Path
        A list of pseudopotential files or a ``dict``/``obj`` that maps
        elements to file paths.
    codes : one or more of {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf", "detect"}, default="detect"
        The name of the code that the pseudos are formatted for, or
        if ``"detect"``, will auto-detect the code name from the
        file extensions.
    Zeff_map : Map of str/Elements to int, optional
        A ``dict`` or ``obj`` mapping elements to their effective nuclear
        charges (Z-valences). If this is supplied, it will override any
        parts of the code that may try to parse the pseudopotential to
        get the Z-valence.
    skip_invalid : bool, default=False (keyword-only)
        If ``True``, then this will emit a warning rather than raise an
        error if a file is not found or if the file does not have a
        valid name.
    """

    file_exts = MappingProxyType({
        # https://www.quantum-espresso.org/Doc/INPUT_PW.html#id268
        "espresso": frozenset({".ncpp", ".upf", ".vdb", ".van", ".rrkj3"}),
        "gamess":   frozenset({".gms", ".gamess"}),
        "vasp":     frozenset({"potcar", ".vasp"}),
        "qmcpack":  frozenset({".xml"}),
        "rmg":      frozenset({".upf", ".xml"}),
        "pyscf":    frozenset({".nwchem", ".gth"})
        })
    known_codes = frozenset(file_exts.keys())

    pseudo_files      : ClassVar[dict[str, str]] = dict()
    labeled_pseudosets: ClassVar[dict[str, dict[str, PseudoSet]]] = dict()


    @staticmethod
    def pseudo_remap(code,pseudos,system):
        if pseudos is None:
            return {}
        if isinstance(pseudos, Mapping) and len(pseudos)>0:
            contains_pseudosets = any(
                isinstance(pseudoset, PseudoSet) for pseudoset in pseudos.values()
                )
            if contains_pseudosets:
                if not all(isinstance(pseudoset, PseudoSet) for pseudoset in pseudos.values()):
                    msg = "A pseudopotential-set mapping must contain only PseudoSet values"
                    raise TypeError(msg)
                code = PseudoSet._check_code_str(code)
                if code not in pseudos:
                    msg = f'pseudopotential set is not available for code "{code}"'
                    raise ValueError(msg)
                if system is None:
                    msg = "A system must be provided with a pseudopotential-set mapping"
                    raise ValueError(msg)
                pps = pseudos[code]
                species = system.structure.species(symbol=True)[1]
                missing = set(species)-pps.pseudos.keys()
                if missing:
                    msg = f'pseudopotential set does not contain species {sorted(missing)}'
                    raise ValueError(msg)
                return {
                    pps.pseudos[element].name: str(pps.pseudos[element])
                    for element in sorted(species)
                    }
        if isinstance(pseudos, Mapping):
            return dict(pseudos)
        if isinstance(pseudos,str):
            label = pseudos
            if label not in PseudoSet.labeled_pseudosets:
                msg = f'pseudopotential set label "{label}" is not registered'
                raise ValueError(msg)
            code = PseudoSet._check_code_str(code)
            if code not in PseudoSet.labeled_pseudosets[label]:
                msg = f'pseudopotential set "{label}" is not available for code "{code}"'
                raise ValueError(msg)
            pps = PseudoSet.labeled_pseudosets[label][code]
            species = system.structure.species(symbol=True)[1]
            missing = set(species)-pps.pseudos.keys()
            if missing:
                msg = f'pseudopotential set "{label}" does not contain species {sorted(missing)}'
                raise ValueError(msg)
            pseudos = [pps.pseudos[e].name for e in sorted(species)]
        missing = set(pseudos)-PseudoSet.pseudo_files.keys()
        if missing:
            msg = f'pseudopotential files {missing} are not present in PseudoSet.pseudo_files'
            raise ValueError(msg)
        return {f:PseudoSet.pseudo_files[f] for f in pseudos}
    #end def pseudo_remap


    def __init__(
        self,
        pseudos     : Iterable[PathLike] | Mapping[Elements | str, PathLike],
        codes       : str | Iterable[str] = "detect",
        Zeff_map    : Mapping[PathLike, int] | None = None,
        *,
        skip_invalid: bool = False,
        ):
        self.pseudos: dict[str, Path] = {}
        if isinstance(pseudos, Mapping | obj):
            for label, psp in pseudos.items():
                psp = Path(psp).resolve()
                if not psp.exists():
                    msg = f"Pseudo file {psp} can not be located"
                    if skip_invalid:
                        warn(msg)
                        continue
                    else:
                        raise FileNotFoundError(msg)
                else:
                    # No need to check if `label` is already defined since
                    # dictionary keys are, by definition, unique.
                    self.pseudos[label] = psp
        else:
            for psp in pseudos:
                psp = Path(psp).resolve()
                if not psp.exists():
                    msg = f"Pseudo file {psp} can not be located"
                    if skip_invalid:
                        warn(msg)
                        continue
                    else:
                        raise FileNotFoundError(msg)

                if psp.name.lower() == "potcar":
                    # POTCARS are stored all with the same name.
                    # The directory they are in is where the actual element is.
                    _, symbol, is_elem = pp_elem_label(psp.parent.name)
                else:
                    _, symbol, is_elem = pp_elem_label(psp.name)

                if not is_elem:
                    msg = (
                       f"Can not determine element for pseudopotential file: {psp}\n"
                        "Pseudopotential file names must be prefixed by an atomic symbol or label\n"
                        "(e.g. Si, Si1, etc)"
                        )
                    if skip_invalid:
                        warn(msg)
                    else:
                        raise ValueError(msg)
                elif symbol in self.pseudos:
                    msg = (
                        "Can not provide multiple pseudos for the same element!\n"
                       f"Duplicate pseudo is at index {pseudos.index(psp)}\n"
                       f"    Existing pseudo file:  {self.pseudos[symbol]}\n"
                       f"    Duplicate pseudo file: {psp}"
                        )
                    raise ValueError(msg)
                else:
                    self.pseudos[symbol] = psp

        if isinstance(codes, str):
            if codes.lower() == "detect":
                self.codes = self._detect_pseudo_code(self.pseudos)
            else:
                self.codes = {PseudoSet._check_code_str(codes)}
        elif isinstance(codes, Iterable):
            if not isinstance(next(iter(codes)), str):
                msg = (
                    "`codes` must be either 'detect', str, or an iterable of str, "
                    "but is an iterable of `{type(next(iter(codes))).__name__}`"
                    )
                raise TypeError(msg)
            self.codes = set([
                PseudoSet._check_code_str(code) for code in codes
                ])
        else:
            msg = f"`codes` must be either 'detect', str, or an iterable of str, but has type `{type(codes).__name__}`"
            raise TypeError(msg)

        self.Zeff_map = dict(Zeff_map) if Zeff_map is not None else {}

        self.pseudo_dirs: set[Path] = set()
        for pseudo in self.pseudos.values():
            if pseudo.name.lower() == "potcar":
                self.pseudo_dirs.add(pseudo.parent.parent) # Parent dir of the POTCAR dirs
            else:
                self.pseudo_dirs.add(pseudo.parent)
    #end def __init__

    @staticmethod
    def _detect_pseudo_code(
        pseudos: Mapping[str, PathLike] | Iterable[PathLike]
        ) -> set[Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"]]:
        """Detect the code based on the suffix of the pseudos."""
        codes = set()
        suffixes = set()
        if isinstance(pseudos, Mapping | obj):
            pseudos = pseudos.values()

        for pseudo in pseudos:
            pseudo = Path(pseudo)
            if pseudo.name.lower() == "potcar":
                suffixes.add(pseudo.name.lower())
            else:
                suffixes.add(pseudo.suffix.lower())

        if len(suffixes) == 0:
            msg = (
                "Can not detect code with no pseudopotentials!\n"
                "If you are initializing an empty PseudoSet, you must provide a code!"
                )
            raise RuntimeError(msg)

        for code_key in PseudoSet.file_exts:
            if suffixes.issubset(PseudoSet.file_exts[code_key]):
                codes.add(code_key)

        if len(codes) == 0:
            msg = (
                "Can not detect code from pseudopotential extensions!\n"
                f"Detected extensions are: {', '.join(suffixes)}"
                )
            raise RuntimeError(msg)

        return codes
    #end def _detect_pseudo_code

    @staticmethod
    def _check_code_str(code: str) -> Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"]:
        """Check to make sure a code string is in the set of known codes.

        Returns
        -------
        str
            Lowercased version of the code.
        """
        clow = code.lower()
        if clow == "pwscf": # Retain alias for backwards compatibility
            warn(
                "Automatically switching code 'pwscf' to 'espresso'.\n"
                "In the future using 'pwscf' will be deprecated, please use 'espresso' instead!"
                )
            clow = "espresso"

        if clow not in PseudoSet.known_codes:
            msg = (
                f"Code '{code}' is not known by Nexus!\n"
                f"Known codes are {list(PseudoSet.known_codes)}"
                )
            raise ValueError(msg)
        else:
            return clow
    #end def _check_code_str


    @classmethod
    def from_dir(
        cls,
        pseudo_dir: PathLike,
        code      : str = "detect",
        Zeff_map  : Mapping[PathLike, int] | None = None,
        ext_filter: str | Iterable[str] | None = None,
        pattern   : str | Pattern | None = None,
        *,
        skip_invalid: bool = False,
        ) -> PseudoSet:
        """Read in pseudopotentials from a directory.

        Parameters
        ----------
        pseudo_dir : PathLike
            The directory from which to read pseudopotentials.
            Does not support nested directories, except for those that
            contain a POTCAR file.
        code : {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf", "detect"}, default="detect"
            The name of the code that the pseudos are formatted for,
            or if ``"detect"``, will auto-detect the code name from the
            file extensions.
        Zeff_map : Map of str/Elements to int, optional
            A ``dict`` or ``obj`` mapping elements to their effective
            nuclear charges (Z-valences). If this is supplied, it will
            override any parts of the code that may try to parse the
            pseudopotential to get the Z-valence.
        ext_filter : str or list of str, optional
            Optionally filter the files in the directory by their
            extension.

            If this is ``None`` it will use the file suffixes
            in ``PseudoSet.file_exts``, unless ``codes="detect"``, in
            which case it will do nothing.

            If this is a string or list of strings, it is assumed the
            string(s) are the file suffixes to filter by. The strings
            should include a leading ``.``, e.g. ``.xml``, not ``xml``.
        pattern : str or Pattern, optional
            A string or regex pattern to use to filter files by name.
        skip_invalid : bool, default=False (keyword-only)
            If ``True``, then this will emit a warning rather than raise an
            error if a file is not found or if the file does not have a
            valid name.

        See Also
        --------
        file_exts : Dictionary mapping codes to file extensions.
        known_codes : Codes known by Nexus.

        Examples
        --------
        Reading pseudos in a directory with only one style of pseudo.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.xml', 'C.ccECP.xml', 'Fe.ccECP.xml']
        >>> psps = PseudoSet.from_dir(pseudo_dir)
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        H: /path/to/pseudo_dir/H.ccECP.xml
        C: /path/to/pseudo_dir/C.ccECP.xml
        Fe: /path/to/pseudo_dir/Fe.ccECP.xml

        Reading in only the UPF pseudos in a directory with UPF and XML
        pseudos.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.xml', 'C.ccECP.xml', 'H.ccECP.upf', 'C.ccECP.upf']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="espresso")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.ccECP.upf
        H: /path/to/pseudo_dir/H.ccECP.upf

        Filtering out two different kinds of pseudos with the same
        extensions.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.upf', 'C.ccECP.upf', 'H.USPP.upf', 'C.USPP.upf']
        >>> psps = PseudoSet.from_dir(pseudo_dir, pattern="USPP")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.USPP.upf
        H: /path/to/pseudo_dir/H.USPP.upf

        Filtering out pseudos by both extension and pattern.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.upf', 'C.ccECP.upf', 'H.USPP.upf', 'C.USPP.upf', 'H.ccECP.xml', 'C.ccECP.xml']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="espresso", pattern="ccECP")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.ccECP.upf
        H: /path/to/pseudo_dir/H.ccECP.upf

        Filtering out VASP pseudos with similar names. Pattern matches
        anything *without* an underscore.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", pattern=r"^((?!_).)*$")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C/POTCAR
        H: /path/to/pseudo_dir/H/POTCAR

        Filtering out VASP pseudos ending with ``_sv``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", pattern=r"_sv$",)
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_sv/POTCAR
        H: /path/to/pseudo_dir/H_sv/POTCAR

        Filtering out VASP pseudos ending with ``_GW``, but do not
        contain ``_sv``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", pattern=r"(?<!sv)_GW",)
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_GW/POTCAR
        H: /path/to/pseudo_dir/H_GW/POTCAR

        Filtering out VASP pseudos ending with ``_sv_GW``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", pattern=r"_sv_GW",)
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_sv_GW/POTCAR
        H: /path/to/pseudo_dir/H_sv_GW/POTCAR
        """
        if code != "detect":
            code = PseudoSet._check_code_str(code)

        if ext_filter is None:
            if code != "detect":
                ext_filter = PseudoSet.file_exts[code]
        elif isinstance(ext_filter, str):
            ext_filter = {ext_filter.lower()}
        elif isinstance(ext_filter, Iterable):
            if not isinstance(next(iter(ext_filter)), str):
                msg = f"`ext_filter` must be either None, str, or an iterable of str, but is {type(ext_filter[0]).__name__}"
                raise TypeError(msg)

            ext_filter = set([ext.lower() for ext in ext_filter])
        else:
            msg = f"`ext_filter` must be either None, str, or an iterable of str, but is {type(ext_filter).__name__}"
            raise TypeError(msg)

        if pattern is not None:
            pattern = re.compile(pattern)

        psp_dir = Path(pseudo_dir).resolve()

        if not psp_dir.exists():
            msg = f"Can not find pseudopotential directory: {psp_dir}"
            raise FileNotFoundError(msg)
        elif not psp_dir.is_dir():
            msg = f"Specified path does not point to a directory: {psp_dir}"
            raise NotADirectoryError(msg)

        pseudos = []
        for pseudo in psp_dir.iterdir():
            if pseudo.is_file() and (ext_filter is None or pseudo.suffix.lower() in ext_filter):
                if pattern is None or pattern.search(pseudo.name) is not None:
                    pseudos.append(pseudo)
            elif pseudo.is_dir() and (ext_filter is None or "potcar" in ext_filter):
                if pattern is None or pattern.search(pseudo.name) is not None:
                    potcar_upper = pseudo / "POTCAR"
                    potcar_lower = pseudo / "potcar"
                    if potcar_upper.exists():
                        pseudos.append(potcar_upper)
                    elif potcar_lower.exists():
                        pseudos.append(potcar_lower)
                else:
                   continue

        if code == "detect":
            code = cls._detect_pseudo_code(pseudos)

        return cls(
            codes        = code,
            pseudos      = pseudos,
            Zeff_map     = Zeff_map,
            skip_invalid = skip_invalid,
            )
    #end def from_dir


    @classmethod
    def from_mixed_dir(
        cls,
        pseudo_dir   : PathLike,
        codes        : str | list[str] | None = None,
        extensions   : Mapping[str, set[str]] | None = None,
        patterns     : Mapping[str, str | Pattern] | None = None,
        code_Zeff_map: Mapping[str, Mapping[str, int]] | None = None,
        *,
        skip_invalid: bool = False,
        ) -> dict[Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"], PseudoSet]:
        """Read in pseudos from a directory with pseudos for more than one code.

        Parameters
        ----------
        pseudo_dir : PathLike
            The directory from which to read pseudopotentials.
            Does not support nested directories, except for those that
            contain a POTCAR file.
        codes : one or more of {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"}, optional
            The code(s) to use to separate the files in the directory
            into their respective groups. If this is set to ``detect``
            and filters is ``None``, then it will use all known codes
            and file extensions to filter the pseudos.
        extensions : Map of str to set of str, optional
            A dictionary mapping codes to the file extensions
            corresponding to those labels. If this is not provided, then
            the filters are automatically populated by the codes in
            ``codes``.
        patterns : Map of str to str or Pattern, optional
            A dictionary mapping codes to strings or regex patterns,
            used to filter out files.
        code_Zeff_map : Map of str to Map of str/Elements to int, optional
            A ``dict`` or ``obj`` for each code that maps elements to
            their effective nuclear charges (Z-valences). If this is
            supplied, it will override any parts of the code that may
            try to parse the pseudopotential to get the Z-valence.
        skip_invalid : bool, default=False (keyword-only)
            If ``True``, then this will emit a warning rather than raise
            an error if a file is not found or if the file does not have
            a valid name.

        Returns
        -------
        pseudos : dict of str: PseudoSet/None
            A map from the labels provided to the function to the
            ``PseudoSet`` objects that were created from the pseudos in
            the directory.

        Notes
        -----
        This function is the most generous in terms of what it is able
        to parse and separate, however it can result in errors if you
        have multiple pseudos with the same file extension for the same
        element. If you have a lot of pseudos with the same extension,
        you should use ``PseudoSet.from_dir()`` and provide a pattern to
        separate the pseudos.

        See Also
        --------
        from_dir : Used to get pseudos after filters have been established.
        file_exts : Dictionary mapping codes to file extensions.
        known_codes : Codes known by Nexus.

        Examples
        --------
        Filter a large collection of pseudos for multiple codes.

        >>> print(contents_of_pseudo_dir)
        pseudo_dir
        ├── C
        │   └── POTCAR
        ├── C.BFD.gms
        ├── C.BFD.gth
        ├── C.ccECP.upf
        ├── C.ccECP.xml
        ├── C.USPP.upf
        ├── H
        │   └── POTCAR
        ├── H.BFD.gms
        ├── H.BFD.gth
        ├── H.ccECP.upf
        ├── H.ccECP.xml
        └── H.USPP.upf
        >>> psps = PseudoSet.from_mixed_dir(
        ...     pseudo_dir=pseudo_dir,
        ...     patterns={"espresso": "ccECP", "rmg": "USPP"}
        ... )
        >>> for code, ps_set in psps.items():
        ...     print(f"{code} pseudos:")
        ...     for lbl, psp in ps_set.pseudos.items():
        ...         print(f"  {lbl}: {psp}")
        espresso pseudos:
        H: /path/to/pseudo_dir/H.ccECP.upf
        C: /path/to/pseudo_dir/C.ccECP.upf
        gamess pseudos:
        C: /path/to/pseudo_dir/C.BFD.gms
        H: /path/to/pseudo_dir/H.BFD.gms
        vasp pseudos:
        C: /path/to/pseudo_dir/C/POTCAR
        H: /path/to/pseudo_dir/H/POTCAR
        qmcpack pseudos:
        C: /path/to/pseudo_dir/C.ccECP.xml
        H: /path/to/pseudo_dir/H.ccECP.xml
        rmg pseudos:
        C: /path/to/pseudo_dir/C.USPP.upf
        H: /path/to/pseudo_dir/H.USPP.upf
        pyscf pseudos:
        C: /path/to/pseudo_dir/C.BFD.gth
        H: /path/to/pseudo_dir/H.BFD.gth
        """
        psp_dir = Path(pseudo_dir).resolve()

        if not psp_dir.exists():
            msg = f"Can not find pseudopotential directory: {psp_dir}"
            raise FileNotFoundError(msg)
        elif not psp_dir.is_dir():
            msg = f"Specified path does not point to a directory: {psp_dir}"
            raise NotADirectoryError(msg)

        if extensions is None:
            if codes is None:
                extensions = PseudoSet.file_exts
            else:
                if isinstance(codes, str):
                    codes = {codes}
                extensions = {}
                for c in codes:
                    code = PseudoSet._check_code_str(c)
                    extensions[code] = PseudoSet.file_exts[code]
        else:
            if codes is None:
                if len(extensions) < len(PseudoSet.known_codes):
                    # codes is None, so we want all codes. Make sure any codes with
                    # unspecified filters are added with their defaults.
                    for code in PseudoSet.known_codes - set(extensions.keys()):
                        extensions[code] = PseudoSet.file_exts[code]

                checked_filters = {}
                for code, suffixes in extensions.items():
                    code = PseudoSet._check_code_str(code)
                    checked_filters[code] = set(suffixes) if suffixes is not None else None

                extensions = checked_filters
            else:
                if isinstance(codes, str):
                    codes = {codes}
                # More filters than codes, or filters provided for unselected codes
                if not set(codes) >= set(extensions.keys()):
                    msg = (
                        "Mismatch between provided filters and codes!\n"
                        f"Provided codes: {tuple(codes)}\n"
                        f"Filter keys:    {tuple(extensions.keys())}"
                        )
                    raise ValueError(msg)
                else:
                    checked_filters = {}
                    for c in codes:
                        code = PseudoSet._check_code_str(c)
                        checked_filters[code] = extensions.get(c)

        if code_Zeff_map is None:
            code_Zeff_map = {}
        elif codes is not None and not set(codes) >= set(code_Zeff_map.keys()):
            msg = (
                "Mismatch between provided code Zeff map and codes!\n"
                f"Provided codes: {tuple(codes)}\n"
                f"Zeff keys:      {tuple(code_Zeff_map.keys())}"
                )
            raise ValueError(msg)

        if patterns is None:
            patterns = {}
        elif codes is not None and not set(codes) >= set(patterns.keys()):
            msg = (
                "Mismatch between provided patterns and codes!\n"
                f"Provided codes: {tuple(codes)}\n"
                f"Pattern keys:   {tuple(patterns.keys())}"
                )
            raise ValueError(msg)

        pseudos = {}
        for code, suffixes in extensions.items():
            try:
                pseudos[code] = PseudoSet.from_dir(
                    pseudo_dir   = psp_dir,
                    code         = code,
                    Zeff_map     = code_Zeff_map.get(code),
                    ext_filter   = suffixes,
                    pattern      = patterns.get(code),
                    skip_invalid = skip_invalid,
                    )
            except ValueError as err:
                msg = format(err) + (
                    f"\n\nDuplicate element detected for code '{code}'\n"
                    f"Either remove '{code}' from the selected codes, or specify "
                    "`filters` and/or `patterns` to ensure the collision does not happen.\n"
                    )
                # Raise from None to prevent exception chain.
                raise RuntimeError(msg) from None

        return pseudos
    #end def from_mixed_dir


    def get_pseudos(
        self,
        system: PhysicalSystem | Iterable[str],
        code  : Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"] | None = None,
        ) -> set[Path] | None:
        """Get the pseudopotential files for the elements in a physical system.

        Parameters
        ----------
        system : PhysicalSystem or list of str
            The system to get pseudopotentials for, or a list of element
            labels.
        code : {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"}, optional
            The name of the code requesting the pseudopotentials.
            If supplied, it will raise an error if the code requested
            does not match the code of the ``PseudoSet``. This provides
            a way to ensure the user does not provide the wrong pseudos
            to a ``generate`` call.

        Returns
        -------
        pseudos : set of Path
            The pseudopotential paths for the given system of elements.
        """
        if code is not None:
            clow = PseudoSet._check_code_str(code)
            if clow not in self.codes:
                msg = f"Tried to get pseudopotentials for {code} from a set of {'/'.join(self.codes)} pseudos!"
                raise ValueError(msg)

        pps = set()
        if isinstance(system, PhysicalSystem):
            if not system.pseudized:
                return None
            else:
                elements = system.ion_labels
        else:
            elements = system

        for label in elements:
            if label not in self.pseudos:
                msg = f"No pseudopotential found for label {label}!"
                raise ValueError(msg)
            pps.add(self.pseudos[label])

        return pps
    #end def get_pseudos


    def get_Zeff(
        self,
        elem_labels: Iterable[Elements | str] | PhysicalSystem,
        *,
        missing_as_ae: bool = False,
        ) -> dict[str, int]:
        """Get the Z-valences for each element in the list of elements.

        Parameters
        ----------
        elem_labels : list of Elements or list of str or PhysicalSystem
            The elements or system to get Z-valences for.
        missing_as_ae : bool, default=False (keyword-only)
            Assume any elements for which a pseudopotential can not be
            found are all-electron, and use their atomic number as the
            value for the Z-valence. If this is not supplied, and if the
            Z-valence for an element is not in ``self.Zeff``, this will
            attempt to extract the Z-valence from the pseudopotential
            file.

        Returns
        -------
        elem_Zeff : dict of str: int
            A dictionary mapping element labels to their effective
            nuclear charges.

        See Also
        --------
        read_upf_z_valence : Used to extract Z-valences from UPF files.
        read_xml_z_valence : Used to extract Z-valences from XML files.
        read_potcar_z_valence : Used to extract Z-valences from POTCAR files.
        """
        if isinstance(elem_labels, PhysicalSystem):
            elem_labels = elem_labels.ion_labels

        elem_labels = set(elem_labels) # Unique only, saves on iteration.
        Z_eff_map = {}
        for label in elem_labels:
            if label in self.Zeff_map:
                Z_eff_map[label] = self.Zeff_map[label]
            elif label in self.pseudos:
                f_ext = self.pseudos[label].suffix.lower()
                if f_ext == ".upf":
                    Z_eff_map[label] = read_upf_z_valence(self.pseudos[label])
                elif f_ext in (".gms", ".gamess"):
                    msg = (
                        "Z-valence parsing not implemented for GAMESS pseudopotentials!\n"
                        "You must supply Z-valences manually until this feature is added."
                        )
                    raise NotImplementedError(msg)
                elif f_ext in ["potcar", ".vasp"]:
                    Z_eff_map[label] = read_potcar_z_valence(self.pseudos[label])
                elif f_ext == ".xml":
                    Z_eff_map[label] = read_qmcpack_xml_z_valence(self.pseudos[label])
                else:
                    msg = f"File extension '{f_ext}' is not parseable by Nexus, can not extract Z-valence!"
                    raise NotImplementedError(msg)
            elif missing_as_ae:
                is_elem, element = Elements.is_element(label, return_element=True)
                if not is_elem:
                    msg = f"Can not determine element for label '{label}'"
                    raise ValueError(msg)
                else:
                    Z_eff_map[label] = element.atomic_number
            else:
                msg = f"No pseudopotential found for label {label}!"
                raise ValueError(msg)

        return Z_eff_map
    #end def get_Zeff


    def __repr__(self) -> str:
        rep = (
            "PseudoSet(\n"
           f"    codes = {self.codes},\n"
            "    pseudos = {"
            )
        if len(self.pseudos) > 0:
            rep += "\n"
            lbl_len = max(map(len, self.pseudos.keys()))+2 # Single quotes around label add 2
            for lbl, pth in self.pseudos.items():
                lbl = f"'{lbl}'"
                rep += f"{' '*8}{lbl:<{lbl_len}}: {pth!r},\n"
            rep += "    },\n"
        else:
            rep += "},\n"

        rep += "    Zeff_map = {"
        if len(self.Zeff_map) > 0:
            lbl_len = max(map(len, self.Zeff_map.keys()))+2
            rep += "\n"
            for lbl, zeff in self.Zeff_map.items():
                lbl = f"'{lbl}'"
                rep += f"{' '*8}{lbl:<{lbl_len}}: {zeff!s},\n"
            rep += "    },\n"
        else:
            rep += "},\n"
        rep += ")\n"
        return rep
    #end def __repr__
#end class PseudoSet