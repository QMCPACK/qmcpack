"""Functionality for gathering, parsing, and validating collections of pseudopotentials."""


from __future__ import annotations

import re
from collections.abc import Collection, Mapping
from copy import deepcopy
from fnmatch import fnmatchcase
from os import PathLike
from pathlib import Path
from types import MappingProxyType
from typing import ClassVar, Literal

from .developer import DevBase, warn
from .generic import nxs_deprecate
from .periodic_table import Elements
from .physical_system import PhysicalSystem
from .utilities import is_valid_filename
from .nexus_base import nexus_core


def pp_elem_label(
    filename: PathLike, *, guard=False,
    ) -> tuple[str, str] | tuple[str, str, bool]:
    """Get the label and atomic symbol of an element from a pseudopotential file name.

    Parameters
    ----------
    filename : PathLike
        Either a :class:`Path` object or the name of a file as a string.
    guard : bool, default=False
        Optionally raise an error if the file name is invalid or if the filename
        does not begin with a valid element prefix.

    Returns
    -------
    elem_label : str
        The label of the pseudopotential file.
    symbol : str
        The atomic symbol of the element that the pseudopotential file describes.
    is_elem : bool
        A bool to indicate whether or not the file has a valid element label.
        If ``guard=True`` this is not returned, and an error is raised if the
        file does not have a valid prefix.

    Examples
    --------
    >>> pp_elem_label("Si.ccECP.xml")
    ('Si', 'Si', True)
    >>> pp_elem_label("Si-USPP.upf")
    ('Si', 'Si', True)
    >>> pp_elem_label("Si_NCPP.upf")
    ('Si', 'Si', True)
    >>> pp_elem_label("Si1.ccECP.xml")
    ('Si1', 'Si', True)
    >>> pp_elem_label("Not_An_Element.ccECP.xml")
    ('Not', 'Not', False)

    Passing a ``Path`` object only gets the file name.

    >>> from pathlib import Path
    >>> pp_elem_label(Path("/path/to/pseudo_dir/Si.ccECP.xml"))
    ('Si', 'Si', True)

    Use ``guard=True`` if you want to raise an error and not deal with an
    invalid output.

    >>> pp_elem_label("Not_An_Element.ccECP.xml", guard=True)
    Traceback (most recent call last):
        ...
    RuntimeError: cannot determine element for pseudopotential file: Not_An_Element.ccECP.xml
    pseudopotential file names must be prefixed by an atomic symbol or label
    (e.g. Si, Si1, etc)
    """
    if isinstance(filename, Path):
        filename = filename.name

    if guard and not is_valid_filename(filename):
        msg = f"Pseudopotential file name {filename} is invalid!"
        raise RuntimeError(msg)

    el = ''
    for c in filename:
        if c in {".", "_", "-"}:
            break
        #end if
        el+=c
    #end for
    elem_label = el
    is_elem, element = Elements.is_element(el, return_element=True)
    if guard:
        if not is_elem:
            msg = (
               f"Cannot determine element for pseudopotential file: {filename}\n"
                "Pseudopotential file names must be prefixed by an atomic symbol or label\n"
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
        msg = (
           f"Could not find Z valence in file: {file!s}\n"
            "You may need to provide the Z valence manually!"
            )
        raise RuntimeError(msg)
    else:
        zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        msg = (
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
        raise RuntimeError(msg)
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
        msg = (
           f"Could not find Z valence in file: {file!s}\n"
            "You may need to provide the Z valence manually!"
            )
        raise RuntimeError(msg)
    else:
        zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        msg = (
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
        raise RuntimeError(msg)
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
                msg = (
                   f"Could not find Z valence in file: {file!s}\n"
                    "You may need to provide the Z valence manually!"
                    )
                raise RuntimeError(msg)
            else:
                zval = float(zval.group(1).lower().replace("d", "e"))

    if zval <= 0 or zval > 118:
        msg = (
            f"Invalid Z-valence found in file, must be in range (0, 118], but is {zval}!"
            )
        raise RuntimeError(msg)
    # Round to 8 digits
    if round(zval, 8).is_integer():
        return int(zval)
    else:
        return zval
#end def read_potcar_z_valence



@nxs_deprecate(since="2.4.0", replacement="generate_pseudoset (https://nexus-workflows.readthedocs.io/en/latest/user_guide/pseudo-handling.html#migrating-from-ppset)")
def ppset(label: str, **codes_pps: Collection[str]):
    """Register pseudopotentials for codes with a label.

    .. deprecated:: 2.4.0
        :func:`ppset` has been replaced by :func:`generate_pseudoset` because
        the labeling system that :func:`ppset` requires checks for the existence
        of the label at runtime, whereas :func:`generate_pseudoset` returns an
        object with a name, and if that name is misspelled then Python will not
        execute and can provide a better diagnostic than Nexus can.

    This is intended as a backwards-compatible interface to not break existing
    user code. Users are suggested to migrate to :func:`generate_pseudoset` or
    calling :class:`PseudoSet` and its methods for pseudopotential handling.

    Parameters
    ----------
    label : str
        Label of the pseudopotential set.
    **codes_pps
        Map of codes to collections of names of the pseudopotential files.

    Notes
    -----
    Using this function is generally discouraged for users as the use of a
    string label for identifying pseudopotentials prevents Python from raising a
    :exc:`NameError` when the label supplied to later functions does not match
    the label used when calling this function.

    For this reason, it is strongly encouraged to use :func:`generate_pseudoset`
    or :class:`PseudoSet` and/or its member functions :meth:`PseudoSet.from_dir`
    and :meth:`PseudoSet.from_mixed_dir`.

    See Also
    --------
    PseudoSet.file_exts : Known file extensions
    PseudoSet.from_dir : Automatically reading pseudos from a directory
    PseudoSet.from_mixed_dir : Automatically reading pseudos from a directory

    Examples
    --------
    After passing ``pseudo_dir`` to :data:`~nexus.settings`, call this with the
    file names for each pseudopotential.

    >>> ppset(
    ...     label   = 'bfd',
    ...     pwscf   = ["C.BFD.upf", "H.BFD.upf", "O.BFD.upf"],
    ...     qmcpack = ["C.BFD.xml", "H.BFD.xml", "O.BFD.xml"],
    ...     )
    """
    if label in PseudoSet.labeled_pseudosets:
        msg = f'Pseudopotential set label "{label}" is already registered!'
        raise ValueError(msg)

    pps_coll = {}
    ref_elements = None
    for code, ppfiles in codes_pps.items():
        missing = set(ppfiles) - PseudoSet.pseudo_files.keys()
        if len(missing)>0:
            msg = f'Pseudopotential files "{missing}" are not present in PseudoSet.pseudo_files!'
            raise ValueError(msg)

        pps = PseudoSet(ppfiles)
        code = PseudoSet._check_code_str(code)
        if code not in pps.codes:
            msg = f'Pseudopotential files provided for code "{code}" are not compatible with that code!'
            raise ValueError(msg)

        elements = set(pps.pseudos)
        if ref_elements is None:
            ref_elements = elements
        elif elements!=ref_elements:
            msg = f'Pseudopotential set "{label}" must contain potentials for the same elements for each code!'
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
        A mapping of elements to their effective nuclear charges (Z-valences).
    pseudo_dirs : set of Path
        The directories that the pseudopotentials are stored in.
    pseudo_files : dict of str: str (class attribute)
        Dictionary mapping pseudopotential file names to their full paths, set
        by passing ``pseudo_dir`` to :data:`~nexus.settings`.
    labeled_pseudosets : dict of str: dict of str: PseudoSet (class attribute)
        Pseudopotential sets registered by label, then compatible code.

    Parameters
    ----------
    pseudos : collection of str/Path or map of str/Elements to Path
        A collection of pseudopotential files or a map of elements to file paths.
    codes : one or more of {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf", "detect"}, default="detect"
        The name of the code that the pseudos are formatted for, or
        if ``"detect"``, will auto-detect the code name from the
        file extensions.
    Zeff_map : Map of str/Elements to int, optional
        A mapping of elements to their effective nuclear charges (Z-valences).
        If this is supplied, it will override any parts of the code that may try
        to parse the pseudopotential file to get the Z-valence.
    skip_invalid : bool, default=False (keyword-only)
        If ``True``, then this will emit a warning rather than raise an
        error if a file is not found or if the file does not have a
        valid name.
    """

    file_exts = MappingProxyType({
        # https://www.quantum-espresso.org/Doc/INPUT_PW.html#id268
        "espresso": frozenset({".ncpp", ".upf", ".vdb", ".van", ".rrkj3"}),
        "gamess":   frozenset({".gms", ".gamess"}),
        "vasp":     frozenset({"potcar", ".potcar", ".vasp"}),
        "qmcpack":  frozenset({".xml", ".data"}), # .data is for CASINO format
        "rmg":      frozenset({".upf", ".xml"}),
        "pyscf":    frozenset({".nwchem", ".gth"})
        })
    known_codes = frozenset(file_exts.keys())

    code_aliases = MappingProxyType({
        "quantum_espresso": "espresso",
        "quantum espresso": "espresso",
        "pwscf": "espresso",
        "qe": "espresso",
        "gms": "gamess",
    })

    pseudo_files      : ClassVar[dict[str, str]] = {}
    labeled_pseudosets: ClassVar[dict[str, dict[str, PseudoSet]]] = {}

    def __init__(
        self,
        pseudos : Collection[PathLike] | Mapping[Elements | str, PathLike],
        codes   : str | Collection[str] = "detect",
        Zeff_map: Mapping[PathLike, int] | None = None,
        *,
        skip_invalid: bool = False,
        ):
        self.pseudos: dict[str, Path] = {}
        if isinstance(pseudos, Mapping):
            for label, psp in pseudos.items():
                psp = Path(psp).resolve()
                if psp.exists():
                    # No need to check if `label` is already defined since
                    # dictionary keys are, by definition, unique.
                    self.pseudos[label] = psp
                elif psp.name in PseudoSet.pseudo_files:
                    self.pseudos[label] = Path(PseudoSet.pseudo_files[psp.name]).resolve()
                else:
                    msg = f"Pseudo file {psp} can not be located"
                    if skip_invalid:
                        warn(msg)
                        continue
                    else:
                        raise FileNotFoundError(msg)
                #end if
        else:
            for psp in pseudos:
                psp = Path(psp).resolve()
                if psp.exists():
                    pass
                elif psp.name in PseudoSet.pseudo_files:
                    psp = Path(PseudoSet.pseudo_files[psp.name]).resolve()
                else:
                    msg = f"Pseudo file {psp} can not be located"
                    if skip_invalid:
                        warn(msg)
                        continue
                    else:
                        raise FileNotFoundError(msg)
                #end if

                if psp.name.lower() == "potcar":
                    if not psp.is_file():
                        msg = (
                            "POTCARs can not be directories!\n"
                            f"{psp}"
                            )
                        raise NotADirectoryError(msg)
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
        elif isinstance(codes, Collection):
            if not isinstance(next(iter(codes)), str):
                msg = (
                    "`codes` must be either 'detect', str, or an collection of str, "
                    f"but is a collection of `{type(next(iter(codes))).__name__}`"
                    )
                raise TypeError(msg)
            self.codes = {PseudoSet._check_code_str(code) for code in codes}
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
        pseudos: Mapping[str, PathLike] | Collection[PathLike]
        ) -> set[Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"]]:
        """Detect the code based on the suffix of the pseudos."""
        codes = set()
        suffixes = set()
        if isinstance(pseudos, Mapping):
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
        if clow in PseudoSet.code_aliases:
            # TODO: Make this a logging call in the future.
            # new_code = PseudoSet.code_aliases[clow]
            # warn(f"Automatically switching code '{clow}' to '{new_code}'.")
            # clow = new_code
            clow = PseudoSet.code_aliases[clow]

        if clow not in PseudoSet.known_codes:
            msg = (
                f"Code '{code}' is not known by Nexus!\n"
                f"Known codes are {list(PseudoSet.known_codes)}"
                )
            raise ValueError(msg)
        else:
            return clow
    #end def _check_code_str


    @staticmethod
    def _normalize_code_map_keys(mapping: Mapping) -> dict:
        """Take a dict with any code keys and normalize them.

        Normalizing in this case means checking the keys and making sure they
        match those in :attr:`PseudoSet.file_exts`.
        """
        mapping = deepcopy(dict(mapping))
        normalized_mapping = {}
        for code in mapping:
            # Normalize all codes in extensions
            checked_code = PseudoSet._check_code_str(code)
            if checked_code in normalized_mapping:
                msg = f"Dictionary supplied two aliases for the same code, found duplicate: '{code}'"
                raise ValueError(msg)

            normalized_mapping[checked_code] = mapping[code]

        return normalized_mapping
    #end def _normalize_code_map_keys


    @classmethod
    def from_dir(
        cls,
        pseudo_dir: PathLike,
        code      : str = "detect",
        extension : str | Collection[str] | None = None,
        include   : str | None = None,
        exclude   : str | None = None,
        Zeff_map  : Mapping[PathLike, int] | None = None,
        *,
        skip_invalid: bool = False,
        ) -> PseudoSet:
        """Read in pseudopotentials from a directory.

        Parameters
        ----------
        pseudo_dir : PathLike
            The directory from which to read pseudopotentials. Does not support
            nested directories, except for those that contain a POTCAR file.
        code : {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf", "detect"}, default="detect"
            The name of the code that the pseudos are formatted for, or if
            ``"detect"``, will auto-detect the code name from the file
            extensions.
        extension : str or list of str, optional
            Optionally filter the files in the directory by their extension.

            If this is ``None`` it will use the file suffixes in
            :attr:`PseudoSet.file_exts`, unless ``codes="detect"``, in which
            case it will do nothing.

            If this is a string or list of strings, it is assumed the string(s)
            are the file suffixes to filter by. The strings should include a
            leading ``.``, e.g. ``.xml``, not ``xml``.
        include : str, optional
            A Unix shell-style wildcard. All files matched by this are included
            in the final :class:`PseudoSet`. By default this matches all files.
        exclude : str, optional
            A Unix shell-style wildcard. All files matched by this are excluded
            in the final :class:`PseudoSet`. By default this matches nothing.
        Zeff_map : Map of str/Elements to int, optional
            A mapping of elements to their effective nuclear charges
            (Z-valences). If this is supplied, it will override any parts of the
            code that may try to parse the pseudopotential to get the Z-valence.
        skip_invalid : bool, default=False (keyword-only)
            If ``True``, then this will emit a warning rather than raise an
            error if a file is not found or if the file does not have a valid
            name.

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

        Reading in only the UPF pseudos in a directory with UPF and XML pseudos.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.xml', 'C.ccECP.xml', 'H.ccECP.upf', 'C.ccECP.upf']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="espresso")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.ccECP.upf
        H: /path/to/pseudo_dir/H.ccECP.upf

        Filtering two different kinds of pseudos with the same extensions.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.upf', 'C.ccECP.upf', 'H.USPP.upf', 'C.USPP.upf']
        >>> psps = PseudoSet.from_dir(pseudo_dir, include="*USPP*")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.USPP.upf
        H: /path/to/pseudo_dir/H.USPP.upf

        Filtering pseudos by both extension and pattern.

        >>> os.listdir(pseudo_dir)
        ['H.ccECP.upf', 'C.ccECP.upf', 'H.USPP.upf', 'C.USPP.upf', 'H.ccECP.xml', 'C.ccECP.xml']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="espresso", include="*ccECP*")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C.ccECP.upf
        H: /path/to/pseudo_dir/H.ccECP.upf

        Filtering out VASP pseudos with similar names. Pattern matches anything
        *without* an underscore.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", exclude="*_*")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C/POTCAR
        H: /path/to/pseudo_dir/H/POTCAR

        Including VASP pseudos ending with ``_sv``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", include="*_sv")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_sv/POTCAR
        H: /path/to/pseudo_dir/H_sv/POTCAR

        Including VASP pseudos ending with ``_GW``, but do not
        contain ``_sv``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", include="*_GW", exclude="*_sv*")
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_GW/POTCAR
        H: /path/to/pseudo_dir/H_GW/POTCAR

        Including only VASP pseudos ending with ``_sv_GW``.

        >>> os.listdir(pseudo_dir)
        ['H_sv_GW', 'C', 'C_GW', 'H_sv', 'H_GW', 'C_sv_GW', 'C_sv', 'H']
        >>> psps = PseudoSet.from_dir(pseudo_dir, code="vasp", include="*_sv_GW",)
        >>> for lbl, pth in psps.pseudos.items(): print(f"{lbl}: {pth}")
        C: /path/to/pseudo_dir/C_sv_GW/POTCAR
        H: /path/to/pseudo_dir/H_sv_GW/POTCAR
        """
        if code.lower() != "detect":
            code = PseudoSet._check_code_str(code)

        if extension is None:
            if code.lower() != "detect":
                extension = PseudoSet.file_exts[code]
        elif isinstance(extension, str):
            extension = {extension.lower()}
        elif isinstance(extension, Collection):
            if not isinstance(next(iter(extension)), str):
                msg = f"`extension` must be either None, str, or a collection of str, but is {type(next(iter(extension))).__name__}"
                raise TypeError(msg)

            extension = set([ext.lower() for ext in extension])
        else:
            msg = f"`extension` must be either None, str, or an iterable of str, but is {type(extension).__name__}"
            raise TypeError(msg)

        psp_dir = Path(pseudo_dir).resolve()

        if not psp_dir.exists():
            msg = f"Can not find pseudopotential directory: {psp_dir}"
            raise FileNotFoundError(msg)
        elif not psp_dir.is_dir():
            msg = f"Specified path does not point to a directory: {psp_dir}"
            raise NotADirectoryError(msg)

        include = include if include is not None else "*"  # "*" matches everything
        exclude = exclude if exclude is not None else "[]" # "[]" matches nothing

        pseudos = []
        for pseudo in sorted(psp_dir.iterdir()):
            if (
                pseudo.is_file()
                and ( # file extension is valid
                    extension is None
                    or pseudo.suffix.lower() in extension
                    )
                and ( # matched by the include pattern, not matched by exclude
                    fnmatchcase(pseudo.stem, include)
                    and not fnmatchcase(pseudo.stem, exclude)
                    )
                ):
                    pseudos.append(pseudo)
            elif (
                pseudo.is_dir()
                and ( # if we are looking for POTCAR files
                    extension is None
                    or "potcar" in extension
                    )
                and ( # checks directory name
                    fnmatchcase(pseudo.name, include)
                    and not fnmatchcase(pseudo.name, exclude)
                    )
                ):
                    potcar_upper = pseudo / "POTCAR"
                    potcar_lower = pseudo / "potcar"
                    if potcar_upper.exists():
                        pseudos.append(potcar_upper)
                    elif potcar_lower.exists():
                        pseudos.append(potcar_lower)

        if code.lower() == "detect":
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
        codes        : str | Collection[str] | None = None,
        extensions   : Mapping[str, set[str]] | str | None = None,
        include      : Mapping[str, str] | str | None = None,
        exclude      : Mapping[str, str] | str | None = None,
        code_Zeff_map: Mapping[str, Mapping[str, int]] | Mapping[str, int] | None = None,
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
            into their respective groups.
        extensions : Map of codes to set of str or a single str, optional
            A dictionary mapping codes to the file extensions corresponding to
            those labels, or a single string that will be applied to all codes.
            If this is not provided, then the filters are automatically
            populated by the codes in ``codes``.
        include : Map of codes to str or single str, optional
            A dictionary mapping codes to a Unix shell-style wildcard.
            All files matched by this are **included** in the final
            :class:`PseudoSet`. By default this matches all files.
        exclude : Map of codes to str or single str, optional
            A dictionary mapping codes to a Unix shell-style wildcard.
            All files matched by this are **excluded** in the final
            :class:`PseudoSet`. By default this matches all files.
        code_Zeff_map : Map of codes to Map of str/Elements to int or Map of str/Elements to int, optional
            A mapping for each code that maps elements to their effective
            nuclear charges (Z-valences), or a single mapping that will be
            applied to all codes. If this is supplied, it will override any
            parts of the code that may try to parse the pseudopotential to get
            the Z-valence.
        skip_invalid : bool, default=False (keyword-only)
            If ``True``, then this will emit a warning rather than raise
            an error if a file is not found or if the file does not have
            a valid name.

        Returns
        -------
        pseudos : dict of str: PseudoSet
            A map from the labels provided to the function to the
            :class:`PseudoSet` objects that were created from the pseudos in the
            directory.

        Notes
        -----
        This function is the most generous in terms of what it is able to parse
        and separate, however it can result in errors if you have multiple
        pseudos with the same file extension for the same element. If you have a
        lot of pseudos with the same extension, you should use
        :meth:`PseudoSet.from_dir()` and provide a pattern to separate the
        pseudos.

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

        if isinstance(codes, str):
            codes = {PseudoSet._check_code_str(codes)}
        elif isinstance(codes, Collection):
            codes = {PseudoSet._check_code_str(code) for code in codes}
        elif codes is not None:
            msg = f"`codes` must be either a string or a collection of strings or None, but is {type(codes).__name__}"
            raise TypeError(msg)

        if extensions is None:
            if codes is None:
                extensions = PseudoSet.file_exts
            else:
                extensions = {}
                for c in codes:
                    code = PseudoSet._check_code_str(c)
                    extensions[code] = PseudoSet.file_exts[code]
        elif isinstance(extensions, str):
            # Single extension for all codes
            if codes is None:
                extensions = {code: extensions for code in PseudoSet.known_codes}
            else:
                extensions = {code: extensions for code in codes}
        else:
            extensions = PseudoSet._normalize_code_map_keys(extensions)
            if codes is None:
                if len(extensions) < len(PseudoSet.known_codes):
                    # codes is None, so we want all codes. Make sure any codes with
                    # unspecified filters are added with their defaults.
                    for code in PseudoSet.known_codes - extensions.keys():
                        extensions[code] = PseudoSet.file_exts[code]

                checked_filters = {}
                for code, suffixes in extensions.items():
                    code = PseudoSet._check_code_str(code)
                    if isinstance(suffixes, str):
                        suffixes = {suffixes}
                    elif isinstance(suffixes, Collection):
                        suffixes = set(suffixes)
                    elif suffixes is not None:
                        msg = (
                            f"File extensions must be a str, a collection of str, or None!\n"
                            f"Received: {type(suffixes).__name__}"
                            )
                        raise TypeError(msg)

                    checked_filters[code] = suffixes

                extensions = checked_filters
            else:
                # More filters than codes, or filters provided for unselected codes
                if not codes >= extensions.keys():
                    msg = (
                        "Mismatch between provided extensions and codes!\n"
                        f"Provided codes: {tuple(codes)}\n"
                        f"Filter keys:    {tuple(extensions.keys())}"
                        )
                    raise ValueError(msg)
                else:
                    checked_filters = {}
                    for c in codes:
                        code = PseudoSet._check_code_str(c)
                        checked_filters[code] = extensions.get(c)
                    extensions = checked_filters
        # extensions is now the source of truth for what codes will be included

        if code_Zeff_map is None:
            code_Zeff_map = {}
        elif all(map(Elements.is_element, code_Zeff_map)):
            # User gave one set of Z valences to apply to all codes
            Zeff_map = deepcopy(code_Zeff_map)
            code_Zeff_map = {code: Zeff_map for code in extensions}
        else:
            code_Zeff_map = PseudoSet._normalize_code_map_keys(code_Zeff_map)

        if not extensions.keys() >= code_Zeff_map.keys():
            msg = (
                "Mismatch between provided code Zeff map and codes!\n"
                f"Provided codes: {tuple(extensions.keys())}\n"
                f"Zeff keys:      {tuple(code_Zeff_map.keys())}"
                )
            raise ValueError(msg)

        if include is None:
            include = {}
        elif isinstance(include, str):
            include = {code: include for code in extensions}
        else:
            include = PseudoSet._normalize_code_map_keys(include)

        if not extensions.keys() >= include.keys():
            msg = (
                "Mismatch between provided include patterns and codes!\n"
                f"Provided codes: {tuple(codes)}\n"
                f"Pattern keys:   {tuple(include.keys())}"
                )
            raise ValueError(msg)

        if exclude is None:
            exclude = {}
        elif isinstance(exclude, str):
            exclude = {code: exclude for code in extensions}
        else:
            exclude = PseudoSet._normalize_code_map_keys(exclude)

        if not extensions.keys() >= exclude.keys():
            msg = (
                "Mismatch between provided exclude patterns and codes!\n"
                f"Provided codes: {tuple(codes)}\n"
                f"Pattern keys:   {tuple(exclude.keys())}"
                )
            raise ValueError(msg)

        pseudos = {}
        for code, suffixes in extensions.items():
            try:
                pseudos[code] = cls.from_dir(
                    pseudo_dir   = psp_dir,
                    code         = code,
                    Zeff_map     = code_Zeff_map.get(code),
                    extension    = suffixes,
                    include      = include.get(code),
                    exclude      = exclude.get(code),
                    skip_invalid = skip_invalid,
                    )
            except ValueError as err:
                msg = str(err) + (
                    f"\n\nProblem detected for code '{code}'\n"
                    f"Either remove '{code}' from the selected codes or specify\n"
                    "`extensions` and/or `include`/`exclude` to ensure the problem does not occur.\n"
                    )
                # Raise from None to prevent exception chain.
                raise ValueError(msg) from None

        sorted_pseudos = {k: v for k, v in sorted(pseudos.items(), key=lambda x: x[0])}
        return sorted_pseudos
    #end def from_mixed_dir


    def _get_pseudos(
        self,
        system: PhysicalSystem | Collection[str],
        code: Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"],
        ) -> dict[str, str]:
        """Private helper function for getting the pseudo files for a given system."""
        code = PseudoSet._check_code_str(code)
        if code not in self.codes:
            msg = f'Pseudopotential set is not available for code "{code}"'
            raise ValueError(msg)

        if system is None:
            # Return all pseudos
            return {psp.name: str(psp) for psp in self.pseudos.values()}
        elif isinstance(system, PhysicalSystem):
            species = set(system.structure.species(symbol=True)[1])
        else:
            species = set()
            for elem in system:
                is_elem, element = Elements.is_element(elem, return_element=True)
                if not is_elem:
                    msg = f"Non-element in provided system: '{element}'"
                else:
                    species.add(element.symbol)
            #end for
        #end if
        missing = species.difference(set(self.pseudos.keys()))
        if missing:
            msg = (
                'Pseudopotential set does not contain the following species:\n'
                f'{sorted(missing)}'
                )
            raise ValueError(msg)

        return {
            self.pseudos[element].name: str(self.pseudos[element])
            for element in sorted(species)
            }
    #end def _get_pseudos


    @staticmethod
    def get_pseudos(
        pseudos: PseudoSet | str | Collection[str] | dict[str, Path] | None,
        system: PhysicalSystem | Collection[str],
        code: Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"],
        ) -> dict[str, str]:
        """Get the pseudopotential files for the elements in a physical system.

        Parameters
        ----------
        pseudos : PseudoSet or str or Collection of str or dict[str: PathLike] or None
            Either a :class:`PseudoSet` object, a string label from
            :func:`ppset`, a list of pseudopotential file names, a dictionary
            mapping file names to file paths, or ``None``.
        system : PhysicalSystem or collection of str or None
            The system to get pseudopotentials for, or a collection of element
            labels, or ``None``.
        code : {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"}
            The name of the code requesting the pseudopotentials. This provides
            a way to ensure the user does not provide the wrong pseudos to a
            ``generate`` call.

        Returns
        -------
        pseudos : dict[str: str]
            A dictionary mapping pseudopotential file names to their full paths.

        Examples
        --------
        Call with an instance of :class:`PseudoSet`.

        >>> psps = PseudoSet(["C.ccECP.xml", "H.ccECP.xml", "O.ccECP.xml"])
        >>> PseudoSet.get_pseudos(
        ...     pseudos=psps,
        ...     system=["C", "H"],
        ...     code="qmcpack",
        ...     )
        {'C.ccECP.xml': '/path/to/pseudo_dir/C.ccECP.xml',
         'H.ccECP.xml': '/path/to/pseudo_dir/H.ccECP.xml'}

        Or call with a label after using :func:`ppset`. (Deprecated)

        >>> ppset(
        ...     label   = 'my_pseudos',
        ...     pwscf   = ["C.ccECP.upf", "H.ccECP.upf", "O.ccECP.upf"],
        ...     qmcpack = ["C.ccECP.xml", "H.ccECP.xml", "O.ccECP.xml"],
        ...     )
        >>> PseudoSet.get_pseudos(
        ...     pseudos="my_pseudos",
        ...     system=["C", "H"],
        ...     code="qmcpack",
        ...     )
        {'C.ccECP.xml': '/path/to/pseudo_dir/C.ccECP.xml',
         'H.ccECP.xml': '/path/to/pseudo_dir/H.ccECP.xml'}

        You can also call with a list of pseudopotential names.

        >>> PseudoSet.get_pseudos(
        ...     pseudos=["C.ccECP.xml", "H.ccECP.xml"],
        ...     system=["C", "H"],
        ...     code="qmcpack",
        ...     )
        {'C.ccECP.xml': '/path/to/pseudo_dir/C.ccECP.xml',
         'H.ccECP.xml': '/path/to/pseudo_dir/H.ccECP.xml'}

        Most commonly this would be called with a :class:`PhysicalSystem` object,
        created here with :func:`~.physical_system.generate_physical_system`.

        >>> sys = generate_physical_system(
        ...     elem = ['C','H'],
        ...     pos  = np.empty((2,3)),
        ...     C    = 4,
        ...     H    = 1,
        ...     )
        >>> PseudoSet.get_pseudos(
        ...     pseudos=["C.ccECP.xml", "H.ccECP.xml"],
        ...     system=sys,
        ...     code="qmcpack",
        ...     )
        {'C.ccECP.xml': '/path/to/pseudo_dir/C.ccECP.xml',
         'H.ccECP.xml': '/path/to/pseudo_dir/H.ccECP.xml'}

        Calling with a dict of :class:`PseudoSet`.

        >>> psps = PseudoSet.from_mixed_dir(
        ...     pseudo_dir="/path/to/pseudo_dir",
        ...     codes={"espresso", "qmcpack"},
        ... )
        >>> for code, ps_set in psps.items():
        ...     print(f"{code} pseudos:")
        ...     for lbl, psp in ps_set.pseudos.items():
        ...         print(f"  {lbl}: {psp}")
        espresso pseudos:
        H: /path/to/pseudo_dir/H.ccECP.upf
        C: /path/to/pseudo_dir/C.ccECP.upf
        qmcpack pseudos:
        C: /path/to/pseudo_dir/C.ccECP.xml
        H: /path/to/pseudo_dir/H.ccECP.xml
        >>> PseudoSet.get_pseudos(
        ...     pseudos=psps,
        ...     system=["C", "H"],
        ...     code="qmcpack",
        ...     )
        {'C.ccECP.xml': '/path/to/pseudo_dir/C.ccECP.xml',
         'H.ccECP.xml': '/path/to/pseudo_dir/H.ccECP.xml'}
        >>> PseudoSet.get_pseudos(
        ...     pseudos=psps,
        ...     system=["C", "H"],
        ...     code="espresso",
        ...     )
        {'C.ccECP.upf': '/path/to/pseudo_dir/C.ccECP.upf',
         'H.ccECP.upf': '/path/to/pseudo_dir/H.ccECP.upf'}
        """
        if pseudos is None:
            return {}

        elif isinstance(pseudos, PseudoSet):
            return pseudos._get_pseudos(system=system, code=code)
        elif isinstance(pseudos, Mapping) and len(pseudos) > 0:
            contains_pseudosets = any(
                isinstance(pseudoset, PseudoSet) for pseudoset in pseudos.values()
                )
            if contains_pseudosets:
                if not all(isinstance(pseudoset, PseudoSet) for pseudoset in pseudos.values()):
                    msg = "A pseudopotential-set mapping must contain only PseudoSet values"
                    raise TypeError(msg)

                code = PseudoSet._check_code_str(code)
                if code not in pseudos:
                    msg = f'Pseudopotential set is not available for code "{code}"'
                    raise ValueError(msg)

                if system is None:
                    msg = "Either provide a list of pseudo names or use a `PhysicalSystem` object for PseudoSet mapping"
                    raise ValueError(msg)

                psps = pseudos[code]
                return psps._get_pseudos(system=system, code=code)

            else: # Non-pseudoset objects as values, likely to occur when re-called.
                return dict(pseudos)
        elif isinstance(pseudos, str):
            label = pseudos
            if label not in PseudoSet.labeled_pseudosets:
                msg = (
                    f'Pseudopotential set label "{label}" is not registered\n'
                    f'Registered labels are {[*PseudoSet.labeled_pseudosets]}'
                    )
                raise KeyError(msg)

            if system is None:
                msg = "Either provide a list of pseudo names or use a `PhysicalSystem` object for a ppset label!"
                raise ValueError(msg)

            code = PseudoSet._check_code_str(code)
            if code not in PseudoSet.labeled_pseudosets[label]:
                msg = (
                    f'Pseudopotential set "{label}" is not available for code "{code}"\n'
                    f'Available codes for {label} are {[*PseudoSet.labeled_pseudosets[label]]}'
                    )
                raise KeyError(msg)

            psps = PseudoSet.labeled_pseudosets[label][code]
            return psps._get_pseudos(system=system, code=code)
        elif isinstance(pseudos, Collection):
            if len(pseudos) == 0:
                return {}
            else:
                missing = set(pseudos) - PseudoSet.pseudo_files.keys()
                if len(missing) > 0:
                    msg = f'The following pseudopotential files are not present in PseudoSet.pseudo_files:\n'
                    for psp in sorted(missing):
                        msg += f"  - {str(psp)}\n"
                    raise FileNotFoundError(msg)
                psps = PseudoSet([PseudoSet.pseudo_files[f] for f in pseudos])
                return psps._get_pseudos(system=system, code=code)
        else:
            msg = f"Expected PseudoSet or str or Collection of str or dict[str: PseudoSet] or None, but got {type(pseudos).__name__}!"
            raise TypeError(msg)
    #end def get_pseudos


    def get_Zeff(
        self,
        elem_labels: Collection[Elements | str] | PhysicalSystem,
        *,
        missing_as_ae: bool = False,
        ) -> dict[str, int]:
        """Get the Z-valences for each element in the list of elements.

        Parameters
        ----------
        elem_labels : list of Elements or list of str or PhysicalSystem
            The elements or system to get Z-valences for.
        missing_as_ae : bool, default=False (keyword-only)
            Assume any elements for which a pseudopotential can not be found are
            all-electron, and use their atomic number as the value for the
            Z-valence. If this is not supplied, and if the Z-valence for an
            element is not in ``self.Zeff``, this will attempt to extract the
            Z-valence from the pseudopotential file.

        Returns
        -------
        elem_Zeff : dict of str: int
            A dictionary mapping element labels to their effective nuclear
            charges.

        See Also
        --------
        read_upf_z_valence : Used to extract Z-valences from UPF files.
        read_qmcpack_xml_z_valence : Used to extract Z-valences from XML files.
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
                elif f_ext in PseudoSet.file_exts["gamess"]:
                    msg = (
                        "Z-valence parsing not implemented for GAMESS pseudopotentials!\n"
                        "You must supply Z-valences manually until this feature is added."
                        )
                    raise NotImplementedError(msg)
                elif f_ext in PseudoSet.file_exts["vasp"]:
                    Z_eff_map[label] = read_potcar_z_valence(self.pseudos[label])
                elif f_ext in PseudoSet.file_exts["qmcpack"]:
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
        rep += ")"
        return rep
    #end def __repr__
#end class PseudoSet


def generate_pseudoset(
    pseudo_dir: PathLike | None = None,
    *,
    code     : str | Collection[str] | None = None,
    extension: Mapping[str, str | Collection[str]] | None = None,
    include  : Mapping[str, str] | str | None = None,
    exclude  : Mapping[str, str] | str | None = None,
    Zeff_map : Mapping[str, Mapping[str, int]] | Mapping[str, int] | None = None,
    **codes_psps: Collection[PathLike] | PathLike,
    ) -> dict[str, PseudoSet]:
    """Generate a dictionary of :class:`PseudoSet`.

    The ideal use of this function is to create a single collection of
    :class:`PseudoSet`'s that represent a single kind or type of pseudopotential.
    See the Notes and Examples for more information.

    Parameters
    ----------
    pseudo_dir : PathLike
        The directory from which to read pseudopotentials. Does not support
        nested directories, except for those that contain a POTCAR file.
        If you do not provide this or ``codes_psps`` it will try to use the
        ``pseudo_dir`` supplied in :data:`~nexus.settings`.
    code : one or more of {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"}, optional
        The code(s) to use to separate the files in the directory
        into their respective groups.
        Incompatible with ``codes_psps``.
    extension : Map of codes to a str/set of str, optional
        A dictionary mapping codes to the file extension(s) corresponding to
        those labels. If this is not provided, then the filters are
        automatically populated by the codes in ``code``.
    include : Map of codes to str or single str, optional
        A dictionary mapping codes to a Unix shell-style wildcard.
        All files matched by this are **included** in the final
        :class:`PseudoSet`. By default this matches all files.
    exclude : Map of codes to str or single str, optional
        A dictionary mapping codes to a Unix shell-style wildcard.
        All files matched by this are **excluded** in the final
        :class:`PseudoSet`. By default this matches all files.
    Zeff_map : Map of codes to Map of str/Elements to int or Map of str/Elements to int, optional
        A mapping for each code that maps elements to their effective nuclear
        charges (Z-valences), or a single mapping that will be applied to all
        codes. If this is supplied, it will override any parts of the code that
        may try to parse the pseudopotential to get the Z-valence.
    **codes_psps
        A mapping of code names to either a collection of pseudopotential file
        names/paths or a path to a directory with pseudopotential files in it.

    Notes
    -----
    This function is designed around the idea that the returned dict represents
    a single type of pseudos, e.g. ccECPs, USPPs, NCPPs, PAWs, with the same
    kind of pseudopotential for each code. That means that, for the same
    element, all of the Z-valences in the pseudos should match, all of the
    radial grids are the same, all pseudos are made with the same functional,
    and so on. Then, the object returned can be passed freely to any function
    that may accept it, e.g. as the ``pseudos`` keyword in
    :func:`~.pwscf.generate_pwscf` or :func:`~.qmcpack.generate_qmcpack`.

    If you are only working with one code, or prefer a class-based interface,
    you may want to use :class:`PseudoSet` directly, along with its methods
    :meth:`PseudoSet.from_dir` and/or :meth:`PseudoSet.from_mixed_dir`.

    Examples
    --------
    Basic usage, all files in one directory.

    >>> print(contents_of_pseudo_dir)
    pseudo_dir
    ├── C.ccECP.gamess
    ├── C.ccECP.nwchem
    ├── C.ccECP.upf
    ├── C.ccECP.xml
    ├── H.ccECP.gamess
    ├── H.ccECP.nwchem
    ├── H.ccECP.upf
    └── H.ccECP.xml
    >>> psps = generate_pseudoset(
    ...     pseudo_dir=pseudo_dir,
    ...     code={"qmcpack", "espresso", "rmg", "pyscf", "gamess"},
    ...     extension={"rmg": ".xml"},
    ... )
    >>> for code, ps_set in psps.items():
    ...     print(f"{code} pseudos:")
    ...     for lbl, psp in ps_set.pseudos.items():
    ...         print(f"  {lbl}: {psp}")
    espresso pseudos:
        C: /path/to/pseudo_dir/C.ccECP.upf
        H: /path/to/pseudo_dir/H.ccECP.upf
    pyscf pseudos:
        H: /path/to/pseudo_dir/H.ccECP.nwchem
        C: /path/to/pseudo_dir/C.ccECP.nwchem
    rmg pseudos:
        H: /path/to/pseudo_dir/H.ccECP.xml
        C: /path/to/pseudo_dir/C.ccECP.xml
    qmcpack pseudos:
        H: /path/to/pseudo_dir/H.ccECP.xml
        C: /path/to/pseudo_dir/C.ccECP.xml
    gamess pseudos:
        C: /path/to/pseudo_dir/C.ccECP.gamess
        H: /path/to/pseudo_dir/H.ccECP.gamess
    """
    if pseudo_dir is None and len(codes_psps) == 0:
        if nexus_core.pseudo_dir is not None:
            pseudo_dir = Path(nexus_core.pseudo_dir).resolve()
        else:
            msg = "Must supply `pseudo_dir` and/or `codes_psps`!"
            raise ValueError(msg)

    if pseudo_dir is not None:
        pseudo_dir = Path(pseudo_dir).resolve()
        if not pseudo_dir.exists():
            msg = "`pseudo_dir` must exist!"
            raise FileNotFoundError(msg)
        elif not pseudo_dir.is_dir():
            msg = "`pseudo_dir` must be a directory!"
            raise NotADirectoryError(msg)

    if len(codes_psps) == 0:
        return PseudoSet.from_mixed_dir(
            pseudo_dir    = pseudo_dir,
            codes         = code,
            extensions    = extension,
            include       = include,
            exclude       = exclude,
            code_Zeff_map = Zeff_map,
            )

    # codes_psps was supplied
    if code is not None:
        msg = "When supplying a direct map of codes to pseudos you cannot pass `code`!"
        raise ValueError(msg)

    if not all([isinstance(psps, Collection | Path) for psps in codes_psps.values()]):
        msg = "Must supply a directory or collection of file paths for direct map!"
        raise TypeError(msg)

    codes_psps = PseudoSet._normalize_code_map_keys(codes_psps)
    pseudosets = {}
    for code, psps in codes_psps.items():
        if pseudo_dir is not None:
            # Can make supplied values relative to the directory for direct mapping
            if isinstance(psps, str | Path):
                psps = (pseudo_dir / psps).resolve()
            elif isinstance(psps, Collection):
                psps = [(pseudo_dir / psp).resolve() for psp in psps]

        if isinstance(psps, str | Path):
            psp_path = Path(psps).resolve()
            if not psp_path.exists():
                msg = f"The path for code {code} does not exist!"
                raise FileNotFoundError(msg)
            elif psp_path.is_dir():
                if isinstance(extension, Mapping):
                    code_ext = extension.get(code)
                elif extension is None:
                    code_ext = extension
                else:
                    msg = f"`extension` must be either a mapping or None, but is {type(extension).__name__}!"
                    raise TypeError(msg)

                if isinstance(include, Mapping):
                    code_inc = include.get(code)
                else:
                    code_inc = include

                if isinstance(exclude, Mapping):
                    code_exc = exclude.get(code)
                else:
                    code_exc = exclude

                if (
                    isinstance(Zeff_map, Mapping)
                    and isinstance(next(iter(Zeff_map.values())), Mapping)
                    ):
                    code_zeff_map = Zeff_map.get(code)
                else:
                    code_zeff_map = Zeff_map

                try:
                    pseudosets[code] = PseudoSet.from_dir(
                        pseudo_dir = psp_path,
                        extension  = code_ext,
                        include    = code_inc,
                        exclude    = code_exc,
                        Zeff_map   = code_zeff_map,
                        )
                except Exception as exc:  # noqa: BLE001
                    msg = (
                        str(exc) + "\n\n"
                        f"Error when processing pseudo directory for code '{code}'"
                        )
                    raise type(exc)(msg)
            else:
                msg = "If you are providing a single path it must be to a directory!"
                raise NotADirectoryError(msg)
        else:
            pseudosets[code] = PseudoSet(pseudos=psps)
    return pseudosets
#end def generate_pseudoset
