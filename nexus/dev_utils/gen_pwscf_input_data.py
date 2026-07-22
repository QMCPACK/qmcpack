#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "packaging>=26.2",
#     "xmltodict>=1.0.4",
# ]
# ///

"""Script for generating the input parameter enums for ``nexus/pwscf_input.py``.

Designed to be run with ``uv``, which will ensure the correct Python
version and get all the required packages.

Requires the ``xmltodict`` and ``packaging`` packages.

Written by Brock Dyer, last updated on July 15, 2026.
"""

from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import textwrap
from typing import Literal
import pprint


try:
    import xmltodict
except ImportError as err:
    raise err("The Python package xmltodict is required for this script!")

try:
    from packaging.version import Version
except ImportError as err:
    raise err("The Python package 'packaging' is required for this script!")

# These are the dates of the commits attached to the INPUT_PW.def
# files for each version tag in the QE GitLab. Update as needed.
QE_DOC_VERSION_DATES = {
    #         DD-MM-YYYY
    Version("4.0.4"): "18-11-2008",
    Version("4.1.0"): "16-07-2009",
    Version("4.1.1"): "01-10-2009",
    Version("4.1.3"): "24-03-2010",
    Version("4.2.1"): "12-06-2010",
    Version("4.3.0"): "13-02-2011",
    Version("4.3.1"): "20-05-2011",
    Version("5.0.0"): "11-05-2012",
    Version("5.0.1"): "04-07-2012",
    Version("5.0.2"): "04-07-2012",
    Version("5.1.0"): "25-04-2014",
    Version("5.1.1"): "16-10-2014",
    Version("5.1.2"): "02-03-2015",
    Version("5.2.0"): "20-06-2015",
    Version("5.2.1"): "30-07-2015",
    Version("5.3.0"): "07-01-2016",
    Version("5.4.0"): "20-04-2016",
    Version("6.0.0"): "19-08-2016",
    Version("6.1.0"): "24-02-2017",
    Version("6.2.0"): "03-09-2017",
    Version("6.2.1"): "25-10-2017",
    Version("6.3.0"): "15-06-2018",
    Version("6.4.0"): "01-03-2019",
    Version("6.4.1"): "03-04-2019",
    Version("6.5.0"): "21-11-2019",
    Version("6.6.0"): "19-07-2020",
    Version("6.7.0"): "30-11-2020",
    Version("6.8.0"): "22-04-2021",
    Version("7.0.0"): "18-12-2021",
    Version("7.1.0"): "08-06-2022",
    Version("7.2.0"): "18-03-2023",
    Version("7.3.0"): "29-11-2023",
    Version("7.3.1"): "21-02-2024",
    Version("7.4.0"): "16-10-2024",
    Version("7.4.1"): "03-02-2025",
    Version("7.5.0"): "30-06-2025",
    Version("7.6.0"): "02-06-2026",
}
EARLIEST = list(QE_DOC_VERSION_DATES.keys())[0]
LATEST = list(QE_DOC_VERSION_DATES.keys())[-1]

def get_total_dict(
    xml_dir: Path = Path("./xml_files"),
    save_json: bool = False,
) -> dict:
    """Take a directory of files labeled ``INPUT_PW.xml`` and parse them into dictionaries.

    Needs the ``xmltodict`` package.

    Notes
    -----
    You can turn the ``INPUT_PW.def`` files into XML by running the
    Quantum ESPRESSO dev tool ``helpdoc`` found at ``dev-tools/helpdoc``.

    If you get an error about a variable ``$basedir`` not existing, you can
    replace the reference to it on line 186 in ``dev-tools/helpdoc.d/helpdoc.tcl``
    with the absolute path to ``dev-tools``, like so:

        namespace eval schema { ::source [file join /path/to/q-e/dev-tools helpdoc.schema] }
    """
    xml_dir = xml_dir.resolve()
    all_inputs_dict = {}
    for file in xml_dir.iterdir():
        if file.suffix != ".xml":
            continue

        qe_version = Version(file.stem.removeprefix("INPUT_PW_"))

        with open(file, "rb") as inp:
            # These are labeled as being ISO-8859-1 encoded, but as best I can tell they
            # are actually UTF-8 encoded.
            inp_dict = xmltodict.parse(
                inp, cdata_key="info", encoding="UTF-8", attr_prefix=""
            )["input_description"]
            # inp_dict = xmltodict.parse(
                # inp, cdata_key="info", encoding="ISO-8859-1", attr_prefix=""
            # )["input_description"]

        temp_dict = {}
        for key, value in inp_dict.items():
            if key not in ["program", "namelist", "card"]:
                continue

            if key == "program":
                temp_dict[key] = f"Quantum ESPRESSO {qe_version!s}"
            else:
                temp_dict[key] = {}

            temp_dict["date"] = QE_DOC_VERSION_DATES[qe_version]

            if isinstance(value, list):
                for item in value:
                    name = item.pop("name")
                    temp_dict[key][name] = item

        all_inputs_dict[qe_version] = temp_dict

    all_inputs_dict = {
        key: val for key, val in sorted(all_inputs_dict.items(), key=lambda x: x[0])
    }

    if save_json:
        with open(f"INPUT_PW_{EARLIEST!s}_to_{LATEST!s}.json", "w") as out:
            out.write(json.dumps(all_inputs_dict, indent=4))

    return all_inputs_dict
#end def get_total_dict


type PwscfInputType = str | bool | int | float | list

@dataclass
class NamelistParamDefinition:
    """Base class for all namelist variables."""

    datatype: tuple[PwscfInputType]
    required: bool
    shape: tuple | None = None
    allowed_values: tuple[PwscfInputType] | None = None
    # dependencies: tuple[NamelistParamDefinition] | None = None
    version_added: Version | None = None
    version_removed: Version | None = None
# end class NamelistParamDefinition


class PWDef:
    """Class for parsing, comparing, and writing information about the namelists in INPUT_PW.def files.

    Designed to work with the generated dict from ``get_total_dict``.
    """

    #: Fortran to Python datatype mapping
    f2p_datatype = {
        "INTEGER": "int",
        "REAL": "float",
        "LOGICAL": "bool",
        "CHARACTER": "str",
    }

    def __init__(self, input_data: dict, version: str | Version):
        self.version = Version(version)

        self.namelists = list(input_data["namelist"].keys())
        self.cards = list(input_data["card"].keys())
        self.namelist_params = self.parse_namelist_params(input_data)

    @staticmethod
    def parse_namelist_params(
        input_data,
    ) -> dict[str, dict[str, NamelistParamDefinition]]:
        """Automatically parse the namelists in the input data.

        This function mostly serves to detect the form of the parameters,
        e.g. ``var``, ``group``, ``vargroup``, ``dimension``, and
        ``multidimension``, then call up their respective parsing
        functions. It then collects all of the parsed parameters and
        organizes them into a dictionary with the top level keys being
        the namelist name, then the next level being the parameter form,
        which is then a dictionary of parameter names as keys and their
        definitions as instances of the ``NamelistParamDefinition``
        dataclass.

        Returns
        -------
        namelist_params : dict[str, dict[str, dict[str, NamelistParamDefinition]]]
            Dictionary of the form

            .. code-block:: python

                namelist_params = {
                    "CONTROL": {
                        "var": {
                            "calculation": NamelistParamDefinition(...),
                            ...
                        },
                    },
                    "SYSTEM": {
                        "var": {
                            "ibrav": NamelistParamDefinition(...),
                            ...
                        },
                        "dimension": { ... },
                        ...
                    },
                    ...
                }
        """
        namelist_params = {}
        # namelist is ("CONTROL", "SYSTEM", ...)
        # param_forms is {"var": list[dict], "dimension": list[dict], ...}
        for namelist, param_forms in input_data["namelist"].items():
            print(f"Parsing namelist {namelist}...")
            namelist_params[namelist] = {}
            for param_form, param_defs in param_forms.items():
                # Apparently sometimes these aren't lists, even
                # though there are ample situations where there
                # are length-one lists of "dimension" in other
                # places.
                # Why not just make them all lists? Please?
                if isinstance(param_defs, dict):
                    param_defs = [param_defs]
                match param_form:
                    case "var":
                        parsed_params = PWDef.parse_namelist_var(param_defs)
                    case "group":
                        parsed_params = PWDef.parse_namelist_group(param_defs)
                    case "vargroup":
                        parsed_params = PWDef.parse_namelist_vargroup(param_defs)
                    case "dimension":
                        parsed_params = PWDef.parse_namelist_dimension(param_defs)
                    case "multidimension":
                        parsed_params = PWDef.parse_namelist_multidimension(param_defs)
                    case _:
                        continue
                # end match
                namelist_params[namelist].update(parsed_params)

        return namelist_params

    @staticmethod
    def parse_namelist_var(var_list: list[dict]) -> dict[str, NamelistParamDefinition]:
        """Parse all of the namelist parameters with form ``var``.

        The output dictionary has their name as the keys and their
        definitions as the values.
        """
        parsed_vars = {}
        for var in var_list:
            name = var["name"]
            # QE 5.0.0 to 6.3.0 doesn't have a type defined for this,
            # but its default is a float, so I guess this should be REAL.
            if name == "Hubbard_J(i,ityp)" and "type" not in var:
                var["type"] = "REAL"
            # I have to use upper() because whoever wrote the doc
            # for dftd3_version forgot their shift key and caps lock
            # key at home.
            datatype = PWDef.f2p_datatype[var["type"].upper()]
            required = var.get("status", "")
            if isinstance(required, dict):
                required = False
            else:
                if required.lower() == "required":
                    required = True
                else:
                    required = False

            # This is only really going to work for things with clearly
            # defined options. Something like rism3d_conv_level doesn't
            # really have options per-se, it just uses the options tag
            # because apparently the info tag was out sick when that was
            # written.
            if "options" in var.keys() and datatype == "str":
                allowed_values = []
                for option in var["options"]["opt"]:
                    # Remove extra single quote marks
                    value = option["val"].replace("'", "")
                    if "," in value:
                        value = [i.strip() for i in value.split(",")]
                    else:
                        value = [value]

                    # Someone decided it'd be a good idea to write
                    # "dftd3_version = <value>" in the dftd3_version var.
                    # So now it's my problem...
                    value = [v.split()[-1] if name in v else v for v in value]
                    allowed_values.extend(value)
                allowed_values = tuple(allowed_values)
            else:
                allowed_values = None

            if "(" in name:
                name = name.split("(")[0]
            parsed_vars[name] = NamelistParamDefinition(
                datatype,
                required,
                None,
                allowed_values,
            )
        return parsed_vars

    @staticmethod
    def parse_namelist_dimension(
        dimension_list: list[dict],
    ) -> dict[str, NamelistParamDefinition]:
        """Parse all of the namelist parameters with form ``dimension``.

        The output dictionary has their name as the keys and their
        definitions as the values.
        """
        parsed_dimensions = {}
        for dim in dimension_list:
            name = dim["name"]

            required = dim.get("status", "")
            if isinstance(required, dict):
                required = False
            else:
                if required.lower() == "required":
                    required = True
                else:
                    required = False

            start = (1,)
            if dim["end"].isdigit():
                end = (int(dim["end"]),)
            else:
                end = (dim["end"],)
            parsed_dimensions[name] = NamelistParamDefinition(
                datatype=PWDef.f2p_datatype[dim["type"].upper()],
                required=required,
                shape=(start, end),
                allowed_values=None,
            )
        return parsed_dimensions

    @staticmethod
    def parse_namelist_multidimension(
        multidimension_list: list[dict],
    ) -> dict[str, NamelistParamDefinition]:
        """Parse all of the namelist parameters with form ``multidimension``.

        The output dictionary has their name as the keys and their
        definitions as the values.
        """
        parsed_multidimensions = {}
        for multidim in multidimension_list:
            name = multidim["name"]
            start = (1,)
            if multidim["end"].isdigit():
                end = (int(multidim["end"]),)
            else:
                end = (multidim["end"],)

            required = multidim.get("status", "")
            if isinstance(required, dict):
                required = False
            else:
                if required.lower() == "required":
                    required = True
                else:
                    required = False

            parsed_multidimensions[name] = NamelistParamDefinition(
                datatype=PWDef.f2p_datatype[multidim["type"].upper()],
                required=required,
                shape=(start, end),
                allowed_values=None,
            )
        return parsed_multidimensions

    @staticmethod
    def parse_namelist_vargroup(
        vargroup_list: list[dict],
    ) -> dict[str, NamelistParamDefinition]:
        """Parse all of the namelist parameters with form ``vargroup``.

        The output dictionary has their name as the keys and their
        definitions as the values.
        """
        parsed_vargroups = {}
        for vargroup in vargroup_list:
            datatype = vargroup["type"].upper()
            datatype_checked = False
            for var in vargroup["var"]:
                name = var["name"]
                if not datatype_checked and "," in datatype:
                    raise NotImplementedError(
                        "Namelist vargroups with heterogeneous types are not implemented!\n"
                        f"Either contact a developer or enter the data for {name} manually!"
                    )
                elif not datatype_checked:
                    datatype = PWDef.f2p_datatype[datatype]
                    datatype_checked = True

                required = vargroup.get("status", "")
                if isinstance(required, dict):
                    required = False
                else:
                    if required.lower() == "required":
                        required = True
                    else:
                        required = False

                if "start" in vargroup.keys() or "end" in vargroup.keys():
                    raise NotImplementedError(
                        "Vargroups of dimensions/multidimensions are not implemented!\n"
                        f"Either contact a developer or enter the data for {name} manually!"
                    )
                else:
                    shape = None

                parsed_vargroups[name] = NamelistParamDefinition(
                    datatype,
                    required,
                    shape,
                    allowed_values=None,
                )
            # end for var in vargroup
        # end for vargroup in vargroup_list
        return parsed_vargroups

    @staticmethod
    def parse_namelist_group(
        group_list: list[dict],
    ) -> dict[str, NamelistParamDefinition]:
        """Parse all of the namelist parameters within a ``group``
        and flatten them out into one dictionary.

        .. important::

            This does not retain the identity of the ``group``'s label.
        """
        parsed_groups = {}
        for group in group_list:
            for group_form, group_defs in group.items():
                # See `parse_namelist_params` for whining about this.
                if isinstance(group_defs, dict):
                    group_defs = [group_defs]

                match group_form:
                    case "var":
                        parsed_params = PWDef.parse_namelist_var(group_defs)
                    case "vargroup":
                        parsed_params = PWDef.parse_namelist_vargroup(group_defs)
                    case "dimension":
                        parsed_params = PWDef.parse_namelist_dimension(group_defs)
                    case "multidimension":
                        parsed_params = PWDef.parse_namelist_multidimension(group_defs)
                    case "group":
                        raise NotImplementedError(
                            "Groups containing form 'group' are not implemented yet!"
                        )
                    case _:
                        continue
                # end match
                parsed_groups.update(parsed_params)
        return parsed_groups

    @staticmethod
    def read_version_json(json_path: Path) -> dict:
        with open(json_path, "r") as json_inp:
            output = json.load(json_inp)
        return output


def get_version_info(input_dict: dict[Version, PWDef]) -> dict[str, dict[str, NamelistParamDefinition]]:
    nmlist_param_vers = defaultdict(dict)
    for version, data in input_dict.items():
        for namelist, params in data.namelist_params.items():
            for param, param_def in params.items():
                if param not in nmlist_param_vers[namelist]:
                    nmlist_param_vers[namelist][param] = param_def
                    nmlist_param_vers[namelist][param].version_added = version
                    nmlist_param_vers[namelist][param].version_removed = version
                else:
                    nmlist_param_vers[namelist][param].version_removed = version

    return dict(nmlist_param_vers)


def write_namelist_param_enum(input_pw: dict[str, dict[str, NamelistParamDefinition]]):
    namelist_enums = {}
    for namelist, param_dict in input_pw.items():
        namelist_string = ""
        name_max = max(map(len, param_dict.keys()))
        for name, param_def in param_dict.items():
            v_added_str = f"'{param_def.version_added!s}'"
            if param_def.version_removed != LATEST:
                v_removed_str = f"'{param_def.version_removed!s}'"
            else:
                v_removed_str = "None"
            namelist_string += (
                f"{name:<{name_max}} = "
                f"{f'{param_def.datatype},':<7}"
                f"{f'{param_def.required},':<7}"
                f"{param_def.shape}, "
                f"{param_def.allowed_values}, "
                f"{v_added_str}, "
                f"{v_removed_str}\n"
            )
        namelist_enums[namelist] = namelist_string
    return namelist_enums


# json_path = Path("./INPUT_PW_{EARLIEST!s}_to_{LATEST!s}.json").resolve()
# json_data = PWDef.read_version_json(json_path)
if __name__ == "__main__":
    if len(sys.argv) > 1:
        xml_path = Path(sys.argv[1]).resolve()
    else:
        xml_path = Path("./input_pw_data.tar.bz2").resolve()

    if not xml_path.exists():
        raise FileNotFoundError(f"Can not find XML data at location {xml_path}")

    if xml_path.is_file():
        if xml_path.suffixes[0] != ".tar":
            raise RuntimeError(
                "Must supply either the XML file directory or a tar file with the XML files"
            )
        with tempfile.TemporaryDirectory() as tmp_dir:
            with tarfile.open(xml_path, mode='r:*') as tar:
                tar.extractall(path=tmp_dir, filter='data')

            pw_data = get_total_dict(Path(tmp_dir) / "xml_files")
    else:
        pw_data = get_total_dict(xml_path)

    print("Parsing namelist data...\n")
    all_data = {}
    for version in QE_DOC_VERSION_DATES:
        print(version)
        all_data[version] = PWDef(pw_data[version], str(version))


    print("\nDone parsing namelist data!\n")
    input_pw_complete = get_version_info(all_data)

    output_file = Path("../nexus/pwscf_input_defs.py").resolve()

    print(f"Writing output to {output_file}\n")
    output_text = f'''\
"""Module for storing input parameters for Quantum ESPRESSO.

The code in this module was auto-generated by a Python script that
uses the intermediate XML representation of QE's ``pw.x`` input
description files.

The script is at ``qmcpack/nexus/dev_utils/gen_pwscf_input_data.py``.

A description of how to turn QE's ``.def`` documentation files into
XML is given in QE's repository, at ``q-e/dev-tools/README.helpdoc``.

This module's code was autogenerated on {datetime.now().strftime("%a %d %b %Y, %I:%M%p")}

The supported QE versions are v{EARLIEST!s}-v{LATEST!s}
"""

from __future__ import annotations
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

PwscfInputType: TypeAlias = str | bool | int | float | Sequence
# For structured array inputs (k-points, weights, reciprocal lattice, etc.)
PwscfArrayInput: TypeAlias = Sequence[float | Sequence[float]]


@dataclass(frozen=True) # We generally don't want any of these to be modified
class NamelistParamDefinition:
    """Base class for all namelist variables."""

    datatype:        tuple[PwscfInputType]
    required:        bool
    shape:           tuple | None = None
    allowed_values:  tuple[PwscfInputType] | None = None
    version_added:   tuple[float] | None = None
    version_removed: tuple[float] | None = None
#end class NamelistDefinition


class NamelistEnumBase(Enum):
    """Abstract base class for the namelist enumerations.

    Provides the ``__new__`` method for the enums.
    """
    def __new__(
        cls,
        datatype:        PwscfInputType,
        required:        bool,
        allowed_values:  tuple[PwscfInputType] | None = None,
        version_added:   tuple[float] | None = None,
        version_removed: tuple[float] | None = None,
    ):
        definition = NamelistParamDefinition.__new__(cls)
        definition._value_ = NamelistParamDefinition(
            datatype,
            required,
            allowed_values,
            version_added,
            version_removed,
        )
        return definition
#end class NamelistEnumBase


'''
    namelist_enums = write_namelist_param_enum(input_pw_complete)
    for key, value in namelist_enums.items():
        value = value.replace("lambda ", "lambda_")
        output_text += (
            f"class {key.title()}Definitions(NamelistParamDefinition, NamelistEnumBase):\n"
            f'    """All variables belonging to the &{key} input namelist for ``pw.x``."""\n\n'
        )
        output_text += textwrap.indent(value, prefix=" "*4)
        output_text += f"#end class {key.title()}Definitions\n\n\n"

    output_text = output_text.rstrip("\n") + "\n"

    with open(output_file, "w") as test_out:
        test_out.write(output_text)
