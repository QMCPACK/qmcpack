from __future__ import annotations
import os
import numpy as np
import numpy.typing as npt
import builtins
from dataclasses import dataclass
from abc import ABC, abstractmethod
from enum import Enum


type PwscfInputType = (
    str
    | bool
    | int
    | float
    | list
)

@dataclass
class NamelistDefinition:
    """Base class for all namelist variables."""

    name: str
    datatype: PwscfInputType
    default: PwscfInputType | None
    allowed_values: tuple[PwscfInputType] | None = None
    dependencies: tuple[NamelistDefinition] | None = None
#end class NamelistDefinition


    # @property
    # def value(self) -> PwscfInputType:
    #     return self._value


    # @value.setter
    # def value(self, val: PwscfInputType):
    #     if val is None and self.required:
    #         raise ValueError(
    #             f"Can not set {self.name} to None, it is required!"
    #         )
    #     elif not isinstance(val, self.datatype):
    #         return TypeError(
    #             f"Expected a value of type {self.datatype} for {self.name} but got {type(val)} instead!"
    #         )
    #     else:
    #         self._value = val


    # def write_value(self) -> str:
    #     match self.datatype:
    #         case builtins.bool:
    #             if self._value:
    #                 return f"   {self.name:<15} = .true."
    #             else:
    #                 return f"   {self.name:<15} = .false."

    #         case builtins.str:
    #             return f"   {self.name:<15} = '{self.value}'"

    #         case builtins.int:
    #             return f"   {self.name:<15} = {self.value}"

    #         case builtins.float:
    #             return f"   {self.name:<15} = {self.value}"

    #         case builtins.list:
    #             string = ""
    #             if isinstance(self.value[0], str):
    #                 for index, val in enumerate(self.value):
    #                     string += f"   {self.name+f'({index+1})':<15} = '{val}'"
    #             else:
    #                 for index, val in enumerate(self.value):
    #                     string += f"   {self.name+f'({index+1})':<15} = {val}"
    #             return string
    #         case _:
    #             raise TypeError("Value is not a valid Pwscf Input Type!")
    # #end def write_value
#end class InputCard


class ControlDefinitions(NamelistDefinition, Enum):
    """Class representing all variables that belongs to the
    &CONTROL input namelist for ``pw.x``.
    """

    def __new__(cls, datatype, default, allowed_values, dependencies): ...

    calculation     = str,   "scf",          ("scf", "nscf", "bands", "relax", "md", "vc-md", "vc-relax"), None
    title           = str,   "",             None, None
    verbosity       = str,   "low",          ("high", "low"), None
    restart_mode    = str,   "from_scratch", ("from_scratch", "restart"), None
    nstep           = int,   (1, 50),        None, (calculation, ("scf, nscf, bands"))
    iprint          = int,   10,             None, None
    tstress         = bool,  (True, False),  None, (calculation, ("vc-md", "vc-relax"))
    tprnfor         = bool,  (True, False),  None, (calculation, ("relax", "md", "vc-md"))
    dt              = float, 20.0,           None, (calculation, ("md"))
    outdir          = str,   "./",           None, None




class SystemVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &SYSTEM input namelist for ``pw.x``.
    """
    ...


class ElectronsVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &ELECTRONS input namelist for ``pw.x``.
    """
    ...


class IonsVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &IONS input namelist for ``pw.x``.
    """
    ...


class CellVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &CELL input namelist for ``pw.x``.
    """
    ...


class FCPVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &FCP input namelist for ``pw.x``.
    """
    ...


class RismVariables(NamelistDefinition, Enum):
    """Class representing a single variable that belongs to the
    &RISM input namelist for ``pw.x``.
    """
    ...


class PWscfNamelist(ABC):
    """Abstract base class for all PWscf namelists."""

    @abstractmethod
    def write_card(self) -> str:
        pass

    @abstractmethod
    def meets_requirements(self) -> bool:
        pass

    @abstractmethod
    def __setattr__(self, name, value):
        pass

