from __future__ import annotations
import os
import numpy as np
import numpy.typing as npt
import builtins
from abc import ABC
from enum import Enum


type PwscfInputType = (
    str
    | bool
    | int
    | float
    | list
)


class NamelistVariableBase(ABC):
    """Abstract base class for all namelist variables."""

    def __init__(
        self,
        name: str,
        datatype: PwscfInputType,
        default: PwscfInputType | None,
        allowed_values: tuple[PwscfInputType] | None = None,
        dependencies: tuple[NamelistVariableBase] | None = None,
    ):
        self.name = name
        self.datatype = datatype
        self.default = default
        self.required = self.default is None
        self.allowed_values = allowed_values
        self.dependencies = dependencies


    @property
    def value(self) -> PwscfInputType:
        return self._value


    @value.setter
    def value(self, val: PwscfInputType):
        if val is None and self.required:
            raise ValueError(
                f"Can not set {self.name} to None, it is required!"
            )
        elif not isinstance(val, self.datatype):
            return TypeError(
                f"Expected a value of type {self.datatype} for {self.name} but got {type(val)} instead!"
            )
        else:
            self._value = val


    def write_value(self) -> str:
        match self.datatype:
            case builtins.bool:
                if self._value:
                    return f"   {self.name:<15} = .true."
                else:
                    return f"   {self.name:<15} = .false."

            case builtins.str:
                return f"   {self.name:<15} = '{self.value}'"

            case builtins.int:
                return f"   {self.name:<15} = {self.value}"

            case builtins.float:
                return f"   {self.name:<15} = {self.value}"

            case builtins.list:
                string = ""
                if isinstance(self.value[0], str):
                    for index, val in enumerate(self.value):
                        string += f"   {self.name+f'({index+1})':<15} = '{val}'"
                else:
                    for index, val in enumerate(self.value):
                        string += f"   {self.name+f'({index+1})':<15} = {val}"
                return string
            case _:
                raise TypeError("Value is not a valid Pwscf Input Type!")
    #end def write_value
#end class InputCard


class ControlVariable(NamelistVariableBase, Enum):
    """Class representing a single variable that belongs to the
    &CONTROL input namelist for ``pw.x``.
    """

    def __new__(cls, datatype, default, allowed_values, dependencies): ...

    calculation = str, "scf", ["scf", "nscf", "bands", "relax", "md"], None
    


class SystemVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &SYSTEM input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class ElectronsVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &ELECTRONS input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class IonsVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &IONS input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class CellVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &CELL input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class FCPVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &FCP input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class RismVariable(NamelistVariableBase):
    """Class representing a single variable that belongs to the
    &RISM input namelist for ``pw.x``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)



