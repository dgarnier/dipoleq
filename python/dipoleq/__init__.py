"""DipolEq: Levitated Dipole Equilibrium Solver

DipolEq is a Python package for solving the equilibrium of a levitated dipole
magnetized plasma. It is a Python wrapper around the C code that was originally
written by Michael Mauel and his group at Columbia University in the early 1999s
called TokaMac. The C code was later modified by Darren Garnier at Columbia
for use in the LDX experiment. The Python wrapper was written by Darren Garnier
at OpenStar Technologies, LTD.

"""

from pathlib import Path
from typing import Any

from typing_extensions import Self

from . import core, input, read, solve, util
from .input import MachineIn


class Machine(core.Machine):
    """Machine class for the Dipole Equilibrium Solver"""

    # don't even think about messing with init
    # and pybind11 superclass

    @classmethod
    def from_fileinput(cls, filename: str | Path) -> Self:
        """Generate a Machine object from a .in file
        Using the Python code to read the file

        Args:
            filename (str): The path to the .in file

            Returns:
                Machine: The Machine object
        """
        input_data = read.input_from_dotin(str(filename))
        return cls.from_input_data(input_data)

    @classmethod
    def from_yaml(cls, filename: str) -> Self:
        """Generate a Machine object from a .yaml file

        Args:
            filename (str): The path to the .yaml file

        Returns:
            Machine: The Machine object
        """
        input_data = read.input_from_yaml(str(filename))
        return cls.from_input_data(input_data)

    @classmethod
    def from_dict(cls, input_data: dict[str, Any]) -> Self:
        """Generate a Machine object from a dictionary

        Args:
            input_data (dict): The data to create the Machine object

        Returns:
            Machine: The Machine object
        """
        verified_data = MachineIn(**input_data)
        return cls.from_input_data(verified_data)

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """Generate a Machine object from a file

        Args:
            filename (str): The path to the file

        Returns:
            Machine: The Machine object
        """

        p = Path(filename)
        filename = str(filename)
        if p.suffix is ".in":
            return cls.from_fileinput(filename)
        elif p.suffix is ".yaml":
            return cls.from_yaml(filename)
        else:
            raise ValueError(f"Unknown file type: {filename}")

    @classmethod
    def from_input_data(cls, input: MachineIn) -> Self:
        m = cls()
        input.initalize_machine(m)
        return m

    def solve(self) -> None:
        solve.solve(self)

    def diff(self, other: Any, verbose: bool = False) -> bool:
        return util.machine_diff(self, other, verbose=verbose)

    def __eq__(self, other: Any) -> bool:
        return not self.diff(other)


__all__ = ["Machine", "MachineIn", "core", "solve", "input", "read", "util"]
