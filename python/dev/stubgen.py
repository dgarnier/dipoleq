import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any, ClassVar

if importlib.util.find_spec("pybind11_stubgen") is None:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pybind11-stubgen"])

from pybind11_stubgen import CLIArgs, arg_parser, stub_parser_from_args
from pybind11_stubgen.printer import Printer
from pybind11_stubgen.structs import Argument, QualifiedName
from pybind11_stubgen.writer import Writer


class MyPrinter(Printer):  # type: ignore[misc]
    """Replace "tokamak" with "Machine" in the argument types"""

    SUBS: ClassVar = {"tokamak": "Machine"}

    def print_argument(self, arg: Argument) -> Any:
        "overload to replace 'tokamak' with 'Machine'"
        result = super().print_argument(arg)
        for old, new in self.SUBS.items():
            result = result.replace(old, new)
        return result


def generate_stubs(module_fp: str, stub_dir: str) -> None:
    """generate stubs for the given module"""

    module_dir = Path(module_fp).parent.absolute()
    module_name = Path(module_fp).stem.split(".")[0]
    out_dir = Path(stub_dir).absolute()
    old_syspath = sys.path
    sys.path.append(str(module_dir))

    sys.argv = [
        "pybind11-stubgen",
        module_name,
        "-o",
        str(stub_dir),
        "--numpy-array-use-type-var",
    ]
    args = arg_parser().parse_args(namespace=CLIArgs())

    parser = stub_parser_from_args(args)
    module = parser.handle_module(
        QualifiedName.from_str(module_name), importlib.import_module(module_name)
    )
    parser.finalize()
    sys.path = old_syspath

    writer = Writer()
    printer = MyPrinter(invalid_expr_as_ellipses=False)
    writer.write_module(module, printer, out_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python stubgen.py <module_filepath> <output_dir>")
        sys.exit(1)
    module_filepath, stub_output_dir = sys.argv[1:]
    generate_stubs(module_filepath, stub_output_dir)
