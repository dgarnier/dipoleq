# Ensure the typing artifacts actually ship in the installed package.
# CI runs pytest inside every built wheel (cibuildwheel test-command),
# so these assertions guard the wheel contents on all released platforms.

import sys
from pathlib import Path

import dipoleq

pkg_dir = Path(dipoleq.__file__).parent


def test_py_typed_marker_installed() -> None:
    """PEP 561 marker must ship so type checkers use our annotations."""
    assert (pkg_dir / "py.typed").is_file()


def test_core_stub_installed() -> None:
    """The generated core.pyi must ship next to the extension module."""
    if sys.platform == "win32":
        # stub generation is skipped on Windows (see python/CMakeLists.txt)
        return
    stub = pkg_dir / "core.pyi"
    assert stub.is_file()
    assert "class Machine" in stub.read_text()
