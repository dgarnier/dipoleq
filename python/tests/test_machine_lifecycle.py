# Regression tests for Machine construction/destruction.
# A bare (uninitialized) Machine used to SIGABRT on teardown because
# free_PsiGrid() freed the NR arrays that only init_PsiGrid() allocates.
# The crash kills the interpreter, so these must run in a subprocess.

import subprocess
import sys


def run_snippet(code: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )


def test_bare_machine_teardown() -> None:
    """Constructing a Machine without init() must not crash on del."""
    result = run_snippet("import dipoleq; m = dipoleq.Machine(); del m")
    assert result.returncode == 0, result.stderr


def test_bare_core_machine_teardown() -> None:
    """Same for the raw core binding."""
    result = run_snippet("import dipoleq; m = dipoleq.core.Machine(); del m")
    assert result.returncode == 0, result.stderr


def test_initialized_machine_teardown() -> None:
    """init()ed Machine teardown must also stay clean."""
    result = run_snippet("import dipoleq; m = dipoleq.core.Machine(); m.init(); del m")
    assert result.returncode == 0, result.stderr
