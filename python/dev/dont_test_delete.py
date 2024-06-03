# test object memory management
# hard to test and find the bug
# so don't actually do this.. instead will
# just leak memory.

import os
from pathlib import Path

from dipoleq import Machine
from dipoleq.core import Machine as _Machine

data_dir = Path(os.path.realpath(__file__)).parent


def test_machine_delete() -> None:
    # see if we get the right current

    m = _Machine(str(data_dir / "beta1.in"))
    m._free()  # free the memory, this is crashing
    del m
    assert True


def test_machine_delete2():
    # see if we get the right current
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    m._free()
    del m
    assert True


if __name__ == "__main__":
    test_machine_delete()
    test_machine_delete2()
