# Test dipoleq with basic LDX beta = 1 equilibrium

from pytest import approx

from dipoleq import Machine, MachineIn
from dipoleq.core import Machine as coreMachine
from dipoleq.read import input_from_dotin
from dipoleq.solve import solve


def props(x):
    return {
        key: getattr(x, key)
        for key in dir(x)
        if not callable(getattr(x, key)) and not key.startswith("__")
    }


def test_FileInput():
    # see if we can read the file and get the right number of coils
    m = coreMachine("beta1.in")
    assert m.NumCoils == 7
    i = 0
    for coil in m.Coils:
        i += 1
    assert i == 7


def test_solve_old():
    # see if we get the right current
    m = coreMachine("beta1.in")
    solve(m)
    assert m.Plasma.Ip == approx(32984, rel=1e-4)


def test_solve_new():
    # see if we get the right current
    m = Machine.from_file("beta1.in")
    m.solve()
    assert m.Plasma.Ip == approx(32984, rel=1e-4)


def test_solve_yaml():
    # see if we get the right current
    m = Machine.from_file("beta1.yaml")
    m.solve()
    assert m.Plasma.Ip == approx(32984, rel=1e-4)


def test_read_dotin():
    mIn = input_from_dotin("beta1.in")


def test_old_new():
    m2 = coreMachine("beta1.in")
    m1 = Machine.from_fileinput("beta1.in")
    m1.Iname = m2.Iname = "test"
    assert m1 == m2


def test_yaml_new():
    m1 = Machine.from_yaml("beta1.yaml")
    m2 = Machine.from_fileinput("beta1.in")
    m1.Iname = m2.Iname = "test"
    assert m1 == m2
