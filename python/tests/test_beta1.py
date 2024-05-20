# Test dipoleq with basic LDX beta = 1 equilibrium

from pytest import approx

from dipoleq import Machine, MachineIn, solve, read_dotin



def test_FileInput():
    # see if we can read the file and get the right number of coils
    m = Machine('beta1.in')
    assert m.NumCoils == 7
    i = 0
    for coil in m.Coils:
        i += 1
    assert i == 7

def test_solve():
    # see if we get the right current
    m = Machine('beta1.in')
    solve(m)
    assert m.Plasma.Ip == approx(32984, rel=1e-4)

def test_read_dotin():
    data = read_dotin('beta1.in')
    
    
    assert data['PsiGrid']['GridType'] == 'Cartesian'
    assert data['Plasma']['NumPlasma'] == 1
    assert data['Coils'][0]['SubCoils'][0]['NumSubCoils'] == 1
    assert data['Shells'][0]['SubShells'][0]['NumSubShells'] == 1
    assert data['Separatricies'][0]['NumSeparatrices'] == 1
    assert data['Measures'][0]['NumMeasures'] == 1
    assert data['Limiters'][0]['NumLimiters'] == 1