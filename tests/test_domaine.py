import numpy as np
from physlab.models.Domaine import Domaine
from physlab.models.Atome import Atome

def test_volume():
    assert Domaine(2, 2, 2).Volume() == 8
    assert Domaine(2, 3, 1).Volume() == 6
    assert np.isclose(Domaine(1e-9, 1e-9, 1e-9).Volume(), 1e-27)
    
def test_mur_droit():
    d = Domaine(2, 2, 2)
    a = Atome("H", 1.0, 0.1, (1.2, 0, 0), (5.0, 1.0, 0))
    d.reflection(a)
    assert np.isclose(a.position[0], 0.8)
    assert a.vitesse[0] == -5.0
    assert a.vitesse[1] == 1.0

def test_coin():
    d = Domaine(2, 2, 2)
    a = Atome("H", 1.0, 0.1, (1.1, -1.1, 0), (3.0, -4.0, 2.0))
    d.reflection(a)
    assert np.allclose(a.position, (0.9, -0.9, 0))
    assert np.allclose(a.vitesse, (-3.0, 4.0, 2.0))

def test_pas_de_collision():
    d = Domaine(2, 2, 2)
    a = Atome("H", 1.0, 0.1, (0.5, 0.5, 0.5), (1.0, 1.0, 1.0))
    d.reflection(a)
    assert np.allclose(a.position, (0.5, 0.5, 0.5))
    assert np.allclose(a.vitesse, (1.0, 1.0, 1.0))