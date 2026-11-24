import numpy as np

from physlab.models.Atome import Atome


def atome(masse, position, vitesse, rayon=1.0):
    return Atome("X", masse, rayon, position, vitesse)


def quantite_mouvement(a, b):
    return a.masse * a.vitesse + b.masse * b.vitesse


def energie_cinetique(a, b):
    return a.kinetic_energy() + b.kinetic_energy()


# ---------------------------------------------------------------- 9.3 : détection
def test_detection_contact_exact():
    a = atome(1.0, (0, 0, 0), (0, 0, 0))
    b = atome(1.0, (2, 0, 0), (0, 0, 0))      # distance = 2 = somme des rayons
    assert a.en_collision(b)


def test_detection_chevauchement():
    a = atome(1.0, (0, 0, 0), (0, 0, 0))
    b = atome(1.0, (1.5, 0, 0), (0, 0, 0))
    assert a.en_collision(b)


def test_detection_trop_loin():
    a = atome(1.0, (0, 0, 0), (0, 0, 0))
    b = atome(1.0, (2.5, 0, 0), (0, 0, 0))
    assert not a.en_collision(b)


# ---------------------------------------------------------------- 9.2 / 9.4 : scénarios
def test_scenario1_choc_frontal_masses_egales():
    # Les deux vitesses doivent simplement s'échanger.
    a = atome(1.0, (-1, 0, 0), (10, 0, 0))
    b = atome(1.0, (1, 0, 0), (-10, 0, 0))
    a.collision(b)
    assert np.allclose(a.vitesse, (-10, 0, 0))
    assert np.allclose(b.vitesse, (10, 0, 0))


def test_scenario2_choc_frontal_masses_differentes():
    # m2 = 4*m1, cible au repos : v1' = (m1-m2)/(m1+m2)*v1 = -6, v2' = 2*m1/(m1+m2)*v1 = 4
    a = atome(1.0, (-1, 0, 0), (10, 0, 0))
    b = atome(4.0, (1, 0, 0), (0, 0, 0))
    p0, e0 = quantite_mouvement(a, b), energie_cinetique(a, b)
    a.collision(b)
    assert np.allclose(a.vitesse, (-6, 0, 0))
    assert np.allclose(b.vitesse, (4, 0, 0))
    assert np.allclose(quantite_mouvement(a, b), p0)
    assert np.isclose(energie_cinetique(a, b), e0)


def test_scenario3_choc_excentre_cible_au_repos():
    # Masses égales, cible au repos, contact à 30° de l'axe x :
    # les deux vitesses finales doivent être perpendiculaires.
    a = atome(1.0, (-np.sqrt(3), 1, 0), (10, 0, 0))
    b = atome(1.0, (0, 0, 0), (0, 0, 0))
    p0, e0 = quantite_mouvement(a, b), energie_cinetique(a, b)
    a.collision(b)
    assert np.isclose(np.dot(a.vitesse, b.vitesse), 0.0, atol=1e-9)
    assert np.allclose(a.vitesse, (2.5, 5 * np.sqrt(3) / 2, 0))
    assert np.allclose(b.vitesse, (7.5, -5 * np.sqrt(3) / 2, 0))
    assert np.allclose(quantite_mouvement(a, b), p0)
    assert np.isclose(energie_cinetique(a, b), e0)


def test_scenario4_atomes_qui_s_eloignent():
    # Déjà en contact mais qui s'éloignent : aucune vitesse ne doit changer.
    a = atome(1.0, (-1, 0, 0), (-5, 0, 0))
    b = atome(1.0, (1, 0, 0), (5, 0, 0))
    a.collision(b)
    assert np.allclose(a.vitesse, (-5, 0, 0))
    assert np.allclose(b.vitesse, (5, 0, 0))


def test_scenario5_conservation_3d_aleatoire():
    # Cas général en 3D : masses et vitesses quelconques.
    rng = np.random.default_rng(0)
    for _ in range(100):
        m1, m2 = rng.uniform(0.5, 5.0, size=2)
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        a = atome(m1, (0, 0, 0), rng.normal(size=3) * 10)
        b = atome(m2, 2.0 * n, rng.normal(size=3) * 10)    # distance = 2 = contact
        p0, e0 = quantite_mouvement(a, b), energie_cinetique(a, b)
        a.collision(b)
        assert np.allclose(quantite_mouvement(a, b), p0)
        assert np.isclose(energie_cinetique(a, b), e0)