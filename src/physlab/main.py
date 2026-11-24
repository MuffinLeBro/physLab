import numpy as np

from physlab.models.Atome import Atome
from physlab.models.Domaine import Domaine
from physlab.models.Simulation import Simulation

# --- Paramètres de la simulation ---
LONGUEURS = (20.0, 20.0, 20.0)   # m : dimensions du domaine (Lx, Ly, Lz)
V0 = (10.0, 20.0, 5.0)           # m/s : vitesse initiale de l'atome de test
DT = 1e-3                        # s : pas de temps
T_FINAL = 10.0                   # s : durée simulée
EXPORT_EVERY = 10                # on exporte un pas sur 10
EXPORT_PATH = "export/simulation"


def creer_atome() -> Atome:
    """Atome de test, au centre du domaine."""
    return Atome(
        formuleChimique="H",
        masse=1.67e-27,
        rayon=0.5,               # m, grand pour être visible dans ParaView
        position=(0.0, 0.0, 0.0),
        vitesse=V0,
    )


def simuler(sim: Simulation, atome: Atome, domaine: Domaine) -> None:
    """Avance la simulation, exporte les images et vérifie que l'atome reste dans le domaine."""
    demi = np.array([domaine.longX, domaine.longY, domaine.longZ]) / 2
    nb_pas = int(round(T_FINAL / DT))

    sim.open_export(EXPORT_PATH)
    try:
        sim.write_step_into_export(time=0.0)
        for n in range(1, nb_pas + 1):
            sim.step(DT)
            assert np.all(np.abs(atome.position) <= demi), f"Atome hors du domaine au pas {n}"
            if n % EXPORT_EVERY == 0:
                sim.write_step_into_export(time=n * DT)
    finally:
        sim.close_export()


def main() -> None:
    domaine = Domaine(*LONGUEURS)
    atome = creer_atome()
    sim = Simulation(np.array([atome]), domaine)
    print(sim)

    e0 = atome.mechanical_energy()
    simuler(sim, atome, domaine)

    print(f"Position finale : {atome.position} m")
    print(f"Vitesse finale  : {atome.vitesse} m/s")
    print(f"Dérive de l'énergie mécanique : {atome.mechanical_energy() - e0:.3e} J")


if __name__ == "__main__":
    main()