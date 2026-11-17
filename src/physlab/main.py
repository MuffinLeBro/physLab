from physlab.models.Atome import Atome
from physlab.models.Simulation import Simulation
import numpy as np

if __name__ == "__main__":
    hydrogene = Atome(
        formuleChimique="H",
        masse=1.67e-27,
        rayon=5.3e-11,
        position=(0.0, 0.0, 0.0),
        vitesse=(100.0, 0.0, 0.0),
    )

    oxygene = Atome(
        formuleChimique="O",
        masse=2.66e-26,
        rayon=6.0e-11,
        position=(2.5e-10, 0.0, 0.0),
        vitesse=(-50.0, 20.0, 0.0),
    )

    helium = Atome(
        formuleChimique="He",
        masse=6.64e-27,
        rayon=3.1e-11,
        position=(-2.5e-10, 0.0, 0.0),
        vitesse=(-80.0, 30.0, 0.0),
    )

    carbone = Atome(
        formuleChimique="C",
        masse=1.99e-26,
        rayon=7.0e-11,
        position=(0.0, 2.5e-10, 0.0),
        vitesse=(40.0, -60.0, 10.0),
    )

    azote = Atome(
        formuleChimique="N",
        masse=2.33e-26,
        rayon=6.5e-11,
        position=(0.0, -2.5e-10, 0.0),
        vitesse=(-30.0, 50.0, -20.0),
    )

    fluor = Atome(
        formuleChimique="F",
        masse=3.15e-26,
        rayon=5.0e-11,
        position=(2.5e-10, 2.5e-10, 0.0),
        vitesse=(70.0, 0.0, 40.0),
    )

    neon = Atome(
        formuleChimique="Ne",
        masse=3.35e-26,
        rayon=3.8e-11,
        position=(-2.5e-10, 2.5e-10, 0.0),
        vitesse=(-20.0, -70.0, 0.0),
    )

    silicium = Atome(
        formuleChimique="Si",
        masse=4.66e-26,
        rayon=1.1e-10,
        position=(2.5e-10, -2.5e-10, 0.0),
        vitesse=(10.0, 10.0, -30.0),
    )

    soufre = Atome(
        formuleChimique="S",
        masse=5.32e-26,
        rayon=1.0e-10,
        position=(-2.5e-10, -2.5e-10, 0.0),
        vitesse=(-60.0, 0.0, 20.0),
    )

    argon = Atome(
        formuleChimique="Ar",
        masse=6.63e-26,
        rayon=7.1e-11,
        position=(0.0, 0.0, 2.5e-10),
        vitesse=(25.0, 45.0, 35.0),
    )

    atomes = np.array([
        hydrogene, oxygene, helium, carbone, azote,
        fluor, neon, silicium, soufre, argon,
    ])
    domaine = (1e-9, 1e-9, 1e-9)

    sim = Simulation(atomes=atomes, domaine=domaine)
    print(sim)

    # Export
    sim.open_export("export/simulation")
    try:
        sim.write_step_into_export(time=0.0)
    finally:
        sim.close_export() 