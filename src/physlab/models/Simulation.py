import os

import h5py as h5
import numpy as np

from physlab.models.Atome import Atome


class Simulation:
    def __init__(self, atomes: np.ndarray, domaine: tuple[float, float, float]):
        self.domaine: np.ndarray = np.array(domaine, dtype=float)
        self.atomes: np.ndarray = np.array(atomes)

        # Export HDF5 / XDMF
        self.h5_file = None
        self.xdmf_file = None
        self.h5_name = ""
        self.step_count = 0

    def __str__(self):
        result = ""
        for atom in self.atomes:
            result += str(atom) + "\n"
        return result

    # ------------------------------------------------------------------
    # Export HDF5 / XDMF
    # ------------------------------------------------------------------
    def _domain_corners(self) -> np.ndarray:
        """
        Les 8 sommets du domaine, centré sur l'origine :
        [-lx/2, lx/2] x [-ly/2, ly/2] x [-lz/2, lz/2], dans l'ordre
        attendu par la topologie Hexahedron de XDMF (face du bas dans le
        sens trigonométrique, puis face du haut).
        """
        hx, hy, hz = self.domaine / 2
        return np.array(
            [
                (-hx, -hy, -hz),
                (hx, -hy, -hz),
                (hx, hy, -hz),
                (-hx, hy, -hz),
                (-hx, -hy, hz),
                (hx, -hy, hz),
                (hx, hy, hz),
                (-hx, hy, hz),
            ],
            dtype=np.float64,
        )

    def open_export(self, filename: str) -> None:
        """
        Crée les fichiers <filename>.h5 et <filename>.xdmf (ainsi que les
        dossiers manquants) et écrit l'en-tête XDMF de la série temporelle,
        y compris le domaine de simulation.

        Parameters:
            filename (str): chemin sans extension, ex: "export/simulation"
        """
        if self.h5_file is not None or self.xdmf_file is not None:
            raise RuntimeError("Un export est déjà ouvert, appelez close_export().")

        base, _ = os.path.splitext(filename)
        directory = os.path.dirname(base)
        if directory:
            os.makedirs(directory, exist_ok=True)

        h5_path = base + ".h5"
        xdmf_path = base + ".xdmf"
        # Le .xdmf référence le .h5 par son nom seul (même dossier)
        self.h5_name = os.path.basename(h5_path)
        self.step_count = 0

        corners = self._domain_corners()

        # --- Fichier .h5 : informations générales + domaine ---
        self.h5_file = h5.File(h5_path, "w")
        self.h5_file.attrs["description"] = "Simulation de mouvement d'atomes"
        self.h5_file.attrs["number_of_atoms"] = len(self.atomes)

        domain_group = self.h5_file.create_group("domain")
        domain_group.attrs["lx"] = self.domaine[0]
        domain_group.attrs["ly"] = self.domaine[1]
        domain_group.attrs["lz"] = self.domaine[2]
        domain_group.attrs["volume"] = float(np.prod(self.domaine))
        domain_group.create_dataset("corners", data=corners)
        self.h5_file.flush()

        # --- Fichier .xdmf : en-tête ---
        # Le domaine est une grille statique, lue depuis le .h5
        h5n = self.h5_name
        self.xdmf_file = open(xdmf_path, "w", encoding="utf-8")
        self.xdmf_file.write(
            f"""<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="boundary" GridType="Uniform">
      <Topology TopologyType="Hexahedron" NumberOfElements="1">
        <DataItem Dimensions="1 8" NumberType="Int" Format="XML">
          0 1 2 3 4 5 6 7
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="8 3" NumberType="Float" Precision="8" Format="HDF">{h5n}:/domain/corners</DataItem>
      </Geometry>
    </Grid>
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
"""
        )
        self.xdmf_file.flush()

    def write_step_into_export(self, time: float = 0.0) -> None:
        """
        Sauvegarde l'état courant (positions, vitesses, masses, rayons) dans
        un nouveau groupe step_XXXX du .h5 et ajoute le bloc Grid
        correspondant dans le .xdmf.

        Parameters:
            time (float): instant physique de cette image
        """
        if self.h5_file is None or self.xdmf_file is None:
            raise RuntimeError("Aucun export ouvert, appelez open_export() d'abord.")

        n = len(self.atomes)
        name = f"step_{self.step_count:04d}"

        positions = np.array(
            [a.position for a in self.atomes], dtype=np.float64
        ).reshape(n, 3)
        velocities = np.array(
            [a.vitesse for a in self.atomes], dtype=np.float64
        ).reshape(n, 3)
        masses = np.array([a.masse for a in self.atomes], dtype=np.float64)
        radii = np.array([a.rayon for a in self.atomes], dtype=np.float64)

        # --- .h5 ---
        group = self.h5_file.create_group(name)
        group.attrs["time"] = time
        group.create_dataset("positions", data=positions)
        group.create_dataset("velocity", data=velocities)
        group.create_dataset("mass", data=masses)
        group.create_dataset("radius", data=radii)
        self.h5_file.flush()

        # --- .xdmf ---
        h5n = self.h5_name
        self.xdmf_file.write(
            f"""      <Grid Name="{name}" GridType="Uniform">
        <Time Value="{time}"/>
        <Topology TopologyType="Polyvertex" NumberOfElements="{n}" NodesPerElement="1"/>
        <Geometry GeometryType="XYZ">
          <DataItem Dimensions="{n} 3" NumberType="Float" Precision="8" Format="HDF">{h5n}:/{name}/positions</DataItem>
        </Geometry>
        <Attribute Name="mass" AttributeType="Scalar" Center="Node">
          <DataItem Dimensions="{n}" NumberType="Float" Precision="8" Format="HDF">{h5n}:/{name}/mass</DataItem>
        </Attribute>
        <Attribute Name="radius" AttributeType="Scalar" Center="Node">
          <DataItem Dimensions="{n}" NumberType="Float" Precision="8" Format="HDF">{h5n}:/{name}/radius</DataItem>
        </Attribute>
        <Attribute Name="velocity" AttributeType="Vector" Center="Node">
          <DataItem Dimensions="{n} 3" NumberType="Float" Precision="8" Format="HDF">{h5n}:/{name}/velocity</DataItem>
        </Attribute>
      </Grid>
"""
        )
        self.xdmf_file.flush()
        self.step_count += 1

    def close_export(self) -> None:
        """Ferme les balises XML ouvertes puis les fichiers .xdmf et .h5."""
        if self.xdmf_file is not None:
            self.xdmf_file.write("    </Grid>\n  </Domain>\n</Xdmf>\n")
            self.xdmf_file.close()
            self.xdmf_file = None
        if self.h5_file is not None:
            self.h5_file.close()
            self.h5_file = None