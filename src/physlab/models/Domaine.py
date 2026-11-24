import numpy as np
class Domaine:
    def __init__(self, longX : float,longY : float, longZ : float):
        self.longX = longX
        self.longY = longY
        self.longZ = longZ

    def __str__(self):
        return f"Volume total de la boite : {self.Volume()} m^3"
    
    def Volume(self) -> float :
        vTot = self.longX * self.longY * self.longZ
        return vTot
    def reflection(self, atome) -> None:
    #Réflexion spéculaire : corrige position et vitesse, axe par axe."""
        demi = (self.longX / 2, self.longY / 2, self.longZ / 2)
        for i in range(3):
            if atome.position[i] > demi[i]:
                atome.position[i] = 2 * demi[i] - atome.position[i]   # symétrie par rapport au mur
                atome.vitesse[i] = -atome.vitesse[i]
            elif atome.position[i] < -demi[i]:
                atome.position[i] = -2 * demi[i] - atome.position[i]
                atome.vitesse[i] = -atome.vitesse[i]