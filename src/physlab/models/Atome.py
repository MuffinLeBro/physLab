import numpy as np
import scipy.constants as cst
class Atome:
    def __init__(self, 
                 formuleChimique : str, 
                 masse : float, #kg
                 rayon : float, #m
                 position : tuple[float,float,float],   
                 vitesse : tuple[float,float,float]):
    
        self.formuleChimique : str = formuleChimique
        self.masse : float = masse
        self.rayon : float = rayon
        self.position : np.ndarray = np.array(position, dtype=float)
        self.vitesse : np.ndarray = np.array(vitesse, dtype=float)


    def __str__(self):
        return (
            f"Formule chimique : {self.formuleChimique}\n"
            f"Masse            : {self.masse} kg\n"
            f"Rayon            : {self.rayon} m\n"
            f"Position         : {self.position} m\n"
            f"Vitesse          : {self.vitesse} m/s"
        )

    def force(self) -> np.ndarray:
        #principe fondamental de la dynamique : F = m*a
        #l'atome est un projectile sans interaction : la seule force est le poids
        #F = m*g  (en N)
   
        g = np.array([0.0, -cst.g, 0.0])
        return self.masse * g
    def acceleration(self)-> np.ndarray:
        #L'accélération est la dérviée de la vitesse et la vitesse et la dérivée de la position
        # a = F/m = m*g/m = g --> la masse se simplifie donc tous les atomes ont la même accélération
        a = self.force()/self.masse
        return a
    def move(self,dt : float,) -> None: 
        #dv/dt = g = a
        #dr/dt = v
         #En isolant v(t+dt) et r(t+dt) on obtient la méthode d'Euler explicite :
        #r(t+dt) = r(t) + dt*v(t)   --> utilise l'ANCIENNE vitesse
        #v(t+dt) = v(t) + dt*a
        #
        #L'ordre des deux lignes compte : on met la position à jour AVANT la vitesse.
        #Erreur globale proportionnelle à dt (méthode d'ordre 1).
        self.position = self.position + dt * self.vitesse
        self.vitesse = self.vitesse + dt * self.acceleration()  

    def kinetic_energy(self) -> float:
        #E_cin = 1/2*m*v²  (v² = produit scalaire de v avec lui-même)
        return 0.5 * self.masse * float(np.dot(self.vitesse, self.vitesse))
 
    def potential_energy(self) -> float:
        #E_pot = m*g*y, le sol (y = 0) est la référence
        gravity : float = cst.G #constante 
        return self.masse * cst.g * self.position[1]
 
    def mechanical_energy(self) -> float:
        #E_mec = E_cin + E_pot
        #La pesanteur est une force conservative --> E_mec devrait rester constante.
        #Avec Euler explicite elle dérive : +1/2*m*g²*dt² par pas, donc
        #+1/2*m*g²*T*dt sur toute la simulation (proportionnel à dt : ordre 1).
        return self.kinetic_energy() + self.potential_energy()

    """
    9.1 :
    ================================================================
    Interaction : Collision élastique entre les 2 sphères. La force est dirigée selon la ligne qui joint les centres lors du contact
    Conservations : la quantité de movuement totale et l'énergie cinétique totale sont conservée
    Hypothèses : Gaz dilué, atomes assimilés à des sphères lisse (hard sphere, sans frottement), colllision binaires et instantanées, aucune force extérieure pendant le choc pas de degré de liberté internes
    =================================================================
    9.2 :
    ================================================================
    1. Choco frontal, masses égales : v1 = (+v,0,0), v2(-v,0,0). les vitesses s'échangent 
    2. Choco frontal, masses différentes : m2 = 4m1 On vérifie la quantité de mouvement et l'énergie cinétique.
    3. Cibles au repos, choc excentré : v2 = 0 avec un paramètre d'impact non nul. Pour des masses égales, les deux vitesses finales sont perpendiculaires. Cela contrôle l'influence de l'angle de contact
    4. Atomes qui s'éloignent, déja en contact : aucun changement de vitesses. Cela évite qu'une collision soit appliquée deux fois
    ================================================================
    """
    def en_collision(self, autre: "Atome") -> bool:
        distance = np.linalg.norm(self.position - autre.position)
        return distance <= self.rayon + autre.rayon
    def collision(self, autre: "Atome") -> None:
        n = self.position - autre.position
        n = n / np.linalg.norm(n)
        u = np.dot(self.vitesse - autre.vitesse, n)
        if u >= 0:          # ils s'éloignent déjà : rien à faire
            return
        M = self.masse + autre.masse
        self.vitesse = self.vitesse - 2 * autre.masse / M * u * n
        autre.vitesse = autre.vitesse + 2 * self.masse / M * u * n