import numpy as np
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
    def moove(masse : float,) -> float: 
        #principe fondamental de la dynamique : F = m*a (a=g=9.81m/s²) --> F = m*g
        
        
        pass
