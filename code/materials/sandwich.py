import numpy as np
#         --------

from ..materials.isotropic import IsotropicMaterial
from ..materials.composite import Laminate, Lamina, LaminaMix

class Sandwich:
    def __init__(self, Composite: Laminate, Core: IsotropicMaterial, tc = None, tf = None):
        """Initializes the Sandwich class with the given composite and core materials.

        Parameters:
        ----------
        Composite: Material
            The composite material of the sandwich.
        Core: Material
            The core material of the sandwich.
        tc: float, optional
            Thickness of the core (default is None).
        tf: float, optional
            Thickness of the foam (default is None).
        """
        self.Composite = Composite
        self.Core = Core
        if tc is None:
            self.tc = Core.thickness
        else:
            self.tc = tc
        if tf is None:
            self.tf = Composite.thickness
        else:
            self.tf = tf
        if self.Gc is None:
            self.Gc = Core.E / (2 * (1 + Core.nu))  # Shear modulus of the core
        else:
            self.Gc = Core.G
        self.d = self.tf + self.tc  # distance between the core and the foam
        self.Ef = Composite.E  # Young's modulus of the foam
        self.Ec = Core.E  # Young's modulus of the core
        self.D = self.get_D()  # D value for the sandwich with a soft core


    def get_D(self):
        """get the D value for the sandwich
        
        Parameters:
        ----------
        self.Ef: float
            Young's modulus of the foam
        self.tf: float
            thickness of the foam
        self.d: float
            distance between the core and the foam
        self.Ec: float
            Young's modulus of the core
        self.tc: float
            thickness of the core
            
        Returns:
        -------
        D: float
            D value for the sandwich with a soft core
        """
        Df = (self.Ef * self.tf ** 3) / 6
        D0 = (self.Ef * self.tf * self.d ** 2) / 2
        Dc = (self.Ec * self.tc ** 3) / 12
        D = 2*Df + D0 + Dc
        return D
    
    def get_S(self):
        """get the S value for the sandwich for a soft core and thin skin
        
        Parameters:
        ----------
        self.Gc: float
            Shear modulus of the core
        self.d: float
            distance between the core and the foam
        self.tc: float
            thickness of the core
        
        Returns:
        -------
        S: float
            S value for the sandwich with a soft core
        """
        S = (self.Gc * self.d ** 2) / self.tc
        return S

    
