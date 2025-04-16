import numpy as np


class Sandwich:
    def __init__(self, Composite, Core, tc = None, tf = None):
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
        self.d = self.tf + self.tc  # distance between the core and the foam
        self.Ef = Composite.E  # Young's modulus of the foam
        self.Ec = Core.E  # Young's modulus of the core
        self.D = self.get_D(self.Ef, self.tf, self.d, self.Ec, self.tc)  # D value for the sandwich with a soft core



    def get_D(Ef, tf, d, Ec, tc):
        """get the D value for the sandwich
        
        Parameters:
        ----------
        Ef: float
            Young's modulus of the foam
        tf: float
            thickness of the foam
        d: float
            distance between the core and the foam
        Ec: float
            Young's modulus of the core
        tc: float
            thickness of the core
            
        Returns:
        -------
        D: float
            D value for the sandwich with a soft core
        """
        Df = (Ef * tf ** 3) / 6
        D0 = (Ef * tf * d ** 2) / 2
        Dc = (Ec * tc ** 3) / 12
        D = 2*Df + D0 + Dc
        return D

    def get_D_soft(Ef, tf, d):
        """get the D value for the sandwich with a soft core
        
        Parameters:
        ----------
        Ef: float
            Young's modulus of the foam
        tf: float
            thickness of the foam
        d: float
            distance between the core and the foam
        
        Returns:
        -------
        D: float
            D value for the sandwich with a soft core
        """
        D0 = (Ef * tf * d ** 2) / 2
        D = D0
        return D

    def get_D_no_core(Ef, tf):
        """get the D value for the sandwich without a core
        
        Parameters:
        ----------
        Ef: float
            Young's modulus of the foam
        tf: float
            thickness of the foam
        
        Returns:
        -------
        D: float
            D value for the sandwich with a soft core
        """
        Df = (Ef * tf ** 3) / 6
        D = 2*Df
        return D

    def get_D_small(Ef, tf , d):
        """get the D value for the sandwich with a small core
        
        Parameters:
        ----------
        Ef: float
            Young's modulus of the foam
        tf: float
            thickness of the foam
        d: float
            distance between the core and the foam
        
        Returns:
        -------
        D: float
            D value for the sandwich with a small core
        """
        Df = (Ef * tf ** 3) / 48
        D0 = (Ef * tf * d ** 2) / 2
        D = Df + D0
        return D