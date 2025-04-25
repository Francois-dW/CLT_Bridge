import numpy as np
#         --------

from ..materials.isotropic import IsotropicMaterial
from ..materials.composite import Laminate, Lamina, LaminaMix

class Sandwich:
    def __init__(self, composite_material: Laminate, core_material: IsotropicMaterial, tc = None, tf = None):
        """Initializes the Sandwich class with the given composite and core materials.

        Parameters:
        ----------
        composite_material: Material
            The composite material of the sandwich.
        core_material: Material
            The core material of the sandwich.
        tc: float, optional (mm)
            Thickness of the core in mm(default is None).
        tf: float, optional (mm)
            Thickness of the foam in mm (default is None).
        """
        self.composite_material = composite_material
        self.core_material = core_material
        if tc is None:
            self.tc = core_material.thickness
        else:
            self.tc = tc
        if tf is None:
            self.tf = composite_material.total_thickness
            print(f"tf is None, using composite_material thickness: {self.tf} mm")
        else:
            self.tf = tf
        if not self.core_material.G:
            self.Gc = core_material.E / (2 * (1 + core_material.nu))  # Shear modulus of the core
        else:
            self.Gc = core_material.G
        self.d = self.tf + self.tc  # distance between the core and the foam
        self.Ef = core_material.E  # Young's modulus of the foam
        self.Ec = composite_material.E  # Young's modulus of the core
        self.D = self.get_D()  # D value for the sandwich with a soft core
        self.S = self.get_S()  # S value for the sandwich with a soft core


    def get_D(self):
        """get the D value for the sandwich
        
        Parameters:
        ----------
        self.Ef: float
            Young's modulus of the foam in Pa
        self.tf: float
            thickness of the foam in mm
        self.d: float
            distance between the core and the foam in mm
        self.Ec: float
            Young's modulus of the core in Pa
        self.tc: float
            thickness of the core in mm
            
        Returns:
        -------
        D: float
            D value for the sandwich with a soft core in Nm^2
        """
        Df = (self.Ef * (self.tf/1000) ** 3) / 6
        D0 = (self.Ef * (self.tf/1000) * (self.d/1000) ** 2) / 2
        Dc = (self.Ec * (self.tc/1000) ** 3) / 12
        D = 2*Df + D0 + Dc
        
        return D
    
    def get_S(self):
        """get the S value for the sandwich for a soft core and thin skin
        
        Parameters:
        ----------
        self.Gc: float
            Shear modulus of the core in Pa
        self.d: float
            distance between the core and the foam in mm
        self.tc: float
            thickness of the core in mm
        
        Returns:
        -------
        S: float
            S value for the sandwich with a soft core in Nm
        """
        S = (self.Gc * (self.d/1000) ** 2) / (self.tc/1000)
        return S

    
