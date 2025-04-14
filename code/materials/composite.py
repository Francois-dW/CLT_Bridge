from typing import List, Dict, Tuple
from .base import BaseMaterial
from .isotropic import IsotropicMaterial
import numpy as np


class CompositeLamina(BaseMaterial):
    """Single ply/unidirectional composite layer."""
    
    def __init__(
        self,
        fiber_material: IsotropicMaterial,
        matrix_material: IsotropicMaterial,
        fiber_volume_ratio: float,
        orientation_angle: float = 0.0  # in degrees
    ) -> None:
        """
        Args:
            fiber_material: Fiber material properties
            matrix_material: Matrix material properties
            fiber_volume_ratio: V_f (fiber volume/total volume)
            orientation_angle: Ply angle relative to principal axis (degrees)
        """
        self.fiber_material = fiber_material
        self.matrix_material = matrix_material
        self.fiber_volume_ratio = fiber_volume_ratio
        self.orientation_angle = orientation_angle

    @property
    def calculated_properties(self) -> Dict[str, float]:
        """Calculate and return material properties."""
        E1 = self.fiber_material.youngs_modulus * self.fiber_volume_ratio + \
            self.matrix_material.youngs_modulus * (1 - self.fiber_volume_ratio)
        E2 = 1 / ( 1 /self.matrix_material.youngs_modulus * (1 - self.fiber_volume_ratio) + \
            1/ (self.fiber_material.youngs_modulus * self.fiber_volume_ratio))
        G12 = 1/ (1/(self.fiber_material.shear_modulus * self.fiber_volume_ratio) + \
            1/(self.matrix_material.shear_modulus * (1 - self.fiber_volume_ratio)))
        poissons_ratio = self.fiber_material.poissons_ratio * self.fiber_volume_ratio + \
            self.matrix_material.poissons_ratio * (1 - self.fiber_volume_ratio)
        return {
            'E1': E1,
            'E2': E2,
            'G12': G12,
            'nu12': poissons_ratio
        }

    @property
    def stiffness_matrix(self) -> List[List[float]]:
        """Get reduced stiffness matrix [Q] for the ply."""
        Q = self.calculated_properties
        E1 = Q['E1']
        E2 = Q['E2']
        G12 = Q['G12']
        nu12 = Q['nu12']
        
        return np.array([
            [E1, nu12 * E1, 0],
            [nu12 * E2, E2, 0],
            [0, 0, G12]
        ])

    def transformed_stiffness_matrix(self) -> List[List[float]]:
        """Get transformed reduced stiffness matrix [Q] for the ply."""
        Q = self.calculated_properties
        E1 = Q['E1']
        E2 = Q['E2']
        G12 = Q['G12']
        nu12 = Q['nu12']
        
        # Transformation matrix calculation
        c = np.cos(self.orientation_angle)
        s = np.sin(self.orientation_angle)
        
        Q11 = E1 / (1 - nu12**2)
        Q22 = E2 / (1 - nu12**2)
        Q12 = nu12 * E2 / (1 - nu12**2)
        
        return np.array([
            [Q11, Q12, 0],
            [Q12, Q22, 0],
            [0, 0, G12]
        ])

    def to_dict(self) -> Dict[str, Any]:
        pass


class CompositeLaminate(BaseMaterial):
    """Stacked composite laminate consisting of multiple plies."""
    
    def __init__(
        self,
        plies: List[CompositeLamina],
        ply_thicknesses: List[float]
    ) -> None:
        """
        Args:
            plies: List of CompositeLamina objects in stacking order
            ply_thicknesses: Corresponding thicknesses for each ply
        """
        self.plies = plies
        self.ply_thicknesses = ply_thicknesses

    @property
    def total_thickness(self) -> float:
        """Total thickness of the laminate."""
        pass

    @property
    def stacking_sequence(self) -> List[float]:
        """Get orientation angles of all plies in order."""
        pass

    def abcd_matrix(self) -> Tuple[
        List[List[float]],  # A matrix
        List[List[float]],  # B matrix
        List[List[float]]   # D matrix
    ]:
        """Calculate ABD matrix for the laminate."""
        pass

    @property
    def flexural_rigidity(self) -> float:
        """Calculate equivalent flexural rigidity of the laminate."""
        pass

    def to_dict(self) -> Dict[str, Any]:
        pass