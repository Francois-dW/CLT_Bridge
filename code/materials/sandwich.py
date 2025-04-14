from typing import Dict, Any
from .base import BaseMaterial
from .isotropic import IsotropicMaterial


class SandwichPanel(BaseMaterial):
    """Sandwich composite material."""
    
    def __init__(
        self,
        face_material: IsotropicMaterial,
        core_material: IsotropicMaterial,
        face_thickness: float,
        core_thickness: float
    ) -> None:
        self.face_material = face_material
        self.core_material = core_material
        self.face_thickness = face_thickness
        self.core_thickness = core_thickness

    @property
    def total_thickness(self) -> float:
        pass

    @property
    def flexural_rigidity(self) -> float:
        pass

    def to_dict(self) -> Dict[str, Any]:
        pass