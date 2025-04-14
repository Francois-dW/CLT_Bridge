from typing import Dict, Any, Optional
from .base import BaseMaterial


class IsotropicMaterial(BaseMaterial):
    """Isotropic material implementation."""
    
    def __init__(
        self,
        youngs_modulus: float,
        poissons_ratio: float,
        density: Optional[float] = None
    ) -> None:
        self.youngs_modulus = youngs_modulus
        self.poissons_ratio = poissons_ratio
        self.density = density

    def to_dict(self) -> Dict[str, Any]:
        pass