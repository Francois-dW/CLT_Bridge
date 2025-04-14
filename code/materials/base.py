from abc import ABC, abstractmethod
from typing import Dict, Any, List


class BaseMaterial(ABC):
    """Abstract base class for all materials."""
    
    @property
    @abstractmethod
    def youngs_modulus_x(self) -> float:
        """Young's modulus in principal x-direction."""
        pass

    @property
    @abstractmethod
    def youngs_modulus_y(self) -> float:
        """Young's modulus in principal y-direction."""
        pass

    @property
    @abstractmethod
    def shear_modulus_xy(self) -> float:
        """In-plane shear modulus."""
        pass

    @property
    @abstractmethod
    def poissons_ratio_xy(self) -> float:
        """Major Poisson's ratio."""
        pass


    @abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        """Serialize material properties to dictionary."""
        pass