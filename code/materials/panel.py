import numpy as np
from ..materials.sandwich import Sandwich

class Panel:
    def __init__(self, sandwich: Sandwich, width: float, length: float, point_load: float = 0.0, distributed_load: float = 0.0):
        """Initialize a panel made of sandwich material.
        
        Parameters:
        -----------
        sandwich : Sandwich
            The sandwich material that makes up the panel
        width : float (m)
            Width of the panel in meters
        length : float (m)
            Length of the panel in meters
        """
        self.sandwich = sandwich
        self.width = width  # m
        self.length = length  # m
        self.point_load = point_load  # N
        self.distributed_load = distributed_load  # N/m²
        
        # Calculate basic properties
        self.area = self.width * self.length  # m²
        self.total_thickness = ((2 * self.sandwich.tf) + self.sandwich.tc)* 1e-3  # m
        self.volume = self.area * self.total_thickness  # m³
        
        # Calculate panel properties
        self.D = self.sandwich.D  # Flexural rigidity
        self.S = self.sandwich.S
        self.g = 9.80665  # m/s² (acceleration due to gravity)
        self.weight = self.calculate_weight()  # N

        # Calculate and store distributed load response
        delta_max, max_shear_force, max_bending_moment = self.calculate_beam_response_distributed_load()
        self.max_bending_moment_distributed_load = max_bending_moment  # Nm
        self.max_shear_force_distributed_load = max_shear_force  # N
        self.delta_max_distributed_load = delta_max  # m

        # Calculate and store point load response
        delta_max, max_shear_force, max_bending_moment = self.calculate_beam_response_point_load()
        self.max_bending_moment_point_load = max_bending_moment  # Nm
        self.max_shear_force_point_load = max_shear_force # N
        self.delta_max_point_load = delta_max  # m


    def calculate_weight(self) -> float:
        """Calculate the total weight of the panel in Newtons.
        
        Returns:
        --------
        float
            Weight of the panel in Newtons
        """
        g = self.g  # m/s² (acceleration due to gravity)
        
        # Calculate volumes
        face_volume = self.area * self.sandwich.tf * 1e-3  # m³
        core_volume = self.area * self.sandwich.tc * 1e-3  # m³
        
        # Calculate masses kg/m³
        face_mass = 2 * face_volume * self.sandwich.composite_material.rho   # kg
        core_mass = core_volume * self.sandwich.core_material.rho  # kg
        
        total_weight = (face_mass + core_mass) * g  # N 
        return total_weight
    
    def calculate_beam_response_distributed_load(self) -> tuple:
        """Calculate midlength deflection, maximum shear force, and maximum bending moment 
        by treating the panel as a beam with its own distributed load.
        
        Returns:
        --------
        tuple
            (midlength_deflection, max_shear_force, max_bending_moment) in meters, Newtons, and Newton-meters
        """
        # Use the panel's distributed load
        q_line = self.distributed_load
        L = self.length  # Beam length (m)
        D = self.D  # Flexural rigidity (Nm²)
        S = self.S  # Shear rigidity (Nm²)
        
        # Maximum deflection at midlength for simply supported beam with uniform line load
        delta_max_bend = (5 * q_line * L**4) / (384 * D)
        delta_max_shear = (q_line * L**2) / (8 * S)  # Shear deflection

        delta_max = delta_max_bend + delta_max_shear  # Total deflection (m)
        
        # Maximum shear force (occurs at the supports)
        max_shear_force = q_line * L / 2  # N
        
        # Maximum bending moment (occurs at the center)
        max_bending_moment = q_line * L**2 / 8  # Nm
        
        return delta_max, max_shear_force, max_bending_moment

    def calculate_beam_response_point_load(self) -> tuple:
        """Calculate midlength deflection, maximum shear force, and maximum bending moment 
        by treating the panel as a beam with a point load at the center.

        Returns:
        --------
        tuple
            (midlength_deflection, max_shear_force, max_bending_moment) in meters, Newtons, and Newton-meters
        """
        L = self.length  # Beam length (m)
        D = self.D  # Flexural rigidity (Nm²)
        S = self.S  # Shear rigidity (Nm²)
        load = self.point_load  # Point load in Newtons (N)

        load_per_width = load / self.width if self.width != 0 else 0  # Convert point load to line load (N/m)

        # Maximum deflection at midlength for simply supported beam with point load at center
        delta_max_bend = (load_per_width * L**3) / (48 * D)
        delta_max_shear = (load_per_width * L) / (4 * S)

        delta_max = delta_max_bend + delta_max_shear  # Total deflection (m)

        # Maximum shear force (occurs at the supports)
        max_shear_force = load_per_width / 2  # N

        # Maximum bending moment (occurs at the center)
        max_bending_moment = load_per_width * L / 4  # Nm

        return delta_max, max_shear_force, max_bending_moment

    def check_against_face_failure(self) -> float:
        """
        Calculate the safety factor for face failure due to normal stress.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        
        Returns:
        --------
        float
            Safety factor for face failure
        """
        max_bending_moment = max(self.max_bending_moment_distributed_load, self.max_bending_moment_point_load)  # Nm
    
        d = self.sandwich.d  # m
        tf = self.sandwich.tf  # Face thickness (m)
        sigma_f_max = self.sandwich.composite_material.sigma_f_max  # Maximum allowable face stress (Pa) TODO: Check if this is correct
        
        # Maximum normal stress in faces
        sigma_max = (max_bending_moment * d) / (2 * self.D)
        
        # Safety factor for face failure
        safety_factor = sigma_f_max / abs(sigma_max)
        return safety_factor

    def check_against_core_shear_failure(self,) -> float:
        """
        Calculate the safety factor for core shear failure.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        
        Returns:
        --------
        float
            Safety factor for core shear failure
        """
        max_shear_force = max(self.max_shear_force_distributed_load, self.max_shear_force_point_load)  # N

        d = self.sandwich.d  
        tau_c_max = self.sandwich.core_material.tau_c_max  # Maximum allowable core shear stress (Pa) TODO: Check if this is correct
        
        # Maximum shear stress in the core
        tau_max = max_shear_force / d
        
        # Safety factor for core shear failure
        safety_factor = tau_c_max / abs(tau_max)
        return safety_factor

    def check_against_face_wrinkling(self) -> float:
        """
        Calculate the safety factor for face wrinkling failure.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        
        Returns:
        --------
        float
            Safety factor for face wrinkling failure
        """
        # Calculate the critical wrinkling stress
        Ef = self.sandwich.composite_material.Ef  # Face modulus of elasticity (Pa)
        Ec = self.sandwich.core_material.Ec  # core_material modulus of elasticity (Pa)
        Gc = self.sandwich.core_material.Gc  # core_material shear modulus (Pa)
        
        sigma_cr = 0.5**3 * np.sqrt(Ef * Ec * Gc)  # Critical wrinkling stress (Pa) WRONG FORMULA
        
        # Maximum normal stress in faces due to bending
        max_bending_moment = max(self.max_bending_moment_distributed_load, self.max_bending_moment_point_load)

        d = self.sandwich.d  # m
        sigma_max = (max_bending_moment * d) / (2 * self.D)
        
        # Safety factor for face wrinkling
        safety_factor = sigma_cr / abs(sigma_max)
        return safety_factor
    
    def check_against_core_compression_failure(self) -> float:
        """
        Calculate the safety factor for core compression failure.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        
        Returns:
        --------
        float
            Safety factor for core compression failure
        """
        point_load = self.point_load
        sigma_c_comp = self.sandwich.core_material.sigma_comp
        min_A = point_load / sigma_c_comp

        return min_A
m
    def __str__(self):
        ret = f"Panel:\n"
        ret += f"  Width: {self.width:.2f} m\n"
        ret += f"  Length: {self.length:.2f} m\n"
        ret += f"  Area: {self.area:.2f} m²\n"
        ret += f"  Total Thickness: {self.total_thickness:.2f} m\n"
        ret += f"  Volume: {self.volume:.2f} m³\n"
        ret += f"  Weight: {self.weight:.2f} N\n"
        ret += f"  Flexural Rigidity (D): {self.D:.2f} Nm²\n"
        ret += f"  Shear Rigidity (S): {self.sandwich.S:.2f} Nm²\n"
        ret += f"  Maximum Bending Moment (Distributed Load): {self.max_bending_moment_distributed_load:.2f} Nm\n"
        ret += f"  Maximum Shear Force (Distributed Load): {self.max_shear_force_distributed_load:.2f} N\n"
        ret += f"  Maximum Deflection (Distributed Load): {self.delta_max_distributed_load:.2f} m\n"
        ret += f"  Maximum Bending Moment (Point Load): {self.max_bending_moment_point_load:.2f} Nm\n"
        ret += f"  Maximum Shear Force (Point Load): {self.max_shear_force_point_load:.2f} N\n"
        ret += f"  Maximum Deflection (Point Load): {self.delta_max_point_load:.2f} m\n"
        return ret
