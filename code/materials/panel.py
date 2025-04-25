import numpy as np
from ..materials.sandwich import Sandwich

class Panel:
    def __init__(self, sandwich: Sandwich, width: float, length: float):
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
        
        # Calculate basic properties
        self.area = self.width * self.length  # m²
        self.total_thickness = ((2 * self.sandwich.tf) + self.sandwich.tc)* 1e-3  # m
        self.volume = self.area * self.total_thickness  # m³
        
        # Calculate panel properties
        self.D = self.sandwich.D  # Flexural rigidity
        self.S = self.sandwich.S
        self.g = 9.80665  # m/s² (acceleration due to gravity)
        self.weight = self.calculate_weight()  # N

        self.max_bending_moment_distributedload = 0.0  # Nm
        self.max_shear_force_distributedload = 0.0  # N
        self.delta_max_distributedload = 0.0  # m

        self.max_bending_moment_pointload = 0.0  # Nm
        self.max_shear_force_pointload = 0.0  # N
        self.delta_max_pointload = 0.0


        self.sigma_max = 0.0  # Pa
        


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
    
    def calculate_beam_response_distributed_load(self, load: float) -> tuple:
        """Calculate midlength deflection, maximum shear force, and maximum bending moment 
        by treating the panel as a beam with uniform line load.
        
        Parameters:
        ----------- 
        load : float
            Uniform load in N/m²
            
        Returns:
        --------
        tuple
            (midlength_deflection, max_shear_force, max_bending_moment) in meters, Newtons, and Newton-meters
        """
        # Treat the panel as a beam with line load
        q_line = load 
        L = self.length  # Beam length (m)
        D = self.D  # Flexural rigidity (Nm²)
        S = self.S  # Shear rigidity (Nm)
        
        # Maximum deflection at midlength for simply supported beam with uniform line load
        delta_max_bend = (5 * q_line * L**4) / (384 * D)
        delta_max_shear = (q_line * L**2) / (8 * S)  # Shear deflection

        delta_max = delta_max_bend + delta_max_shear  # Total deflection (m)
        
        # Maximum shear force (occurs at the supports)
        max_shear_force = q_line * L / 2  # N
        
        # Maximum bending moment (occurs at the center)
        max_bending_moment = q_line * L**2 / 8  # Nm
        
        return delta_max, max_shear_force, max_bending_moment

    def calculate_beam_response_point_load(self, load: float) -> tuple:
        """Calculate midlength deflection, maximum shear force, and maximum bending moment 
        by treating the panel as a beam with a point load at the center.
        
        Parameters:
        ----------- 
        load : float
        Point load in Newtons (N)
        
        Returns:
        --------
        tuple
        (midlength_deflection, max_shear_force, max_bending_moment) in meters, Newtons, and Newton-meters
        """
        L = self.length  # Beam length (m)
        D = self.D  # Flexural rigidity (Nm²)
        S = self.S  # Shear rigidity (Nm²)

        load /= self.width  # Convert point load to line load (N/m)
        
        # Maximum deflection at midlength for simply supported beam with point load at center
        delta_max_bend = (load * L**3) / (48 * D)
        delta_max_shear = (load * L) / (4*S)

        delta_max = delta_max_bend + delta_max_shear  # Total deflection (m)
        
        # Maximum shear force (occurs at the supports)
        max_shear_force = load / 2  # N
        
        # Maximum bending moment (occurs at the center)
        max_bending_moment = load * L / 4  # Nm
        
        return delta_max, max_shear_force, max_bending_moment

    def check_against_face_failure(self,  point_load: float, distributed_load: float) -> float:
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
        max_bending_moment = max(self.max_bending_moment_distributedload, self.max_bending_moment_pointload)  # Nm
    
        d = self.sandwich.d  # m
        tf = self.sandwich.tf  # Face thickness (m)
        sigma_f_max = self.sandwich.composite_material.sigma_f_max  # Maximum allowable face stress (Pa)
        
        # Maximum normal stress in faces
        sigma_max = (max_bending_moment * d) / (2 * self.D)
        
        # Safety factor for face failure
        safety_factor = sigma_f_max / abs(sigma_max)
        return safety_factor

    def check_against_core_shear_failure(self, point_load: float, distributed_load: float) -> float:
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
        max_shear_force = max(self.max_shear_force_distributedload, self.max_shear_force_pointload)  # N

        d = self.sandwich.d  
        tau_c_max = self.sandwich.core_material.tau_c_max  # Maximum allowable core shear stress (Pa)
        
        # Maximum shear stress in the core
        tau_max = max_shear_force / d
        
        # Safety factor for core shear failure
        safety_factor = tau_c_max / abs(tau_max)
        return safety_factor

    def check_against_face_wrinkling(self, point_load: float, distributed_load: float) -> float:
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
        
        # sigma_cr = 0.5 * (Ef * Ec * Gc)**(1/3)  # Critical wrinkling stress (Pa) WRONG FORMULA
        
        # Maximum normal stress in faces due to bending
        max_bending_moment = max(self.max_bending_moment_distributedload, self.max_bending_moment_pointload)

        d = self.sandwich.d  # m
        sigma_max = (max_bending_moment * d) / (2 * self.D)
        
        # Safety factor for face wrinkling
        safety_factor = sigma_cr / abs(sigma_max)
        return safety_factor

    def calculate_max_global_stress(self, load: float) -> tuple:
        """Calculate maximum stresses in the faces.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        
        Returns:
        --------
        tuple
            (max_normal_stress, max_shear_stress) in Pa
        """
        q = load  # N/m²
        a = self.length  # m
        d = self.sandwich.d  # m
        
        # Maximum bending moment for simply supported plate
        max_bending_moment = max(self.max_bending_moment_distributedload, self.max_bending_moment_pointload)

        # Maximum normal stress in faces
        sigma_max = (M_max * d) / (2 * self.D)
        
        # Maximum shear stress in core
        tau_max = (q * a) / (2 * d)
        
        return sigma_max, tau_max

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
        return ret
