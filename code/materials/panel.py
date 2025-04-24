import numpy as np
from ..materials.sandwich import Sandwich

class Panel:
    def __init__(self, sandwich: Sandwich, width: float, length: float):
        """Initialize a panel made of sandwich material.
        
        Parameters:
        -----------
        sandwich : Sandwich
            The sandwich material that makes up the panel
        width : float
            Width of the panel in meters
        length : float
            Length of the panel in meters
        """
        self.sandwich = sandwich
        self.width = width  # m
        self.length = length  # m
        
        # Calculate basic properties
        self.area = self.width * self.length  # m²
        self.total_thickness = (2 * self.sandwich.tf) + self.sandwich.tc  # m
        self.volume = self.area * self.total_thickness  # m³
        
        # Calculate panel properties
        self.D = self.sandwich.D  # Flexural rigidity
        self.weight = self.calculate_weight()  # N
        
    def calculate_weight(self) -> float:
        """Calculate the total weight of the panel in Newtons.
        
        Returns:
        --------
        float
            Weight of the panel in Newtons
        """
        g = 9.81  # m/s²
        
        # Calculate volumes
        face_volume = self.area * self.sandwich.tf  # m³
        core_volume = self.area * self.sandwich.tc  # m³
        
        # Calculate masses (convert density from g/cm³ to kg/m³)
        face_mass = 2 * face_volume * (self.sandwich.Composite.rho * 1000)  # kg
        core_mass = core_volume * (self.sandwich.Core.rho * 1000)  # kg
        
        total_weight = (face_mass + core_mass) * g  # N
        return total_weight
    
    def calculate_beam_response_uniform_load(self, load: float) -> tuple:
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
        q_line = load * self.width  # Convert surface load (N/m²) to line load (N/m)
        L = self.length  # Beam length (m)
        D = self.D  # Flexural rigidity (Nm²)
        S = self.sandwich.S  # Shear rigidity (Nm²)
        
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
            S = self.sandwich.S  # Shear rigidity (Nm²)
            
            # Maximum deflection at midlength for simply supported beam with point load at center
            delta_max_bend = (load * L**3) / (48 * D)
            delta_max_shear = (load * L) / 4*S

            delta_max = delta_max_bend + delta_max_shear  # Total deflection (m)
            
            # Maximum shear force (occurs at the supports)
            max_shear_force = load / 2  # N
            
            # Maximum bending moment (occurs at the center)
            max_bending_moment = load * L / 4  # Nm
            
            return max_bending_moment, max_shear_force, delta_max

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
        max_bending_moment, _, _ = self.calculate_beam_response_uniform_load(point_load)

        d = self.sandwich.d  # m
        tf = self.sandwich.tf  # Face thickness (m)
        sigma_f_max = self.sandwich.Composite.sigma_f_max  # Maximum allowable face stress (Pa)
        
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
        _, max_shear_force, _ = self.calculate_beam_response_uniform_load(load)

        d = self.sandwich.d  
        tau_c_max = self.sandwich.Core.tau_c_max  # Maximum allowable core shear stress (Pa)
        
        # Maximum shear stress in the core
        tau_max = max_shear_force / d
        
        # Safety factor for core shear failure
        safety_factor = tau_c_max / abs(tau_max)
        return safety_factor

        def check_against_face_wrinkling(self, load: float) -> float:
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
            Ef = self.sandwich.Composite.Ef  # Face modulus of elasticity (Pa)
            Ec = self.sandwich.Core.Ec  # Core modulus of elasticity (Pa)
            Gc = self.sandwich.Core.Gc  # Core shear modulus (Pa)
            
            #sigma_cr = 0.5 * (Ef * Ec * Gc)**(1/3)  # Critical wrinkling stress (Pa) WRONG FORMULA
            
            # Maximum normal stress in faces due to bending
            max_bending_moment, _, _ = self.calculate_beam_response_uniform_load(load)
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
        M_max = (q * a**2) / 8
        
        # Maximum normal stress in faces
        sigma_max = (M_max * d) / (2 * self.D)
        
        # Maximum shear stress in core
        tau_max = (q * a) / (2 * d)
        
        return sigma_max, tau_max

    def __str__(self):
        return (f"Sandwich Panel:\n"
                f"Dimensions: {self.length}m x {self.width}m x {self.total_thickness}m\n"
                f"Area: {self.area:.2f} m²\n"
                f"Volume: {self.volume:.3f} m³\n"
                f"Weight: {self.weight:.1f} N\n"
                f"Flexural Rigidity (D): {self.D:.2e} Nm²")