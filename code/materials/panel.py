import numpy as np
from ..materials.sandwich import Sandwich
from code.materials.load import Load
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
        self.volume = self.area * self.total_thickness * 1e-3 # m³
        
        # Calculate panel properties
        self.D = self.sandwich.D  # Flexural rigidity
        self.S = self.sandwich.S
        self.g = 9.80665  # m/s² (acceleration due to gravity)
        self.core_mass, self.face_mass, self.mass, self.weight = self.calculate_weight()  # N

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

        total_mass = face_mass + core_mass  # kg
        total_weight = (face_mass + core_mass)  * g # N

        print (f"Panel mass: {total_mass} kg")
        
        return core_mass, face_mass, total_mass, total_weight
    
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

    def check_against_face_failure(self, print_values=False, print_only_max=False) -> float:
        """
        Calculate the safety factor for face failure due to normal stress.
        
        Parameters:
        -----------
        load : float
            Uniform load in N/m²
        print_values : bool
            If True, print the values for each ply
         

        Returns:
        --------
        float
            Safety factor for face failure
        """
        # Calculate the maximum normal stress in the faces due to bending
        max_bending_moment = max(self.max_bending_moment_distributed_load, self.max_bending_moment_point_load)

        d = self.sandwich.d  # m
        loading_laminate_bottom = Load(
            Nx=max_bending_moment / (d*1e-3),  # Normal force in x-direction (N/m²)
            Ny=0,  # Normal force in y-direction (N/m²)
            Nxy=0,  # Shear force in xy-plane (N/m²)
            Mx=0,  # Bending moment about x-axis (Nm/m²)
            My=0,  # Bending moment about y-axis (Nm/m²)
            Mxy=0   # Twisting moment (Nm/m²)
        )
        loading_laminate_top = Load(
            Nx=-(max_bending_moment / (d*1e-3)),  # Normal force in x-direction (N/m²)
            Ny=0,  # Normal force in y-direction (N/m²)
            Nxy=0,  # Shear force in xy-plane (N/m²)
            Mx=0,  # Bending moment about x-axis (Nm/m²)
            My=0,  # Bending moment about y-axis (Nm/m²)
            Mxy=0   # Twisting moment (Nm/m²)
        )

        # Calculate the Tsai-Wu, Tsai-Hill, and Maximum Stress criteria for both bottom and top faces

        # tsai-wu
        values_tsai_wu, failure_tsai_wu= self.sandwich.composite_material.tsai_wu_laminate(loading_laminate_bottom)

        # Bottom face (tension) tsai-hill
        values_tsai_hill_bottom, failure_tsai_hill_bottom = self.sandwich.composite_material.tsai_hill_laminate(loading_laminate_bottom)
        failure_max_bottom = self.sandwich.composite_material.maximum_stress_laminate(loading_laminate_bottom)

        # Top face (compression) tsai-hill
        
        values_tsai_hill_top, failure_tsai_hill_top = self.sandwich.composite_material.tsai_hill_laminate(loading_laminate_top)
        failure_max_top = self.sandwich.composite_material.maximum_stress_laminate(loading_laminate_top)

        if failure_tsai_wu:
            print("Failure in the bottom laminate according to Tsai-Wu criterion")
            print(f"Values: {values_tsai_wu}")

        if failure_tsai_hill_bottom:
            print("Failure in the bottom laminate according to Tsai-Hill criterion")
            print(f"Values: {values_tsai_hill_bottom}")

        if failure_max_bottom:
            print("Failure in the bottom laminate according to Maximum Stress criterion")

        if failure_tsai_hill_top:
            print("Failure in the top laminate according to Tsai-Hill criterion")
            print(f"Values: {values_tsai_hill_top}")

        if failure_max_top:
            print("Failure in the top laminate according to Maximum Stress criterion")
        
        if print_values:
            for i, ply in enumerate(self.sandwich.composite_material.plies):
                print(f"Ply {i+1} with angle {ply.angle}°:")
                print(f"  Tsai-Wu: {values_tsai_wu[i]} ")
                print(f"  Tsai-Hill top: {values_tsai_hill_top[i]}")
                print(f"  Tsai-Hill bottom: {values_tsai_hill_bottom[i]}")
                print()

        if print_only_max:
            print(f"Max Tsai-Wu: {max(values_tsai_wu)}")
            print(f"Max Tsai-Hill bottom: {max(values_tsai_hill_bottom)}")
            print(f"Max Tsai-Hill top: {max(values_tsai_hill_top)}")

        return values_tsai_wu, failure_tsai_wu, values_tsai_hill_bottom, failure_tsai_hill_bottom, values_tsai_hill_top, failure_tsai_hill_top, failure_max_bottom, failure_max_top

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

        d = self.sandwich.d/1000  # m (convert from mm to m)  
        tau_c_max = self.sandwich.core_material.sigma_shear  # Maximum allowable core shear stress (Pa) 
        
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
        Ef = self.sandwich.composite_material.E  # Face modulus of elasticity (Pa)
        Ec = self.sandwich.core_material.E  # core_material modulus of elasticity (Pa)
        Gc = self.sandwich.core_material.G  # core_material shear modulus (Pa)
        
        sigma_cr = 0.5**3 * np.sqrt(Ef * Ec * Gc)  # Critical wrinkling stress (Pa)
        print(f"Critical wrinkling stress: {sigma_cr} Pa")
        
        # Maximum normal stress in faces due to bending
        max_bending_moment = max(self.max_bending_moment_distributed_load, self.max_bending_moment_point_load)

        d = self.sandwich.d/1000  # m
        sigma_max = (max_bending_moment / d) 
        
        # Safety factor for face wrinkling
        safety_factor = abs( sigma_cr) / sigma_max 
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
        sigma_c_comp = self.sandwich.core_material.Rm
        min_A = point_load / sigma_c_comp 
        safety_factor = (0.12*0.12)/ min_A


        return safety_factor
    
    def calculate_load_vector_point_load_laminate(self) -> np.ndarray:
        """Calculate the loads in each ply of the sandwich panel.
        
        Returns:
        --------
        np.ndarray
            Array of loads in each ply (N)
        """
        # Initialize an array to hold the loads in each ply
        
        Nx = self.max_bending_moment_point_load / self.sandwich.d
        Ny = 0
        Nxy = 0
        Mx = 0
        My = 0
        Mxy = 0

        return np.array([Nx, Ny, Nxy, Mx, My, Mxy])
    
    def calculate_load_vector_distributed_load_laminate(self) -> np.ndarray:
        """Calculate the loads in each ply of the sandwich panel.
        
        Returns:
        --------
        np.ndarray
            Array of loads in each ply (N)
        """
        # Initialize an array to hold the loads in each ply
        
        Nx = self.max_bending_moment_distributed_load / self.sandwich.d
        Ny = 0
        Nxy = 0
        Mx = 0
        My = 0
        Mxy = 0

        return np.array([Nx, Ny, Nxy, Mx, My, Mxy])
    

    def __str__(self):
        ret = f"Panel:\n"
        ret += f"  Width: {self.width} m\n"
        ret += f"  Length: {self.length} m\n"
        ret += f"  Area: {self.area} m²\n"
        ret += f"  Total Thickness: {self.total_thickness} m\n"
        ret += f"  Volume: {self.volume} m³\n"
        ret += f"  Mass: {self.mass} kg\n"
        ret += f"  Core Mass: {self.core_mass} kg\n"
        ret += f"  Face Mass: {self.face_mass} kg\n"
        ret += f"  Total Weight: {self.weight} N\n"
        ret += f"  Flexural Rigidity (D): {self.D} Nm²\n"
        ret += f"  Shear Rigidity (S): {self.sandwich.S} Nm²\n"
        ret += f"  Maximum Bending Moment (Distributed Load): {self.max_bending_moment_distributed_load} Nm\n"
        ret += f"  Maximum Shear Force (Distributed Load): {self.max_shear_force_distributed_load} N\n"
        ret += f"  Maximum Deflection (Distributed Load): {self.delta_max_distributed_load} m\n"
        ret += f"  Maximum Bending Moment (Point Load): {self.max_bending_moment_point_load} Nm\n"
        ret += f"  Maximum Shear Force (Point Load): {self.max_shear_force_point_load} N\n"
        ret += f"  Maximum Deflection (Point Load): {self.delta_max_point_load} m\n"
        return ret
