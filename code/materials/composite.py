import numpy as np
from code.materials.isotropic import Matrix, Reinforcement
from code.materials.load import Load

class Lamina:
    def __init__(self, E1, E2, G12, G23, nu12,  sigma_1t, sigma_1c, sigma_2t, sigma_2c, sigma_shear, nu13=None, nu23=None, thickness=None, G13=None,nu21= None, rho=None, angle=0):
        """
        Initializes the Lamina class with the given properties.

        parameters:
        ------------
        E1 : float
            Young's modulus in the fibre direction (Pa)
        E2 : float
            Young's modulus in the transverse direction (Pa)
        G12 : float
            In-plane shear modulus (Pa)
        G23 : float
            Out-of-plane shear modulus (Pa)
        nu12 : float
            Poisson's ratio for the 1-2 plane
        nu23 : float, optional
            Poisson's ratio for the 2-3 plane (default is None, set to nu12 if not provided)
        nu13 : float, optional
            Poisson's ratio for the 1-3 plane (default is None, set to nu12 if not provided)
        G13 : float, optional
            Out-of-plane shear modulus (Pa) (default is None, set to G12 if not provided)
        nu21 : float, optional
            Poisson's ratio for the 2-1 plane (default is None, set to nu12 if not provided)
        rho : float, optional
            Density (kg/m^3) (default is None)
        sigma_1t : float
            Longitudinal tensile strength (Pa)
        sigma_1c : float
            Longitudinal compressive strength (Pa)
        sigma_2t : float
            Transverse tensile strength (Pa)
        sigma_2c : float
            Transverse compressive strength (Pa)
        sigma_shear : float
            Shear strength (Pa)
        angle : float, optional
            Fiber orientation in degrees (default is 0)
        thickness : float, optional (mm)
            Thickness of the ply in mm (default is None)
        """
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        if G13 is None:
            self.G13 = G12
        else:
            self.G13 = G13
        self.G23 = G23
        self.nu12 = nu12
        if nu23 is None:
            self.nu23 = nu12 
        else:
            self.nu23 = nu23
        if nu13 is None:
            self.nu13 = nu12 
        else:
            self.nu13 = nu13
        if nu21 is None:
            self.nu21 = (E2/E1)*nu12
        else:
            self.nu21 = nu21
        self.rho = rho
        self.sigma_1t = sigma_1t # Longitudinal tensile strength in Pa
        self.sigma_1c = sigma_1c # Longitudinal compressive strength in Pa
        self.sigma_2t = sigma_2t # Transverse tensile strength in Pa
        self.sigma_2c = sigma_2c # Transverse compressive strength in Pa
        self.sigma_shear = sigma_shear # Shear strength in Pa
        self.angle = angle  # Fiber orientation in degrees
        self.thickness = thickness  # Thickness of the ply in mm
        self.Q12 = self.get_Q12_matrix()  # Stiffness matrix in local coordinate system
        self.Qxy = self.transform_Qxy_matrix()  # Stiffness matrix in global coordinate system
        self.T = self.get_T_matrix()  # Transformation matrix for strain transformation
        
    def get_Q12_matrix(self):
        """
        Build the Q and S matrices that will serve to calculate the stress and strain in the laminate
        in the local coordinate system

        Parameters:
        ------------
        E1 : float
            Young's modulus in the fibre direction
        E2 : float
            Young's modulus in the transverse direction
        nu12 : float
            Poisson's ratio for the 1-2 plane
        nu21 : float
            Poisson's ratio for the 2-1 plane
        G12 : float
            In-plane shear modulus

        Returns:
        ---------
        Q12 : np.array (size 3x3)
            The Q matrix in local coordinate system
        """
        Q12 = np.array([
            [self.E1 / (1 - self.nu12 * self.nu21), self.E2 * self.nu12 / (1 - self.nu12 * self.nu21), 0],
            [self.E1 * self.nu21 / (1 - self.nu12 * self.nu21), self.E2 / (1 - self.nu12 * self.nu21), 0],
            [0, 0, self.G12]
        ])

        return Q12
    
    def transform_Qxy_matrix(self):
        """
        Transform the Q matrix to the global coordinate system of the ply

        Parameters:
        -----------
        self: Lamina
            The Lamina object

        Returns:
        --------
        Q_xy : np.array (size 3x3)
            The Q matrix in the global coordinate system
        """
        Q12 = self.get_Q12_matrix()
        ply_angle_x1 = self.angle  # Angle in degrees
        c = np.cos(np.radians(ply_angle_x1))
        s = np.sin(np.radians(ply_angle_x1))
        T = np.array([
            [c**2, s**2, -2*c*s],
            [s**2, c**2, 2*c*s],
            [c*s, -c*s, c**2 - s**2]
        ])
        Qxy = T@Q12@T.T
        return Qxy
    
    def get_T_matrix(self):
        """
        Build the transformation matrix T that will serve to calculate transform the global strain to local strain

        Parameters:
        ------------
        angle : float
            The angle of the ply in degrees

        Returns:
        ---------
        T : np.array (size 3x3)
            The transformation matrix
        """
        c = np.cos(np.radians(self.angle))
        s = np.sin(np.radians(self.angle))
        T = np.array([
            [c**2, s**2, c*s],
            [s**2, c**2, -c*s],
            [-2*c*s, 2*c*s, c**2 - s**2]
        ])
        return T
    
    def calculate_local_stress(self, global_strain):
        """
        Calculate the local stress in the ply using the global strain and the transformation matrix.

        Parameters:
        -----------
        global_strain : np.array (size 3)
            The strain vector for the ply in the global coordinate system [epsilon_x, epsilon_y, gamma_xy]

        Returns:
        --------
        local_stress : np.array (size 3)
            Stress vector in loca coordinate system [sigma_x, sigma_y, tau_xy]
        """
        local_strain = np.dot(self.T, global_strain)
        local_stress = np.dot(self.Q12, local_strain)
        return local_stress
    
    def calculate_global_stress(self, global_strain):
        """
        Calculate the global stress in the ply using the local strain and the transformation matrix.

        Parameters:
        -----------
        global_strain : np.array (size 3)
            The strain vector for the ply in the global coordinate system [epsilon_x, epsilon_y, gamma_xy]

        Returns:
        --------
        global_stress : np.array (size 3)
            The global stress vector [sigma_x, sigma_y, tau_xy]
        """
        global_stress = self.calculate_local_stress(global_strain)
        global_stress = np.dot(self.T.T, global_stress)
        return global_stress
    
    def failure_Tsai_Wu(self, global_strain):
        """
        Calculate the Tsai-Wu failure criterion for all plies in the laminate.

        Parameters:
        -----------
        global_strain : np.array (size 3)
            The strain vector for the ply in the global coordinate system [epsilon_x, epsilon_y, gamma_xy]

        Returns:
        --------
        values : np.array (size n_nodes x n_plies)
            The value of the Tsai-Wu failure criterion for each ply in the laminate at each node
        """
        
        stresses = self.calculate_local_stress(global_strain)
        # Gather strengths from this ply
        F1t = self.sigma_1t
        F1c = self.sigma_1c
        F2t = self.sigma_2t
        F2c = self.sigma_2c
        F6 = self.sigma_shear
        f1 = 1 / F1t - 1 / F1c
        f2 = 1 / F2t - 1 / F2c
        f11 = 1 / (F1t * F1c)
        f22 = 1 / (F2t * F2c)
        f66 = 1 / (F6 ** 2)
        f12 = -0.5 * np.sqrt(f11 * f22)
        sigma1 = stresses[0]
        sigma2 = stresses[1]
        sigma6 = stresses[2]
        value = abs(f1 * sigma1 + f2 * sigma2 + f11 * sigma1 ** 2 + f22 * sigma2 ** 2 + f66 * sigma6 ** 2 + 2 * f12 * sigma1 * sigma2)
        return value
        
    def failure_Tsai_Hill(self, global_strain):
        """
        Calculate the Tsai-Hill failure criterion for this ply using the global strain.

        Parameters:
        -----------
        global_strain : np.array (size 3)
            The strain vector for the ply in the global coordinate system [epsilon_x, epsilon_y, gamma_xy]

        Returns:
        --------
        value : float
            The value of the Tsai-Hill failure criterion
        """
        # Calculate local stress from global strain
        stress = self.calculate_local_stress(global_strain)
        # Use ply strengths
        F1t = self.sigma_1t
        F1c = self.sigma_1c
        F2t = self.sigma_2t
        F2c = self.sigma_2c
        F6 = self.sigma_shear

        F1 = F1t if stress[0] >= 0 else F1c
        F2 = F2t if stress[1] >= 0 else F2c

        value = (stress[0] ** 2) / (F1 ** 2) + (stress[1] ** 2) / (F2 ** 2) + (stress[2] ** 2) / (F6 ** 2) - (stress[0] * stress[1]) / (F1 ** 2)
        return value
    
    
    
    def __str__(self):
        result = f"Lamina Properties:\n"
        result += f"Angle: {self.angle} degrees\n"
        result += f"Thickness: {self.thickness} mm\n"
        result += f"Density: {self.rho} kg/m^3\n"
        result += f"local Stiffness matrix Q12:\n{self.Q12}\n"
        result += f"global Stiffness Matrix Qxy:\n{self.Qxy}\n"
        return result


class Laminate:
    def __init__(self, plies, total_thickness=None, ply_thicknesses=None, density=None):
        """
        Initializes the Laminate class with the given plies.

        Parameters:
        ----------
        plies: array of Lamina objects
            The plies that make up the laminate.
        total_thickness: float, optional (mm)
            The total total_thickness of the laminate in mm (default is None).
        ply_thicknesses: array of floats, optional (mm)
            The thicknesses of each ply in mm (default is None).
        density: float, optional (kg/m^3)
            The density of the laminate in kg/m^3 (default is None).
        
        
        """
        self.plies = plies
        self.n_plies = len(plies)       
        self._initialize_thicknesses(total_thickness, ply_thicknesses)
        A, B , D = self.calculate_A_B_D_matrices()
        self.A = A
        self.B = B
        self.D = D
        self.ABD = np.block([[A, B], [B, D]])
        self.inv_ABD = np.linalg.inv(self.ABD)
        self.E = self.get_equivalent_elasticity()
        if density is None:
            self.rho = self.get_equivalent_density()
        else:
            self.rho = density
        
    def _initialize_thicknesses(self, total_thickness, ply_thicknesses):
            """
            Initializes the thicknesses of the plies and the total thickness of the laminate.

            Parameters:
            ----------
            total_thickness: float, optional (mm)
                The total thickness of the laminate in mm.
            ply_thicknesses: array of floats, optional (mm)
                The thicknesses of each ply in mm.

            Raises:
            -------
            ValueError: If thickness information is insufficient or inconsistent.
            """
            # Check if individual ply thicknesses are defined in Lamina objects
            if all(ply.thickness is not None for ply in self.plies):
                # Use the thicknesses defined in each Lamina object
                self.ply_thicknesses = np.array([ply.thickness for ply in self.plies])
                self.total_thickness = sum(self.ply_thicknesses)
                # Optional: Check consistency if total_thickness was also provided
                if total_thickness is not None and not np.isclose(self.total_thickness, total_thickness):
                     print(f"Warning: Provided total_thickness ({total_thickness}) does not match the sum of individual ply thicknesses ({self.total_thickness}). Using the sum. Thickness used: {self.total_thickness}")
                if ply_thicknesses is not None and not np.array_equal(self.ply_thicknesses, ply_thicknesses):
                     print(f"Warning: Provided ply_thicknesses do not match the thicknesses defined in Lamina objects. Using Lamina object thicknesses.Thicknesses used: {self.ply_thicknesses}")

            elif total_thickness is not None:
                # Use the provided total_thickness
                self.total_thickness = total_thickness
                if ply_thicknesses is None:
                    # If not given, calculate the thickness of each ply equally
                    ply_thickness = total_thickness / self.n_plies
                    self.ply_thicknesses = np.full(self.n_plies, ply_thickness)
                    for i, ply in enumerate(self.plies):
                        ply.thickness = self.ply_thicknesses[i]
                else:
                    # If given, assign the thickness to each ply
                    if len(ply_thicknesses) != self.n_plies:
                        raise ValueError("The length of ply_thicknesses must be equal to the number of plies.")
                    if not np.isclose(sum(ply_thicknesses), total_thickness):
                         raise ValueError("The sum of ply_thicknesses must be equal to the total_thickness.")
                    self.ply_thicknesses = np.array(ply_thicknesses)
                    for i, ply in enumerate(self.plies):
                        ply.thickness = self.ply_thicknesses[i]
            else:
                raise ValueError("Either individual ply thicknesses must be defined in Lamina objects, or total_thickness must be provided.")

    def calculate_height_list(self):
        total_thickness = 0
        ply_thicknesses = np.array([])
        for ply in self.plies:
            total_thickness += ply.thickness
            np.append(ply_thicknesses, ply.thickness)
        
        mid_plane = total_thickness / 2
        height_list = np.zeros(self.n_plies + 1)
        height_list[0] = mid_plane
        for i in range(1, self.n_plies + 1):
            height_list[i] = height_list[i - 1] - (self.plies[i - 1].thickness)

        return height_list
    
    def get_equivalent_elasticity(self):
        """
        Calculate the equivalent elasticity modulus of a laminate at each node.

        Parameters:
        -----------
        A : np.array (size x3x3)
            The ABD stiffness matrix for each node in the laminate.
        self.tot_height : np.array (size N)
            The laminate total_thickness at each node.

        Returns:
        --------
        E_laminate : np.array (size N)
            The equivalent elasticity modulus of the laminate at each node.
        """
        A = self.A
        E_laminate = (1/(self.total_thickness * 1e-3)) * (( 2* A[0, 2] * A[0, 1] * A[1, 2] - A[0, 2]**2 * A[1, 1] - A[0, 1]**2 * A[2, 2] + A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2]**2)) / (A[1, 1] * A[2, 2] - A[1, 2]**2))

        return E_laminate
    
    def get_equivalent_density(self):
        """
        Calculate the equivalent density of a laminate at each node.

        Parameters:
        -----------
        plies : list of Lamina objects
            The plies that make up the laminate.

        self.tot_height : np.array (size N)
            The laminate total_thickness at each node.

        Returns:
        --------
        rho_laminate : float
        """
        rho_laminate = 0
        for ply in self.plies:
            ply_density = ply.rho
            ply_thickness = ply.thickness
            rho_laminate += ply_density * ply_thickness
        rho_laminate /= self.total_thickness
        return rho_laminate

    def calculate_A_B_D_matrices(self):
        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))

        # Calculate height list (assuming it's in mm)
        height_list_mm = self.calculate_height_list()
        # Convert heights to meters for this specific calculation
        height_list_m = height_list_mm /1000  # Convert mm to m

        z_index = 1
        for ply in self.plies:
            # Assuming ply.Qxy is in Pascals (N/m^2)
            Qxy = ply.Qxy 
            
            # Get heights in meters for the current ply
            h_k_m = height_list_m[z_index - 1]
            h_km1_m = height_list_m[z_index]

            # Calculate contributions using meters
            A += (h_k_m - h_km1_m) * Qxy  # [A] = (m) * (N/m^2) = N/m
            B += ((h_k_m**2 - h_km1_m**2)/2) * Qxy  # [B] = (m^2) * (N/m^2) = N
            D += ((h_k_m**3 - h_km1_m**3)/3) * Qxy  # [D] = (m^3) * (N/m^2) = N*m
            z_index += 1
         # At the end of the calculate_A_B_D_matrices function, add:
            
            tolerance = 1e-9
            A[np.abs(A) < tolerance] = 0.0
            B[np.abs(B) < tolerance] = 0.0
            D[np.abs(D) < tolerance] = 0.0
        return A, B, D

    def calculate_strain_mid(self, loads: Load):
        """
        Calculate the strain at the mid-plane of the laminate using the inverse ABD matrix and loads.

        parameters:
        -----------
        loads : Load object
            The applied loads and moments.

        returns:
        --------
        strain_middle_plane : np.array (size 6)
            The strain at the mid-plane of the laminate with first 3 components being strain and last 3 components being curvature.
        """

        strain_middle_plane = np.dot(self.inv_ABD, loads.array)
        return strain_middle_plane
    
    def calculate_strain(self, loads: Load):
        """
        Calculate the strain in each ply of the laminate based on the applied loads.
        
        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.

        Returns:
        --------
        strain : np.array (size n_plies x 3)
            The strain in each ply of the laminate.
        """
        strain_middle_plane = self.calculate_strain_mid(loads)
        strain = strain_middle_plane[0:3]
        curvature = strain_middle_plane[3:6]
        height_list = self.calculate_height_list()
        strain = np.zeros((self.n_plies, 3))  # Initialize strain array for each ply
        for i in range(1, self.n_plies+1):
            z_k = height_list[i - 1]
            z_km1 = height_list[i]
            strain[i - 1] = strain_middle_plane[0:3] + (curvature * (z_k - (z_km1 + z_k) / 2))
        return strain

    def tsai_wu_laminate(self, loads: Load):
        """
        Calculate the Tsai-Wu failure criterion for all plies in the laminate.

        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.
        
        returns:
        --------
        values : np.array (size n_plies x 3)
            The value of the Tsai-Wu failure criterion for each ply in the laminate at each node
        failure : bool
            True if any ply fails, False otherwise.
        """
        #initialize the failure variable to False
        failure = False
        values = np.zeros(self.n_plies)  # Initialize values array for each ply
        # Calculate the strain of the laminate at each ply
        strain_global = self.calculate_strain(loads)
        for i, ply in enumerate(self.plies):
            # Calculate Tsai-Wu failure criterion for each ply
            value = ply.failure_Tsai_Wu(strain_global[i])
            values[i] = value
            # Check if the Tsai-Wu criterion value exceeds 1.0 (failure condition)
            if value >= 1.0:
                print(f"Failure in Lamina {i+1}: Tsai-Wu criterion value = {value:.4e}")
                failure = True
        return values, failure
    
    def tsai_hill_laminate(self, loads: Load):
        """
        Calculate the Tsai-Hill failure criterion for all plies in the laminate.

        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.
        
        returns:
        --------
        values : np.array (size n_plies x 3)
            The value of the Tsai-Hill failure criterion for each ply in the laminate at each node
        failure : bool
            True if any ply fails, False otherwise.
        """
        #initialize the failure variable to False
        failure = False
        values = np.zeros(self.n_plies)

        # Calculate the strain of the laminate at each ply
        strain_global = self.calculate_strain(loads)
        for i, ply in enumerate(self.plies):
            # Calculate Tsai-Hill failure criterion for each ply
            value = ply.failure_Tsai_Hill(strain_global[i])
            values[i] = value
            # Check if the Tsai-Hill criterion value exceeds 1.0 (failure condition)
            if value >= 1.0:
                print(f"Failure in Lamina {i+1}: Tsai-Hill criterion value = {value:.4e}")
                failure = True
        return values, failure
    
    def maximum_stress_laminate(self, loads: Load):
        """
        Calculate the maximum stress failure criterion for all plies in the laminate.

        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.
        
        returns:
        --------
        values : np.array (size n_plies x 3)
            The value of the maximum stress failure criterion for each ply in the laminate at each node
        failure : bool
            True if any ply fails, False otherwise.
        """
        #initialize the failure variable to False
        failure = False
        # Calculate the strain of the laminate at each ply
        strain_global = self.calculate_strain(loads)

        for i, ply in enumerate(self.plies):
            # Calculate local stress from global strain
            local_stress = ply.calculate_local_stress(strain_global[i])

            #Compare local stresses with ply strengths
            if local_stress[0] > 0:
                if local_stress[0] > ply.sigma_1t:
                    print(f"Failure in Lamina {i+1}: Maximum stress criterion value = {local_stress[0]:.4e} compared to {ply.sigma_1t:.4e} the tensile strength in the fiber direction")
                    failure = True
            else:
                if abs(local_stress[0]) > ply.sigma_1c:
                    print(f"Failure in Lamina {i+1}: Maximum stress criterion value = {local_stress[0]:.4e} compared to {ply.sigma_1c:.4e} the compressive strength in the fiber direction")
                    failure = True
            if local_stress[1] > 0:
                if local_stress[1] > ply.sigma_2t:
                    print(f"Failure in Lamina {i+1}: Maximum stress criterion value = {local_stress[1]:.4e} compared to {ply.sigma_2t:.4e} the tensile strength in the transverse direction")
                    failure = True
            else:
                if abs(local_stress[1]) > ply.sigma_2c:
                    print(f"Failure in Lamina {i+1}: Maximum stress criterion value = {local_stress[1]:.4e} compared to {ply.sigma_2c:.4e} the compressive strength in the transverse direction")
                    failure = True
            if abs(local_stress[2]) > ply.sigma_shear:
                print(f"Failure in Lamina {i+1}: Maximum stress criterion value = {local_stress[2]:.4e} compared to {ply.sigma_shear:.4e} the shear strength")
                failure = True
            
        return failure
    
    def print_local_stresses(self, loads: Load):
        """
        Print the stresses in the each lamina.

        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.
        """
        strain = self.calculate_strain(loads)
        for i, ply in enumerate(self.plies):
            local_stress = ply.calculate_local_stress(strain[i])
            print(f"Lamina {i+1} local stress:\n{local_stress}\n")
    
    def print_global_stresses(self, loads: Load):
        """
        Print the stresses and moments in each lamina.

        Parameters:
        -----------
        loads : Load object
            The applied loads and moments.
        """
        strain = self.calculate_strain(loads)
        for i, ply in enumerate(self.plies):
            global_stress = ply.calculate_global_stress(strain[i])
            print(f"Lamina {i+1} global stress:\n{global_stress}\n")
    


    def __str__(self):
        np.set_printoptions(precision=4, suppress=False, formatter={'float_kind':lambda x: f"{x:.4e}"})
        result = "Laminate Characteristics:\n"
        result += f"Total Thickness: {self.total_thickness} mm\n"
        result += f"Number of Plies: {self.n_plies}\n"
        result += f"Density: {self.rho} kg/m^3\n"
        result += f"Equivalent Elasticity Modulus: {self.E:.4e} Pa\n"
        result += f"Plies:\n"
        result += f"thicknesses: {np.array2string(self.ply_thicknesses, separator=', ')} mm\n"
        result += f"plies list thicknesses: {np.array2string(self.calculate_height_list(), separator=', ')} mm\n"
        result += f"\nMatrix A:\n{self.A}\n"
        result += f"\nMatrix B:\n{self.B}\n"
        result += f"\nMatrix D:\n{self.D}\n"
        result += f"\nABD Matrix:\n{self.ABD}\n"
        result += f"\nInverse ABD Matrix:\n{self.inv_ABD}\n"
        return result
    
    

# Example usage of the library
if __name__ == "__main__":
    # Define matrix and reinforcement properties
    matrix = Matrix(E=4.3e9, nu=0.4, rho=1.2, G=1640e6, Rm=135e6)
    reinforcement = Reinforcement(E_l=390e9, E_t=12e9, G=20000e6, nu_lt=0.35, nu_tl=0.01, rho=1.8, Rm=2540e6)

    # Define a single ply
    ply = Lamina(matrix, reinforcement, fiber_volume_ratio=0.4, fiber_areal_weight=200, orientation=45,
              sigma_lt=2978e6, sigma_tt=36e6, sigma_cl=2100e6, sigma_ct=200e6, tau_lt=70e6)

    # Define applied loads and moments
    loads = Load(Nx=80, My=10)

    # Define the laminate with a single ply
    laminate = Laminate([ply])

    # Print stiffness matrix for verification
    print("Stiffness Matrix Q and Transformed Q:")
    S, Sxy, Q, Qxy = ply.stiffness_matrix()
    print("Q:\n", Q)
    print("Qxy:\n", Qxy)

    # Calculate and print loads
    laminate.print_stresses(loads)

    # Calculate and print Tsai-Wu failure criterion
    laminate.tsai_wu(loads)



#being worked on the LaminaMix class need to have the strength calculated form values of the matrix and the reinforcement
class LaminaMix:
    def __init__(self, matrix: Matrix, reinforcement: Reinforcement, fiber_volume_ratio, fiber_areal_weight, orientation=0,
                 sigma_lt=0, sigma_tt=0, sigma_cl=0, sigma_ct=0, tau_lt=0):
        self.matrix = matrix
        self.reinforcement = reinforcement
        self.fiber_volume_ratio = fiber_volume_ratio
        self.orientation = orientation  # Fiber orientation in degrees
        self.fiber_areal_weight = fiber_areal_weight  # Fiber areal weight in g/m2
        self.sigma_lt = sigma_lt  # Longitudinal tensile strength in Pa
        self.sigma_tt = sigma_tt  # Transverse tensile strength in Pa
        self.sigma_ct = sigma_ct  # Transverse compressive strength in Pa
        self.sigma_cl = sigma_cl  # Longitudinal compressive strength in Pa
        self.tau_lt = tau_lt  # Longitudinal-transverse shear strength in Pa
        self.thickness = self.get_ply_thickness()  # Calculate ply thickness
        self.Rm = self.calculate_Rm()  # Calculate tensile strength
        S12, Sxy, Q12, Qxy = self.stiffness_matrix()  # Calculate stiffness matrix
        self.S12 = S12  # Compliance matrix
        self.Sxy = Sxy  # Transformed compliance matrix
        self.Q12 = Q12  # Stiffness matrix
        self.Qxy = Qxy  # Transformed stiffness matrix


    def get_ply_thickness(self):
        # Calculate ply thickness from fiber areal weight
        matrix_volume_ratio = 1 - self.fiber_volume_ratio
        thickness = (self.fiber_areal_weight * 1e-3) / (matrix_volume_ratio * self.reinforcement.rho)
        return thickness

    def stiffness_matrix(self):
        # Calculate the stiffness matrix of the ply
        matrix_volume_ratio = 1 - self.fiber_volume_ratio
        E_l = (self.fiber_volume_ratio * self.matrix.E) + (matrix_volume_ratio * self.reinforcement.E_l)
        E_t = self.matrix.E / (self.fiber_volume_ratio + ((self.matrix.E / self.reinforcement.E_t) * matrix_volume_ratio))
        nu_lt = self.reinforcement.nu_lt * matrix_volume_ratio + self.matrix.nu * self.fiber_volume_ratio
        nu_tl = nu_lt * (E_t / E_l)
        G_lt = self.matrix.G / (self.fiber_volume_ratio + ((self.matrix.G / self.reinforcement.G) * matrix_volume_ratio))

        theta = np.radians(self.orientation)
        c = np.cos(theta)
        s = np.sin(theta)

        S = np.array([[1 / (E_l * 1e-6), -nu_tl / (E_t * 1e-6), 0],
                       [-nu_lt / (E_l * 1e-6), 1 / (E_t * 1e-6), 0],
                       [0, 0, 1 / (G_lt * 1e-6)]])
        Q = np.linalg.inv(S)

        T = np.array([[c**2, s**2, 2 * c * s],
                      [s**2, c**2, -2 * c * s],
                      [-c * s, c * s, c**2 - s**2]])
        T1 = np.array([[c**2, s**2, -2 * c * s],
                       [s**2, c**2, 2 * c * s],
                       [c * s, -c * s, c**2 - s**2]])

        TS = np.dot(T.T, S)
        T1Q = np.dot(T1, Q)
        Sxy = np.dot(TS, T)
        Qxy = np.dot(T1Q, T1.T)

        return S, Sxy, Q, Qxy

    def calculate_Rm(self):
        Rm = self.reinforcement.Rm * ((1 - self.fiber_volume_ratio) + ((self.matrix.E / self.reinforcement.E_l) * self.fiber_volume_ratio))
        return Rm

    def __str__(self):
        thickness = self.thickness
        S, Sxy, Q, Qxy = self.stiffness_matrix()
        Rm = self.calculate_Rm()
        matrix_volume_ratio = 1 - self.fiber_volume_ratio
        E_l = (self.fiber_volume_ratio * self.matrix.E) + (matrix_volume_ratio * self.reinforcement.E_l)
        E_t = self.matrix.E / (self.fiber_volume_ratio + ((self.matrix.E / self.reinforcement.E_t) * matrix_volume_ratio))
        nu_lt = self.reinforcement.nu_lt * matrix_volume_ratio + self.matrix.nu * self.fiber_volume_ratio
        nu_tl = nu_lt * (E_t / E_l)
        G_lt = self.matrix.G / (self.fiber_volume_ratio + ((self.matrix.G / self.reinforcement.G) * matrix_volume_ratio))
        rho = matrix_volume_ratio * self.reinforcement.rho + self.fiber_volume_ratio * self.matrix.rho

        result = "Lamina Characteristics:\n"
        result += f"Rm: {Rm * 1e-6:.3f} MPa\n"
        result += f"Density: {rho:.3f} g/cm3\n"
        result += f"Elastic Moduli: E_l = {E_l * 1e-9:.3f} GPa, E_t = {E_t * 1e-9:.3f} GPa\n"
        result += f"Shear Modulus: G_lt = {G_lt * 1e-9:.3f} GPa\n"
        result += f"Poisson's Ratios: nu_tl = {nu_tl:.3f}, nu_lt = {nu_lt:.3f}\n"
        result += f"Thickness: {thickness:.3f} mm\n"
        result += f"\nStiffness Matrix Q:\n{Q}"
        result += f"\nTransformed Stiffness Matrix Qxy:\n{Qxy}"
        result += f"\nCompliance Matrix S:\n{S}"
        result += f"\nTransformed Compliance Matrix Sxy:\n{Sxy}"

        return result