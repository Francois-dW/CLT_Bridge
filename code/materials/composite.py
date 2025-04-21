import numpy as np
from isotropic import Matrix, Reinforcement

class Lamina:
    def __init__(self, E1, E2, G12, G23, nu12, nu13, nu23, sigma_1l, sigma_1c, sigma_2t, sigma_2c, sigma_12, thickness= None, G13=None,nu21= None, rho=None, angle=0):
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        if G13 is None:
            self.G13 = G12
        else:
            self.G13 = G13
        self.G23 = G23
        self.nu12 = nu12
        if nu21 is None:
            self.nu21 = nu12 * E2 / E1
        else:
            self.nu21 = nu21
        self.nu13 = nu13
        self.nu23 = nu23
        self.rho = rho
        self.sigma_1l = sigma_1l
        self.sigma_1c = sigma_1c
        self.sigma_2t = sigma_2t
        self.sigma_2c = sigma_2c
        self.sigma_12 = sigma_12
        self.angle = angle  # Fiber orientation in degrees
        self.thickness = thickness  # Thickness of the ply in mm
    

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
            [self.E1 / (1 - self.nu12 * self.nu21), self.E1 * self.nu21 / (1 - self.nu12 * self.nu21), 0],
            [self.E2 * self.nu12 / (1 - self.nu12 * self.nu21), self.E2 / (1 - self.nu12 * self.nu21), 0],
            [0, 0, self.G12]
        ])

        return Q12
    

    def get_S12_matrix(self):
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
        S12 : np.array (size 3x3)
            The Q matrix in local coordinate system
        """
        S12 = np.array([
            [1 / self.E1, -self.nu12 / self.E1, 0],
            [-self.nu21 / self.E2, 1 / self.E2, 0],
            [0, 0, 1 / self.G12]
        ])

        return S12
        
    
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
            [-c*s, c*s, c**2 - s**2]
        ])
        Qxy = np.dot(np.dot(np.linalg.inv(T), Q12), T)
        return Qxy
    
    def transform_Sxy_matrix(self):
        """
        Transform the S matrix to the global coordinate system of the ply

        Parameters:
        -----------
        self: Lamina
            The Lamina object

        Returns:
        --------
        S_xy : np.array (size 3x3)
            The S matrix in the global coordinate system
        """
        S12 = self.get_S12_matrix()
        ply_angle_x1 = self.angle
        c = np.cos(np.radians(ply_angle_x1))
        s = np.sin(np.radians(ply_angle_x1))
        T = np.array([
            [c**2, s**2, 2*c*s],
            [s**2, c**2, -2*c*s],
            [-c*s, c*s, c**2 - s**2]
        ])
        Sxy = np.dot(np.dot(np.linalg.inv(T), S12), T)
        return Sxy


class Laminate:
    def __init__(self, plies, total_thickness=None, ply_thicknesses=None):
        """
        Initializes the Laminate class with the given plies.

        Parameters:
        ----------
        plies: array of Lamina objects
            The plies that make up the laminate.
        total_thickness: float, optional
            The total thickness of the laminate (default is None).
        ply_thicknesses: array of floats, optional
            The thicknesses of each ply (default is None).
        
        """
        self.plies = plies
        self.n_plies = len(plies)

        # Check if individual ply thicknesses are defined in Lamina objects
        if all(ply.thickness is not None for ply in plies):
            # Use the thicknesses defined in each Lamina object
            self.thickness = sum(ply.thickness for ply in plies)
            self.thicknesses = np.array([ply.thickness for ply in self.plies])
        elif total_thickness is not None:
            # Use the provided total_thickness and divide it equally among the plies
            self.thickness = total_thickness
            if ply_thicknesses is None:
                #If not given, calculate the thickness of each ply
                ply_thickness = total_thickness / self.n_plies
                for ply in self.plies:
                    ply.thickness = ply_thickness
            else:
                #If given, assign the thickness to each ply
                if len(ply_thicknesses) != self.n_plies:
                    raise ValueError("The length of ply_thicknesses must be equal to the number of plies.")
                for i, ply in enumerate(self.plies):
                    ply.thickness = ply_thicknesses[i]
        else:
            raise ValueError("Either individual ply thicknesses must be defined in Lamina objects or total_thickness must be provided.")

    def calculate_height_list(self):
        total_thickness = 0
        ply_thicknesses = []
        for ply in self.plies:
            thickness = ply.thickness()
            ply_thicknesses.append(thickness)
            total_thickness += thickness

        half_thickness = total_thickness / 2

        height_list = [-half_thickness]  # Add the height of the first layer
        z_prev = -half_thickness

        for thickness in ply_thicknesses:
            z = z_prev + thickness
            height_list.append(z)  # Add the height of the current layer
            z_prev = z

        height_list.sort()

        return height_list

    def calculate_A_B_D_matrices(self):
        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))

        height_list = self.calculate_height_list()

        z_index = 1
        for ply in self.plies:
            Q_transformed = ply.stiffness_matrix()[3]
            A += (height_list[z_index] - height_list[z_index - 1]) * Q_transformed  # [A] = N/mm
            B += 0.5 * ((height_list[z_index]**2) - (height_list[z_index - 1]**2)) * Q_transformed  # [B] = N
            D += (1/3) * ((height_list[z_index]**3) - (height_list[z_index - 1]**3)) * Q_transformed  # [D] = Nmm
            z_index += 1

        return A, B, D

    def calculate_ABD_matrix(self):
        A, B, D = self.calculate_A_B_D_matrices()
        ABD = np.vstack((np.hstack((A, B)), np.hstack((B, D))))
        inv_ABD = np.linalg.inv(ABD)
        return ABD, inv_ABD

    def calculate_strain(self, stresses):
        ABD, inv_ABD = self.calculate_ABD_matrix()
        strain = np.dot(inv_ABD, stresses.to_array())
        return strain

    def calculate_stresses(self, stresses):
        strain = self.calculate_strain(stresses)
        height_list = self.calculate_height_list()

        stresses_mid = np.zeros((len(self.plies), 3))  # Stress at the mid-plane of each ply
        stresses_top = np.zeros((len(self.plies), 3))  # Stress at the top of each ply
        stresses_bot = np.zeros((len(self.plies), 3))  # Stress at the bottom of each ply

        z_index = 1
        for ply in self.plies:
            theta = np.radians(ply.orientation)
            c = np.cos(theta)
            s = np.sin(theta)
            T = np.array([[c**2, s**2, 2 * c * s],
                          [s**2, c**2, -2 * c * s],
                          [-c * s, c * s, c**2 - s**2]])
            Q_transformed = ply.stiffness_matrix()[3]
            strain_top = strain[0:3] + strain[3:6] * height_list[z_index - 1]
            strain_mid = strain[0:3] + strain[3:6] * ((height_list[z_index - 1] + height_list[z_index]) / 2)
            strain_bot = strain[0:3] + strain[3:6] * height_list[z_index]

            stress_top_max = np.dot(Q_transformed, strain_top)
            stress_mid_max = np.dot(Q_transformed, strain_mid)
            stress_bot_max = np.dot(Q_transformed, strain_bot)

            stress_top = np.dot(T, stress_top_max)
            stress_mid = np.dot(T, stress_mid_max)
            stress_bot = np.dot(T, stress_bot_max)

            stresses_mid[z_index - 1] = stress_mid
            stresses_top[z_index - 1] = stress_top
            stresses_bot[z_index - 1] = stress_bot

            z_index += 1

        return stresses_mid, stresses_top, stresses_bot

    def print_stresses(self, stresses):
        stresses_mid, stresses_top, stresses_bot = self.calculate_stresses(stresses)
        result = "Stress:\n"
        for i, ply in enumerate(self.plies):
            result += f"Lamina {i+1}:\n"
            result += f"  Stress at mid-plane: {stresses_mid[i]} MPa\n"
            result += f"  Stress at top: {stresses_top[i]} MPa\n"
            result += f"  Stress at bottom: {stresses_bot[i]} MPa\n"
        print(result)

    def calculate_stresses_principal(self, stresses):
        strain = self.calculate_strain(stresses)
        height_list = self.calculate_height_list()

        stresses_mid = np.zeros((len(self.plies), 3))  # Stress at the mid-plane of each ply
        stresses_top = np.zeros((len(self.plies), 3))  # Stress at the top of each ply
        stresses_bot = np.zeros((len(self.plies), 3))  # Stress at the bottom of each ply

        z_index = 1
        for ply in self.plies:
            Q_transformed = ply.stiffness_matrix()[3]

            strain_top = strain[0:3] + (strain[3:6] * height_list[z_index - 1])
            strain_mid = strain[0:3] + (strain[3:6] * ((height_list[z_index - 1] + height_list[z_index]) / 2))
            strain_bot = strain[0:3] + (strain[3:6] * height_list[z_index])

            stress_top_max = np.dot(Q_transformed, strain_top)
            stress_mid_max = np.dot(Q_transformed, strain_mid)
            stress_bot_max = np.dot(Q_transformed, strain_bot)

            stresses_mid[z_index - 1] = stress_mid_max
            stresses_top[z_index - 1] = stress_top_max
            stresses_bot[z_index - 1] = stress_bot_max

            z_index += 1

        return stresses_mid, stresses_top, stresses_bot

    def print_stresses_principal(self, stresses):
        stresses_mid, stresses_top, stresses_bot = self.calculate_stresses_principal(stresses)
        result = "Principal Stress:\n"
        for i, ply in enumerate(self.plies):
            result += f"Lamina {i+1}:\n"
            result += f"  Stress at mid-plane: {stresses_mid[i]} MPa\n"
            result += f"  Stress at top: {stresses_top[i]} MPa\n"
            result += f"  Stress at bottom: {stresses_bot[i]} MPa\n"
        print(result)

    def tsai_wu(self, stresses):
        stresses_mid, stresses_top, stresses_bot = self.calculate_stresses(stresses)
        tsai_wu_total = 0
        for i, ply in enumerate(self.plies):
            F1 = (1 / (ply.sigma_lt * 1e-6)) - (1 / (ply.sigma_cl * 1e-6))
            F11 = 1 / ((ply.sigma_lt * 1e-6) * (ply.sigma_cl * 1e-6))
            F2 = (1 / (ply.sigma_tt * 1e-6)) - (1 / (ply.sigma_ct * 1e-6))
            F22 = 1 / ((ply.sigma_tt * 1e-6) * (ply.sigma_ct * 1e-6))
            F66 = 1 / (ply.tau_lt * 1e-6)**2
            F12 = -np.sqrt(F11 * F22) / 2
            sigma_x, sigma_y, tau_xy = stresses_mid[i]
            tsai_wu_value = F1 * sigma_x + F2 * sigma_y + F11 * sigma_x**2 + F22 * sigma_y**2 + F66 * tau_xy**2 + 2 * F12 * sigma_x * sigma_y
            tsai_wu_total += tsai_wu_value
            if tsai_wu_value >= 1:
                print(f"Lamina {i+1} has failed by Tsai-Wu criterion")
            print(f"Tsai-Wu Lamina {i+1}: {tsai_wu_value}")
        if tsai_wu_total >= 1:
            print(f"Laminate has failed by Tsai-Wu criterion, total: {tsai_wu_total}")
        else:
            print(f"Tsai-Wu total {tsai_wu_total}")

    def print_tsai_wu_factors(self):
        for i, ply in enumerate(self.plies):
            F1 = (1 / (ply.sigma_lt * 1e-6)) - (1 / (ply.sigma_cl * 1e-6))
            F11 = 1 / ((ply.sigma_lt * 1e-6) * (ply.sigma_cl * 1e-6))
            F2 = (1 / (ply.sigma_tt * 1e-6)) - (1 / (ply.sigma_ct * 1e-6))
            F22 = 1 / ((ply.sigma_tt * 1e-6) * (ply.sigma_ct * 1e-6))
            F66 = 1 / (ply.tau_lt * 1e-6)**2
            F12 = -np.sqrt(F11 * F22) / 2
            print(f"Lamina {i+1}:")
            print(f"F1: {F1}")
            print(f"F11: {F11}")
            print(f"F12: {F12}")
            print(f"F2: {F2}")
            print(f"F22: {F22}")
            print(f"F66: {F66}")
            print(f"total: {F1 + F2 + F11 + F22 + F66 + 2 * F12}")

    def __str__(self):
        result = "Laminate Characteristics:\n"
        result += f"Total Thickness: {self.thickness():.3f} mm\n"
        A, B, D = self.calculate_A_B_D_matrices()
        ABD, inv_ABD = self.calculate_ABD_matrix()
        result += f"\nMatrix A:\n{A}"
        result += f"\nMatrix B:\n{B}"
        result += f"\nMatrix D:\n{D}"
        result += f"\nABD Matrix:\n{ABD}"
        result += f"\nInverse ABD Matrix:\n{inv_ABD}"
        return result

class Stress:
    def __init__(self, Nx=0, Ny=0, Nxy=0, Mx=0, My=0, Mxy=0):
        self.Nx = Nx  # Normal force in x-direction (N/mm)
        self.Ny = Ny  # Normal force in y-direction (N/mm)
        self.Nxy = Nxy  # Shear force in xy-plane (N/mm)
        self.Mx = Mx  # Bending moment about x-axis (Nmm)
        self.My = My  # Bending moment about y-axis (Nmm)
        self.Mxy = Mxy  # Twisting moment (Nmm)

    def __str__(self):
        result = "Stresses and Moments:\n"
        result += f"Nx = {self.Nx:.4e} N/mm\n"
        result += f"Ny = {self.Ny:.4e} N/mm\n"
        result += f"Nxy = {self.Nxy:.4e} N/mm\n"
        result += f"Mx = {self.Mx:.4e} Nmm\n"
        result += f"My = {self.My:.4e} Nmm\n"
        result += f"Mxy = {self.Mxy:.4e} Nmm\n"
        return result

    def to_array(self):
        return np.array([self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy])

# Example usage of the library
if __name__ == "__main__":
    # Define matrix and reinforcement properties
    matrix = Matrix(E=4.3e9, nu=0.4, rho=1.2, G=1640e6, Rm=135e6)
    reinforcement = Reinforcement(E_l=390e9, E_t=12e9, G=20000e6, nu_lt=0.35, nu_tl=0.01, rho=1.8, Rm=2540e6)

    # Define a single ply
    ply = Lamina(matrix, reinforcement, fiber_volume_ratio=0.4, fiber_areal_weight=200, orientation=45,
              sigma_lt=2978e6, sigma_tt=36e6, sigma_cl=2100e6, sigma_ct=200e6, tau_lt=70e6)

    # Define applied stresses and moments
    stresses = Stress(Nx=80, My=10)

    # Define the laminate with a single ply
    laminate = Laminate([ply])

    # Print stiffness matrix for verification
    print("Stiffness Matrix Q and Transformed Q:")
    S, S_transformed, Q, Q_transformed = ply.stiffness_matrix()
    print("Q:\n", Q)
    print("Q_transformed:\n", Q_transformed)

    # Calculate and print stresses
    laminate.print_stresses(stresses)

    # Calculate and print Tsai-Wu failure criterion
    laminate.tsai_wu(stresses)



#being worked on the LaminaMix class need to have the strength calculated form values of the matrix and the reinforcement
# class LaminaMix:
#     def __init__(self, matrix, reinforcement, fiber_volume_ratio, fiber_areal_weight, orientation=0,
#                  sigma_lt=0, sigma_tt=0, sigma_cl=0, sigma_ct=0, tau_lt=0):
#         self.matrix = matrix
#         self.reinforcement = reinforcement
#         self.fiber_volume_ratio = fiber_volume_ratio
#         self.orientation = orientation  # Fiber orientation in degrees
#         self.fiber_areal_weight = fiber_areal_weight  # Fiber areal weight in g/m2
#         self.sigma_lt = sigma_lt  # Longitudinal tensile strength in Pa
#         self.sigma_tt = sigma_tt  # Transverse tensile strength in Pa
#         self.sigma_ct = sigma_ct  # Transverse compressive strength in Pa
#         self.sigma_cl = sigma_cl  # Longitudinal compressive strength in Pa
#         self.tau_lt = tau_lt  # Longitudinal-transverse shear strength in Pa

#     def thickness(self):
#         # Calculate ply thickness from fiber areal weight
#         matrix_volume_ratio = 1 - self.fiber_volume_ratio
#         thickness = (self.fiber_areal_weight * 1e-3) / (matrix_volume_ratio * self.reinforcement.rho)
#         return thickness

#     def stiffness_matrix(self):
#         # Calculate the stiffness matrix of the ply
#         matrix_volume_ratio = 1 - self.fiber_volume_ratio
#         E_l = (self.fiber_volume_ratio * self.matrix.E) + (matrix_volume_ratio * self.reinforcement.E_l)
#         E_t = self.matrix.E / (self.fiber_volume_ratio + ((self.matrix.E / self.reinforcement.E_t) * matrix_volume_ratio))
#         nu_lt = self.reinforcement.nu_lt * matrix_volume_ratio + self.matrix.nu * self.fiber_volume_ratio
#         nu_tl = nu_lt * (E_t / E_l)
#         G_lt = self.matrix.G / (self.fiber_volume_ratio + ((self.matrix.G / self.reinforcement.G) * matrix_volume_ratio))

#         theta = np.radians(self.orientation)
#         c = np.cos(theta)
#         s = np.sin(theta)

#         S = np.array([[1 / (E_l * 1e-6), -nu_tl / (E_t * 1e-6), 0],
#                        [-nu_lt / (E_l * 1e-6), 1 / (E_t * 1e-6), 0],
#                        [0, 0, 1 / (G_lt * 1e-6)]])
#         Q = np.linalg.inv(S)

#         T = np.array([[c**2, s**2, 2 * c * s],
#                       [s**2, c**2, -2 * c * s],
#                       [-c * s, c * s, c**2 - s**2]])
#         T1 = np.array([[c**2, s**2, -2 * c * s],
#                        [s**2, c**2, 2 * c * s],
#                        [c * s, -c * s, c**2 - s**2]])

#         TS = np.dot(T.T, S)
#         T1Q = np.dot(T1, Q)
#         S_transformed = np.dot(TS, T)
#         Q_transformed = np.dot(T1Q, T1.T)

#         return S, S_transformed, Q, Q_transformed

#     def calculate_Rm(self):
#         Rm = self.reinforcement.Rm * ((1 - self.fiber_volume_ratio) + ((self.matrix.E / self.reinforcement.E_l) * self.fiber_volume_ratio))
#         return Rm

#     def __str__(self):
#         thickness = self.thickness()
#         S, S_transformed, Q, Q_transformed = self.stiffness_matrix()
#         Rm = self.calculate_Rm()
#         matrix_volume_ratio = 1 - self.fiber_volume_ratio
#         E_l = (self.fiber_volume_ratio * self.matrix.E) + (matrix_volume_ratio * self.reinforcement.E_l)
#         E_t = self.matrix.E / (self.fiber_volume_ratio + ((self.matrix.E / self.reinforcement.E_t) * matrix_volume_ratio))
#         nu_lt = self.reinforcement.nu_lt * matrix_volume_ratio + self.matrix.nu * self.fiber_volume_ratio
#         nu_tl = nu_lt * (E_t / E_l)
#         G_lt = self.matrix.G / (self.fiber_volume_ratio + ((self.matrix.G / self.reinforcement.G) * matrix_volume_ratio))
#         rho = matrix_volume_ratio * self.reinforcement.rho + self.fiber_volume_ratio * self.matrix.rho

#         result = "Lamina Characteristics:\n"
#         result += f"Rm: {Rm * 1e-6:.3f} MPa\n"
#         result += f"Density: {rho:.3f} g/cm3\n"
#         result += f"Elastic Moduli: E_l = {E_l * 1e-9:.3f} GPa, E_t = {E_t * 1e-9:.3f} GPa\n"
#         result += f"Shear Modulus: G_lt = {G_lt * 1e-9:.3f} GPa\n"
#         result += f"Poisson's Ratios: nu_tl = {nu_tl:.3f}, nu_lt = {nu_lt:.3f}\n"
#         result += f"Thickness: {thickness:.3f} mm\n"
#         result += f"\nStiffness Matrix Q:\n{Q}"
#         result += f"\nTransformed Stiffness Matrix Q_transformed:\n{Q_transformed}"
#         result += f"\nCompliance Matrix S:\n{S}"
#         result += f"\nTransformed Compliance Matrix S_transformed:\n{S_transformed}"

#         return result