class Matrix:
    def __init__(self, E=0, G=0, nu=0, rho=0, Rm=0):
        self.E = E  # Young's modulus in Pa
        self.nu = nu  # Poisson's ratio
        self.rho = rho  # Density in g/cm3
        self.G = G  # Shear modulus in Pa
        self.Rm = Rm  # Tensile strength in Pa
    
    def __repr__(self):
        return f"Matrix(E={self.E}, nu={self.nu}, G={self.G}, rho={self.rho}, Rm={self.Rm})"

class Reinforcement:
    def __init__(self, E_l=0, E_t=0, G=0, nu_lt=0, nu_tl=0, rho=0, Rm=0):
        self.E_l = E_l  # Longitudinal Young's modulus in Pa
        self.E_t = E_t  # Transverse Young's modulus in Pa
        self.nu_lt = nu_lt  # Longitudinal Poisson's ratio
        self.nu_tl = nu_tl  # Transverse Poisson's ratio
        self.rho = rho  # Density in g/cm3
        self.G = G  # Shear modulus in Pa
        self.Rm = Rm  # Tensile strength in Pa

    def __repr__(self):
        return f"Reinforcement(E_l={self.E_l}, E_t={self.E_t}, G={self.G}, nu_lt={self.nu_lt}, nu_tl={self.nu_tl}, rho={self.rho}, Rm={self.Rm})"
    
class IsotropicMaterial:
    def __init__(self, E, nu, G, rho, Rm, sigma_t= None, sigma_c= None, sigma_shear= None):
        self.E = E  # Young's modulus in Pa
        self.nu = nu  # Poisson's ratio
        self.G = G  # Shear modulus in Pa
        self.rho = rho  # Density in g/cm3
        self.Rm = Rm  # Tensile strength in Pa
        self.sigma_t = sigma_t # Tensile strength in Pa
        self.sigma_c = sigma_c # Compressive strength in Pa
        self.sigma_shear = sigma_shear # Shear strength in Pa

    def __repr__(self):
        return f"IsotropicMaterial(E={self.E}, nu={self.nu}, G={self.G}, rho={self.rho}, Rm={self.Rm})"