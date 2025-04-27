import numpy as np

class Load:
    def __init__(self, Nx=0, Ny=0, Nxy=0, Mx=0, My=0, Mxy=0, array=None):
        if array is not None:
            array = np.asarray(array)
            if array.shape != (6,):
                raise ValueError("Array must have 6 elements: [Nx, Ny, Nxy, Mx, My, Mxy]")
            self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy = array
        else:
            self.Nx = Nx  # Normal force in x-direction (N/mm)
            self.Ny = Ny  # Normal force in y-direction (N/mm)
            self.Nxy = Nxy  # Shear force in xy-plane (N/mm)
            self.Mx = Mx  # Bending moment about x-axis (Nmm)
            self.My = My  # Bending moment about y-axis (Nmm)
            self.Mxy = Mxy  # Twisting moment (Nmm)
        self.array = np.array([self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy])

    def __str__(self):
        result = "Stresses and Moments:\n"
        result += f"Nx = {self.Nx:.4e} N/mm\n"
        result += f"Ny = {self.Ny:.4e} N/mm\n"
        result += f"Nxy = {self.Nxy:.4e} N/mm\n"
        result += f"Mx = {self.Mx:.4e} Nmm\n"
        result += f"My = {self.My:.4e} Nmm\n"
        result += f"Mxy = {self.Mxy:.4e} Nmm\n"
        return result