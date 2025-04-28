import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add the parent directory of the script's directory to the Python path
script_dir = os.path.dirname(os.path.abspath(__file__)) # c:\LSC\CLT_Bridge
parent_dir = os.path.dirname(script_dir) # c:\LSC
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir) # Add c:\LSC to sys.path

# Modify imports to use the full path from the perspective of c:\LSC
from code.materials.isotropic import Matrix, Reinforcement, IsotropicMaterial
from code.materials.composite import Laminate, Lamina, LaminaMix
from code.materials.sandwich import Sandwich
from code.materials.panel import Panel

def main():
    # 1. Create materials

    divinycell_H45= IsotropicMaterial(
        E= 50e6,    # Pa # Compressive modulus
        nu=0.4,
        G=15e6,    # Pa
        rho=48 ,   # kg/m^3
        Rm=0.6e6,      # Pa Compressive Strength 
        sigma_shear=0.56e6 # Pa Shear Strength
    )

    divinycell_H100= IsotropicMaterial(
        E= 135 * 1e6,    # Pa # Compressive modulus
        nu=0.4,
        G=35 * 1e6,    # Pa
        rho=100 ,   # kg/m^3
        Rm=2* 1e6,     # Pa Compressive Strength 
        sigma_shear=1.6e6 # Pa Shear Strength
    )

    divinycell_H200= IsotropicMaterial(
        E= 310 * 1e6,    # Pa # Compressive modulus
        nu=0.4,
        G=73 * 1e6,    # Pa
        rho=200 ,   # kg/m^3
        Rm=5.4 * 1e6,      # Pa Compressive Strength 
        sigma_shear=3.5e6 # Pa Shear Strength
    )

    # 2. Create composite laminate
    lamina0 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        rho=1600,    # kg/m^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=0,  # degrees
        thickness=6 #mm
    )

    lamina90 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1600,    # kg/m^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=90,  # degrees 
        thickness=2 #mm
    )

    lamina45 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1600,    # kg/m^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=45,  # degrees
        thickness=1.
    )

    lamina_45 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1600,    # kg/m^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=-45,  # degrees
        thickness=1. #mm  
    )

    laminate = Laminate(
        plies=[lamina45, lamina_45, lamina90, lamina0, lamina0, lamina90, lamina_45, lamina45],  # 8 layers
    )

    # 3. Define requirements and loads
    point_load = 2000 * 9.81 # N
    uniform_load = 5000 * 9.81 # N/m²
    panel_width = 4.200 # m
    panel_length = 4.0 # m

    # 4. Parametric Study Setup
    tc_values = np.linspace(1, 200, 50) # Core thickness from 1mm to 200mm (50 steps)
    sf_shear_list = []
    sf_wrinkling_list = []
    sf_compression_list = []
    # sf_face_failure_list = [] # Assuming check_against_face_failure returns a safety factor

    # 5. Run Parametric Study
    for tc in tc_values:
        # Create sandwich with current core thickness (convert mm to m)
        sandwich = Sandwich(
            composite_material=laminate,
            core_material=divinycell_H200,
            tc=tc  # Convert mm 
        )

        # Create panel
        panel = Panel(
            sandwich=sandwich,
            width=panel_width,
            length=panel_length,
            point_load=point_load,
            distributed_load=uniform_load
        )

        # Calculate and store safety factors
        sf_shear = panel.check_against_core_shear_failure()
        sf_wrinkling = panel.check_against_face_wrinkling()
        sf_compression = panel.check_against_core_compression_failure()
        # sf_face = panel.check_against_face_failure(False) # Assuming it returns SF

        sf_shear_list.append(sf_shear)
        sf_wrinkling_list.append(sf_wrinkling)
        sf_compression_list.append(sf_compression)
        # sf_face_failure_list.append(sf_face)

    # 6. Plotting Results
    # Create separate figures for each safety factor

    # Plot Core Shear Failure SF
    plt.figure(figsize=(10, 6))
    plt.plot(tc_values, sf_shear_list, label='Core Shear Failure SF', marker='o')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Core Shear Failure Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot Face Wrinkling SF
    plt.figure(figsize=(10, 6))
    plt.plot(tc_values, sf_wrinkling_list, label='Face Wrinkling SF', marker='s', color='orange')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Face Wrinkling Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot Core Compression SF
    plt.figure(figsize=(10, 6))
    plt.plot(tc_values, sf_compression_list, label='Core Compression SF', marker='^', color='green')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Core Compression Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot Face Failure SF (if calculated)
    # if sf_face_failure_list:
    #     plt.figure(figsize=(10, 6))
    #     plt.plot(tc_values, sf_face_failure_list, label='Face Failure SF', marker='d', color='red')
    #     plt.xlabel('Core Thickness (tc) [mm]')
    #     plt.ylabel('Safety Factor')
    #     plt.title('Face Failure Safety Factor vs. Core Thickness')
    #     plt.legend()
    #     plt.grid(True)
    #     plt.ylim(bottom=0)

    plt.show() # Show all plots

    # --- Remove or comment out the old single panel analysis --- 
    # print(f"panel:\n{panel}")
    # panel.check_against_face_failure(True)
    # safety_factor_shear = panel.check_against_core_shear_failure()
    # print(f"Safety factor against core shear failure: {safety_factor_shear:.2f}")
    # safety_factor_wrinkiling = panel.check_against_face_wrinkling()
    # print(f"Safety factor against face wrinkling: {safety_factor_wrinkiling}")
    # safety_factor_compression = panel.check_against_core_compression_failure()
    # print(f"Safety factor against core compression : {safety_factor_compression}")

    # Remove the panel_fail block as sandwich_failure is not defined
    # panel_fail = Panel(
    #     sandwich=sandwich_failure,
    #     width=4.200,  # m
    #     length=4,   # m
    #     point_load=point_load,  # N
    #     distributed_load=uniform_load  # N/m²
    # )

    # Remove bridge section
    # # 3. Create bridge
    # ...
    # print(f"Max deflection: {results['max_deflection']:.4f} m")

if __name__ == "__main__":
    main()