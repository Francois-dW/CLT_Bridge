import sys
import os

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

    divinycell = IsotropicMaterial(
        E=2.5e9,    # Pa
        nu=0.3,
        G=1.2e9,    # Pa
        rho=0.05,   # g/cm^3 # Note: Ensure units are consistent (e.g., kg/m^3)
        Rm=1e6      # Pa
    )

    
    # 2. Create composite laminate
    lamina0 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1.6,    # g/cm^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=0,  # degrees
        thickness=8 
    )

    lamina90 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1.6,    # g/cm^3 # Corrected unit comment consistency
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=90,  # degrees 
        thickness=8. #mm
    )

    lamina45 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1.6,    # g/cm^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=45,  # degrees
        thickness=4.
    )

    lamina_45 = Lamina(
        E1=39e9,   # Pa
        E2=9.8e9,    # Pa
        nu12=0.3,
        G12=2.8e9,    # Pa
        G23=2e9,    # Pa
        G13=2.8e9,    # Pa
        rho=1.6,    # g/cm^3
        sigma_1t=1100e6,  # Pa
        sigma_1c=600e6,   # Pa
        sigma_2t=20e6,    # Pa
        sigma_2c=140e6,    # Pa
        sigma_shear=70e6,  # Pa
        angle=-45,  # degrees
        thickness=4. #mm  
    )

    laminate = Laminate(
        plies=[lamina45, lamina_45, lamina0, lamina90, lamina90, lamina0, lamina_45, lamina45],  # 4 layers of carbon fiber
        density=1900,  # kg/m^3
    )

    test_laminate = Laminate(
        plies=[lamina90],  # 4 layers of carbon fiber
        density=1900,  # kg/m^3
    )
    # print(laminate)
    # 3. Create sandwich & panel
    sandwich = Sandwich(
        composite_material=laminate,
        core_material=divinycell,
        tc=200 # mm # Assuming core thickness unit is meters based on panel dimensions
    )
    # print(f"Composite equivalent elastic modulus: {sandwich.composite_material.E} Pa")
    # print(sandwich)
    point_load = 2000  * 9.81 # N
    uniform_load = 5000 * 9.81 #N/m²
    panel = Panel(
        sandwich=sandwich,
        width=4.200,  # m
        length=4,   # m
        point_load=point_load,  # N
        distributed_load=uniform_load  # N/m²
    )

    print(f"panel:\n{panel}")
    
     

    panel.check_against_face_failure()
    panel.check_against_core_shear_failure()
    panel.check_against_face_wrinkling()
    panel.check_against_core_compression_failure()
    # print(panel)


    # # 3. Create bridge
    # my_bridge = Bridge(
    #     spans=3,
    #     span_lengths=[30.0, 40.0, 30.0],  # meters
    #     cross_section=laminate,
    #     supports=["fixed", "pinned", "pinned", "fixed"]
    # )
    
    # # 4. Perform analysis
    # results = my_bridge.analyze(
    #     loading={"type": "uniform", "magnitude": 5000}  # N/m
    # )
    
    # print(f"Max deflection: {results['max_deflection']:.4f} m")

if __name__ == "__main__":
    main()