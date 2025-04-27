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

    divinycell_H45= IsotropicMaterial(
        E= 50 * 1e6,    # Pa # Compressive modulus
        nu=0.4,
        G=15 * 1e6,    # Pa
        rho=48 ,   # kg/m^3
        Rm=0.6 * 1e6      # Pa Compressive Strength 
    )

    divinycell_H100= IsotropicMaterial(
        E= 135 * 1e6,    # Pa # Compressive modulus
        nu=0.4,
        G=35 * 1e6,    # Pa
        rho=100 ,   # kg/m^3
        Rm=2* 1e6      # Pa Compressive Strength 
    )

    divinycell_H200= IsotropicMaterial(
        E= 310 * 1e6,    # Pa # Compressive modulus
        nu=0.4,
        G=73 * 1e6,    # Pa
        rho=200 ,   # kg/m^3
        Rm=5.4 * 1e6      # Pa Compressive Strength 
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
    #print(lamina90)

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
        plies=[lamina45, lamina_45, lamina0, lamina90, lamina90, lamina0, lamina_45, lamina45],  # 4 layers of carbon fiber
    )
    laminate_break = Laminate(
        plies=[lamina0]  # 4 layers of carbon fiber
    )

    # 3. Create sandwich & panel
    sandwich = Sandwich(
        composite_material=laminate,
        core_material=divinycell_H200,
        tc=200 # mm # Assuming core thickness unit is meters based on panel dimensions
    )
    # print(f"Composite equivalent elastic modulus: {sandwich.composite_material.E} Pa")
    # print(sandwich)

    #requirements:
    # total thickness of the panel (tf + tc) <= 250 mm
    # able to withstand a point load of 2000 kg * 9.81 m/s² 
    # able to wistand distributed load of 5000 kg * 9.81 m/s²
    #max deflection of panel <= 40 mm
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

    panel.check_against_face_failure(True)
    # panel.check_against_core_shear_failure()
    # panel.check_against_face_wrinkling()
    # panel.check_against_core_compression_failure()
    # panel.check_against_face_failure()
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