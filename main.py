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
        thickness=6. #mm
    )
    #print(lamina90)

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



    # 3. Create sandwich & panel
    sandwich = Sandwich(
        composite_material=laminate,
        core_material=divinycell_H100,
        tc=0.200  # m # Assuming core thickness unit is meters based on panel dimensions
    )

    print(sandwich)

    panel = Panel(
        sandwich=sandwich,
        width=4.200,  # m
        length=4.000   # m
    )


    # 4. Perform calculations
    panel_weight = panel.calculate_weight()
    print(f"Panel weight: {panel_weight:.4f} N")
    
    # Calculate deflection with uniform load, assuming a uniform load of 5000 N/m²
    uniform_load = 5000 * 9.81 #N/m²
    delta_max, max_shear_force, max_bending_moment = panel.calculate_beam_response_distributed_load(load=uniform_load)  # m
    print(f"Panel midlength deflection with uniform load: {delta_max:.4f} m")
    print(f"Panel max shear force with uniform load: {max_shear_force:.4f} N")
    print(f"Panel max bending moment with uniform load: {max_bending_moment:.4f} Nm")

    # Calculate deflection with point load at midlength, assuming a point load of 2000 N without areal distribution
    point_load = 2000  * 9.81 # N
    delta_max, max_shear_force, max_bending_moment = panel.calculate_beam_response_point_load(load=point_load)  # m
    print(f"Panel midlength deflection with point load: {delta_max:.4f} m")
    print(f"Panel max shear force with point load: {max_shear_force:.4f} N")
    print(f"Panel max bending moment with point load: {max_bending_moment:.4f} Nm")   

    
    print(panel)


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