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
        rho=0.05,   # g/cm^3
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
        sigma_shear=50e6,  # Pa
        angle=0,  # degrees    
    )

    lamina90 = Lamina(
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
        sigma_shear=50e6,  # Pa
        angle=90,  # degrees    
    )
    laminate = Laminate(
        plies=[lamina0, lamina90, lamina90, lamina0],  # 3 layers of carbon fiber
        ply_thicknesses=0.003  # 3mm total thickness
    )
    

    # 3. Create sandwich & panel
    sandwich = Sandwich(
        face_material=laminate,
        core_material=divinycell,
        tf=0.001,  # 1mm face thickness
        tc=0.02    # 20mm core thickness
    )

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
    panel_deflection = panel.calculate_midlength_deflection_as_beam_uniform_load(load=uniform_load)  # m
    print(f"Panel midlength deflection: {panel_deflection:.4f} m")

    # Calculate deflection with point load at midlength, assuming a point load of 2000 N without areal distribution
    point_load = 2000  * 9.81 # N
    panel_deflection_point = panel.calculate_midlength_deflection_as_beam_point_load(load=point_load)  # m
    print(f"Panel midlength deflection with point load: {panel_deflection_point:.4f} m")    


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