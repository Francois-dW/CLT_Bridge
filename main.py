from code.materials.isotropic import Matrix, Reinforcement, IsotropicMaterial
from code.materials.composite import Laminate, Lamina, LaminaMix
from code.materials.sandwich import Sandwich
from code.materials.panel import Panel

def main():
    # 1. Create materials
    carbon_fiber = Reinforcement(
        E_l=230e9,  # Pa
        E_t=15e9,   # Pa
        G=30e9,     # Pa
        nu_lt=0.2,
        nu_tl=0.3,
        rho=1.6,    # g/cm^3
        Rm=4e9      # Pa
    )
    epoxy = Matrix(
        E=3e9,      # Pa
        G=1.2e9,    # Pa
        nu=0.35,
        rho=1.2,    # g/cm^3
        Rm=80e6     # Pa
    )

    divinycell = IsotropicMaterial(
        E=2.5e9,    # Pa
        nu=0.3,
        G=1.2e9,    # Pa
        rho=0.05,   # g/cm^3
        Rm=1e6      # Pa
    )

    
    # 2. Create composite laminate
    lamina = LaminaMix(
        reinforcement=carbon_fiber,
        matrix=epoxy,
        fiber_volume_ratio=0.6,  # 60% fiber volume fraction
        fiber_areal_weight=200,  # 200 g/m²
        # thickness=0.001,  # 1mm
        # angle=0  # 0 degrees
    )
    laminate = Laminate(
        layers=[lamina, lamina, lamina],  # 3 layers of carbon fiber
        thickness=0.003  # 3mm total thickness
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