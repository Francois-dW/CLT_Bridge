from code.materials.isotropic import Matrix, Reinforcement, IsotropicMaterial
from code.materials.composite import Laminate, Lamina
from code.core.beam import EulerBeam

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
    lamina = Lamina(
        reinforcement=carbon_fiber,
        matrix=epoxy,
        thickness=0.001,  # 1mm
        angle=0  # 0 degrees
    )
    laminate = Laminate(
        layers=[lamina, lamina, lamina],  # 3 layers of carbon fiber
        thickness=0.003  # 3mm total thickness
    )
    
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