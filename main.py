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

    divinycell_H45 = IsotropicMaterial(
        E=50e6,          # Pa - Compressive Modulus (Nominal)
        nu=0.4,
        G=15e6,          # Pa - Shear Modulus (Nominal)
        rho=48,          # kg/m³
        Rm=0.6e6,        # Pa - Compressive Strength (Nominal)
        sigma_shear=0.56e6  # Pa - Shear Strength (Nominal)
    )

    divinycell_H60 = IsotropicMaterial(
        E=70e6,
        nu=0.4,
        G=20e6,
        rho=60,
        Rm=0.9e6,
        sigma_shear=0.76e6
    )

    divinycell_H80 = IsotropicMaterial(
        E=90e6,
        nu=0.4,
        G=27e6,
        rho=80,
        Rm=1.4e6,
        sigma_shear=1.15e6
    )

    divinycell_H100 = IsotropicMaterial(
        E=135e6,
        nu=0.4,
        G=35e6,
        rho=100,
        Rm=2.0e6,
        sigma_shear=1.6e6
    )

    divinycell_H130 = IsotropicMaterial(
        E=170e6,
        nu=0.4,
        G=50e6,
        rho=130,
        Rm=3.0e6,
        sigma_shear=2.2e6
    )

    divinycell_H160 = IsotropicMaterial(
        E=200e6,
        nu=0.4,
        G=60e6,
        rho=160,
        Rm=3.4e6,
        sigma_shear=2.6e6
    )

    divinycell_H200 = IsotropicMaterial(
        E=310e6,
        nu=0.4,
        G=73e6,
        rho=200,
        Rm=5.4e6,
        sigma_shear=3.5e6
    )

    divinycell_H250 = IsotropicMaterial(
        E=400e6,
        nu=0.4,
        G=97e6,
        rho=250,
        Rm=7.2e6,
        sigma_shear=4.5e6
    )

    # Group all foams for iteration
    foams = {
        'H45': divinycell_H45,
        'H60': divinycell_H60,
        'H80': divinycell_H80,
        'H100': divinycell_H100,
        'H130': divinycell_H130,
        'H160': divinycell_H160,
        'H200': divinycell_H200,
        'H250': divinycell_H250
    }
    # Assign colors and markers for all foams (extend as needed)
    foam_colors = {
        'H45': 'blue',
        'H60': 'cyan',
        'H80': 'purple',
        'H100': 'orange',
        'H130': 'brown',
        'H160': 'magenta',
        'H200': 'green',
        'H250': 'black'
    }
    foam_markers = {
        'H45': 'o',
        'H60': 'v',
        'H80': 'D',
        'H100': 's',
        'H130': 'P',
        'H160': 'X',
        'H200': '^',
        'H250': '*'
    }


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
        thickness=5 #mm
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
        thickness=0.25    )

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
        thickness=0.25 #mm
    )

    laminate = Laminate(
        plies=[lamina45, lamina_45, lamina90, lamina0, lamina0, lamina90, lamina_45, lamina45],  # 8 layers
    )

    # Additional design with specified ply thicknesses
    lamina0_6mm_new = Lamina(
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
        angle=0,  # degrees
        thickness=6 # mm
    )

    lamina90_4mm_new = Lamina(
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
        thickness=4 # mm
    )

    lamina45_1mm_new = Lamina(
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
        thickness=1 # mm
    )

    lamina_45_1mm_new = Lamina(
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
        thickness=1 # mm
    )

    laminate_custom_2 = Laminate(
        plies=[lamina45_1mm_new, lamina_45_1mm_new, lamina90_4mm_new, lamina0_6mm_new, lamina0_6mm_new, lamina90_4mm_new, lamina_45_1mm_new, lamina45_1mm_new],  # 8 layers
    )

    # Additional design with all ply thicknesses set to 4mm
    lamina0_4mm = Lamina(
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
        angle=0,  # degrees
        thickness=4 # mm
    )

    lamina90_4mm = Lamina(
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
        thickness=4 # mm
    )

    lamina45_4mm = Lamina(
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
        thickness=4 # mm
    )

    lamina_45_4mm = Lamina(
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
        thickness=4 # mm
    )

    laminate_4mm = Laminate(
        plies=[lamina45_4mm, lamina_45_4mm, lamina90_4mm, lamina0_4mm, lamina0_4mm, lamina90_4mm, lamina_45_4mm, lamina45_4mm],  # 8 layers
    )

    # 3. Define requirements and loads
    point_load = 2000 * 9.81 # N
    uniform_load = 5000 * 9.81 # N/m²
    panel_width = 4.200 # m
    panel_length = 4.0 # m

    # 4. Parametric Study Setup
    tc_values = np.linspace(1, 200, 50) # Core thickness from 1mm to 200mm (50 steps)
    # Dictionary to store results for each foam type
    results = {foam_name: {'shear': [], 'wrinkling': [], 'compression': []} for foam_name in foams}


    # 5. Run Parametric Study for each foam type
    failures = {foam_name: {'tsai_wu': [], 'tsai_hill': [], 'max_stress': [], 'deflection': []} for foam_name in foams}
    selected_laminate = laminate_custom_2  # Change this variable to switch laminate designs
    for foam_name, core_material in foams.items():
        print(f"Running analysis for foam: {foam_name}")
        for tc in tc_values:
            sandwich = Sandwich(
                composite_material=selected_laminate,  # Use the global variable here
                core_material=core_material,
                tc=tc
            )
            panel = Panel(
                sandwich=sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
            sf_shear = panel.check_against_core_shear_failure()
            sf_wrinkling = panel.check_against_face_wrinkling()
            sf_compression = panel.check_against_core_compression_failure()
            results[foam_name]['shear'].append(sf_shear)
            results[foam_name]['wrinkling'].append(sf_wrinkling)
            results[foam_name]['compression'].append(sf_compression)

            # Use panel.check_against_face_failure for all failure checks
            _, fail_tsai_wu, _, fail_tsai_hill, _, _, fail_max_bottom, fail_max_top = panel.check_against_face_failure()
            failures[foam_name]['tsai_wu'].append(fail_tsai_wu)
            failures[foam_name]['tsai_hill'].append(fail_tsai_hill)
            failures[foam_name]['max_stress'].append(fail_max_bottom or fail_max_top)
            # Deflection check (40mm = 0.04m)
            deflection_fail = (panel.delta_max_distributed_load > 0.04) or (panel.delta_max_point_load > 0.04)
            failures[foam_name]['deflection'].append(deflection_fail)

    # 6. Plotting Results
    # Helper to get color for each point
    def get_point_color(failures_list, idx, default_color):
        return 'red' if failures_list[idx] else default_color

    # Plot Core Shear Failure SF for all foams
    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        last_fail_idx = None
        for i in range(len(tc_values)):
            fail = failures[foam_name]['tsai_wu'][i] or failures[foam_name]['deflection'][i]
            if fail:
                last_fail_idx = i
        for i, tc in enumerate(tc_values):
            fail = failures[foam_name]['tsai_wu'][i] or failures[foam_name]['deflection'][i]
            if not fail:
                plt.scatter(tc, results[foam_name]['shear'][i], color=foam_colors[foam_name], marker='+', s=30)
        if last_fail_idx is not None:
            plt.scatter(tc_values[last_fail_idx], results[foam_name]['shear'][last_fail_idx], color='red', marker='+', s=30)
        plt.plot(tc_values, results[foam_name]['shear'], color=foam_colors[foam_name], alpha=0.5, label=f'{foam_name}')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Core Shear Failure Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot Face Wrinkling SF for all foams
    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        last_fail_idx = None
        for i in range(len(tc_values)):
            fail = failures[foam_name]['tsai_hill'][i] or failures[foam_name]['deflection'][i]
            if fail:
                last_fail_idx = i
        for i, tc in enumerate(tc_values):
            fail = failures[foam_name]['tsai_hill'][i] or failures[foam_name]['deflection'][i]
            if not fail:
                plt.scatter(tc, results[foam_name]['wrinkling'][i], color=foam_colors[foam_name], marker='+', s=30)
        if last_fail_idx is not None:
            plt.scatter(tc_values[last_fail_idx], results[foam_name]['wrinkling'][last_fail_idx], color='red', marker='+', s=30)
        plt.plot(tc_values, results[foam_name]['wrinkling'], color=foam_colors[foam_name], alpha=0.5, label=f'{foam_name}')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Face Wrinkling Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot Core Compression SF for all foams
    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        last_fail_idx = None
        for i in range(len(tc_values)):
            fail = failures[foam_name]['max_stress'][i] or failures[foam_name]['deflection'][i]
            if fail:
                last_fail_idx = i
        for i, tc in enumerate(tc_values):
            fail = failures[foam_name]['max_stress'][i] or failures[foam_name]['deflection'][i]
            if not fail:
                plt.scatter(tc, results[foam_name]['compression'][i], color=foam_colors[foam_name], marker='+', s=30)
        if last_fail_idx is not None:
            plt.scatter(tc_values[last_fail_idx], results[foam_name]['compression'][last_fail_idx], color='red', marker='+', s=30)
        plt.plot(tc_values, results[foam_name]['compression'], color=foam_colors[foam_name], alpha=0.5, label=f'{foam_name}')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Safety Factor')
    plt.title('Core Compression Safety Factor vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # 7. Calculate and plot weight vs. core thickness for each foam
    weights = {foam_name: [] for foam_name in foams}
    for foam_name, core_material in foams.items():
        for tc in tc_values:
            sandwich = Sandwich(
                composite_material=selected_laminate,
                core_material=core_material,
                tc=tc
            )
            panel = Panel(
                sandwich=sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
            weights[foam_name].append(panel.mass)

    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        last_fail_idx = None
        for i in range(len(tc_values)):
            fail = failures[foam_name]['tsai_wu'][i] or failures[foam_name]['deflection'][i]
            if fail:
                last_fail_idx = i
        for i, tc in enumerate(tc_values):
            fail = failures[foam_name]['tsai_wu'][i] or failures[foam_name]['deflection'][i]
            if not fail:
                plt.scatter(tc, weights[foam_name][i], color=foam_colors[foam_name], marker='+', s=30)
        if last_fail_idx is not None:
            plt.scatter(tc_values[last_fail_idx], weights[foam_name][last_fail_idx], color='red', marker='+', s=30)
        plt.plot(tc_values, weights[foam_name], color=foam_colors[foam_name], alpha=0.5, label=f'{foam_name}')
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Panel Weight [kg]')
    plt.title('Panel Weight vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    

    # 8. Plot max deflection vs. core thickness for point load and distributed load
    max_deflections_point = {foam_name: [] for foam_name in foams}
    max_deflections_distributed = {foam_name: [] for foam_name in foams}
    for foam_name, core_material in foams.items():
        for tc in tc_values:
            sandwich = Sandwich(
                composite_material=selected_laminate,
                core_material=core_material,
                tc=tc
            )
            panel = Panel(
                sandwich=sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
            max_deflections_point[foam_name].append(panel.delta_max_point_load)
            max_deflections_distributed[foam_name].append(panel.delta_max_distributed_load)

    # Plot max deflection for point load
    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        plt.plot(tc_values, max_deflections_point[foam_name], label=f'{foam_name}', color=foam_colors[foam_name])
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Max Deflection (m)')
    plt.title('Max Deflection (Point Load) vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)

    # Plot max deflection for distributed load
    plt.figure(figsize=(10, 6))
    for foam_name in foams:
        plt.plot(tc_values, max_deflections_distributed[foam_name], label=f'{foam_name}', color=foam_colors[foam_name])
    plt.xlabel('Core Thickness (tc) [mm]')
    plt.ylabel('Max Deflection (m)')
    plt.title('Max Deflection (Distributed Load) vs. Core Thickness')
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    
    plt.show()

    good_sandwich = Sandwich(
        composite_material=laminate,
        core_material=divinycell_H130,
        tc=170 # mm
    )

    design1_sandwich = Sandwich(
        composite_material=laminate_4mm,
        core_material=divinycell_H250,
        tc=120 # mm
    )

    design2_sandwich = Sandwich(
        composite_material=laminate_custom_2,
        core_material=divinycell_H100,
        tc=140 # mm
    )

    good_panel = Panel(
                sandwich=good_sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
    
    design1_panel = Panel(
                sandwich=design1_sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
    design2_panel = Panel(
                sandwich=design2_sandwich,
                width=panel_width,
                length=panel_length,
                point_load=point_load,
                distributed_load=uniform_load
            )
    
    security_compression = design1_panel.check_against_core_compression_failure()
    design1_panel.check_against_face_failure(print_only_max=True)
    security_wrinkling = design1_panel.check_against_face_wrinkling()
    security_shear = design1_panel.check_against_core_shear_failure()
    print(f"Security Compression: {security_compression}")
    print(f"Security Wrinkling: {security_wrinkling}")
    print(f"Security Shear: {security_shear}")
    print(design1_panel)   
if __name__ == "__main__":
    main()