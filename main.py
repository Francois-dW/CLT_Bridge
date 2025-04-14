"""
Main module for bridge analysis program.
"""

from .bridge import Bridge
from .utils.visualization import plot_deflection, plot_stress_distribution

def analyze_bridge(bridge_config):
    """
    Main analysis function for the bridge.
    
    Args:
        bridge_config (dict): Configuration dictionary for the bridge
        
    Returns:
        dict: Analysis results including deflections, stresses, etc.
    """
    # Create bridge instance
    bridge = Bridge.from_config(bridge_config)
    
    # Perform analysis
    results = bridge.analyze()
    
    # Visualize results
    if bridge_config.get('visualize', True):
        plot_deflection(results['deflection'])
        plot_stress_distribution(results['stress'])
    
    return results

if __name__ == "__main__":
    # Example configuration
    config = {
        "length": 20.0,  # meters
        "supports": ["pinned", "pinned"],
        "loading": {"type": "uniform", "magnitude": 5000},  # N/m
        "cross_section": {
            "type": "sandwich",
            "face_thickness": 0.01,
            "core_thickness": 0.1,
            "face_material": {"E": 70e9, "nu": 0.3},
            "core_material": {"E": 0.5e9, "nu": 0.35}
        }
    }
    
    results = analyze_bridge(config)
    print("Maximum deflection:", results['max_deflection'])