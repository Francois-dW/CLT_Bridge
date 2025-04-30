# CLT_Bridge

## Overview
The `CLT_Bridge` project is a computational tool designed to analyze and simulate composite laminate and sandwich structures. It provides a framework for evaluating the mechanical properties, failure criteria, and performance of composite materials under various loading conditions.

# Abstract
The CLT_Bridge project is a computational framework designed for the analysis and optimization of composite laminate and sandwich structures. Leveraging classical laminate theory and advanced failure criteria, the framework enables engineers and researchers to evaluate the mechanical performance of composite materials under various loading conditions. Key features include parametric studies, visualization tools, and support for multiple failure criteria such as Tsai-Wu, Tsai-Hill, and Maximum Stress. The project is modular, allowing for easy integration and extension, and is implemented in Python for accessibility and flexibility. Applications range from aerospace to civil engineering, where lightweight and high-strength materials are critical.

## Features
- **Composite Laminate Analysis**: Supports the creation and analysis of composite laminates with multiple plies, each with customizable material properties and orientations.
- **Sandwich Panel Analysis**: Models sandwich structures with composite skins and isotropic cores, calculating properties like flexural rigidity and shear rigidity.
- **Failure Criteria**: Implements Tsai-Wu, Tsai-Hill, and Maximum Stress failure criteria for composite laminates.
- **Parametric Studies**: Enables parametric studies to evaluate the performance of materials under varying conditions, such as core thickness and applied loads.
- **Visualization**: Provides tools for visualizing results, including safety factors, deflections, and weight comparisons.

## File Structure
- `main.py`: Entry point for running analyses and parametric studies.
- `code/`: Contains the core implementation of the project.
  - `materials/`: Defines material models for isotropic and composite materials, as well as sandwich structures.
    - `composite.py`: Defines the `Lamina` and `Laminate` classes for composite materials.
    - `isotropic.py`: Defines isotropic material properties.
    - `sandwich.py`: Defines the `Sandwich` class for sandwich panel analysis.
    - `panel.py`: Defines the `Panel` class for analyzing sandwich panels.
  - `utils/`: Contains utility functions for visualization and other tasks.
- `README.md`: Documentation for the project.
- `LICENSE`: Licensing information.

## Getting Started
1. Clone the repository:
   ```bash
   git clone <repository-url>
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the main script:
   ```bash
   python main.py
   ```

## Usage
### Defining Materials
Materials can be defined using the classes in `materials/`. For example:
```python
from materials.composite import Lamina, Laminate

lamina = Lamina(
    E1=39e9, E2=9.8e9, G12=2.8e9, G23=2e9, nu12=0.3,
    sigma_1t=1100e6, sigma_1c=600e6, sigma_2t=20e6, sigma_2c=140e6, sigma_shear=70e6,
    angle=0, thickness=5
)
laminate = Laminate(plies=[lamina])
```

### Running Analyses
Use the `Panel` and `Sandwich` classes to analyze structures:
```python
from materials.sandwich import Sandwich
from materials.panel import Panel

sandwich = Sandwich(composite_material=laminate, core_material=core, tc=50)
panel = Panel(sandwich=sandwich, width=2.0, length=4.0, point_load=1000)
```

### Visualizing Results
Results can be visualized using the `utils.visualization` module.

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
