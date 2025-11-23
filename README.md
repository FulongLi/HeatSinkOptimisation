# Heat Sink Optimization

![license - CC BY 4.0](https://img.shields.io/badge/license-CC--BY-green)
![type - hardware](https://img.shields.io/badge/type-hardware-blue)
![category - power electronics](https://img.shields.io/badge/category-power%20electronics-lightgrey)
![status - archived](https://img.shields.io/badge/status-archived-red)

A Python-based optimization tool for plate fin heat sink design using differential evolution algorithms. This project minimizes both thermal resistance and heat sink mass to find optimal geometric parameters.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Mathematical Models](#mathematical-models)
- [Optimization Parameters](#optimization-parameters)
- [Author](#author)
- [Citation](#citation)
- [License](#license)
- [Project Status](#project-status)

## Overview

This project implements a multi-objective optimization approach for plate fin heat sinks, balancing thermal performance (minimizing thermal resistance) and material efficiency (minimizing mass). The optimization uses differential evolution to explore the design space and generate a Pareto front of optimal solutions.

![Heat Sink Geometry](heat_sink_geometry.png)

## Features

- **Thermal Analysis**: Calculates thermal resistance using Elenbaas number, heat transfer coefficients, and fin efficiency models
- **Multi-objective Optimization**: Minimizes both thermal resistance and heat sink mass simultaneously
- **Differential Evolution**: Uses scipy's differential evolution algorithm for robust global optimization
- **Visualization**: Generates Pareto front plots and 3D heat sink geometry visualizations
- **Material Support**: Configurable material properties (copper, aluminum, air)

## Requirements

### Software
- Python 3.6 or higher

### Python Packages
- **NumPy** (>= 1.19.0): Numerical computations
- **SciPy** (>= 1.5.0): Optimization algorithms (differential evolution)
- **Matplotlib** (>= 3.3.0): Plotting and visualization

## Getting Started

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/HeatSinkOptimisation.git
cd HeatSinkOptimisation
```

2. Install required dependencies:
```bash
pip install numpy scipy matplotlib
```

Or use a requirements file (if available):
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

Run the main optimization script:
```bash
python HS_main.py
```

### What the Script Does

The script performs the following steps:

1. **Optimization**: Runs 100 optimization iterations using differential evolution
   - Population size: 15
   - Maximum iterations: 1000 per run
   - Strategy: 'best1bin'

2. **Results Analysis**: 
   - Generates a Pareto front showing the trade-off between thermal resistance and mass
   - Identifies the optimal solution from the final run

3. **Output**:
   - Prints optimal geometric parameters to console
   - Displays thermal resistance and mass values
   - Generates two visualizations:
     - Pareto front scatter plot (thermal resistance vs. mass)
     - 3D visualization of the optimal heat sink geometry

### Example Output

```
Optimal Heat Sink Parameters:
Length (L): 0.XXXX m
Width (W): 0.XXXX m
Height of fins (H): 0.XXXX m
Thickness of fins (t): 0.XXXX m
Spacing between fins (b): 0.XXXX m
Base plate thickness (t_bp): 0.XXXX m
Optimal Thermal Resistance (Rha): X.XXXX K/W
Optimal Heat Sink Mass (m_h): X.XXXX kg
```

## Project Structure

```
HeatSinkOptimisation/
├── HS_main.py                          # Main optimization script
├── heat_sink_optimization_summary.html # Detailed mathematical documentation
├── heat_sink_geometry.png             # Reference geometry image
└── README.md                           # This file
```

## Mathematical Models

The optimization uses several key thermal models:

1. **Elenbaas Number**: Characterizes natural convection heat transfer
   - Relates buoyancy forces to viscous forces in natural convection
   - Used to determine the heat transfer regime

2. **Heat Transfer Coefficient**: Calculated from Elenbaas number
   - Determines the convective heat transfer rate
   - Accounts for natural convection between fins

3. **Fin Efficiency**: Accounts for temperature variation along fins
   - Measures the effectiveness of fins in transferring heat
   - Depends on fin geometry and material properties

4. **Thermal Resistance**: Overall heat sink thermal performance
   - Inverse of the product of heat transfer coefficient and effective surface area
   - Lower values indicate better thermal performance

5. **Mass Calculation**: Total heat sink mass based on geometry
   - Includes base plate and all fins
   - Uses material density (copper by default)

For detailed mathematical formulations and derivations, see [`heat_sink_optimization_summary.html`](heat_sink_optimization_summary.html).

## Optimization Parameters

The optimization searches for optimal values of six geometric parameters:

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Length | L | 0.05 - 0.1 m | Heat sink length in the flow direction |
| Width | W | 0.05 - 0.1 m | Heat sink width (perpendicular to flow) |
| Fin Height | H | 0.02 - 0.1 m | Height of individual fins |
| Fin Thickness | t | 0.001 - 0.005 m | Thickness of each fin |
| Fin Spacing | b | 0.002 - 0.01 m | Gap between adjacent fins |
| Base Plate Thickness | t_bp | 0.002 - 0.01 m | Thickness of the base plate |

The optimization uses a weighted objective function (50% thermal resistance, 50% mass) to find the best compromise between thermal performance and material usage.

## Author

**Dr. Fulong Li**

This project was created for heat sink optimization research and development in power electronics applications.

## Citation

If you use this work in your research, please cite:

```
Heat Sink Optimization - A Python-based optimization tool for plate fin heat sink design
Created by Dr. Fulong Li
```

## License

This project is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

For more details, see the [LICENSE](LICENSE) file or visit [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/).

## Project Status

⚠️ **Note**: This project is currently archived and is no longer actively maintained.
