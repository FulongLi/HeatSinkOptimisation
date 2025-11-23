# Heat Sink Optimization

A Python-based optimization tool for plate fin heat sink design using differential evolution algorithms. This project minimizes both thermal resistance and heat sink mass to find optimal geometric parameters.

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

- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/HeatSinkOptimisation.git
cd HeatSinkOptimisation
```

2. Install required dependencies:
```bash
pip install numpy scipy matplotlib
```

## Usage

Run the main optimization script:
```bash
python HS_main.py
```

The script will:
1. Perform 100 optimization runs using differential evolution
2. Generate a Pareto front showing the trade-off between thermal resistance and mass
3. Display optimal heat sink parameters
4. Visualize the optimal heat sink geometry in 3D

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
2. **Heat Transfer Coefficient**: Calculated from Elenbaas number
3. **Fin Efficiency**: Accounts for temperature variation along fins
4. **Thermal Resistance**: Overall heat sink thermal performance
5. **Mass Calculation**: Total heat sink mass based on geometry

For detailed mathematical formulations, see `heat_sink_optimization_summary.html`.

## Optimization Parameters

The optimization searches for optimal values of:
- **L**: Length (0.05 - 0.1 m)
- **W**: Width (0.05 - 0.1 m)
- **H**: Fin height (0.02 - 0.1 m)
- **t**: Fin thickness (0.001 - 0.005 m)
- **b**: Fin spacing (0.002 - 0.01 m)
- **t_bp**: Base plate thickness (0.002 - 0.01 m)

## Author

Created by Dr. Fulong Li

## License

[Specify your license here]
