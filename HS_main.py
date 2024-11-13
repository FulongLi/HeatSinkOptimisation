import numpy as np
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants for air and material properties
RHO_AIR = 1.095  # Density of air (kg/m^3)
BETA = 1 / 300  # Thermal expansion coefficient for air (approximation at 300 K) (1/K)
G = 9.81  # Acceleration due to gravity (m/s^2)
C_P_AIR = 1007  # Specific heat capacity of air (J/kg·K)
MU_AIR = 1.95e-5  # Absolute viscosity of air (kg/s·m)
K_F_AIR = 0.0274  # Thermal conductivity of air (W/m·K)
K_ALUMINUM = 196  # Thermal conductivity of aluminum (W/m·K)
RHO_ALUMINUM = 2700  # Density of aluminum (kg/m^3)
K_COPPER = 380  # Thermal conductivity of copper (W/m·K)
RHO_COPPER = 8960*0.97  # Density of copper (kg/m^3)

# Step 1: Calculate Elenbaas number
def elenbaas_number(fin_spacing, fin_length, delta_temp):
    return (RHO_AIR ** 2 * BETA * G * C_P_AIR * fin_spacing ** 4 * delta_temp) / (MU_AIR * K_F_AIR * fin_length)

# Step 2: Calculate heat transfer coefficient using the updated formula
def heat_transfer_coefficient(elenbaas, fin_spacing, k_f=K_F_AIR):
    return (k_f / fin_spacing) * ((576 / elenbaas ** 2) + (2.873 / np.sqrt(elenbaas))) ** (-1 / 2)

# Step 3: Calculate fin efficiency (η)
def fin_efficiency(h, fin_height, fin_thickness, k=K_COPPER):
    m = np.sqrt(2 * h / (k * fin_thickness))
    return np.tanh(m * fin_height) / (m * fin_height)

# Step 4: Calculate the total surface area based on the geometry definition from the image
def total_surface_area(W, L, H, t, b, eta):
    n = int((W - t) / (t + b)) + 1  # Number of fins
    base_area = W * L - n * L * t  # Base area adjusted for fins
    fin_area = 2 * H * L  # Fin area considering both sides of each fin
    return base_area + n * fin_area * eta

# Step 5: Calculate thermal resistance (Rha)
def thermal_resistance(L, W, H, t, b, t_bp, delta_temp):
    # Step 1: Calculate Elenbaas number
    el = elenbaas_number(b, L, delta_temp)

    # Step 2: Calculate heat transfer coefficient
    h = heat_transfer_coefficient(el, b, k_f=K_F_AIR)

    # Step 3: Calculate fin efficiency
    eta = fin_efficiency(h, H, t, k=K_COPPER)

    # Step 4: Calculate total surface area
    total_area = total_surface_area(W, L, H, t, b, eta)

    # Step 5: Calculate thermal resistance
    return 1 / (h * total_area)

# Step 6: Calculate heat sink mass (m_h)
def heat_sink_mass(W, L, H, t, t_bp, b, rho=RHO_COPPER):
    n = int((W - t) / (t + b)) + 1  # Number of fins
    return rho * (W * L * t_bp + n * L * H * t)

# Define the fitness function for optimization
def fitness(params):
    L, W, H, t, b, t_bp = params
    Rha = thermal_resistance(L, W, H, t, b, t_bp, delta_temp=40)
    mh = heat_sink_mass(W, L, H, t, t_bp, b)
    # Normalize Rha and mh by dividing by their respective maximum values to scale them
    Rha_norm = Rha # Assume a reasonable maximum value for normalization
    mh_norm = mh    # Assume a reasonable maximum value for normalization
    # 50% weight on minimizing mass and 50% on minimizing thermal resistance
    return 0.5 * Rha_norm + 0.5 * mh_norm

# Define bounds for each parameter based on Table 7.1
bounds = [
    (0.05, 0.1),    # Length (L) in meters
    (0.05, 0.1),    # Width (W) in meters
    (0.02, 0.1),    # Height of fins (H) in meters
    (0.001, 0.005), # Thickness of fins (t) in meters
    (0.002, 0.01),  # Spacing between fins (b) in meters
    (0.002, 0.01)   # Base plate thickness (t_bp) in meters
]

# Run the differential evolution optimization
results = []
for _ in range(100):  # Running multiple times to create a Pareto front
    result = differential_evolution(fitness, bounds, strategy='best1bin', maxiter=1000, popsize=15)
    Rha = thermal_resistance(*result.x, delta_temp=40)
    mh = heat_sink_mass(result.x[1], result.x[0], result.x[2], result.x[3], result.x[5], result.x[4])
    results.append((Rha, mh))

# Extract results for Pareto front
Rha_values, mh_values = zip(*results)

# Create a combined figure for Pareto front and 3D heat sink geometry
fig = plt.figure(figsize=(14, 8))

# Plot Pareto front
ax1 = fig.add_subplot(121)
ax1.scatter(mh_values, Rha_values, c='b', marker='o')
ax1.set_xlabel('Heat Sink Mass (kg)')
ax1.set_ylabel('Thermal Resistance (K/W)')
ax1.set_title('Pareto Front for Heat Sink Optimization')
ax1.grid(True)

# Extract optimal parameters and calculate final values
optimal_params = result.x
optimal_L, optimal_W, optimal_H, optimal_t, optimal_b, optimal_t_bp = optimal_params
optimal_Rha = thermal_resistance(optimal_L, optimal_W, optimal_H, optimal_t, optimal_b, optimal_t_bp, delta_temp=40)
optimal_mh = heat_sink_mass(optimal_W, optimal_L, optimal_H, optimal_t, optimal_t_bp, optimal_b)

print("Optimal Heat Sink Parameters:")
print(f"Length (L): {optimal_L:.4f} m")
print(f"Width (W): {optimal_W:.4f} m")
print(f"Height of fins (H): {optimal_H:.4f} m")
print(f"Thickness of fins (t): {optimal_t:.4f} m")
print(f"Spacing between fins (b): {optimal_b:.4f} m")
print(f"Base plate thickness (t_bp): {optimal_t_bp:.4f} m")
print(f"Optimal Thermal Resistance (Rha): {optimal_Rha:.4f} K/W")
print(f"Optimal Heat Sink Mass (m_h): {optimal_mh:.4f} kg")

# Plot the heat sink shape in 3D based on optimal parameters
ax2 = fig.add_subplot(122, projection='3d')
n = int((optimal_W - optimal_t) / (optimal_t + optimal_b)) + 1  # Number of fins
fin_positions = [i * (optimal_t + optimal_b) for i in range(n)]

# Plot base plate
ax2.bar3d(0, 0, 0, optimal_W, optimal_L, optimal_t_bp, color='gray', alpha=0.6, label='Base Plate')

# Plot fins
for pos in fin_positions:
    ax2.bar3d(pos, 0, optimal_t_bp, optimal_t, optimal_L, optimal_H, color='blue', alpha=0.7)

ax2.set_xlabel('Width (m)')
ax2.set_ylabel('Length (m)')
ax2.set_zlabel('Height (m)')
ax2.set_title('3D Heat Sink Geometry')
plt.legend()

plt.tight_layout()
plt.show()
