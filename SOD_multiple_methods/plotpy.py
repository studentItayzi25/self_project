import numpy as np
import matplotlib.pyplot as plt
import os

# Define the file paths
base_path = r"\\wsl.localhost\Ubuntu\home\itay_22\fortran\self learning\self_project\self_project\SOD_multiple_methods"
conservative_and_flux_file = os.path.join(base_path, "conservative_and_flux.dat")
fluxes_glf_file = os.path.join(base_path, "fluxes_GLF.dat")
fluxes_hll_file = os.path.join(base_path, "fluxes_HLL.dat")
fluxes_hllc_file = os.path.join(base_path, "fluxes_HLLC.dat")
final_file = os.path.join(base_path, "final.dat")

# Function to read data from a file
def read_data(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} not found.")
    try:
        data = np.loadtxt(file_path, skiprows=1)
    except ValueError as e:
        print(f"Error reading {file_path}: {e}")
        data = np.genfromtxt(file_path, skip_header=1, filling_values=np.nan)
    return data

# Read the data
try:
    conservative_and_flux_data = read_data(conservative_and_flux_file)
    fluxes_glf_data = read_data(fluxes_glf_file) if os.path.exists(fluxes_glf_file) else None
    fluxes_hll_data = read_data(fluxes_hll_file) if os.path.exists(fluxes_hll_file) else None
    fluxes_hllc_data = read_data(fluxes_hllc_file) if os.path.exists(fluxes_hllc_file) else None
    final_data = read_data(final_file)
except FileNotFoundError as e:
    print(e)
    exit(1)

# Plot the conservative and flux data
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
plt.plot(conservative_and_flux_data[:, 0], conservative_and_flux_data[:, 1], label='rho')
plt.plot(conservative_and_flux_data[:, 0], conservative_and_flux_data[:, 2], label='rho*u')
plt.plot(conservative_and_flux_data[:, 0], conservative_and_flux_data[:, 3], label='E')
plt.xlabel('x')
plt.ylabel('Conservative Variables')
plt.legend()
plt.title('Conservative Variables vs x')

# Plot the GLF fluxes data if it exists
if fluxes_glf_data is not None:
    plt.subplot(2, 2, 2)
    plt.plot(fluxes_glf_data[:, 0], fluxes_glf_data[:, 1], label='F_glf')
    plt.plot(fluxes_glf_data[:, 0], fluxes_glf_data[:, 2], label='F_left')
    plt.plot(fluxes_glf_data[:, 0], fluxes_glf_data[:, 3], label='F_right')
    plt.xlabel('x')
    plt.ylabel('GLF Fluxes')
    plt.legend()
    plt.title('GLF Fluxes vs x')

# Plot the HLL fluxes data if it exists
if fluxes_hll_data is not None:
    plt.subplot(2, 2, 3)
    plt.plot(fluxes_hll_data[:, 0], fluxes_hll_data[:, 1], label='F_hll')
    plt.plot(fluxes_hll_data[:, 0], fluxes_hll_data[:, 2], label='F_left')
    plt.plot(fluxes_hll_data[:, 0], fluxes_hll_data[:, 3], label='F_right')
    plt.xlabel('x')
    plt.ylabel('HLL Fluxes')
    plt.legend()
    plt.title('HLL Fluxes vs x')

# Plot the HLLC fluxes data if it exists
if fluxes_hllc_data is not None:
    plt.subplot(2, 2, 4)
    plt.plot(fluxes_hllc_data[:, 0], fluxes_hllc_data[:, 1], label='F_hllc')
    plt.plot(fluxes_hllc_data[:, 0], fluxes_hllc_data[:, 2], label='F_left')
    plt.plot(fluxes_hllc_data[:, 0], fluxes_hllc_data[:, 3], label='F_right')
    plt.xlabel('x')
    plt.ylabel('HLLC Fluxes')
    plt.legend()
    plt.title('HLLC Fluxes vs x')

plt.tight_layout()
plt.show()

# Plot the final data
plt.figure(figsize=(12, 6))
plt.plot(final_data[:, 0], final_data[:, 1], label='rho')
plt.plot(final_data[:, 0], final_data[:, 2], label='rho*u')
plt.plot(final_data[:, 0], final_data[:, 3], label='E')
plt.xlabel('x')
plt.ylabel('Final Variables')
plt.legend()
plt.title('Final Variables vs x')
plt.show()