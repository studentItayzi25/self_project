import numpy as np
import matplotlib.pyplot as plt

# Function to read data from files
def read_file(filename, skip_header=1):
    return np.loadtxt(filename, skiprows=skip_header, delimiter=None)

# Plot initial and final states
def plot_initial_final(initial_file, final_file, output_file):
    initial_data = read_file(initial_file)
    final_data = read_file(final_file)

    x_init, rho_init, p_init, u_init = initial_data[:, 0], initial_data[:, 1], initial_data[:, 2], initial_data[:, 3]
    x_final, rho_final, p_final, u_final = final_data[:, 0], final_data[:, 1], final_data[:, 2], final_data[:, 3]

    plt.figure(figsize=(12, 8))
    
    # Density plot
    plt.subplot(3, 1, 1)
    plt.plot(x_init, rho_init, label="Initial", linestyle="--")
    plt.plot(x_final, rho_final, label="Final")
    plt.ylabel("Density (rho)")
    plt.legend()
    plt.grid()

    # Velocity plot
    plt.subplot(3, 1, 2)
    plt.plot(x_init, u_init, label="Initial", linestyle="--")
    plt.plot(x_final, u_final, label="Final")
    plt.ylabel("Velocity (u)")
    plt.legend()
    plt.grid()

    # Pressure plot
    plt.subplot(3, 1, 3)
    plt.plot(x_init, p_init, label="Initial", linestyle="--")
    plt.plot(x_final, p_final, label="Final")
    plt.xlabel("x")
    plt.ylabel("Pressure (p)")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

# Plot fluxes
def plot_flux(flux_file, output_file):
    flux_data = read_file(flux_file)
    x, F_conv, F_left, F_right = flux_data[:, 0], flux_data[:, 1], flux_data[:, 2], flux_data[:, 3]

    plt.figure(figsize=(10, 6))
    plt.plot(x, F_conv, label="Convective Flux", color="blue")
    plt.plot(x, F_left, label="Left Flux", linestyle="--", color="green")
    plt.plot(x, F_right, label="Right Flux", linestyle="--", color="red")
    plt.xlabel("x")
    plt.ylabel("Flux")
    plt.title("Fluxes Across the Domain")
    plt.legend()
    plt.grid()
    plt.savefig(output_file)
    plt.show()

# Plot residuals
def plot_residuals(residual_file, output_file):
    residual_data = read_file(residual_file)
    x, r_rho, r_rhou, r_E = residual_data[:, 0], residual_data[:, 1], residual_data[:, 2], residual_data[:, 3]

    plt.figure(figsize=(10, 6))
    plt.plot(x, r_rho, label="Residual of rho", color="blue")
    plt.plot(x, r_rhou, label="Residual of rho*u", color="orange")
    plt.plot(x, r_E, label="Residual of E", color="green")
    plt.xlabel("x")
    plt.ylabel("Residuals")
    plt.title("Residuals Across the Domain")
    plt.legend()
    plt.grid()
    plt.savefig(output_file)
    plt.show()

# Main
if __name__ == "__main__":
    base_path = "/home/itay/CFD_class/self_project/"
    
    plot_initial_final(
        initial_file=base_path + "shock_tube_output.dat",
        final_file=base_path + "final.dat",
        output_file=base_path + "initial_final_plot.png"
    )

    plot_flux(
        flux_file=base_path + "fluxes.dat",
        output_file=base_path + "flux_plot.png"
    )

    plot_residuals(
        residual_file=base_path + "residuals.dat",
        output_file=base_path + "residuals_plot.png"
    )
