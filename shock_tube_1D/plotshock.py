import numpy as np
import matplotlib.pyplot as plt
import os

def read_file(filename, skip_header=1):
    if not os.path.isfile(filename):
        print(f"Error: {filename} not found.")
        return np.array([])

    try:
        data = np.loadtxt(filename, skiprows=skip_header, delimiter=None)
    except ValueError as e:
        print(f"Error reading {filename}: {e}")
        data = []
        with open(filename, 'r') as file:
            for line in file:
                try:
                    row = np.fromstring(line, sep=' ')
                    if len(row) == 4:  # Adjust this based on expected number of columns
                        data.append(row)
                except ValueError:
                    continue
        data = np.array(data)
    return data

def plot_initial_final(initial_file, final_file, output_file):
    initial_data = read_file(initial_file)
    final_data = read_file(final_file)

    if initial_data.size == 0 or final_data.size == 0:
        print("Error: No data to plot.")
        return

    x_init, rho_init, rho_u_init, E_init = initial_data[:, 0], initial_data[:, 1], initial_data[:, 2], initial_data[:, 3]
    x_final, rho_final, rho_u_final, E_final = final_data[:, 0], final_data[:, 1], final_data[:, 2], final_data[:, 3]

    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.plot(x_init, rho_init, '--', label="Initial")
    plt.plot(x_final, rho_final, label="Final")
    plt.ylabel(r"$\rho$ (Density)")
    plt.legend()
    plt.grid()

    plt.subplot(3, 1, 2)
    plt.plot(x_init, rho_u_init, '--', label="Initial")
    plt.plot(x_final, rho_u_final, label="Final")
    plt.ylabel(r"$\rho u$ (Momentum)")
    plt.legend()
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(x_init, E_init, '--', label="Initial")
    plt.plot(x_final, E_final, label="Final")
    plt.xlabel("x")
    plt.ylabel("Energy")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

initial_file = r"\\wsl.localhost\Ubuntu\home\itay_22\fortran\self learning\self_project\self_project\shock_tube_1D\time_step_50.dat"
final_file = r"\\wsl.localhost\Ubuntu\home\itay_22\fortran\self learning\self_project\self_project\shock_tube_1D\time_step_1000.dat"
output_file_initial_final = r"\\wsl.localhost\Ubuntu\home\itay_22\fortran\self learning\self_project\self_project\shock_tube_1D\comparison_plot.png"

plot_initial_final(initial_file, final_file, output_file_initial_final)
