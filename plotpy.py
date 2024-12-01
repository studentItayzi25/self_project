import numpy as np
import matplotlib.pyplot as plt

def plot_data(filename, title, labels):
    data = np.genfromtxt(filename, skip_header=1, filling_values=np.nan)
    x = np.arange(data.shape[0])
    
    plt.figure()
    for i in range(min(data.shape[1], len(labels))):
        plt.plot(x, data[:, i], label=labels[i])
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{filename.split('.')[0]}.png")
    plt.show()

# Plot grid.dat
plot_data('grid.dat', 'Grid Data', ['x center', 'x face'])

# Plot initial_condition.dat
plot_data('initial_condition.dat', 'Initial Condition', ['rho', 'u', 'p', 'c'])

# Plot left_right_states.dat
plot_data('left_right_states.dat', 'Left and Right States', ['rho_left', 'u_left', 'p_left', 'c_left', 'rho_right', 'u_right', 'p_right', 'c_right'])

# Plot U.dat
plot_data('U.dat', 'U Data', ['rho', 'u', 'p'])

# Plot F.dat
plot_data('F.dat', 'F Data', ['rho left', 'u left', 'E left', 'rho right', 'u right', 'E right'])

# Plot F_flux.dat
plot_data('F_flux.dat', 'F Flux Data', ['rho flux', 'u flux', 'E flux'])

# Plot residuals.dat
plot_data('residuals.dat', 'Residuals Data', ['rho residual', 'u residual', 'E residual'])

# Plot final_condition.dat
plot_data('final_condition.dat', 'Final Condition', ['rho', 'u', 'p', 'c'])