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

# Plot final_condition.dat
plot_data('final_condition.dat', 'Final Condition', ['rho', 'u', 'p', 'c'])