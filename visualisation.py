import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.animation import FuncAnimation

df = pd.read_csv("displacement.txt", sep=";", header=None, names=["x", "y", "u", "v"])

w = 2420
x = df["x"].to_numpy()
y = df["y"].to_numpy()
u = df["u"].to_numpy()
v = df["v"].to_numpy()

x = np.concatenate((x, -x))
y = np.concatenate((y, y))

u = np.concatenate((u, -u))
v = np.concatenate((v, v))

fig, ax = plt.subplots()

# Tracé de la géométrie initiale
ax.plot(x, y, "bo")
x_limits = np.array([-0.1, 0.1])
y_limits = np.array([0, 0.15])
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.autoscale(False)
ax.set_aspect('equal', 'box')

def update(frame):
    # Convert frequency to radians per second
    w_rad = 2 * np.pi * w

    # Compute coordinates of deformed nodes with a vibrating frequency
    xdef = x + u * np.sin(frame * w_rad) * 1e9
    ydef = y + v * np.sin(frame * w_rad) * 1e9

    # Compute norm of eigenvectors
    norm = np.sqrt(u**2 + v**2)

    # Clear previous plot
    ax.cla()

    # Plot heat map
    vmin, vmax = np.min(norm), np.max(norm)
    ax.scatter(xdef, ydef, c=norm, cmap="coolwarm", vmin=vmin, vmax=vmax)
    #ax.plot(xdef, ydef,"bo")
    # Adjust x and y axis
    
    # Add title and legend
    ax.set_title(f"Eigenvectors at time {frame}")

    # Fix the axis limits
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_aspect('equal', 'box')


# Call update with frame=0 to initialize the plot
update(0)

# Création de l'animation
animation = FuncAnimation(fig, update, frames=range(50), interval=50, repeat=True)
plt.show()


