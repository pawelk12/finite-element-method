import os.path
import threading

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.animation import FuncAnimation

def simulate(global_matrix_h, global_matrix_c, vector_p, time_step, simulation_time, init_temp, n, nodes):
    t0 = np.full((n, 1), init_temp)

    points = np.zeros((n, 2))
    for i,node in enumerate(nodes):
        points[i][0] = node.x
        points[i][1] = node.y

    min_x = np.min(points[:, 0])
    max_x = np.max(points[:, 0])
    min_y = np.min(points[:, 1])
    max_y = np.max(points[:, 1])

    fig = plt.figure(figsize=(6, 6))
    fig.patch.set_facecolor('black')

    grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100))

    def update(curr_time):
        nonlocal t0
        c_div_t = global_matrix_c / time_step
        t1 = np.linalg.inv(global_matrix_h + c_div_t) @ (c_div_t @ t0 + vector_p)
        print(np.min(t1), np.max(t1))
        grid_z = sp.interpolate.griddata(points, t0.flatten(), (grid_x, grid_y), method='cubic')
        contour = plt.contourf(grid_x, grid_y, grid_z, levels=1000, cmap='coolwarm')
        t0 = t1
        curr_time += time_step
        return contour

    gif_path = os.path.join('..', 'assets', 'animation.gif')
    animation = FuncAnimation(fig, func=update, frames=np.arange(0, simulation_time, time_step))
    animation.save(gif_path, writer='pillow', fps=5)
    plt.close(fig)


    #final plot after simulation
    grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100))
    grid_z = sp.interpolate.griddata(points, t0.flatten(), (grid_x, grid_y), method='cubic')
    plt.figure(figsize=(8, 6))
    plt.contourf(grid_x, grid_y, grid_z, levels=1000, cmap='coolwarm')
    plt.colorbar(label='temperature')
    plt.show()