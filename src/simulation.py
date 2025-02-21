import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def simulate(global_matrix_h, global_matrix_c, vector_p, time_step, simulation_time, init_temp, n, nodes):

    curr_time = 0
    t0 = matrix = np.full((n, 1), init_temp)

    points = np.zeros((n, 2))
    for i,node in enumerate(nodes):
        points[i][0] = node.x
        points[i][1] = node.y

    min_x = np.min(points[:, 0])
    max_x = np.max(points[:, 0])
    min_y = np.min(points[:, 1])
    max_y = np.max(points[:, 1])
    print(points)

    while curr_time < simulation_time:
        c_div_t = global_matrix_c / time_step
        t1 = np.linalg.inv(global_matrix_h + c_div_t) @ (c_div_t @ t0 + vector_p)
        print(np.min(t1), np.max(t1))
        t0 = t1
        curr_time+=time_step

    grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100))
    grid_z = sp.interpolate.griddata(points, t0.flatten(), (grid_x, grid_y), method='cubic')
    plt.figure(figsize=(8, 6))
    plt.contourf(grid_x, grid_y, grid_z, levels=1000, cmap='coolwarm')
    plt.colorbar(label='temperature')
    plt.show()
