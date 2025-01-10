import numpy as np

def simulate(global_matrix_h, global_matrix_c, vector_p, time_step, simulation_time, init_temp, n):

    curr_time = 0
    t0 = matrix = np.full((n, 1), init_temp)
    while curr_time < simulation_time:
        C_div_T = global_matrix_c / time_step
        t1 = np.linalg.inv(global_matrix_h + C_div_T) @ (C_div_T @ t0 + vector_p)
        print(np.min(t1), np.max(t1))
        t0 = t1
        curr_time+=time_step