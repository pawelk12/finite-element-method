import numpy as np
import os
from models import Grid, ElemUniv, GlobalMatrixH, GlobalMatrixC, GlobalVectorP
from parse_data import read_global_data, read_coords, read_elements, read_bc
from simulation import simulate
from utils import init_integration_points_and_weights
from utils import get_args

def main():
    filename, n_integration_points = get_args()

    filepath = os.path.join('..', 'test_grid_data', filename)

    global_data = read_global_data(filepath)
    nodes = read_coords(filepath)
    bc_list = read_bc(filepath)
    elements = read_elements(filepath)

    for i,node in enumerate(nodes, start=1):
        if i in bc_list:
            node.bc = True

    grid = Grid(global_data.nN, global_data.nE, nodes, elements)


    integration_points, weights = init_integration_points_and_weights(n_integration_points)
    elem_univ = ElemUniv(integration_points, weights, global_data.alfa, global_data.tot)

    global_matrix_h = GlobalMatrixH(global_data.nN)
    global_matrix_c = GlobalMatrixC(global_data.nN)
    global_vector_p = GlobalVectorP(global_data.nN)

    for element in grid.elements:
        element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))
        element.initialize_matrixH(elem_univ, npc=len(integration_points),
                                   conductivity=global_data.conductivity,
                                   weights=weights)

        element.initialize_matrixC(integration_points, npc=len(integration_points),
                                   density=global_data.density,
                                   specific_heat=global_data.specific_heat,
                                   weights=weights)


        element.calculate_hbc_from_template(grid, elem_univ)
        element.add_hbc_matrix_to_h()
        element.agregate_matrixes_h(global_matrix_h)
        element.agregate_matrixes_c(global_matrix_c)
        element.calculate_vector_p_from_template(grid, elem_univ)
        element.agregate_vectors_p(global_vector_p)

    #steady
    #{t}=−1 * [H]^−1 * {P}
    vector_t = -1 * (np.linalg.inv(global_matrix_h.matrix_h) @ global_vector_p.vector_P)

    #unsteady
    simulate(global_matrix_h.matrix_h,
             global_matrix_c.matrix_c,
             global_vector_p.vector_P,
             global_data.simulation_step_time,
             global_data.simulation_time,
             global_data.initial_temp,
             global_data.nN,
             grid.nodes)

if __name__ == "__main__":
    main()