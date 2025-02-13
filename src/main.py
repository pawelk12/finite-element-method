import numpy as np
from models import Grid, ElemUniv, GlobalMatrixH, GlobalMatrixC, GlobalVectorP
from parse_data import read_global_data, read_coords, read_elements, read_bc
from simulation import simulate
from utils import get_args


filename, n_integration_points = get_args()

print(filename)
print(n_integration_points)


#filepath = "../test_grid_data/Test1_4_4.txt"
#filepath = "../test_grid_data/Test2_4_4_MixGrid.txt"
filepath = "../test_grid_data/Test3_31_31_kwadrat.txt"


global_data = read_global_data(filepath)
nodes = read_coords(filepath)
bc_list = read_bc(filepath)

for i,node in enumerate(nodes, start=1):
    if i in bc_list:
        node.bc = True

elements = read_elements(filepath)
grid = Grid(global_data.nN, global_data.nE, nodes, elements)

#two-point Gaussian quadrature
'''
integration_points = [(-1/math.sqrt(3), -1/math.sqrt(3)), 
                      (1/math.sqrt(3), -1/math.sqrt(3)),
                      (1/math.sqrt(3), 1/math.sqrt(3)), 
                      (-1/math.sqrt(3), 1/math.sqrt(3))]
weights_for_integration_points = [(1, 1), (1, 1), (1, 1), (1, 1)]

'''

#three-point Gaussian quadrature
'''
integration_points = [(-math.sqrt(3)/math.sqrt(5), -math.sqrt(3)/math.sqrt(5)), 
                      (-math.sqrt(3)/math.sqrt(5), 0),
                      (-math.sqrt(3)/math.sqrt(5), math.sqrt(3)/math.sqrt(5)), 
                      (0, -math.sqrt(3)/math.sqrt(5)),
                      (0,0),
                      (0, math.sqrt(3)/math.sqrt(5)),
                      (math.sqrt(3)/math.sqrt(5), -math.sqrt(3)/math.sqrt(5)),
                      (math.sqrt(3)/math.sqrt(5), 0),
                      (math.sqrt(3)/math.sqrt(5), math.sqrt(3)/math.sqrt(5))]

weights_for_integration_points = [(5/9, 5/9), (5/9, 8/9), (5/9, 5/9),
                                  (8/9, 5/9), (8/9, 8/9), (8/9, 5/9),
                                  (5/9, 5/9), (5/9, 8/9), (5/9, 5/9)]
'''
                                  

#four-point Gaussian quadrature

#   pc4   pc8   pc12  pc16

#   pc3   pc7   pc11  pc15

#   pc2   pc6   pc10  pc14
 
#   pc1   pc5   pc9   pc13
     

integration_points = [(-0.861136, -0.861136),
                      (-0.861136, -0.339981),
                      (-0.861136, 0.339981),
                      (-0.861136, 0.861136),
                      (-0.339981, -0.861136),
                      (-0.339981, -0.339981),
                      (-0.339981, 0.339981),
                      (-0.339981, 0.861136),
                      (0.339981, -0.861136),
                      (0.339981, -0.339981),
                      (0.339981, 0.339981),
                      (0.339981, 0.861136),
                      (0.861136, -0.861136),
                      (0.861136, -0.339981),
                      (0.861136,  0.339981),
                      (0.861136, 0.861136)]

weights_for_integration_points = [(0.347855, 0.347855), (0.347855, 0.652145), (0.347855, 0.652145), (0.347855, 0.347855),
                                  (0.652145, 0.347855), (0.652145, 0.652145), (0.652145, 0.652145), (0.652145, 0.347855),
                                  (0.652145, 0.347855), (0.652145, 0.652145), (0.652145, 0.652145), (0.652145, 0.347855),
                                  (0.347855, 0.347855), (0.347855, 0.652145), (0.347855, 0.652145), (0.347855, 0.347855)]



elem_univ = ElemUniv(integration_points, weights_for_integration_points, global_data.alfa, global_data.tot)

global_matrix_h = GlobalMatrixH(global_data.nN)
global_matrix_c = GlobalMatrixC(global_data.nN)
global_vector_p = GlobalVectorP(global_data.nN)

for element in grid.elements:
    element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))
    element.initialize_matrixH(elem_univ, npc=len(integration_points),
                               conductivity=global_data.conductivity,
                               weights=weights_for_integration_points)

    element.initialize_matrixC(integration_points, npc=len(integration_points),
                               density=global_data.density,
                               specific_heat=global_data.specific_heat,
                               weights=weights_for_integration_points)


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
         global_data.nN)
