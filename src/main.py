import math
from models import Grid, ElemUniv, GlobalMatrixH, GlobalVectorP
from parse_data import read_global_data, read_coords, read_elements, read_bc

'''

### Test1
# Test element, 2 points integration shcema with 2D Gauss
conductivity = 30
nodes = [Node(x=0.0, y=0.0), Node(x=0.025, y=0.0), Node(x=0.025, y=0.025), Node(x=0.0, y=0.025)]
element = Element(id=1, nodes_ids=(1, 2, 3, 4))

grid = Grid(nN=4, nE=1, nodes=nodes, elements=element)

integration_points = [(-1/math.sqrt(3), -1/math.sqrt(3)), 
                      (1/math.sqrt(3), -1/math.sqrt(3)),
                      (-1/math.sqrt(3), 1/math.sqrt(3)), 
                      (1/math.sqrt(3), 1/math.sqrt(3))]
#
#   npc3  npc4  w2
#
#   npc1  npc2  w1
#    w1    w2
#
weights_for_integration_points = [(1, 1), (1, 1), (1, 1), (1, 1)]
elem_univ = ElemUniv(integration_points)

element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))
element.initialize_matrixH(elem_univ, npc=len(integration_points), conductivity=conductivity, weights=weights_for_integration_points) 


print("Matrix H for element")
print(element.final_matrix_H)


'''
##############################################
### Test2 3 integration points 2D Gauss schema
'''
conductivity = 30
nodes = [Node(x=0.0, y=0.0), Node(x=0.025, y=0.0), Node(x=0.025, y=0.025), Node(x=0.0, y=0.025)]
element = Element(id=1, nodes_ids=(1, 2, 3, 4))

grid = Grid(nN=4, nE=1, nodes=nodes, elements=element)


#
#  5/9 npc3  npc6  npc9
#
#  8/9 npc2  npc5  npc8
#
#  5/9 npc1  npc4  npc7
#       5/9   8/9   5/9


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


elem_univ = ElemUniv(integration_points)
element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))
#print(element) #printing jacobians
element.initialize_matrixH(elem_univ, npc=len(integration_points), conductivity=conductivity, weights=weights_for_integration_points)
print('Matrix H for element')
print(element.final_matrix_H)
'''

##################################
# Test 3, 2 points of integration for non-square element shape
'''
conductivity = 30
nodes = [Node(x=0.01, y=-0.01), Node(x=0.025, y=0.0), Node(x=0.025, y=0.025), Node(x=0.0, y=0.025)]
element = Element(id=1, nodes_ids=(1, 2, 3, 4))

grid = Grid(nN=4, nE=1, nodes=nodes, elements=element)
integration_points = [(-1/math.sqrt(3), -1/math.sqrt(3)), 
                      (1/math.sqrt(3), -1/math.sqrt(3)), 
                      (1/math.sqrt(3), 1/math.sqrt(3)),
                      (-1/math.sqrt(3), 1/math.sqrt(3))]
#
#   npc4  npc3  w2
#
#   npc1  npc2  w1
#    w1    w2
#
weights_for_integration_points = [(1, 1), (1, 1), (1, 1), (1, 1)]
elem_univ = ElemUniv(integration_points)

element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))

element.initialize_matrixH(elem_univ, npc=len(integration_points), conductivity=conductivity, weights=weights_for_integration_points)
print("Matrix H for element")
print(element.final_matrix_H)
'''


# filepath = "../test_grid_data/Test1_4_4.txt"
filepath = "../test_grid_data/Test2_4_4_MixGrid.txt"
global_data = read_global_data(filepath)
nodes = read_coords(filepath)
bc_list = read_bc(filepath)
print(bc_list)

for i,node in enumerate(nodes, start=1):
    if i in bc_list:
        node.bc = True

elements = read_elements(filepath)
grid = Grid(global_data.nN, global_data.nE, nodes, elements)


print(global_data)

print("\nNode coords:")
for node in grid.nodes:
    print(node)

for element in elements:
    print(f"element ID: {element.id}, nodes ID: {element.nodes_ids}")




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
global_vector_p = GlobalVectorP(global_data.nN)

for element in grid.elements:
    element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))
    element.initialize_matrixH(elem_univ, npc=len(integration_points),
                               conductivity=global_data.conductivity,
                               weights=weights_for_integration_points)
    
    #element.agregate_matrixes_h(global_matrix_h)
    #element.calculate_hbc_from_template(grid, elem_univ)

    element.calculate_hbc_from_template(grid, elem_univ)
    element.add_hbc_matrix_to_h()
    element.agregate_matrixes_h(global_matrix_h)
    element.calculate_vector_p_from_template(grid, elem_univ)
    element.agregate_vectors_p(global_vector_p)


global_matrix_h.print_matrix()
global_vector_p.print_vector()