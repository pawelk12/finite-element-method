import math
from models import GlobalData, Node, Element, Grid, ElemUniv, Jakobian
from parse_data import read_global_data, read_coords, read_elements

'''
filename = "Test1_4_4.txt"
global_data = read_global_data(filename)
nodes = read_coords(filename)
elements = read_elements(filename)
grid = Grid(global_data.nN, global_data.nE, nodes, elements)


print(global_data)

print("\nWspolrzedne nodow:")
for node in grid.nodes:
    print(node)

for element in elements:
    print(f"element ID: {element.id}, id noda: {element.nodes_ids}")




integration_points = [(-1/math.sqrt(3), -1/math.sqrt(3)), 
                      (1/math.sqrt(3), -1/math.sqrt(3)),
                      (1/math.sqrt(3), 1/math.sqrt(3)), 
                      (-1/math.sqrt(3), 1/math.sqrt(3))]


elem_univ = ElemUniv(integration_points)


grid = Grid(nN=4, nE=1, nodes=nodes, elements=elements)


for element in grid.elements:
    element.initialize_jakobian(elem_univ, grid, npc=len(integration_points))


for element in grid.elements:
    for i, jakobian in enumerate(element.jakobian):
        print(f"Element ID: {element.id}, Punkt całkowania {i}")
        print("J:", jakobian.J)
        print("J1:", jakobian.J1)      
        print("detJ:", jakobian.detJ)  
        print()

'''

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
