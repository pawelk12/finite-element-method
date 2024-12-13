from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np
import math


@dataclass
class Node:
    x: float
    y: float
    bc: bool

@dataclass
class GlobalData:
    simulation_time: int
    simulation_step_time: int
    conductivity: int
    alfa: int
    tot: int
    initial_temp: int
    density: int
    specific_heat: int
    nN: int
    nE: int


class Jakobian:
    def __init__(self, elem_univ, element_coords: List[Tuple[float, float]], npc: int):
        self.J: np.ndarray = np.zeros((2, 2))
        self.J1: np.ndarray = np.zeros((2, 2))
        self.detJ: float = 0.0

        dNdxi = elem_univ.dNdxi[npc]
        dNdeta = elem_univ.dNdeta[npc]

        for i in range(len(element_coords)):
            self.J[0, 0] += dNdxi[i] * element_coords[i][0]  # dxdxi
            self.J[0, 1] += dNdxi[i] * element_coords[i][1]  # dydxi
            self.J[1, 0] += dNdeta[i] * element_coords[i][0]  # dxdeta
            self.J[1, 1] += dNdeta[i] * element_coords[i][1]  # dydeta

        self.J1[0, 0] = self.J[1, 1]
        self.J1[0, 1] = -self.J[0, 1]
        self.J1[1, 0] = -self.J[1, 0]
        self.J1[1, 1] = self.J[0, 0]


        self.detJ = np.linalg.det(self.J)

    def __str__(self):
        J_str = np.array2string(self.J, precision=5, separator=', ')
        J1_str = np.array2string(self.J1, precision=5, separator=', ')
        return f"Jakobian(J={J_str}, \nJ1={J1_str},\n detJ={self.detJ})"



class MatrixH:
    def __init__(self, elem_univ, scaling_matrix: np.ndarray, iteration: int, npc: int, conductivity: int, detJ: float):
        self.scaling_matrix: np.darray = scaling_matrix
        self.matrix_H: np.darray = np.zeros((4,4))
        dNdxi = elem_univ.dNdxi[iteration]
        dNdeta = elem_univ.dNdeta[iteration]
        dNidx_mat: np.ndarray = np.zeros((4,1))
        dNidy_mat: np.ndarray = np.zeros((4,1))

        for i in range (0, 4):
            dNidx = scaling_matrix[0 ,0] * dNdxi[i] + scaling_matrix[0 ,1] * dNdeta[i]
            dNidy = scaling_matrix[1 ,0] * dNdxi[i] + scaling_matrix[1 ,1] * dNdeta[i]
            dNidx_mat[i, 0] = dNidx
            dNidy_mat[i, 0] = dNidy
        
        
        self.matrix_H = conductivity * (dNidx_mat @ dNidx_mat.T + dNidy_mat @ dNidy_mat.T) * detJ
        #print(f"Matrix H for integration point {iteration + 1}")
        #print(np.array2string(self.matrix_H, precision=3, separator=', '))

@dataclass
class Element:
    id: int
    nodes_ids: Tuple[int, int, int, int]
    jakobian: List[Jakobian] = None
    matrixes_H: List[MatrixH] = field(default_factory=list)

    def get_nodes_coords(self, grid: 'Grid') -> List[Tuple[float, float]]:
        return [(grid.nodes[node_id - 1].x, grid.nodes[node_id - 1].y) for node_id in self.nodes_ids]

    def get_nodes_boundary_conditions(self, grid: 'Grid') -> List[bool]:
        return [grid.nodes[node_id - 1].bc for node_id in self.nodes_ids]

    def initialize_jakobian(self, elem_univ, grid: 'Grid', npc: int):
        element_coords = self.get_nodes_coords(grid)
        self.jakobian = [Jakobian(elem_univ, element_coords, i) for i in range(npc)]
    
    def initialize_matrixH(self, elem_univ, npc: int, conductivity: int, weights: List[Tuple]):
        for i in range(npc):
            J1: np.ndarray = self.jakobian[i].J1
            detJ: float = self.jakobian[i].detJ
            rev_detJ: float = 1/detJ
            scaling_matrix: np.ndarray = rev_detJ * J1
            self.matrixes_H.append(MatrixH(elem_univ, scaling_matrix, i, npc, conductivity, detJ))


        self.final_matrix_H: np.ndarray = np.zeros((4,4))
        for i in range(npc):
            self.matrixes_H[i].matrix_H
            #print(self.matrixes_H[i].matrix_H)
            self.final_matrix_H += self.matrixes_H[i].matrix_H * weights[i][0] * weights[i][1]
        print(self.final_matrix_H)

    def agregate_matrixes_h(self, agregated_matrix_h):
        print(self.nodes_ids)
        agregation_formula = []
        for i in self.nodes_ids:
            for j in self.nodes_ids:
                agregation_formula.append((i,j))

        for i,elem in enumerate(self.final_matrix_H.flat):
            agregated_matrix_h.matrix_h[agregation_formula[i][0]-1, agregation_formula[i][1]-1] += elem


    def calculate_hbc_from_template(self, grid: 'Grid', elem_univ: 'ElemUniv'):
        nodes_coords = self.get_nodes_coords(grid)
        nodes_bc = self.get_nodes_boundary_conditions(grid)
        self.hbc = np.zeros((4,4))

        for i in range(4):  
            next_i = (i + 1) % 4  
            if nodes_bc[i] == True and nodes_bc[next_i] == True:
                pitagorean_distance = math.sqrt(pow(nodes_coords[next_i][0] - nodes_coords[i][0], 2) + 
                                        pow(nodes_coords[next_i][1] - nodes_coords[i][1], 2))
                detJ = pitagorean_distance / 2
                self.hbc += elem_univ.hbc_templates[i] * detJ

        print(f"Hbc - element {self.id}")
        print(self.hbc)
        print("----------------")

    def add_hbc_matrix_to_h(self):
        print("Matrix H")
        print(self.final_matrix_H)
        print("Matrix hbc")
        print(self.hbc)
        print("------")
        self.final_matrix_H = self.final_matrix_H + self.hbc

    def calculate_vector_p_from_template(self, grid: 'Grid', elem_univ: 'ElemUniv'):
        nodes_coords = self.get_nodes_coords(grid)
        nodes_bc = self.get_nodes_boundary_conditions(grid)
        self.p = np.zeros((4,1))

        for i in range(4):  
            next_i = (i + 1) % 4  
            if nodes_bc[i] == True and nodes_bc[next_i] == True:
                pitagorean_distance = math.sqrt(pow(nodes_coords[next_i][0] - nodes_coords[i][0], 2) + 
                                        pow(nodes_coords[next_i][1] - nodes_coords[i][1], 2))
                detJ = pitagorean_distance / 2
                self.p += elem_univ.vector_p_templates[i] * detJ

        print(f"P - element {self.id}")
        print(self.p)
        print("----------------")
        
    def agregate_vectors_p(self, agregated_vector_p):
        print(self.nodes_ids)
        agregation_formula = []
        for i in self.nodes_ids:
            agregation_formula.append(i)

        for i,elem in enumerate(self.p.flat):
            agregated_vector_p.vector_P[agregation_formula[i]-1] += elem


    def __str__(self):
        if self.jakobian:
            jakobian_str = "\n  ".join(str(j) for j in self.jakobian)
        else:
            jakobian_str = "None"
        return f"Element(id={self.id}, nodes_ids={self.nodes_ids}, jakobian=[\n  {jakobian_str}\n])"


@dataclass
class Grid:
    nN: int
    nE: int
    nodes: List[Node] = field(default_factory=list)
    elements: List[Element] = field(default_factory=list)

class ElemUniv:
    def __init__(self, integration_points, weights, alfa, tot):
        self.dNdxi = []
        self.dNdeta = []
        npc = math.sqrt(len(integration_points))
        # bottom, right, top, left
        self.surfaces = [np.zeros((int(npc), 4)), np.zeros((int(npc), 4)),
            np.zeros((int(npc), 4)), np.zeros((int(npc), 4))]

        self.hbc_templates = [np.zeros((4, 4)), np.zeros((4, 4)),
            np.zeros((4, 4)), np.zeros((4, 4))]

        self.vector_p_templates = [np.zeros((4, 1)), np.zeros((4, 1)),
            np.zeros((4, 1)), np.zeros((4, 1))]
        
        for point in integration_points:
            dN1dxi = -0.25 * (1 - point[1])
            dN2dxi = 0.25 * (1 - point[1])
            dN3dxi = 0.25 * (1 + point[1])
            dN4dxi = -0.25 * (1 + point[1])
            self.dNdxi.append((dN1dxi, dN2dxi, dN3dxi, dN4dxi))

            dN1deta = -0.25 * (1 - point[0])
            dN2deta = -0.25 * (1 + point[0])
            dN3deta = 0.25 * (1 + point[0])
            dN4deta = 0.25 * (1 - point[0])
            self.dNdeta.append((dN1deta, dN2deta, dN3deta, dN4deta))



        #version for four-point Gaussian quadrature:
        #points_bottom_tuples = integration_points[::4]
        points_bottom_tuples = integration_points[::int(npc)]
        points_bottom = [] # integration poins casted into bottom edge of element
        points_top =[]

        for point in points_bottom_tuples:
            point = list(point)
            point[1] = -1
            points_bottom.append(point.copy())
            point[1] = 1
            points_top.append(point.copy())

        points_right_tuples = integration_points[:int(npc):]
        #points_right_tuples = integration_points[:4:]
        points_right = []
        points_left = []

        for point in points_right_tuples:
            point = list(point)
            point[0] = -1
            points_left.append(point.copy())
            point[0] = 1
            points_right.append(point.copy())
        

        for i,point in enumerate(points_bottom):
            bottom_surface = self.surfaces[0]
            bottom_surface[i][0] = 0.25*(1 - point[0])*(1 - point[1])
            bottom_surface[i][1] = 0.25*(1 + point[0])*(1 - point[1])
            bottom_surface[i][2] = 0.25*(1 + point[0])*(1 + point[1])
            bottom_surface[i][3] = 0.25*(1 - point[0])*(1 + point[1])

        for i,point in enumerate(points_right):
            bottom_surface = self.surfaces[1]
            bottom_surface[i][0] = 0.25*(1 - point[0])*(1 - point[1])
            bottom_surface[i][1] = 0.25*(1 + point[0])*(1 - point[1])
            bottom_surface[i][2] = 0.25*(1 + point[0])*(1 + point[1])
            bottom_surface[i][3] = 0.25*(1 - point[0])*(1 + point[1])


        for i,point in enumerate(points_top):
            bottom_surface = self.surfaces[2]
            bottom_surface[i][0] = 0.25*(1 - point[0])*(1 - point[1])
            bottom_surface[i][1] = 0.25*(1 + point[0])*(1 - point[1])
            bottom_surface[i][2] = 0.25*(1 + point[0])*(1 + point[1])
            bottom_surface[i][3] = 0.25*(1 - point[0])*(1 + point[1])

        for i,point in enumerate(points_left):
            bottom_surface = self.surfaces[3]
            bottom_surface[i][0] = 0.25*(1 - point[0])*(1 - point[1])
            bottom_surface[i][1] = 0.25*(1 + point[0])*(1 - point[1])
            bottom_surface[i][2] = 0.25*(1 + point[0])*(1 + point[1])
            bottom_surface[i][3] = 0.25*(1 - point[0])*(1 + point[1])


        for x in self.surfaces:
            print(np.array2string(x))

        for j,surface in enumerate(self.surfaces):
            for i in range(int(npc)):
                row = np.array(surface[i]).reshape(1, -1)
                multiplied_matrix = row.T @ row
                multiplied_matrix_times_weight = multiplied_matrix * weights[i][1] * alfa
                #print(multiplied_matrix_times_weight)
                self.hbc_templates[j] += multiplied_matrix_times_weight
                multiplied_vector_p = row.T * tot * alfa *weights[i][1]
                self.vector_p_templates[j] += multiplied_vector_p
            


        print("hbc matrixes for every element without jacobi det")
        print("bottom - right - top- left")
        for i in range(4):
            print(self.hbc_templates[i])
        
                


    def __str__(self):
        return f"dNdxi: {self.dNdxi}\ndNdeta: {self.dNdeta}"
    


@dataclass
class GlobalMatrixH:
    nN: int
    matrix_h: np.ndarray = field(init=False)
    def __post_init__(self):
        self.matrix_h = np.zeros((self.nN, self.nN))

    def print_matrix(self):
        print(np.array2string(self.matrix_h, precision=5, separator=', ', max_line_width=200))


@dataclass
class GlobalVectorP:
    nN: int
    vector_P: np.ndarray = field(init=False)
    def __post_init__(self):
        self.vector_P = np.zeros((self.nN, 1))

    def print_vector(self):
        print(np.array2string(self.vector_P, precision=5, separator=', ', max_line_width=200))