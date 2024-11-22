from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np


@dataclass
class Node:
    x: float
    y: float



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

    def get_node_coords(self, grid: 'Grid') -> List[Tuple[float, float]]:
        return [(grid.nodes[node_id - 1].x, grid.nodes[node_id - 1].y) for node_id in self.nodes_ids]

    def initialize_jakobian(self, elem_univ, grid: 'Grid', npc: int):
        element_coords = self.get_node_coords(grid)
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
        #print(self.final_matrix_H)

    def agregate_matrixes_h(self, agregated_matrix_h):
        print(self.nodes_ids)
        agregation_formula = []
        for i in self.nodes_ids:
            for j in self.nodes_ids:
                agregation_formula.append((i,j))

        for i,elem in enumerate(self.final_matrix_H.flat):
            agregated_matrix_h.matrix_h[agregation_formula[i][0]-1, agregation_formula[i][1]-1] += elem

    def __str__(self):
        if self.jakobian:
            jakobian_str = "\n  ".join(str(j) for j in self.jakobian)
        else:
            jakobian_str = "None"
        return f"Element(id={self.id}, nodes_ids={self.nodes_ids}, jakobian=[\n  {jakobian_str}\n])"

@dataclass
class Node:
    x: float
    y: float

@dataclass
class Grid:
    nN: int
    nE: int
    nodes: List[Node] = field(default_factory=list)
    elements: List[Element] = field(default_factory=list)

class ElemUniv:
    def __init__(self, integration_points):
        self.dNdxi = []
        self.dNdeta = []
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