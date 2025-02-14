from typing import List
from models import GlobalData, Node, Element
import sys

def read_global_data(filename: str) -> GlobalData:
    data = []
    try:
        with open(filename, "r") as file:
            for line in file:
                if line.startswith("*Node"):
                    break
                data.append(int(line.split()[-1]))
    except FileNotFoundError:
        print("Failed to open the file.")
        sys.exit(1)

    return GlobalData(*data)

def read_coords(filename: str) -> List[Node]:
    nodes = []
    try:
        with open(filename, "r") as file:
            read = False
            for line in file:
                if line.startswith("*Node"):
                    read = True
                    continue
                if line.startswith("*Element"):
                    break
                if read:
                    parts = line.split(",")
                    x, y = float(parts[1]), float(parts[2].strip())
                    nodes.append(Node(x, y, 0))
    except FileNotFoundError:
        print("Failed to open the file.")
        sys.exit(1)
    return nodes

def read_bc(filename: str) -> List[int]:
    bc_nodes_list = []
    try:
        with open(filename, "r") as file:
            read = False
            for line in file:
                if line.startswith("*BC"):
                    read = True
                    continue
                if read:
                    parts = line.split(",")
                    bc_nodes_list = [int(part.strip()) for part in parts]
    except FileNotFoundError:
        print("Failed to open the file.")
        sys.exit(1)
    return bc_nodes_list

def read_elements(filename: str) -> List[Element]:
    elements = []
    try:
        with open(filename, "r") as file:
            read = False
            for line in file:
                if line.startswith("*Element"):
                    read = True
                    continue
                if line.startswith("*BC"):
                    break
                if read:
                    my_line = line.strip().split(",")
                    element_id = int(my_line[0])
                    ids = tuple(map(int, my_line[1:5]))
                    elements.append(Element(id=element_id, nodes_ids=ids))
    except FileNotFoundError:
        print("Failed to open the file.")
        sys.exit(1)
    return elements
