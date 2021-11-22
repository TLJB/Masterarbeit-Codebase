import sys
import re
import os
import getopt
from typing import TextIO
import numpy as np
import math


class Node:
    def __init__(
        self,
        id: int = 0,
        x_coord: float = 0.0,
        y_coord: float = 0.0,
        z_coord: float = 0.0,
    ):
        self.id = id
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord


class Element:
    def __init__(self):
        id = -1
        type = -1
        n_tags = 2
        tags = [-1 for i in range(n_tags)]


class Connectivity:
    def __init__(self, elements: int = -1):
        self.n_elem = elements
        self.n_inter = -1
        self.n_bulk = -1
        self.list = [Element(i) for i in range(self.n_elem)]


class PhysicalName:
    def __init__(self, geo_type: int = -1, number: int = -1, name: str = ""):
        self.dimension = geo_type
        self.tag = number
        self.name = name


def check_format(fobj: TextIO) -> bool:
    # move forward until section $MeshFormat

    Section = "MeshFormat"
    BeginSection = "$" + Section
    EndSection = "$End" + Section
    stop = False
    while not stop:
        line = fobj.readline()
        if (BeginSection) in line:
            stop = True
            break
    stop = False
    while not stop:
        line = fobj.readline()
        if EndSection in line:
            return False
        line_split = line.split()
        if line_split[0] == "2.2":
            return True
        if not line:
            return False


def read_physical_names(fobj: TextIO) -> list[PhysicalName]:
    Section = "PhysicalNames"
    BeginSection = "$" + Section
    EndSection = "$End" + Section
    list_of_names = []
    stop = False
    while not stop:
        line = fobj.readline()
        if BeginSection in line:
            stop = True
            break
    stop = False
    while not stop:
        line = fobj.readline()
        if EndSection in line:
            stop = True
            break
        line_split = line.split()
        if len(line_split) == 1:
            n_names = int(line_split[0])
        if len(line_split) == 3:
            list_of_names.append(
                PhysicalName(
                    int(line_split[0]), int(line_split[1]), line_split[2]
                )
            )
    return list_of_names


def read_nodes(fobj: TextIO) -> list[Node]:

    Section = "Nodes"
    BeginSection = "$" + Section
    EndSection = "$End" + Section
    stop = False
    while not stop:
        line = fobj.readline()
        if BeginSection in line:
            stop = True
            break
    stop = False
    while not stop:
        line = fobj.readline()
        if EndSection in line:
            stop = True
            break
        line_split = line.split()
        if len(line_split) == 1:
            n_nodes = int(line_split[0])
            nodes = [Node() for i in range(n_nodes)]
            iter = 0
        if len(line_split) == 4:
            nodes[iter].id = int(float(line_split[0]))
            nodes[iter].x_coord = float(line_split[1])
            nodes[iter].y_coord = float(line_split[2])
            nodes[iter].z_coord = float(line_split[3])
            iter += 1
    return nodes


def read_elements(fobj: TextIO, names: list[PhysicalName]) -> list[Element]:

    n_inter_elem = -1

    Section = "Elements"
    BeginSection = "$" + Section
    EndSection = "$End" + Section
    stop = False
    while not stop:
        line = fobj.readline()
        if BeginSection in line:
            stop = True
            break
    stop = False
    while not stop:
        line = fobj.readline()
        if EndSection in line:
            stop = True
            break
        line_split = line.split()
        if len(line_split) == 1:
            n_elements = int(line_split[0])
            elements = [Element() for i in range(n_elements)]
            iter = 0
        if not len(line_split) == 1:
            elements[iter].id = int(line_split[0])
            elements[iter].type = int(line_split[1])
            elements[iter].n_tags = int(line_split[2])
            elements[iter].tags = [
                int(_m)
                for _m in line_split[2 + 1 : 2 + int(line_split[2]) + 1]
            ]
            elements[iter].nodes = [
                int(_m) for _m in line_split[2 + int(line_split[2]) + 1 :]
            ]
            iter += 1

    return elements


def volume(nodal_values: np.ndarray) -> float:
    if nodal_values.shape[0] == 4:
        # 2d quadrilateral
        V = 0
        qp = (
            math.sqrt(3)
            / 3
            * np.array(
                [
                    [-1, -1],
                    [1, -1],
                    [-1, 1],
                    [1, 1],
                ]
            )
        )
        w8 = [1 for i in range(qp.shape[0])]
        for i, xi in enumerate(qp):
            gamma = (
                1
                / 4
                * np.array(
                    [
                        [-(1.0 - xi[1]), -(1.0 - xi[0])],
                        [(1.0 - xi[1]), -(1.0 + xi[0])],
                        [(1.0 + xi[1]), (1.0 + xi[0])],
                        [-(1.0 + xi[1]), (1.0 - xi[0])],
                    ]
                )
            )
            Je = np.transpose(nodal_values[:, :-1:]) @ gamma
            detJe = np.linalg.det(Je)
            V += detJe * w8[i]
    elif nodal_values.shape[0] == 8:
        # 3d hexahedra
        V = 0
        qp = (
            math.sqrt(3)
            / 3
            * np.array(
                [
                    [-1, -1, -1],
                    [1, -1, -1],
                    [-1, -1, 1],
                    [1, -1, 1],
                    [-1, 1, -1],
                    [1, 1, -1],
                    [-1, 1, 1],
                    [1, 1, 1],
                ]
            )
        )
        w8 = [1 for i in range(qp.shape[0])]
        for i, xi in enumerate(qp):
            gamma = (
                1
                / 8
                * np.array(
                    [
                        [
                            -(1.0 - xi[1]) * (1.0 + xi[2]),
                            -(1.0 - xi[0]) * (1.0 + xi[2]),
                            (1.0 - xi[0]) * (1.0 - xi[1]),
                        ],
                        [
                            (1.0 - xi[1]) * (1.0 + xi[2]),
                            -(1.0 + xi[0]) * (1.0 + xi[2]),
                            (1.0 + xi[0]) * (1.0 - xi[1]),
                        ],
                        [
                            (1.0 + xi[1]) * (1.0 + xi[2]),
                            (1.0 + xi[0]) * (1.0 + xi[2]),
                            (1.0 + xi[0]) * (1.0 + xi[1]),
                        ],
                        [
                            -(1.0 + xi[1]) * (1.0 + xi[2]),
                            (1.0 - xi[0]) * (1.0 + xi[2]),
                            (1.0 - xi[0]) * (1.0 + xi[1]),
                        ],
                        [
                            -(1.0 - xi[1]) * (1.0 - xi[2]),
                            -(1.0 - xi[0]) * (1.0 - xi[2]),
                            -(1.0 - xi[0]) * (1.0 - xi[1]),
                        ],
                        [
                            (1.0 - xi[1]) * (1.0 - xi[2]),
                            -(1.0 + xi[0]) * (1.0 - xi[2]),
                            -(1.0 + xi[0]) * (1.0 - xi[1]),
                        ],
                        [
                            (1.0 + xi[1]) * (1.0 - xi[2]),
                            (1.0 + xi[0]) * (1.0 - xi[2]),
                            -(1.0 + xi[0]) * (1.0 + xi[1]),
                        ],
                        [
                            -(1.0 + xi[1]) * (1.0 - xi[2]),
                            (1.0 - xi[0]) * (1.0 - xi[2]),
                            -(1.0 - xi[0]) * (1.0 + xi[1]),
                        ],
                    ]
                )
            )
            Je = np.transpose(nodal_values) @ gamma
            detJe = np.linalg.det(Je)
            V += detJe * w8[i]

    else:
        tb = sys.exc_info()[2]
        raise NotImplementedError(
            "Element type not implemented"
        ).with_traceback(tb)

    return V


def invert_negative_volumes(
    nodes: list, elements: list, names: list
) -> [bool, list]:

    cell_invert = False
    for index, elem in enumerate(elements):
        if set(elem.tags).isdisjoint(
            [phys.tag for phys in names if "interface" in phys.name]
        ):
            cell_invert = True
            nodal_values = np.array(
                [
                    [
                        nodes[i - 1].x_coord,
                        nodes[i - 1].y_coord,
                        nodes[i - 1].z_coord,
                    ]
                    for i in elem.nodes
                ]
            )

            if volume(nodal_values) < 0:
                if len(elem.nodes) == 4:
                    elements[index].nodes = elem.nodes[::-1]
                elif len(elem.nodes) == 8:
                    elements[index].nodes = (
                        elem.nodes[int(len(elem.nodes) / 2) :: 1]
                        + elem.nodes[: int(len(elem.nodes) / 2) : 1]
                    )
            nodal_values = np.array(
                [
                    [
                        nodes[i - 1].x_coord,
                        nodes[i - 1].y_coord,
                        nodes[i - 1].z_coord,
                    ]
                    for i in elem.nodes
                ]
            )
            if volume(nodal_values) < 0:
                tb = sys.exc_info()[2]
                raise RuntimeError(
                    "Element {} could not be inverted to have positive volume".format(
                        index + 1
                    )
                ).with_traceback(tb)

    return cell_invert, elements


def add_interface_elements(
    nodes: list[Node], elements: list[Element], names: list[PhysicalName]
) -> [list[Node], list[Element]]:
    class Duplicates:
        def __init__(self):
            self.first = []
            self.second = []

    duplicates = Duplicates()

    bulk_elem_type = -1
    bulk_elem_nodes = -1

    interface_tag = -1
    bulk_1_tag = -1
    bulk_2_tag = -1
    for phys in names:
        if "interface" in phys.name:
            interface_tag = phys.tag
        elif "bulk_1" in phys.name:
            bulk_1_tag = phys.tag
        elif "bulk_2" in phys.name:
            bulk_2_tag = phys.tag

    for elem in elements:
        if any([int(tag) == interface_tag for tag in elem.tags]):
            for node in elem.nodes:
                duplicates.first.append(node)
        if (
            any([int(tag) == bulk_1_tag for tag in elem.tags])
            and bulk_elem_type == -1
        ):
            bulk_elem_type = elem.type
            bulk_elem_nodes = len(elem.nodes)

    duplicates.first = list(dict.fromkeys(duplicates.first))
    duplicates.first.sort()
    for node in nodes:
        if node.id in duplicates.first:
            new_id = len(nodes) + 1
            nodes.append(
                Node(new_id, node.x_coord, node.y_coord, node.z_coord)
            )
            duplicates.second.append(new_id)
    node_is_offset = [False for i in nodes]

    for elem in elements:
        if not set(elem.tags).isdisjoint(
            [phys.tag for phys in names if "interface" in phys.name]
        ):
            if len(elem.nodes) == 2:
                Q = np.array(
                    [
                        [math.cos(math.pi / 2), -math.sin(math.pi / 2), 1],
                        [math.sin(math.pi / 2), math.cos(math.pi / 2), 1],
                        [0, 0, 1],
                    ]
                )
                nodal_values = np.array(
                    [
                        [
                            nodes[i - 1].x_coord,
                            nodes[i - 1].y_coord,
                            nodes[i - 1].z_coord,
                        ]
                        for i in elem.nodes
                    ]
                )
                vec1 = nodal_values[1] - nodal_values[0]
                normal_vector = Q @ vec1
                normal_vector /= np.linalg.norm(normal_vector)
            elif len(elem.nodes) == 4:
                normal_vector = 0
                nodal_values = np.array(
                    [
                        [
                            nodes[i - 1].x_coord,
                            nodes[i - 1].y_coord,
                            nodes[i - 1].z_coord,
                        ]
                        for i in elem.nodes
                    ]
                )
                vec1 = nodal_values[1] - nodal_values[0]
                vec2 = nodal_values[2] - nodal_values[1]
                vec3 = nodal_values[3] - nodal_values[2]
                vec4 = nodal_values[0] - nodal_values[1]
                normal_vector = (
                    np.cross(vec1, -vec2)
                    + np.cross(vec2, -vec1)
                    + np.cross(vec3, -vec2)
                    + np.cross(vec4, -vec3)
                )
                normal_vector /= np.linalg.norm(normal_vector)

            for node in elem.nodes:
                if not node_is_offset[
                    duplicates.second[duplicates.first.index(node)] - 1
                ]:
                    nodes[
                        duplicates.second[duplicates.first.index(node)] - 1
                    ].x_coord += (normal_vector[0] * 1e-12)*(1)
                    nodes[
                        duplicates.second[duplicates.first.index(node)] - 1
                    ].y_coord += (normal_vector[1] * 1e-12)*(1)
                    nodes[
                        duplicates.second[duplicates.first.index(node)] - 1
                    ].z_coord += (normal_vector[2] * 1e-12)*(1)
                    node_is_offset[
                        duplicates.second[duplicates.first.index(node)] - 1
                    ] = True

    for elem in elements:
        if any([int(tag) == interface_tag for tag in elem.tags]):
            elem.type = bulk_elem_type
            tmp_nodes = elem.nodes[:2:]
            for node in tmp_nodes[2::-1]:
                tmp_nodes.append(
                    duplicates.second[duplicates.first.index(node)]
                )
            if len(elem.nodes) == 4:
                tmp_nodes.append(elem.nodes[3])
                tmp_nodes.append(elem.nodes[2])
                for node in elem.nodes[2::1]:
                    tmp_nodes.append(
                        duplicates.second[duplicates.first.index(node)]
                    )
            elem.nodes = tmp_nodes

        elif any([int(tag) == bulk_2_tag for tag in elem.tags]):
            for index, node in enumerate(elem.nodes):
                if node in duplicates.first:
                    elem.nodes[index] = duplicates.second[
                        duplicates.first.index(node)
                    ]

    return nodes, elements


def invert_interfaces(
    nodes: list[Node], elements: list[Element], names: list[PhysicalName]
) -> list:
    for index, elem in enumerate(elements):
        if not set(elem.tags).isdisjoint(
            [phys.tag for phys in names if "interface" in phys.name]
        ):
            cell_invert = True
            nodal_values = np.array(
                [
                    [
                        nodes[i - 1].x_coord,
                        nodes[i - 1].y_coord,
                        nodes[i - 1].z_coord,
                    ]
                    for i in elem.nodes
                ]
            )

            if volume(nodal_values) < 0:
                elements[index].nodes = (
                    elem.nodes[int(len(elem.nodes) / 2) :: 1]
                    + elem.nodes[: int(len(elem.nodes) / 2) : 1]
                )
            nodal_values = np.array(
                [
                    [
                        nodes[i - 1].x_coord,
                        nodes[i - 1].y_coord,
                        nodes[i - 1].z_coord,
                    ]
                    for i in elem.nodes
                ]
            )
            if volume(nodal_values) < 0:
                tb = sys.exc_info()[2]
                raise RuntimeError(
                    "Element {} could not be inverted to have positive volume".format(
                        index + 1
                    )
                ).with_traceback(tb)

    return elements


def write_mesh(
    fobj: TextIO,
    nodes: list[Node],
    elements: list[Element],
    names: list[PhysicalName],
) -> bool:
    # write mesh info
    fobj.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    # write physical names
    fobj.write("$PhysicalNames\n")
    fobj.write("{}\n".format(len(names)))
    for name in names:
        fobj.write("{} {} {}\n".format(name.dimension, name.tag, name.name))
    fobj.write("$EndPhysicalNames\n")
    # write nodes
    fobj.write("$Nodes\n")
    fobj.write("{}\n".format(len(nodes)))
    for node in nodes:
        fobj.write(
            "{} {} {} {}\n".format(
                node.id, node.x_coord, node.y_coord, node.z_coord
            )
        )
    fobj.write("$EndNodes\n")
    # write elements
    fobj.write("$Elements\n")
    fobj.write("{}\n".format(len(elements)))
    for elem in elements:
        fobj.write("{} ".format(elem.id))
        fobj.write("{} ".format(elem.type))
        fobj.write("{} ".format(elem.n_tags))
        for tag in elem.tags:
            fobj.write("{} ".format(tag))
        for node in elem.nodes:
            fobj.write("{} ".format(node))
        fobj.write("\n")
    fobj.write("$EndElements\n")


def main(argv):

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print("interface_creator.py -i <inputfile> -o <outputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("interface_creator.py -i <inputfile> -o <outputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg.strip()
        elif opt in ("-o", "--ofile"):
            outputfile = arg.strip()
    if not os.path.exists(inputfile):
        tb = sys.exc_info()[2]
        raise OSError(
            "inputfile {} not found".format(inputfile)
        ).with_traceback(tb)

    fin = open(inputfile, "r")
    if not check_format(fin):
        tb = sys.exc_info()[2]
        raise SyntaxError(
            "GMSH Format must be legacy version 2.2"
        ).with_traceback(tb)
    physical_names = read_physical_names(fin)
    nodes = read_nodes(fin)
    elements = read_elements(fin, physical_names)
    fin.close()
    cell_invert, elements = invert_negative_volumes(
        nodes, elements, physical_names
    )

    nodes, elements = add_interface_elements(nodes, elements, physical_names)
    if cell_invert:
        invert_interfaces(nodes, elements, physical_names)

    fout = open(outputfile, "w")
    write_mesh(fout, nodes, elements, physical_names)


if __name__ == "__main__":
    main(sys.argv[1:])
