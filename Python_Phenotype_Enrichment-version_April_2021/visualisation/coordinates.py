"""
;===================================================================================================
; Title:   Visualization of enriched phenotypes: coordinates
; Authors: Alexandra, Ivana, Irene, Thao
;===================================================================================================

"""
import random


def gene_y_position(y, total_genes, count_n):
    """ Creates y coordinate of a gene node
     Y coordinate must be placed lower on the y axis as a previous one.
     As long as the count is lower than the total number of genes, the y position is produced.

     :param y: y coordinate/position of a previous node
     :param total_genes: the total number of genes to be visualised as nodes
     :param count_n: the current round of iteration as integer
     """
    # TODO do something about this; the highest y coordinate is quite low
    if count_n <= total_genes:
        y = (count_n / total_genes) * 100
    else:
        pass
    return y


def check_coordinates(coordinates, coordinate_float, x1, x2, x3, x4):
    """ Checks if node y coordinate (position) is already present in list of positions
    If it is not true, the position is returned and added to the list of positions
    otherwise, the function is called upon itself.

    :param coordinates: list where the non-overlapping positions are stored
    :param coordinate_float: the y coordinate (position) of related gene node
    :param x1: the left border of left part of graph
    :param x2: the right border of left part of graph
    :param x3: the left border of right part of graph
    :param x4: the right border of right part of graph
    """

    # x coordinate - left or right side of gene nodes:
    if bool(random.getrandbits(1)) is True:
        x_coordinate = random.uniform(x1, x2)
    else:
        x_coordinate = random.uniform(x3, x4)

    # since during the first iterration through df, position_float=1 and we don't wnt all to be too close
    # y coordinate
    if coordinate_float == 1:
        y_coordinate = random.uniform(0, 100)
    else:
        y_coordinate = coordinate_float + random.uniform(0, 20)

    coordinate = (x_coordinate, y_coordinate)

    # chek if given coordinates exist already
    if coordinate not in coordinates:
        coordinates.append(coordinate)
    else:
        check_coordinates(coordinates, coordinate_float, x1, x2, x3, x4)  # recursioon

    return x_coordinate, y_coordinate
