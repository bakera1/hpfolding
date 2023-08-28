from termcolor import colored, cprint
import numpy as np

from collections import namedtuple
from typing import List
from enum import Enum

Edge = namedtuple("Edge", "left top right bottom")
Corner = namedtuple("Corner", "top_left top_right bottom_left bottom_right")


class GridPointType(Enum):
    """ simple enum to classify each point in a grid """
    Edge = 'Edge'
    Corner = 'Corner'
    Inner = 'Inner'


class EdgeLink:
    """ used to track edge links """

    def __init__(self, p, q):
        self.p = p
        self.q = q

    def __repr__(self):
        return f'{self.p}->{self.q}'

    def __str__(self):
        return f'{self.p}->{self.q}'

    def has_traversed(self, candidate_p, candidate_q, verbose=False):
        if verbose:
            print(
                f'inside has_traversed check {(self.p == candidate_p and self.q == candidate_q) or (self.q == candidate_p and self.p == candidate_q)}')
        return True if (self.p == candidate_p and self.q == candidate_q) or (
                self.q == candidate_p and self.p == candidate_q) else False


def classify_p_in_grid(p: int, n: int) -> GridPointType:
    first_column = [i for i in range(1, n ** 2, n)]
    last_column = [i for i in range(n, n ** 2, n)]
    top_row = [i for i in range(1, n)]
    bottom_row = [i for i in range(n ** 2 - n, n ** 2 + 1)]
    corners = [1, n, n ** 2 - n + 1, n ** 2]
    if p in corners:
        return GridPointType.Corner, Corner(p == 1, p == n, p == n ** 2 - n + 1, p == n ** 2)
    elif p in first_column or p in last_column or p in top_row or p in bottom_row:
        return GridPointType.Edge, Edge(p in first_column, p in top_row, p in last_column, p in bottom_row)
    else:
        return GridPointType.Inner, None


def generate_neighbours(p: int, n: int) -> list():
    """ generate the N_p neighbours to point p in the grid
    general 4 points for point in the grid
    corner case only 2 neighbours
    first or last row, or first or last column only has """

    # check against bad data points
    if p not in [x for x in range(n ** 2 + 1)] or p == 0:
        raise Exception(f'{p} outside range {0} to {n ** 2}')
    output = list()
    grid_type, edge = classify_p_in_grid(p=p, n=n)
    # print(grid_type)

    # inner easiest case here
    if grid_type in [GridPointType.Inner]:
        output = list([p - 1, p - n, p + n, p + 1])
    # classify the location on the edges
    elif grid_type in [GridPointType.Edge]:
        if edge.left:
            output = [p - n, p + n, p + 1]
        elif edge.top:
            output = [p - 1, p + n, p + 1]
        elif edge.right:
            output = [p - 1, p - n, p + n]
        elif edge.bottom:
            output = [p - 1, p - n, p + 1]
    # handle the corners
    elif grid_type in [GridPointType.Corner]:
        if edge.top_left:
            output = list([p + 1, p + n])
        if edge.top_right:
            output = list([n - 1, n + n])
        if edge.bottom_left:
            output = list([p - n, p + 1])
        if edge.bottom_right:
            output = list([p - 1, n ** 2 - n])
    return sorted([i for i in list(set(output)) if i != 0])


def wrap_generate_neighbours(p: int, n: int) -> List[EdgeLink]:
    candidate_neighbours: list(int) = generate_neighbours(p=p, n=n)
    for candidate in candidate_neighbours:
        yield EdgeLink(p=p, q=candidate)


def generate_unique_edge(n: int) -> List[EdgeLink]:
    ### function to build unique set of edges to traverse from grid of size n ###

    edge_link_link = list()

    # iterate the entire grid
    for p in range(1, n ** 2 + 1):

        # print(f'p::{p}')
        # iterate all candidate neighbours
        for candidate_neighbour in wrap_generate_neighbours(p=p, n=n):
            # print(candidate_neighbour)

            if len(edge_link_link) == 0:
                # print(f'{candidate_neighbour} not traversed aaa {str(candidate_neighbour)}')
                edge_link_link.append(candidate_neighbour)
                continue

            check = False
            for edge in edge_link_link:
                check = edge.has_traversed(candidate_p=candidate_neighbour.p, candidate_q=candidate_neighbour.q)
                # print(f'check edge.has_traversed {edge.p} {edge.q} -->{candidate_neighbour.p} {candidate_neighbour.q} check {check}')

                if check:
                    # break if we have already traversed this link
                    break
                else:
                    check *= True
            if not check:
                # print(f'{candidate_neighbour} not traversed')
                edge_link_link.append(candidate_neighbour)
            # else:
            # print(f'check {check}')

    return edge_link_link


def compute_overcount_pairs(protein: str, match: str, mode=0) -> int:
    ### count the number of distinct matches of match in protein ###
    for position in range(len(protein) - 1):
        sequence = f'{protein[position]}{protein[position + 1]}'
        # print(f'{position}->{sequence} ')
        if mode == 0:
            yield 1 if sequence == match else 0
        else:
            yield position if sequence == match else None


def wrap_compute_count_pairs(protein: str, match='11'):
    return sum([count for count in compute_overcount_pairs(protein=protein, match=match)])


def compute_position_match(protein: str, match='1'):
    return [position + 1 for position in range(len(protein) - 1) if protein[position] == '1']


class SquareGrid:

    def __init__(self, n):
        """ @param n : int size of one side of the square grid """
        self.n = n
        self.grid = np.array(range(1, n ** 2 + 1)).reshape(n, n)
        self.redline = colored(' + ', 'red')
        self.greenline = colored(' + ', 'green')
        self.blackline = colored(' + ', 'black')

        self.cache_lookup = dict()
        dd = np.array(range(self.n ** 2)).reshape((self.n, self.n))
        rows, cols = dd.shape
        for p in range(1, n ** 2 + 1):
            for x in range(0, rows):
                for y in range(0, cols):
                    self.cache_lookup[f'x_i_p_x_i_p_{p}_{dd[x, y] + 1}'] = (x, y, 0)

    def convert_xy_to_dict(self, xy):
        return_dict = dict()
        for item in xy:
            if item[0] in return_dict:
                return_dict[item[0]][0].append(item[1])
                return_dict[item[0]][1].append(item[2])
            else:
                return_dict[item[0]] = [[item[1]], [item[2]]]
        return return_dict

    def parse_problem_variable_list(self, location_list):
        # 'x_i_p_x_i_p_1_13: (1)'
        p0 = [var.name.replace('x_i_p_x_i_p_', '') for var in location_list]
        for var in p0:
            index, p = var.split('_')
            yield index, p

    def render_var(self, prob, protein):
        #print('allocation:',
        #      ["%s: (%d)" % (v.name, v.varValue) for v in prob.variables() if 'x_i_p' in v.name and v.varValue == 1.0])
        location_list = []
        index = 0
        for v in prob.variables():
            if 'x_i_p' in v.name and v.varValue == 1.0:
                (a, b, cdefault) = self.cache_lookup[v.name]
                location_list.append((a, b, int(protein[index])))
                index += 1
        self.render(xy=location_list,protein=protein)

    def render(self, xy=list(), protein=None):
        lookup_dict = self.convert_xy_to_dict(xy)
        #print(lookup_dict)
        rows, cols = self.grid.shape
        for x in range(0, rows):
            row_list = []
            for y in range(0, cols):
                row_list.append('|')
            output = ''
            for i, item in enumerate(row_list):
                if x in lookup_dict:
                    if i in lookup_dict[x][0]:
                        index = lookup_dict[x][0].index(i)
                        if lookup_dict[x][1][index] in [1]:
                            output += self.greenline
                        elif lookup_dict[x][1][index] in [0]:
                            output += self.redline
                    else:
                        output += self.blackline
                else:
                    output += self.blackline
            print(output)
        print(protein)
