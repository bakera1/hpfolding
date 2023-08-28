import unittest

from pulp import LpVariable, LpProblem, LpStatus, LpMaximize, GLPK, value, lpSum

from hp.core.common import generate_neighbours, wrap_compute_count_pairs, generate_unique_edge, compute_position_match, \
    SquareGrid


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.protein = """1001011"""

    def test_protein_solve(self):

        # setup LpProblem with sample Protein number of grid points
        prob = LpProblem('Protein_Folding_HP_Model', LpMaximize)
        n = len(self.protein)

        # we given each point in the n * n grid an identifer n**2
        i_n_by_n = [f'x_i_p_{i}_{p}' for i in range(1, n + 1) for p in range(1, n ** 2 + 1)]

        # i position in the protein string on any point in the n ** n grid
        i_assigned_p_var = LpVariable.dicts("x_i_p", \
                                            i_n_by_n, \
                                            lowBound=0, \
                                            upBound=1, \
                                            cat='Integer')

        for i in range(1, n + 1):
            prob += lpSum([i_assigned_p_var[f'x_i_p_{i}_{p}'] \
                           for p in range(1, n ** 2 + 1)]) == 1

        # ensure that no point on the grid is assigned more than one position in the
        # string; for each pair of positions i,j in the string; no grid point

        for p in range(1, n ** 2 + 1):
            # prob += lpSum([i_assigned_p_var[f'x_i_p_{i}_{p}'] + i_assigned_p_var[f'x_i_p_{j}_{p}'] for (i,j) in zip(range(1, n+1, 2), range(2,n+2, 2))]) <= 1
            prob += lpSum([i_assigned_p_var[f'x_i_p_{i}_{p}'] for i in range(1, n + 1)]) <= 1

        # ensure connectedness
        # need to ensure that adjacent positions on the string are assigned to
        # neighbouring points on the grid

        for i in range(1, n):
            for p in range(1, n ** 2 + 1):
                prob += i_assigned_p_var[f'x_i_p_{i}_{p}'] - lpSum(
                    [i_assigned_p_var[f'x_i_p_{i + 1}_{q}'] for q in generate_neighbours(p=p, n=n)]) <= 0

        # inequalities to detect contacts
        grid = [f'i_p_{p}' for p in range(1, n ** 2 + 1)]
        # i position in the protein string on any point in the n ** n grid
        i_p_var = LpVariable.dicts("i_p", \
                                   grid, \
                                   lowBound=0, \
                                   upBound=1, \
                                   cat='Integer')

        # for each position in the grid
        match_list = compute_position_match(protein=self.protein)

        for p in range(1, n ** 2 + 1):
            prob += lpSum([i_assigned_p_var[f'x_i_p_{i}_{p}'] for i in range(1, n + 1) if i in match_list]) - \
                    i_p_var[f'i_p_{p}'] == 0

        edge_list = [f'c_{edge.p}_{edge.q}' for edge in generate_unique_edge(n=n)]

        # i position in the protein string on any point in the n ** n grid
        c_p_q_var = LpVariable.dicts("c_p_q", \
                                     edge_list, \
                                     lowBound=0, \
                                     upBound=1, \
                                     cat='Integer')

        for edge in generate_unique_edge(n=n):
            prob += i_p_var[f'i_p_{edge.p}'] + i_p_var[f'i_p_{edge.q}'] - 2 * c_p_q_var[f'c_{edge.p}_{edge.q}'] >= 0

        # objective function
        offset = wrap_compute_count_pairs(protein=self.protein, match='11')
        offset_var = LpVariable('offset', lowBound=offset, upBound=offset)
        prob += lpSum([c_p_q_var[f'c_{edge.p}_{edge.q}'] for edge in generate_unique_edge(n=n)]) - offset_var

        prob.writeLP("simple.lp")

        status = prob.solve(GLPK(msg=0))

        if LpStatus[status] == 'Optimal':
            print(LpStatus[status])
            print('allocation:', ["%s: (%.2f)" % (v.name, v.varValue) for v in prob.variables()])

            g = SquareGrid(n=n)
            g.render_var(prob=prob, protein=self.protein)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
