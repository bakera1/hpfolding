import unittest

from hp.core.common import classify_p_in_grid, GridPointType, generate_neighbours, \
    wrap_compute_count_pairs, generate_unique_edge, compute_position_match, SquareGrid


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.results_generate_neighbours = {1: [2, 4],
                                            2: [1, 3, 5],
                                            3: [2, 6],
                                            4: [1, 5, 7],
                                            5: [2, 4, 6, 8],
                                            6: [3, 5, 9],
                                            7: [4, 8],
                                            8: [5, 7, 9],
                                            9: [6, 8]}

        self.local_n = 3
        self.results_generate_unique_edge = [12, 14, 23, 25, 36, 45, 47, 56, 58, 69, 78, 89]
        self.result_square_grid = {3: [[1, 2, 3], [1, 0, 0]], 2: [[2], [0]], 4: [[2], [1]]}

    def test_square_grid(self):
        g = SquareGrid(n=11)
        d = g.convert_xy_to_dict(xy=[(3, 1, 1), (3, 2, 0), (3, 3, 0), (2, 2, 0), (4, 2, 1)])
        self.assertEqual(True, d == self.result_square_grid)

    def test_compute_position_match2(self):
        self.assertEqual(True, compute_position_match(protein='110') == [1, 2])

    def test_generate_unique_edge(self):
        self.assertEqual(True, [int(f'{edge.p}{edge.q}') for edge in
                                generate_unique_edge(n=self.local_n)] == self.results_generate_unique_edge)

    def test_wrap_compute_count_pairs(self):
        self.assertEqual(True, wrap_compute_count_pairs(protein='1110101100', match='11') == 3)

    def test_compute_position_match(self):
        self.assertEqual(True, compute_position_match(protein='1110101100') == [1, 2, 3, 5, 7, 8])

    def test_generate_neighbours(self):
        for p in range(1, self.local_n ** 2 + 1):
            # print(f'{p} {generate_neighbours(p=p, n=local_n)}')
            self.assertEqual(True, self.results_generate_neighbours[p] == generate_neighbours(p=p, n=self.local_n))

    def test_classify_p_in_grid(self):
        self.assertEqual(True, classify_p_in_grid(p=4, n=8)[0] == GridPointType.Edge and \
                         classify_p_in_grid(p=1, n=8)[0] == GridPointType.Corner and \
                         classify_p_in_grid(p=64, n=8)[0] == GridPointType.Corner)  # add assertion here


if __name__ == '__main__':
    unittest.main()
