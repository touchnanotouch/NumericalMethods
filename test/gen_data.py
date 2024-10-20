import re

import numpy as np
import random as rd
import pandas as pd


# n_s = [10, 25, 50, 75, 100]
n_s = [10, 25, 50, 75, 100, 250, 500, 750, 1000]
# n_s = [1000]

for n in n_s:
    min_tens = -int("10" * (len(str(n)) - 1))
    max_tens = int("10" * len(str(n)))

    min_val = -1 * min_tens / 100
    max_val = 1 * max_tens / 100

    matrix = np.random.uniform(min_val, max_val, size=(n, n))

    for i in range(n):
        diagonal_element = np.abs(matrix[i, i])
        other_elements = np.sum(np.abs(matrix[i, :i])) + np.sum(np.abs(matrix[i, i+1:]))
        matrix[i, i] = diagonal_element + other_elements + 1

    b = [rd.uniform(-10 * min_tens, 10 * max_tens) for _ in range(n)]

    sol = np.linalg.solve(matrix, b)

    with open(f"matrix_n{n}.txt", "w") as mf:
        np.savetxt(mf, matrix, fmt="%.2f")
    with open(f"vector_n{n}.txt", "w") as vf:
        np.savetxt(vf, b, fmt="%.2f")
    with open(f"solution_n{n}.txt", "w") as sf:
        np.savetxt(sf, sol, fmt="%.2f")
