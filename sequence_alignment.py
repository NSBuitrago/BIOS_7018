import numpy as np
from Bio.Align import substitution_matrices


def global_align(x, y):
    # implementation of the Needleman-Wunsh global alignment algorithm
    # using the Blosum50 substitution matrix
    blosum50 = substitution_matrices.load('BLOSUM50')
    f_matrix = np.zeros((len(y)+1, len(x)+1))
    f_matrix[0, :] = np.arange(0, (len(x)+1)*-8, -8)
    f_matrix[:, 0] = np.arange(0, (len(y)+1)*-8, -8)

    for i in range(0, len(f_matrix[0, :])-1):
        for j in range(0, len(f_matrix[:, 0])-1):
            sub = blosum50.get((x[i], y[j]))
            f_matrix[j+1, i+1] = max((f_matrix[j, i] + sub),
                                     (f_matrix[j, i+1] - 8), (f_matrix[j+1, i] - 8))

    return f_matrix


x = 'HEAGAWGHEE'
y = 'PAWHEAE'
if __name__ == "__main__":
    f_matrix = global_align(x, y)
