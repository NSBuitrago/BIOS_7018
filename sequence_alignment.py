import numpy as np
from Bio.Align import substitution_matrices


def global_align(x, y):
    """implementation of the Needleman-Wunsh global alignment algorithm:
            input: x, y = protein sequences
            output: x_align, y_align = global alignment of sequences x and y
    """
    # using the Blosum50 substitution matrix
    blosum50 = substitution_matrices.load('BLOSUM50')

    # Initialize matrices
    f_matrix = np.zeros((len(y)+1, len(x)+1))
    f_matrix[0, :] = np.arange(0, (len(x)+1)*-8, -8)
    f_matrix[:, 0] = np.arange(0, (len(y)+1)*-8, -8)
    traceback = {}

    # F matrix computation
    for i in range(0, len(f_matrix[0, :])-1):
        traceback[0, i] = 'h_gap'
        for j in range(0, len(f_matrix[:, 0])-1):
            traceback[j, 0] = 'v_gap'
            align = f_matrix[j, i] + blosum50.get((x[i], y[j]))
            v_gap = f_matrix[j, i+1] - 8
            h_gap = f_matrix[j+1, i] - 8
            f_value = max(align, v_gap, h_gap)
            f_matrix[j+1, i+1] = f_value

            # Assign direction to coordinates
            if f_value == align:
                traceback[j+1, i+1] = 'align'
            elif f_value == v_gap:
                traceback[j+1, i+1] = 'v_gap'
            else:
                traceback[j+1, i+1] = 'h_gap'

    # Recover Traceback
    index = max(traceback.keys())
    x_align, y_align = '', ''
    while True:
        try:
            if index == (0, 0):
                break
            elif traceback[index] == 'align':
                y_align = y_align + y[index[0]-1]
                x_align = x_align + x[index[1]-1]
                index = (index[0] - 1, index[1] - 1)
            elif traceback[index] == 'v_gap':
                y_align = y_align + y[index[0]-1]
                x_align = x_align + '-'
                index = (index[0] - 1, index[1])
            else:
                y_align = y_align + '-'
                x_align = x_align + x[index[1]-1]
                index = (index[0], index[1] - 1)
        except KeyError:
            break

    x_align = "".join(reversed(x_align))
    y_align = "".join(reversed(y_align))

    return x_align, y_align


def local_align(x, y):
    """implementation of the Smith-Watermann local alignment algorithm:
            input: x, y = protein sequences
            output: x_align, y_align = local alignment of sequences x and y
    """
    # using the Blosum50 substitution matrix
    blosum50 = substitution_matrices.load('BLOSUM50')

    # Initialize matrices
    f_matrix = np.zeros((len(y)+1, len(x)+1))
    f_matrix[0, :] = np.arange(0, (len(x)+1)*-8, -8)
    f_matrix[:, 0] = np.arange(0, (len(y)+1)*-8, -8)
    traceback = {}

    # F matrix computation
    for i in range(0, len(f_matrix[0, :])-1):
        traceback[0, i] = 'h_gap'
        for j in range(0, len(f_matrix[:, 0])-1):
            traceback[j, 0] = 'v_gap'
            align = f_matrix[j, i] + blosum50.get((x[i], y[j]))
            v_gap = f_matrix[j, i+1] - 8
            h_gap = f_matrix[j+1, i] - 8
            # local alignment includes 0 in the recurrence relation
            f_value = max(align, v_gap, h_gap, 0)
            f_matrix[j+1, i+1] = f_value

            # Assign direction to coordinates
            if f_value == align:
                traceback[j+1, i+1] = 'align'
            elif f_value == v_gap:
                traceback[j+1, i+1] = 'v_gap'
            elif f_value == h_gap:
                traceback[j+1, i+1] = 'h_gap'
            else:
                traceback[j+1, i+1] = 'stop'

    # Recover Traceback
    max_f_index = np.where(f_matrix == f_matrix.max())
    index = (max_f_index[0][0], max_f_index[1][0])
    x_align, y_align = '', ''
    while True:
        try:
            if index == (0, 0):
                break
            elif traceback[index] == 'align':
                y_align = y_align + y[index[0]-1]
                x_align = x_align + x[index[1]-1]
                index = (index[0] - 1, index[1] - 1)
            elif traceback[index] == 'v_gap':
                y_align = y_align + y[index[0]-1]
                x_align = x_align + '-'
                index = (index[0] - 1, index[1])
            elif traceback[index] == 'h_gap':
                y_align = y_align + '-'
                x_align = x_align + x[index[1]-1]
                index = (index[0], index[1] - 1)
            else:
                break
        except KeyError:
            break

    x_align = "".join(reversed(x_align))
    y_align = "".join(reversed(y_align))

    return x_align, y_align, f_matrix


x = 'HEAGAWGHEE'
y = 'PAWHEAE'
if __name__ == "__main__":
    x_align, y_align, f_matrix = local_align(x, y)
