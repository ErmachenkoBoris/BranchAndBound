import urllib.request
import numpy as np
import cplex


def read_graph(file_path):
    _n = _m = -1
    data = urllib.request.urlopen(file_path)
    for line in data:
        line = line.decode('ascii')
        line = line.strip('\n')
        if line.startswith('p'):
            _n = int(line.split(' ')[2])
            _m = int(line.split(' ')[3])
            break
    _confusion_matrix = np.zeros((_n, _n))
    for line in data:
        line = line.decode('ascii')
        line = line.strip('\n')
        if line.startswith('e'):
            i = int(line.split(' ')[1]) - 1
            j = int(line.split(' ')[2]) - 1
            _confusion_matrix[i, j] = 1
            _confusion_matrix[j, i] = 1
    return _n, _m, _confusion_matrix


if __name__ == '__main__':
    # path = "http://iridia.ulb.ac.be/~fmascia/files/DIMACS/brock200_2.clq"
    path = "http://iridia.ulb.ac.be/~fmascia/files/DIMACS/C125.9.clq"
    # path = "http://iridia.ulb.ac.be/~fmascia/files/DIMACS/keller4.clq"
    n, m, confusion_matrix = read_graph(path)

    max_clique_model = cplex.Cplex()

    max_clique_model.variables.add(names=["y" + str(i) for i in range(n)],
                                   types=[max_clique_model.variables.type.continuous for i in range(n)])
    for i in range(n):
        max_clique_model.variables.set_lower_bounds(i, 0.0)

    constrains = []
    constrains_names = []
    constrains_types = []
    constrains_right_parts = []

    for i in range(confusion_matrix.shape[0]):
        for j in range(i + 1, confusion_matrix.shape[1]):
            if confusion_matrix[i][j] == 0.0:
                constrains.append([["y" + str(i), "y" + str(j)], [1, 1]])
                constrains_names.append("constraint_" + str(i) + "_" + str(j))
                constrains_types.append('L')
                constrains_right_parts.append(1.0)

    for i in range(n):
        max_clique_model.linear_constraints.add(
            lin_expr=constrains,
            rhs=constrains_right_parts,
            names=constrains_names,
            senses=constrains_types
        )

    for i in range(n):
        max_clique_model.objective.set_linear("y" + str(i), 1)
    max_clique_model.objective.set_sense(max_clique_model.objective.sense.maximize)

    max_clique_model.solve()
    values = max_clique_model.solution.get_values()
    result_size = 0
    print('values ', values)
    result = 0
    for v in values:
        result = result + v

    print("Result: ", result)

    print("-------")

    print(values)