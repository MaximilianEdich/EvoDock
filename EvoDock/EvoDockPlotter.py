
import matplotlib.pyplot as plt
import numpy as np
import ast

def test_3D():
    X = [1, 2, 3]
    Y = [2, 5, 8]
    Z = [6, 4, 5]
    colors=["#0000FF", "#00FF00", "#FF0066"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')


    for i in range(len(X)):
        ax.scatter(X[i], Y[i], Z[i], color=colors[i])
    plt.show()
    return


def test_2D(name, title, exp):
    labels_alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    labels_hydrophob = ['F', 'I', 'W', 'L', 'V', 'M', 'Y', 'C', 'A', 'T', 'H', 'G', 'S', 'Q', 'R', 'K', 'N', 'E', 'P', 'D']

    labels = labels_hydrophob
    rev_labels = []
    for i in labels[::-1]:
        rev_labels.append(i)

    # init empy map
    data = np.random.random(size=(20, 20))
    for i in range(20):
        for j in range(20):
            data[i][j] = 0

    # load fitness with AAs
    results = []
    min = -1
    max = -1
    file = open("results_exp_final/" + exp + "/EvoDock_history_" + name + ".txt", 'r')
    lines = file.readlines()
    file.close()
    for line in lines[1:]:
        load_individual = [[], 0]
        content = line.strip().split(';')
        load_individual[1] = float(content[1]) ** 1
        load_individual[0] = ast.literal_eval(content[0])

        results.append(load_individual)
        if min == -1 or load_individual[1] < min:
            min = load_individual[1]
        if max == -1 or load_individual[1] > max:
            max = load_individual[1]

    print(min)
    print(min + (max - min) * 0.5)
    print(max)
    print()

    # map results
    for entry in results:
        x_coords = labels.index(entry[0][0])
        y_coords = rev_labels.index(entry[0][1])
        data[y_coords][x_coords] = entry[1]

    plt.imshow(data, interpolation='nearest')
    plt.xticks(range(20), labels)
    plt.yticks(range(20), rev_labels)
    plt.title(title)
    plt.xlabel("amino acid 1st position")
    plt.ylabel("amino acid 2nd position")
    plt.show()
    return

def skala():
    data = []
    labels = []
    for i in range(5):
        data.append([i])
        labels.append(i)
    plt.imshow(data, interpolation='nearest')
    plt.xticks(range(5), labels)
    plt.yticks(range(5), labels)
    plt.xlabel("amino acid 1st position")
    plt.ylabel("amino acid 2nd position")
    plt.show()
    return

test_2D("1RSZ", "1RSZ", "exp8")
test_2D("1GWR", "1GWR", "exp8")
test_2D("1NJE", "1NJE", "exp8")
test_2D("1BRA_h2o", "1BRA", "exp8")
test_2D("2PQL_h2o", "2PQL", "exp8")
test_2D("1OH0", "1OH0", "exp8")
test_2D("1BVG", "1BVG", "exp8")

#test_2D("1BRA_h2o_struct", "1BRA", "exp5")
#test_2D("1NJE_struct", "1NJE", "exp5")
#test_2D("1OH0_struct", "1OH0", "exp5")
#test_2D("1RSZ_struct", "1RSZ", "exp5")

#test_2D("2PQL_h2o", "2PQL", "exp5")
#test_2D("2PQL_h2o_struct", "2PQL", "exp5")

#test_2D("1GWR_struct", "1GWR", "exp5")
#test_2D("1BVG_struct", "1BVG", "exp5")