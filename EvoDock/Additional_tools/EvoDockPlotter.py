
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


def test_2D(path, title):
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
    file = open(path, 'r')
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

        case = 5
        if case == 1:
            if load_individual[0] == ['S', 'P']:
                load_individual[1] = 185997301.01523507
        elif case == 2:
            pass
        elif case == 3:
            if load_individual[0] == ['E', 'I']:
                load_individual[1] = 185993375.9506488
        elif case == 5:
            if load_individual[0] == ['C', 'N']:
                load_individual[1] = 152

        if load_individual[1] < 0:
            load_individual[1] = 0

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

#test_2D("results_exp_additional/2AIK/2AIK/EvoDock_history_exp1.txt", "2AIK - X 71, X 73")
#test_2D("results_exp_additional/2AIK/2AIK/EvoDock_history_exp2.txt", "2AIK - A 154, X 73")
#test_2D("results_exp_additional/2AIK/2AIK/EvoDock_history_exp3.txt", "2AIK - A 176, X 71")

#test_2D("results_exp_additional/2AIK/struct/EvoDock_history_Exp1_struct.txt", "2AIK - X 71, X 73")
test_2D("results_exp_additional/2AIK/struct/EvoDock_history_Exp2_struct.txt", "2AIK - A 154, X 73")
#test_2D("results_exp_additional/2AIK/struct/EvoDock_history_Exp3_struct.txt", "2AIK - A 176, X 71")