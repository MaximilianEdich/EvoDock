
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


def test_2D():
    labels = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
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
    file = open("results_exp_final/exp5_2D/1NJA_2D_new.txt", 'r')
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
    print(max)
    print(min)



    # map results
    for entry in results:
        x_coords = labels.index(entry[0][0])
        y_coords = rev_labels.index(entry[0][1])
        data[y_coords][x_coords] = entry[1]


    plt.imshow(data, interpolation='nearest')
    plt.xticks(range(20), labels)
    plt.yticks(range(20), rev_labels)
    plt.show()
    return

test_2D()