

file = open("BLOSUM62.txt")
lines = file.readlines()
file.close()

# load head row
labels_cols = lines[0].strip().split(' ')
print(labels_cols)
new_list = []
for item in labels_cols:
    if item != '':
        new_list.append(item)
labels_cols = new_list
print(labels_cols)
print(len(labels_cols))

# load head column
labels_rows = []
for line in lines[1:]:
    labels_rows.append(line[0])
print(labels_rows)
print(len(labels_rows))

# check if identical
if len(labels_rows) != len(labels_cols):
    exit("Matrix is not quadratical!")
index = 0
for item in labels_cols:
    if labels_cols[index] != labels_rows[index]:
        exit("Matrix labels for row and columns are not identical!")
    index += 1

# build matrix
matrix = []
for line in lines[1:]:
    row_load = line[1:].strip().split(' ')
    row = []
    for item in row_load:
        if item != '':
            row.append(item)
    matrix.append(row)


input_AA = ['Y', 'S']
row_i = labels_rows.index(input_AA[1])
col_i = labels_cols.index(input_AA[0])
print(matrix[row_i][col_i])
