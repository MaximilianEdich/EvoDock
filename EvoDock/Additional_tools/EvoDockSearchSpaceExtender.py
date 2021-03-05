import argparse
import ast

from pyrosetta import init, pose_from_file

init()

parser = argparse.ArgumentParser(description="EXtends the search space.")
parser.add_argument("-pf", "--population_file", type=str, required=True, help="Text files containing population.")
parser.add_argument("-aai", "--res_ids_input", type=str, required=True, nargs='+', help="Res ids of input.")
parser.add_argument("-p", "--protein", type=str, required=True, help="Reference Protein.")
parser.add_argument("-aax", "--res_ids_extend", type=str, required=True, nargs='+', help="Res ids used to extend.")
parser.add_argument("-o", "--out_path", type=str, required=True, help="Output path.")
args = parser.parse_args()

file = open(args.population_file)
lines = file.readlines()
file.close()


# insert numbers at right position
input_aa = list(args.res_ids_input).copy()
for aa_pos in args.res_ids_extend:
    index = 0
    added = False
    for aa in input_aa:
        if int(aa_pos) < int(aa):
            input_aa.insert(index, aa_pos)
            added = True
            break
        index += 1
    if not added:
        input_aa.append(aa_pos)

# mark extension positions
extension = []
for aa in input_aa:
    extend = False
    for aa_pos in args.res_ids_extend:
        if aa == aa_pos:
            extend = True
            break
    if extend:
        extension.append('X')
    else:
        extension.append('_')

pose = pose_from_file(args.protein)
extension_mask = []
index = 0
for aa in input_aa:
    if extension[index] == "_":
        extension_mask.append("_")
    else:
        extension_mask.append(pose.residue(int(aa)).name1())
    index += 1
print(extension_mask)

population = []
for line in lines[1:]:
    load_individual = [[], 0]
    content = line.strip().split(';')
    load_individual[1] = float(content[1])
    load_individual[0] = ast.literal_eval(content[0])
    # extend
    new_aas = []
    index = 0
    for pos in extension_mask:
        if pos == "_":
            new_aas.append(load_individual[0][index])
            index += 1
        else:
            new_aas.append(pos)
    load_individual[0] = new_aas
    population.append(load_individual)


population.sort(reverse=False)
file_content = "Entries in extended population (extended from " + str(args.res_ids_input) + " to " + str(input_aa) + "): " + str(len(population)) + "\n"
for individual in population:
    file_content += str(individual[0]) + ";" + str(individual[1]) + "\n"

out_file = open(args.out_path, 'w')
out_file.write(file_content)
out_file.close()



