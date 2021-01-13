import argparse
import ast

parser = argparse.ArgumentParser(description="EXtends the search space.")
parser.add_argument("-pf", "--population_files", type=str, required=True, nargs='+', help="Text files containing population.")
parser.add_argument("-o", "--out_path", type=str, required=True, help="Output path.")
args = parser.parse_args()

combined = []
for file_name in args.population_files:
    file = open(file_name)
    lines = file.readlines()
    file.close()

    for line in lines[1:]:
        load_individual = [[], 0]
        content = line.strip().split(';')
        load_individual[1] = float(content[1])
        load_individual[0] = ast.literal_eval(content[0])
        is_duplicate = False
        for item in combined:
            if item[0] == load_individual[0]:
                is_duplicate = True
                print("merging scores of duplicates...")
                print(load_individual)
                print(item)
                new_score = (load_individual[1] + item[1]) / 2
                print("new score: " + str(new_score))
                item[1] = new_score
        if not is_duplicate:
            combined.append(load_individual)
print(len(combined))


def get_score(input_individual):
    return input_individual[1]


combined.sort(reverse=True,key=get_score)
file_content = "Entries in combined population (" + str(args.population_files) + "): " + str(len(combined)) + "\n"
for individual in combined:
    file_content += str(individual[0]) + ";" + str(individual[1]) + "\n"

out_file = open(args.out_path, 'w')
out_file.write(file_content)
out_file.close()