import argparse

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-o", "--out_name", type=str, required=True,
                    help="Outputpath and name, only a .txt is added.")
args = parser.parse_args()

start_stuffix = 2
start_repeat_until = 26
replace = "#"

lines = []
lines.append("python3 ResultsReader.py -r export3/export3/mutate_1ogx/results#/ -i 5\n"
             "python3 ResultsReader.py -r export3/export3/mutate_1rsz/results#/ -i 6"
             )

commands = ""
for line in lines:
    suffix = start_stuffix
    repeat_until_inclusive = start_repeat_until
    for _ in range(repeat_until_inclusive - start_stuffix + 1):
        commands += line.replace(replace, str(suffix)) + "\n"
        suffix += 1

file = open(args.out_name + ".txt", 'w')
file.write(commands)
file.close()