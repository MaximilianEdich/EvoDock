import argparse

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-o", "--out_name", type=str, required=True,
                    help="Outputpath and name, only a .txt is added.")
args = parser.parse_args()

start_stuffix = 1
start_repeat_until = 5
replace = "#"

lines = []
lines.append("python3 EvoDock.py -s is_pre/is_1ANE_#.txt -r routine.txt -o is_pre/1ANE/results#/\n"
             "python3 EvoDock.py -s is_pre/is_1ANE_#.txt -r routine.txt -o is_pre/1ANE/results#/\n"
             "python3 EvoDock.py -s is_pre/is_1ANE_#.txt -r routine.txt -o is_pre/1ANE/results#/\n"
             "python3 EvoDock.py -s is_pre/is_1ANE_#.txt -r routine.txt -o is_pre/1ANE/results#/\n"
             "python3 EvoDock.py -s is_pre/is_1ANE_#.txt -r routine.txt -o is_pre/1ANE/results#/\n"
             "python3 ResultsReader.py -r is_pre/1ANE/results#/ -i 0"
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