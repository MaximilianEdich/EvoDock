import argparse

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-i", "--input_file", type=str, required=True,
                    help="Base input file.")
parser.add_argument("-o", "--output", type=str, required=True,
                    help="Output folder.")
args = parser.parse_args()


file = open(args.input_file)
content = file.read()
file.close()

repN = "N"
repN_in = [5, 10, 15, 20, 25, 30]

repX = "X"
repX_in = [1, 2]

repS = "S"
repS_in = [5, 10, 15, 20, 25, 30]

repR = "R"
repR_in = [10]

for n in repN_in:
    for s in repS_in:
        # get and edit content
        new_c = str(content)
        new_c = new_c.replace(repN, str(n))
        new_c = new_c.replace(repS, str(s))
        new_c = new_c.replace("Z", str(int(s/2)))
        # get and edit file name
        name = str(args.input_file).split('/')[-1]
        name = name.replace(repN, str(n))
        name = name.replace(repS, str(s))
        file = open(args.output + name, "w")
        file.write(new_c)
        file.close()