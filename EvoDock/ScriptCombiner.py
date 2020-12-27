
import argparse

parser = argparse.ArgumentParser(description="Combines two script halves to one.")
parser.add_argument("-hd", "--head", type=str, required=True, help="Upper half.")
parser.add_argument("-b", "--body", type=str, required=True, help="Lower half.")
parser.add_argument("-o", "--out", type=str, required=True, help="Out-path + name, write .txt !.")
args = parser.parse_args()


head_file = open(args.head)
head = head_file.read()
head_file.close()

body_file = open(args.body)
body = body_file.read()
body_file.close()

result = head + "\n" + body

r_file = open(args.out, 'w')
r_file.write(result)
r_file.close()
