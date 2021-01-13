import argparse

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-o", "--out_name", type=str, required=True,
                    help="Outputpath and name, only a .txt is added.")
args = parser.parse_args()

start_stuffix = 5
step_size = 5
start_repeat_until = 30
replace = "_S_"

lines = []
lines.append("python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 1/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 2/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 3/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 4/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 5/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 6/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 7/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 8/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 9/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 10/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 11/\n"
            "python3 EvoDock.py -s is_exp_final/exp6/is_exp6_2PQL_h2o_N25.txt -r is_exp_final/exp6/C/routineC_5__S_.txt -o results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -fn 12/\n"
            "python3 ResultsReader.py -r results_exp_final/exp6/2PQL_H2O/C/n25_s_S_/ -i 0"
             )

commands = ""
for line in lines:
    suffix = start_stuffix
    repeat_until_inclusive = start_repeat_until
    for _ in range(int((repeat_until_inclusive - start_stuffix) / step_size + 1)):
        commands += line.replace(replace, str(suffix)) + "\n"
        suffix += step_size

file = open(args.out_name + ".txt", 'w')
file.write(commands)
file.close()