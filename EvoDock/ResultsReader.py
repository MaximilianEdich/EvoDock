import argparse
import os
import ast
import math
import xlsxwriter

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-r", "--results", type=str, required=True,
                    help="Results Folder with all results")
args = parser.parse_args()

dirs = os.listdir(args.results)

summary = []
for directory in dirs:
    print(directory)
    if os.path.isfile(directory):
        print(directory)
    history = open(args.results + directory + "/EvoDock_history.txt", 'r')
    lines = history.readlines()
    history.close()
    for line in lines[1:]:
        load_individual = [[], 0]
        content = line.strip().split(';')
        load_individual[1] = float(content[1])
        load_individual[0] = ast.literal_eval(content[0])
        added = False
        for sum in summary:
            if sum[0] == load_individual[0]:
                sum.append(load_individual[1])
                added = True
        if not added:
            summary.append(load_individual)

times = []
for directory in dirs:
    run_info = open(args.results + directory + "/EvoDock_run_information.txt", 'r')
    lines = run_info.readlines()
    run_info.close()
    time = lines[1][10:].strip().split('.')[0]
    times.append(time)

# init content and worksheet
out_content = ""
workbook = xlsxwriter.Workbook(args.results + "ResultsReader_out_table.xlsx")
worksheet = workbook.add_worksheet()


for entry in summary:
    out_content += str(entry) + "\n"
out_content += str(times) + "\n"
sum_time = 0
for time in times:
    hours = int(time[0:1])
    mins = int(time[2:4])
    secs = int(time[5:7])
    sum_time += (60 * 60 * hours) + (60 * mins) + secs
average_time = sum_time / len(times)
average_hours = math.floor(average_time / float(60 * 60))
average_mins = math.floor((average_time - average_hours) / 60)
average_secs = int(average_time % 60)
if len(str(average_hours)) == 1:
    average_hours = "0" + str(average_hours)
if len(str(average_mins)) == 1:
    average_mins = "0" + str(average_mins)
if len(str(average_secs)) == 1:
    average_secs = "0" + str(average_secs)
out_content += str(average_hours) + ":" + str(average_mins) + ":" + str(average_secs)
workbook.close()

out_file = open(args.results + "ResultsReader_out.txt", 'w')
out_file.write(out_content)
out_file.close()