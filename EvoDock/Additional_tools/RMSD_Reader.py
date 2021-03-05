import argparse
import os
import ast
import math
import xlsxwriter

from pyrosetta import init, pose_from_file
from pyrosetta.rosetta.core.scoring import all_atom_rmsd

init()

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-r", "--reference", type=str, required=True,
                    help="Reference structure.")
parser.add_argument("-f", "--folder", type=str, required=True,
                    help="Folder with other structures.")
parser.add_argument("-s", "--suffix", type=str, required=True,
                    help="Suffix of output.")
parser.add_argument("-n", "--number_of_poses", type=str, required=True,
                    help="Number of poses in results.")
args = parser.parse_args()

files = os.listdir(args.folder)
pdbs = []
for file in files:
    if file[-3:] == "pdb":
        pdbs.append(file)

# reference pose and compare
pose_ref = pose_from_file(args.reference)

# calculate RMSDs
results = []
for file in pdbs:
    pose = pose_from_file(args.folder + file)
    try:
        results.append([file, all_atom_rmsd(pose_ref, pose)])
    except RuntimeError:
        print("Error, try vice versa...")
        results.append([file, all_atom_rmsd(pose, pose_ref)])
        print("... success!")

# output
results_text = ""
workbook = xlsxwriter.Workbook(args.folder + "RMDSs_" + args.suffix + ".xlsx")
worksheet = workbook.add_worksheet()
row = 0
for entry in results:
    row += 1
    results_text += entry[0] + "\t" + str(entry[1]) + "\n"
    worksheet.write(row, 1, str(entry[0]))
    worksheet.write(row, 2, float(entry[1]))
for x in range(int(args.number_of_poses)):
    row = 2 + (21 * x)
    row2 = 2 + (21 * (x + 1))
    worksheet.write_formula(row, 4, "=MITTELWERT(C" + str(row + 1) + ":C" + str(row2 - 1) + ")")
workbook.close()

results_file = open(args.folder + "RMSDs_" + args.suffix + ".txt", 'w')
results_file.write(results_text)
results_file.close()