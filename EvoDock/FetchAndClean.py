import argparse

from pyrosetta import init, toolbox, pose_from_pdb
from pyrosetta.toolbox import cleanATOM, pose_from_rcsb

init()

parser = argparse.ArgumentParser(description="...")
parser.add_argument("-p", "--prot", type=str, required=True,
                    help="Protein id")
args = parser.parse_args()

identifier = str(args.prot)
load = pose_from_rcsb(identifier)

pose = pose_from_pdb(identifier.capitalize() + ".pdb")
# init ligand
