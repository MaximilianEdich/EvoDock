from pyrosetta import init, toolbox, pose_from_pdb
from pyrosetta.toolbox import cleanATOM, pose_from_rcsb

init()

identifier = "1OGX"
load = pose_from_rcsb(identifier)

pose = pose_from_pdb(identifier.capitalize() + ".pdb")
# init ligand
