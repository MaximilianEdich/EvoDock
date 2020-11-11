
# try importing PyRosetta modules
try:
    from pyrosetta import init, toolbox, pose_from_file, dump_pdb, get_fa_scorefxn, standard_packer_task
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.core.pack.task import PackerTask
    from pyrosetta.rosetta.core.pack.task import operation

    from pyrosetta.rosetta.protocols.ligand_docking import LigandArea
    from pyrosetta.rosetta.protocols.ligand_docking import InterfaceBuilder
    from pyrosetta.rosetta.protocols.ligand_docking import MoveMapBuilder
    from pyrosetta.rosetta.protocols.ligand_docking import StartFrom
    from pyrosetta.rosetta.protocols.ligand_docking import Translate
    from pyrosetta.rosetta.protocols.ligand_docking import Rotate
    from pyrosetta.rosetta.protocols.ligand_docking import SlideTogether
    from pyrosetta.rosetta.protocols.ligand_docking import HighResDocker
    from pyrosetta.rosetta.protocols.ligand_docking import FinalMinimizer
    from pyrosetta.rosetta.protocols.ligand_docking import InterfaceScoreCalculator
    from pyrosetta.rosetta.protocols.rosetta_scripts import ParsedProtocol

    from pyrosetta.rosetta.utility.tag import Tag

except ImportError as e:
    exit("ImportError in the module \"pyrosetta\": " + str(e))

init()


def is_application_module():
    """
    This function is essential to verify, that it is a is_application_module module and is required, so it cannot be
    used as a module of a different type.
    """
    return True


def perform_docking(docking_input, out_path):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in a new pdb file. During this process a pml file is
    generated, containing the PyMOL script that performs the mutagenesis.
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param out_path: The path leading to the output files specific to the mutation.
    :param mutations: List of mutations relative to the original protein. An empty string represents no mutation while
    any substitution is represented by the given single letter code of the amino acid.
    :return: None. The generated files are of interest.
    """
    print("DOCKING!")
    docking_results = []

    # get score function
    score_fxn = get_fa_scorefxn()

    # ligand areas
    docking_sidechain_x = LigandArea()
    docking_sidechain_x.chain_ = 'X'
    docking_sidechain_x.cutoff_ = 6.0
    docking_sidechain_x.add_nbr_radius_ = True
    docking_sidechain_x.all_atom_mode_ = True
    docking_sidechain_x.minimize_ligand_ = 10

    final_sidechain_x = LigandArea()
    final_sidechain_x.chain_ = 'X'
    final_sidechain_x.cutoff_ = 6.0
    final_sidechain_x.add_nbr_radius_ = True
    final_sidechain_x.all_atom_mode_ = True

    final_backbone_x = LigandArea()
    final_backbone_x.chain_ = 'X'
    final_backbone_x.cutoff_ = 7.0
    final_backbone_x.add_nbr_radius_ = False
    final_backbone_x.all_atom_mode_ = True
    final_backbone_x.Calpha_restraints_ = 0.3

    # interface builders


    # move map builders


    # mover



    return docking_results