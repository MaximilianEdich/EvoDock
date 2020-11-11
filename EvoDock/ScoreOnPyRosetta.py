
# try importing PyRosetta modules
try:
    from pyrosetta import init, get_fa_scorefxn
    from pyrosetta.rosetta.core.pose import Pose

except ImportError as e:
    exit("ImportError in the module \"pyrosetta\": " + str(e))

init()


def is_scoring_module():
    """
    This function is essential to verify, that it is a scoring module and is required, so it cannot be used as
    a module of a different type.
    """
    return True


def calculate_fitness(docking_results):
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

    if docking_results is None:
        return 0

    # get score function
    score_fxn = get_fa_scorefxn()
    print(docking_results)

    score_sum = 0
    for pose in docking_results:
        score_sum += score_fxn(pose)
    score = (score_sum / float(len(docking_results))) * -1

    return score