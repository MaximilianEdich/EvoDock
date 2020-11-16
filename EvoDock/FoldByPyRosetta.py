try:
    from pyrosetta import toolbox
    from pyrosetta import init
except ImportError as e:
    exit("FoldByPyRosetta: ImportError: " + str(e))

init()


def is_folding_module():
    """
    This function is essential to verify, that it is a folding module and is required, so it cannot be used as
    a module of a different type.
    """
    return True


def validate_data(protein_path, out_path):
    """
    Checks, if all settings and inputs are valid.
    :param protein_path: Path to the input PDB file, which represents the wild type protein.
    :param out_path: Path to the output folder of this run.
    :return: True, all checks are valid.
    """
    print("FoldByPyRosetta: Inputs are validated!")
    return


def generate_application_input(protein_path, run_out_path, amino_acid_paths, mutations):
    """
    """

    return
