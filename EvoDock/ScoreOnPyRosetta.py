# region Imports and init
try:
    from pyrosetta import init

    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.core.scoring import ScoreType

    from pyrosetta.rosetta.protocols.ligand_docking import InterfaceScoreCalculator
    from pyrosetta.rosetta.std import vector_std_string

except ImportError as e:
    exit("ScoreOnPyRosetta: ImportError: " + str(e))

init()

# endregion

weight_mutagenesis = 0
weight_application = 1


def is_scoring_module():
    """
    This function is essential to verify, that it is a scoring module and is required, so it cannot be used as
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
    print("ScoreOnPyRosetta: Inputs are validated!")
    return


def calculate_fitness(specific_results, score_function_mutagenesis, score_function_application):
    """
    :return:
    """
    if specific_results is None:
        return 0

    mutagenesis_results = specific_results[0]
    application_results = specific_results[1]

    # get score functions
    inter_face = InterfaceScoreCalculator()
    inter_face.score_fxn(score_function_application)
    vec = vector_std_string("X")
    inter_face.chains(vec)

    lowest_energy = 0
    lowest_energy_pose = 0
    lowest_energy_sub_pose = 0
    outer_pose = 0
    for pose_list in application_results:
        inner_pose = 0
        for pose in pose_list:
            inter_face.apply(pose)
            energy = score_function_application(pose)
            print(energy)
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_pose = outer_pose
                lowest_energy_sub_pose = inner_pose
            inner_pose += 1
        outer_pose += 1

    print ("lowest energy pose: " + str(lowest_energy_pose) + " - " + str(lowest_energy_sub_pose) + " with energy:")
    print(lowest_energy)
    score_mutate = score_function_mutagenesis(mutagenesis_results[lowest_energy_pose]) * -1
    score_apply = score_function_application(application_results[lowest_energy_pose][lowest_energy_sub_pose]) * -1

    print("score mutagenesis: " + str(score_mutate) + " | weight: " + str(weight_mutagenesis))
    print("score application: " + str(score_apply) + " | weight: " + str(weight_application))
    score = (score_mutate * weight_mutagenesis) + (score_apply * weight_application)
    score = score / (weight_mutagenesis + weight_application)
    print("score final: " + str(score))

    return score
