

VERSION = "0.21_03_05"


# region Imports and init
try:
    from pyrosetta import init

    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.core.scoring import ScoreType

    from pyrosetta.rosetta.protocols.ligand_docking import InterfaceScoreCalculator
    from pyrosetta.rosetta.protocols.ligand_docking import StartFrom
    from pyrosetta.rosetta.std import vector_std_string
    from pyrosetta.rosetta.numeric import xyzVector_double_t

except ImportError as e:
    exit("ScoreByPyRosetta: ImportError: " + str(e))

init()

# endregion

# region Fixed values and variables
WEIGHT_MUTAGENESIS = "-weight-mutate"
WEIGHT_APPLICATION = "-weight-apply"

weight_mutagenesis = 1
weight_application = 1

# endregion


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

    print("ScoreByPyRosetta: Inputs are validated!")
    return


def print_documentation():
    """
    Print the documentation of this pipeline module.
    :return: None.
    """
    print("ScoreByPyRosetta:")
    print("Version: " + str(VERSION))
    print("Module-specific parameters:")
    print(WEIGHT_MUTAGENESIS + "<float>")
    print("\tSet the weight on the energy from results of Mutagenesis Module.")
    print(WEIGHT_APPLICATION + "<float>")
    print("\tSet the weight on the energy from results of Application Module.")
    return


def parameter_handling(params):
    """
    Handle parameter which are written in the initial settings file and are specified as score module parameter.
    :param params: list of parameter name and values, where the first element is always the name.
    :return: None.
    """
    if params[0] == WEIGHT_MUTAGENESIS:
        global weight_mutagenesis
        try:
            weight_mutagenesis = float(params[1])
        except IndexError:
            exit("ERROR in ScoreByPyRosetta: Argument(s) missing! " + str(params))
        except ValueError:
            exit("ERROR in ScoreByPyRosetta: Weight value has to be float! " + str(params))
    elif params[0] == WEIGHT_APPLICATION:
        global weight_application
        try:
            weight_application = float(params[1])
        except IndexError:
            exit("ERROR in ScoreByPyRosetta: Argument(s) missing! " + str(params))
        except ValueError:
            exit("ERROR in ScoreByPyRosetta: Weight value has to be float! " + str(params))
    else:
        exit("ERROR in ScoreByPyRosetta: unknown argument(s): " + str(params))


def calculate_fitness_score(specific_results, score_function_mutagenesis, score_function_application):
    """
    Use the results from Mutagenesis and Application modules to calculate the fitness determining final score.
    :param specific_results: A list containing exactly two elements. The first one is a list, containing all results
    from the Mutagenesis module. The second one is a list, containing all results from the Application module.
    :param score_function_mutagenesis: The PyRosetta score function used in Mutagenesis module.
    :param score_function_application: The PyRosetta score function used in Application module.
    :return:
    """
    if specific_results is None:
        return 0

    if specific_results[0] == []:
        # mutagenesis results are empty, nothing to evaluate
        return 0

    # fetch results for easier access
    mutagenesis_results = specific_results[0]
    application_results = specific_results[1]

    # initiate StartFrom Mover
    start_from = StartFrom()
    start_from.chain('X')
    vec = xyzVector_double_t(1000, 1000, 1000)
    start_from.add_coords(vec)

    # if application was skipped
    if application_results == [[]] or application_results == []:
        lowest_energy = 0
        lowest_energy_pose_index = 0
        outer_pose = 0
        for pose in mutagenesis_results:
            energy_complex = score_function_application(pose)
            pose_seperated = Pose()
            pose_seperated.assign(pose)
            start_from.apply(pose_seperated)
            energy_separated = score_function_application(pose_seperated)
            energy = energy_complex - energy_separated

            print(energy)
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_pose_index = outer_pose
            outer_pose += 1
        print("lowest energy pose: " + str(lowest_energy_pose_index) + " with energy:")
        print(lowest_energy)
        score_mutate = score_function_mutagenesis(mutagenesis_results[lowest_energy_pose_index]) * -1
        score_apply = lowest_energy * -1
        score = (score_mutate * weight_mutagenesis) + (score_apply * weight_application)
        print("score final: " + str(score))

        results = [score, mutagenesis_results[lowest_energy_pose_index], None]
        return results

    # application was not skipped, evaluate all results
    lowest_energy = 0
    lowest_energy_pose_index = 0
    lowest_energy_sub_pose_index = 0
    outer_pose = 0
    for pose_list in application_results:
        inner_pose = 0
        for pose in pose_list:
            energy_complex = score_function_application(pose)
            pose_seperated = Pose()
            pose_seperated.assign(pose)
            start_from.apply(pose_seperated)
            energy_separated = score_function_application(pose_seperated)
            energy = energy_complex - energy_separated

            print(energy)
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_pose_index = outer_pose
                lowest_energy_sub_pose_index = inner_pose
            inner_pose += 1
        outer_pose += 1

    print ("lowest energy pose: "
           "" + str(lowest_energy_pose_index) + " - " + str(lowest_energy_sub_pose_index) + " with energy:")
    print(lowest_energy)
    score_mutate = score_function_mutagenesis(mutagenesis_results[lowest_energy_pose_index]) * -1
    score_apply = lowest_energy * -1

    print("score mutagenesis: " + str(score_mutate) + " | weight: " + str(weight_mutagenesis))
    print("score application: " + str(score_apply) + " | weight: " + str(weight_application))
    score = (score_mutate * weight_mutagenesis) + (score_apply * weight_application)
    print("score final: " + str(score))

    results = [score, mutagenesis_results[lowest_energy_pose_index],
               application_results[lowest_energy_pose_index][lowest_energy_sub_pose_index]]
    return results
