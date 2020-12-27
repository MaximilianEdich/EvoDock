

VERSION = "0.20_12"


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
    exit("ScoreOnPyRosetta: ImportError: " + str(e))

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

    print("ScoreOnPyRosetta: Inputs are validated!")
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
            exit("ERROR in ScoreOnPyRosetta: Argument(s) missing! " + str(params))
        except ValueError:
            exit("ERROR in ScoreOnPyRosetta: Weight value has to be float! " + str(params))
    elif params[0] == WEIGHT_APPLICATION:
        global weight_application
        try:
            weight_application = float(params[1])
        except IndexError:
            exit("ERROR in ScoreOnPyRosetta: Argument(s) missing! " + str(params))
        except ValueError:
            exit("ERROR in ScoreOnPyRosetta: Weight value has to be float! " + str(params))
    else:
        exit("ERROR in ScoreOnPyRosetta: unknown argument(s): " + str(params))


def calculate_fitness(specific_results, score_function_mutagenesis, score_function_application):
    """
    :return:
    """
    if specific_results is None:
        return 0

    mutagenesis_results = specific_results[0]
    application_results = specific_results[1]

    if mutagenesis_results == []:
        return 0

    if application_results == [[]] or application_results == []:
        lowest_energy = 0
        lowest_energy_pose_index = 0
        outer_pose = 0
        for pose in mutagenesis_results:
            energy = score_function_mutagenesis(pose)
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_pose_index = outer_pose
            outer_pose += 1
        print("lowest energy pose: " + str(lowest_energy_pose_index) + " with energy:")
        print(lowest_energy)
        score_mutate = score_function_mutagenesis(mutagenesis_results[lowest_energy_pose_index]) * -1
        score = (score_mutate * weight_mutagenesis)
        print("score final: " + str(score))

        results = [score, mutagenesis_results[lowest_energy_pose_index], None]
        return results

    # get score functions
    '''
    inter_face = InterfaceScoreCalculator()
    inter_face.score_fxn(score_function_application)
    svec = vector_std_string()
    svec.append('A')
    svec.append('X')
    inter_face.chains(svec)'''

    start_from = StartFrom()
    start_from.chain('X')
    vec = xyzVector_double_t(1000, 1000, 1000)
    start_from.add_coords(vec)

    lowest_energy = 0
    lowest_energy_pose_index = 0
    lowest_energy_sub_pose_index = 0
    outer_pose = 0
    for pose_list in application_results:
        inner_pose = 0
        for pose in pose_list:
            #inter_face.apply(pose)
            #energy = score_function_application(pose)

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
