

VERSION = "0.21_01_03"


# region Imports and init
try:
    from pyrosetta import init, pose_from_file, get_fa_scorefxn, create_score_function
    from pyrosetta.rosetta.core.pose import Pose

    from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser

    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.core.scoring import ScoreType

except ImportError as e:
    exit("DockByPyRosetta: ImportError: " + str(e))

init()

# endregion

# region Fixed values and variables
NONE = "NONE"
TRUE = "TRUE"
FALSE = "FALSE"

SET_SCORE_FUNCTION = "-score-function"
XML_PROT_PATH = "-xml-protocol-path"
XML_SUBSTITUTION = "-xml-substitution"
SAVE_PDB = "-save-pdb"
MAKE_POSES = "-make-poses"

score_function = None
docking_protocol = None

xml_protocol_path = ""
xml_subst_list = []
save_pdb = True
make_poses = 20
fixed_files_path = None

# endregion


def is_application_module():
    """
    This function is essential to verify, that it is a application module and is required, so it cannot be
    used as a module of a different type.
    """
    return True


def validate_data(protein_path, out_path):
    """
    Checks, if all settings and inputs are valid.
    :param protein_path: Path to the input PDB file, which represents the wild type protein.
    :param out_path: Path to the output folder of this run.
    :return: True, all checks are valid.
    """
    # check score function
    if score_function is None:
        exit("ERROR in DockByPyRosetta: Score function file was not specified. "
             "Use this in the initial settings file:\n"
             "'-module-param-apply " + SET_SCORE_FUNCTION + " <path/name>', where <path/name> is either a path to the "
                                                            ".txt file containing all details or a name of a Rosetta "
                                                            "score function.")
    # check xml protocol
    # check if path is not empty
    global xml_protocol_path
    if xml_protocol_path == "":
        exit("ERROR in DockByPyRosetta: XML file was not specified. "
             "Do this with '-module-param-apply " + XML_PROT_PATH + " <path>'"
             " in the initial settings file.")
    # try opening the file
    xml_content = ""
    try:
        xml = open(xml_protocol_path, 'r')
        xml_content = xml.read()
        xml.close()
    except FileNotFoundError as e:
        print("ERROR in DockByPyRosetta: XML file not found: " + xml_protocol_path)
        exit(e)
    except PermissionError as e:
        print("ERROR in DockByPyRosetta: Permission denied for XML file : " + xml_protocol_path)
        exit(e)
    # do substitutions, if required
    if len(xml_subst_list) > 0:
        for key_value in xml_subst_list:
            xml_content = xml_content.replace(key_value[0], key_value[1])
        name = xml_protocol_path.split('/')
        xml_protocol_path = out_path + "/" + name[len(name) - 1]
        print("Generated new XML protocol at: " + xml_protocol_path)
        xml = open(xml_protocol_path, 'w')
        xml.write(xml_content)
        xml.close()
    # parse (new) protocol
    xml_object = RosettaScriptsParser()
    try:
        global docking_protocol
        docking_protocol = xml_object.generate_mover(xml_protocol_path)
    except RuntimeError as e:
        print("ERROR in DockByPyRosetta: XML protocol cannot be parsed: " + xml_protocol_path)
        exit(e)
    print("DockByPyRosetta: XML-protocol is validated!")

    print("DockByPyRosetta: Inputs are validated!")
    return


def print_documentation():
    # TODO
    print("DockByPyRosetta:")
    print("Version: " + str(VERSION))
    return


def set_score_function(path_name):
    """
    Set the Rosetta score function of this module.
    :param path_name: Name of one of Rosetta's score functions or path to a .txt file with score function base and
    weights.
    :return: None.
    """
    global score_function
    if path_name[-4:] == ".txt":
        # open file and create new score function
        try:
            score_file = open(path_name, 'r')
            score_file_content = score_file.readlines()
            score_file.close()
        except FileNotFoundError:
            score_file_content = None
            exit("ERROR in DockByPyRosetta: Score function file not found! " + path_name)
        except PermissionError:
            score_file_content = None
            exit("ERROR in DockByPyRosetta: Cannot open score function file! Permission denied! " + path_name)
        # create empty score function
        score_function = ScoreFunction()
        if score_file_content[0].strip() != "":
            # first line is score function name, try to load it
            try:
                score_function = create_score_function(score_file_content[0].strip())
            except RuntimeError as e:
                print("ERROR in DockByPyRosetta: Cannot set specified score function from line 1!")
                print("How to use the score function text file:")
                print("1st line: Set base score function from Rosetta. If none is desired, keep 1st line blank")
                print("all other lines: write lines in the format '<weight> <value>', separated by a space.")
                exit(e)
        # recreate or set weights
        for line in score_file_content[1:]:
            line = line.strip()
            key_value = line.split(' ')
            try:
                score_function.set_weight(getattr(ScoreType, key_value[0]), float(key_value[1]))
            except AttributeError as e:
                print("ERROR in DockByPyRosetta: Unknown ScoreType.")
                exit(e)
            except ValueError as e:
                print("ERROR in DockByPyRosetta: Illegal value.")
                exit(e)
    else:
        # arg1 is score function name, try to load respective score function
        try:
            score_function = create_score_function(path_name)
        except RuntimeError as e:
            print("ERROR in MutateByPyRosetta: Cannot set specified score function!")
            exit(e)
    print("DockByPyRosetta: specified score function:")
    print(score_function)
    return


def get_score_function():
    return score_function


def parameter_handling(params):
    """
        Handle parameter which are written in the initial settings file and are specified as application module
        parameter.
        :param params: list of parameter name and values, where the first element is always the name.
        :return: None.
        """
    # TODO explicit error with how to use it right
    if params[0] == SET_SCORE_FUNCTION:
        try:
            set_score_function(params[1])
        except IndexError:
            exit("ERROR in DockByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == SAVE_PDB:
        global save_pdb
        if params[1] == TRUE:
            save_pdb = True
        else:
            save_pdb = False
    elif params[0] == MAKE_POSES:
        try:
            number = int(params[1])
            if number >= 0:
                global make_poses
                make_poses = number
        except ValueError as e:
            print("ERROR in DockByPyRosetta: wrong value type of argument(s): " + str(params))
            exit(e)
    elif params[0] == XML_PROT_PATH:
        global xml_protocol_path
        xml_protocol_path = params[1]
        if xml_protocol_path[-4:] != ".xml":
            print("WARNING in DockByPyRosetta: specified XML file may not is a true XML file! Wrong file-ending!")
    elif params[0] == XML_SUBSTITUTION:
        try:
            xml_subst_list.append([params[1], params[2]])
        except IndexError as e:
            print("ERROR in DockByPyRosetta: '" + XML_SUBSTITUTION + " <k> <v>' requires two arguments, where "
                                                                     "<k> is the key, a string that in the XML that is"
                                                                     "substituted by <v>, the value.")
            exit(e)
    else:
        exit("ERROR in DockByPyRosetta: unknown flag: " + str(params))

    return


def make_xml(out_path):

    return


def prepare_files_for_tool(protein_path, out_path):
    """
    Prepares input files for the main task.
    :return: True, if preparation was successful.
    """

    make_xml(out_path)

    return


def perform_application(application_input, out_path):
    """
    Use the PyRosetta Pose objects from the application_input and perform on each the specified ligand docking protocol.
    Each mutant variant from the input may leads to several sub variants as a result from several ligand docking.
    :param application_input: list with PyRosetta poses.
    :param out_path: path of the mutant, including output folder path and file name of this mutant without suffix.
    :return: list of lists, where each of these lists contains poses that are results from the ligand docking.
    Example of the structure for 2 input poses and 3 results per pose:
    application_input: [p1, p2] -> results: [[p11, p12, p13], [p21, p22, p23]], where each p_ is a pose.
    """
    print("Start DockByPyRosetta")
    docking_results = []
    # get score function
    score_fxn = score_function

    # create copies per input pose
    poses = []
    for pose in application_input:
        pose_list = []
        for x in range(make_poses):
            new_pose = Pose()
            new_pose.assign(pose)
            pose_list.append(new_pose)
        poses.append(pose_list.copy())

    # get protocol, apply it to all poses
    protocol = docking_protocol
    suffix = 0
    for pose_list in poses:
        suffix_2 = 0
        suffix += 1
        pose_results = []
        for pose in pose_list:
            suffix_2 += 1
            print("Score before docking: " + str(score_fxn(pose)))
            protocol.apply(pose)
            print("Score after docking: " + str(score_fxn(pose)))
            if save_pdb:
                pose.dump_pdb(out_path + str(suffix) + "_" + str(suffix_2) + "_docking.pdb")
            pose_results.append(pose)
        docking_results.append(pose_results)
    return docking_results


def perform_application_with_pdb(pdb_input, out_path):
    """
    Reads in paths to PDBs to recreate PyRosetta Pose objects and uses these as an input for the function
    perform_application(application_input, out_path).
    :param pdb_input: list of paths leading to PDB files.
    :param out_path: path of the mutant, including output folder path and file name of this mutant without suffix.
    :return: list of lists, where each of these lists contains poses that are results from the ligand docking.
    Example of the structure for 2 input poses and 3 results per pose:
    application_input: [p1, p2] -> results: [[p11, p12, p13], [p21, p22, p23]], where each p_ is a pose.
    """
    pose_input = []
    for path in pdb_input:
        try:
            pose = pose_from_file(path)
            pose_input.append(pose)
        except RuntimeError:
            exit("ERROR in DockByPyRosetta: Cannot load pose from PDB input, check if output is generated.")

    return perform_application(pose_input, out_path)


def get_compatibility_check(input_obj):
    """
    Check, if the given input is a pose and therefore this module is compatible with pose input/output.
    :param input_obj: Input, a pose is expected.
    :return: True, if input is expected type, otherwise False.
    """
    new_pose = Pose()
    try:
        new_pose.assign(input_obj)
    except ValueError:
        return False

    return True


def save_pdb_file(pose, out_path):
    """
    Save a given pose as pdb. Used to save poses later after evaluation.
    :param pose: Input PyRosetta pose object.
    :param out_path: Outpath.
    :return: None.
    """
    pose.dump_pdb(out_path + ".pdb")
    return
