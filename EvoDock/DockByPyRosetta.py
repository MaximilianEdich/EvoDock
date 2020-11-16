
# try importing PyRosetta modules
try:
    from pyrosetta import init, pose_from_file, get_fa_scorefxn, create_score_function
    from pyrosetta.rosetta.core.pose import Pose

    from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser

    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.core.scoring import ScoreType

except ImportError as e:
    exit("DockByPyRosetta: ImportError: " + str(e))

init()

NONE = "NONE"
TRUE = "TRUE"
FALSE = "FALSE"

XML_PROT_PATH = "-xml-protocol-path"
SAVE_PDB = "-save-pdb"
MAKE_POSES = "-make-poses"

score_function = None
docking_protocol = None

xml_protocol_path = ""
save_pdb = True
make_poses = 1
fixed_files_path = None


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
    if xml_protocol_path == "":
        exit("ERROR in DockByPyRosetta: XML file was not specified. "
             "Do this with '-module-param-apply " + XML_PROT_PATH + " <path>'"
             " in the initial settings file.")
    try:
        xml = open(xml_protocol_path, 'r')
        xml.close()
    except FileNotFoundError as e:
        print("ERROR in DockByPyRosetta: XML file not found: " + xml_protocol_path)
        exit(e)
    except PermissionError as e:
        print("ERROR in DockByPyRosetta: Permission denied for XML file : " + xml_protocol_path)
        exit(e)
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


def set_score_function(path_name, arg1_is_path=True):
    global score_function
    if arg1_is_path:
        # open file and create new score function
        score_function = ScoreFunction()
        score_function.set_weight(ScoreType.fa_intra_rep, 0.004)
    else:
        # arg1 is string, try to load respective score function
        score_function = create_score_function(path_name)

    return


def get_score_function():
    return score_function


def parameter_handling(params):
    if params[0] == SAVE_PDB:
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

    # score function
    global score_function
    score_function = create_score_function("ligand")
    score_function.set_weight(ScoreType.fa_intra_rep, 0.004)
    score_function.set_weight(ScoreType.fa_elec, 0.42)
    score_function.set_weight(ScoreType.hbond_bb_sc, 1.3)
    score_function.set_weight(ScoreType.hbond_sc, 1.3)
    score_function.set_weight(ScoreType.rama, 0.2)


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

    poses = []
    for pose in application_input:
        pose_list = []
        for x in range(make_poses):
            new_pose = Pose()
            new_pose.assign(pose)
            pose_list.append(new_pose)
        poses.append(pose_list.copy())

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
