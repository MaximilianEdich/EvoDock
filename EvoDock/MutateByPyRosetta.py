
VERSION = "0.21_01_10"


# region Imports and init
try:
    from pyrosetta import init, pose_from_file, get_fa_scorefxn, create_score_function, standard_packer_task
    from pyrosetta.toolbox import mutate_residue
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.relax import ClassicRelax
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMinMover
    from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover
    from pyrosetta.rosetta.protocols.simple_moves import SmallMover
    from pyrosetta.rosetta.protocols.simple_moves import ShearMover

    from pyrosetta.rosetta.core.kinematics import MoveMap

    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.core.scoring import ScoreType

except ImportError as e:
    exit("MutateByPyRosetta: ImportError: " + str(e))

init()

# endregion

# region Fixed values and variables
NONE = "NONE"
TRUE = "TRUE"
FALSE = "FALSE"

SET_SCORE_FUNCTION = "-score-function"
RELAX = "-relax"
PREP_ONLY = "-prepare-only"
EXTRA_RELAX = "-extra-relax"
RELAX_FAST = "FAST-RELAX"
RELAX_CLASSIC = "CLASSIC-RELAX"
RELAX_BOTH = "RELAX-BOTH"
RELAX_MAX_ITER = "-relax-max-iter"

SAVE_PDB = "-save-pdb"
PACK_RADIUS = "-pack-radius"
ROTAMER_MOVES_NUMBER = "-rotamer-moves-number"
BACKBONE_MOVES_NUMBER = "-backbone-moves-number"
MAKE_POSES = "-make-poses"
SET_KT = "-set-kT"
SET_N_MOVES = "-set-n-moves"
ROTAMER_MOVER = "-rotamer-mover"
BACKBONE_MOVER = "-backbone-mover"
BACKBONE_MOVER_CHI = "-backbone-mover-modify-chi"

# rotamer movers
ROTAMER_MOVERS = []
MOVER_ID_ROTAMER_TRIALS_MIN_MOVER = "ROTAMER-TRIALS-MIN-MOVER"
ROTAMER_MOVERS.append(MOVER_ID_ROTAMER_TRIALS_MIN_MOVER)
MOVER_ID_PACK_ROTAMERS_MOVER = "PACK-ROTAMERS-MOVER"
ROTAMER_MOVERS.append(MOVER_ID_PACK_ROTAMERS_MOVER)
MOVER_ID_ROTAMER_TRIALS_MOVER = "ROTAMER-TRIALS-MOVER"
ROTAMER_MOVERS.append(MOVER_ID_ROTAMER_TRIALS_MOVER)

BACKBONE_MOVERS = []
MOVER_ID_SMALL_MOVER = "SMALL-MOVER"
BACKBONE_MOVERS.append(MOVER_ID_SMALL_MOVER)
MOVER_ID_SHEAR_MOVER = "SHEAR-MOVER"
BACKBONE_MOVERS.append(MOVER_ID_SHEAR_MOVER)
MOVER_ID_MIN_MOVER = "MIN-MOVER"
BACKBONE_MOVERS.append(MOVER_ID_MIN_MOVER)

score_function = get_fa_scorefxn()

relax_mode = None
relax_max_iter = 0
prep_only = False
save_pdb = False
pack_radius = 8
number_of_rotamer_moves = 1
number_of_backbone_moves = 1
make_poses = 1
rotamer_mover_id = MOVER_ID_ROTAMER_TRIALS_MIN_MOVER
backbone_mover_id = NONE
backbone_mover_chi = True
kT = 1.0
n_moves = 5
tolerance = 0.01
extra_relax = None

# endregion


def is_mutation_module():
    """
    This function is essential to verify, that it is a mutation module and is required, so it cannot be used as
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
    try:
        pose = pose_from_file(protein_path)
    except RuntimeError:
        exit("ERROR in MutateByPyRosetta: Cannot load pose from " + str(protein_path))
    print("MutateByPyRosetta: Pose is validated!")
    print("MutateByPyRosetta: Inputs are validated!")
    return True


def print_documentation():
    # TODO
    print("MutateByPyRosetta:")
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
            exit("ERROR in MutateByPyRosetta: Score function file not found! " + path_name)
        except PermissionError:
            score_file_content = None
            exit("ERROR in MutateByPyRosetta: Cannot open score function file! Permission denied! " + path_name)
        # create empty score function
        score_function = ScoreFunction()
        if score_file_content[0].strip() != "":
            # first line is score function name, try to load it
            try:
                score_function = create_score_function(score_file_content[0].strip())
            except RuntimeError as e:
                print("ERROR in MutateByPyRosetta: Cannot set specified score function from line 1!")
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
            except AttributeError as error:
                print("ERROR in MutateByPyRosetta: Unknown ScoreType.")
                exit(error)
            except ValueError as error:
                print("ERROR in MutateByPyRosetta: Illegal value.")
                exit(error)
    else:
        # arg1 is score function name, try to load respective score function
        try:
            score_function = create_score_function(path_name)
        except RuntimeError as e:
            print("ERROR in MutateByPyRosetta: Cannot set specified score function!")
            exit(e)
    print("MutateByPyRosetta: specified score function:")
    print(score_function)
    return


def get_score_function():
    return score_function


def get_initial_amino_acids(protein_path, amino_acid_paths):
    """
    Get the single-letter amino acid code from the residues of interest.
    :param protein_path: path to a PDB.
    :param amino_acid_paths: list of residue ids.
    :return: list of single-letter aa codes, respectively to the input list amino_amino_acid_paths.
    """
    pose = pose_from_file(protein_path)
    result = []
    for aa_pos in amino_acid_paths:
        result.append(pose.residue(int(aa_pos)).name1())
    return result


def get_reformatted_amino_acids(protein_path, amino_acid_paths):
    pose = pose_from_file(protein_path)
    new_list = []
    for entry in amino_acid_paths:
        new_list.append(pose.pdb_info().pdb2pose(entry[0], int(entry[1])))

    return new_list


def parameter_handling(params):
    """
    Handle parameter which are written in the initial settings file and are specified as mutagenesis module parameter.
    :param params: list of parameter name and values, where the first element is always the name.
    :return: None.
    """
    # TODO explicit error with how to use it right
    if params[0] == BACKBONE_MOVER_CHI:
        try:
            global backbone_mover_chi
            if params[1] == TRUE:
                backbone_mover_chi = True
            else:
                backbone_mover_chi = False
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == SET_SCORE_FUNCTION:
        try:
            set_score_function(params[1])
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == PREP_ONLY:
        global prep_only
        try:
            if params[1] == TRUE:
                prep_only = True
            else:
                prep_only = False
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == RELAX:
        global relax_mode
        try:
            if params[1] == NONE:
                relax_mode = None
            elif params[1] == RELAX_FAST or params[1] == RELAX_CLASSIC or params[1] == RELAX_BOTH:
                relax_mode = params[1]
            else:
                exit("ERROR in MutateByPyRosetta: unknown argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == EXTRA_RELAX:
        global extra_relax
        try:
            if params[1] == NONE:
                extra_relax = None
            elif params[1] == RELAX_FAST or params[1] == RELAX_CLASSIC:
                extra_relax = params[1]
            else:
                exit("ERROR in MutateByPyRosetta: unknown argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == RELAX_MAX_ITER:
        global relax_max_iter
        try:
            number = int(params[1])
            if number < 1:
                exit("ERROR in MutateByPyRosetta: argument must be int >= 1: " + str(params))
            relax_max_iter = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type, int >= 1 expected: " + str(params))
    elif params[0] == SAVE_PDB:
        global save_pdb
        try:
            if params[1] == TRUE:
                save_pdb = True
            else:
                save_pdb = False
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == PACK_RADIUS:
        try:
            radius = float(params[1])
            if radius >= 0:
                global pack_radius
                pack_radius = radius
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type, float expected: " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == ROTAMER_MOVES_NUMBER:
        try:
            number = int(params[1])
            if number >= 0:
                global number_of_rotamer_moves
                number_of_rotamer_moves = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == BACKBONE_MOVES_NUMBER:
        try:
            number = int(params[1])
            if number >= 0:
                global number_of_backbone_moves
                number_of_backbone_moves = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == MAKE_POSES:
        try:
            number = int(params[1])
            if number >= 0:
                global make_poses
                make_poses = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == SET_KT:
        try:
            number = float(params[1])
            if number >= 0:
                global kT
                kT = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == SET_N_MOVES:
        try:
            number = int(params[1])
            if number >= 0:
                global n_moves
                n_moves = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == ROTAMER_MOVER:
        try:
            mover = str(params[1])
            global rotamer_mover_id
            rotamer_mover_id = mover
            # check argument
            valid_mover = False
            for mover in ROTAMER_MOVERS:
                if rotamer_mover_id == mover:
                    valid_mover = True
                    break
            if not valid_mover and not rotamer_mover_id == NONE:
                exit("Error in MutateByPyRosetta: Invalid rotamer mover defined!")
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    elif params[0] == BACKBONE_MOVER:
        try:
            mover = str(params[1])
            global backbone_mover_id
            backbone_mover_id = mover
            # check argument
            valid_mover = False
            for mover in BACKBONE_MOVERS:
                if backbone_mover_id == mover:
                    valid_mover = True
                    break
            if not valid_mover and not backbone_mover_id == NONE:
                exit("Error in MutateByPyRosetta: Invalid backbone mover defined!")
        except IndexError:
            exit("ERROR in MutateByPyRosetta: argument(s) missing at: " + str(params))
    else:
        exit("ERROR in MutateByPyRosetta: unknown flag: " + str(params))

    return


def preparation_result_path(protein_path, out_path):
    if relax_mode is not None:
        prot_split = protein_path.split('/')
        prot_name = prot_split[len(prot_split) - 1]
        return out_path + "/" + prot_name[:(len(prot_name) - 3)] + "relaxed.pdb"
    else:
        return protein_path


def prepare_files_for_tool(protein_path, out_path):
    """
    Prepares input files for the main task.
    :return: True, if preparation was successful.
    """
    pose = pose_from_file(protein_path)

    # set and check score function
    score_fxn = get_fa_scorefxn()

    # relax structure
    if relax_mode is not None:
        test_pose = Pose()
        test_pose.assign(pose)
        score1 = score_fxn(pose)
        content = "Scores:\n"
        content += "before: " + str(score1) + "\n"
        print("\nrelax structure...\n")
        if relax_mode == RELAX_FAST or relax_mode == RELAX_BOTH:
            f_relax = FastRelax()
            f_relax.set_scorefxn(score_fxn)
            if relax_max_iter > 0:
                f_relax.max_iter(relax_max_iter)
            f_relax.apply(pose)
            relax_out_path_name = preparation_result_path(protein_path, out_path)
            pose.dump_pdb(relax_out_path_name)
            print("\nrelax done! Saved relaxed strucutre as " + str(relax_out_path_name) + "\n")
            score2 = score_fxn(pose)
            content += "fast: " + str(score2) + "\n"
        if relax_mode == RELAX_CLASSIC or relax_mode == RELAX_BOTH:
            c_relax = ClassicRelax()
            c_relax.set_scorefxn(score_fxn)
            if relax_max_iter > 0:
                c_relax.max_iter(relax_max_iter)
            c_relax.apply(test_pose)
            relax_out_path_name = preparation_result_path(protein_path, out_path)
            pose.dump_pdb(relax_out_path_name)
            print("\nrelax done! Saved relaxed strucutre as " + str(relax_out_path_name) + "\n")
            score3 = score_fxn(test_pose)
            content += "classic: " + str(score3) + "\n"
        run_info_file = open(out_path + "/RELAX_information.txt", 'w')
        run_info_file.write(content)
        run_info_file.close()
        if prep_only:
            exit("Exit EvoDock. Preparation finished, saved relaxed structure(s).")

    return True


def generate_application_input(protein_path, out_path, amino_acid_paths, mutations):
    """
    """
    print(protein_path)

    # get score function
    score_fxn = score_function

    # load pose
    load_pose = pose_from_file(protein_path)

    # create copies
    poses = []
    for x in range(make_poses):
        new_pose = Pose()
        new_pose.assign(load_pose)
        poses.append(new_pose)

    # apply mutations
    index = 0
    for mut in mutations:
        if mut != '':
            for pose in poses:
                mutate_residue(pose, int(amino_acid_paths[index]), mut)
        index += 1

    # make residues of interest packable and calculate centers of change
    index = 0
    centers = []
    for mut in mutations:
        if mut != '':
            centers.append(load_pose.residue(int(amino_acid_paths[index])).nbr_atom_xyz())
        index += 1
    # make residues in range of these packable
    pack_list = []
    # rosetta is 1-indexed, so do in range(1, n+1)
    for i in range(1, load_pose.total_residue() + 1):
        for c in centers:
            if c.distance_squared(load_pose.residue(i).nbr_atom_xyz()) <= pack_radius * pack_radius:
                # residue is in range, pack residue
                if (i not in pack_list
                        and not load_pose.residue(i).name() == "CYS:disulfide" and not load_pose.residue(
                            i).is_ligand()):
                    pack_list.append(i)
    print(pack_list)
    print(len(pack_list))

    # init move map
    move_map = MoveMap()
    for p in pack_list:
        move_map.set_bb(p, True)
        if backbone_mover_chi:
            move_map.set_chi(p, True)

    print(move_map)
    # for each pose
    mutagenesis_output = []
    suffix = 0
    for pose in poses:
        suffix += 1

        # init PackerTask
        packer_task = standard_packer_task(pose)
        packer_task.restrict_to_repacking()
        packer_task.temporarily_fix_everything()

        # apply restrictions
        for p in pack_list:
            packer_task.temporarily_set_pack_residue(p, True)

        # init Rotamer Mover and apply repacking
        backbone_mover = None

        if backbone_mover_id == MOVER_ID_SMALL_MOVER:
            backbone_mover = SmallMover(move_map, kT, n_moves)
        elif backbone_mover_id == MOVER_ID_SHEAR_MOVER:
            backbone_mover = ShearMover(move_map, kT, n_moves)
        elif backbone_mover_id == MOVER_ID_MIN_MOVER:
            backbone_mover = MinMover()
            backbone_mover.score_function(score_fxn)
            backbone_mover.movemap(move_map)
            backbone_mover.tolerance(tolerance)

        # init Rotamer Mover and apply repacking
        print(packer_task)
        rotamer_mover = None
        if rotamer_mover_id == MOVER_ID_ROTAMER_TRIALS_MIN_MOVER:
            rotamer_mover = RotamerTrialsMinMover(score_fxn, packer_task)
        elif rotamer_mover_id == MOVER_ID_PACK_ROTAMERS_MOVER:
            rotamer_mover = PackRotamersMover(score_fxn, packer_task)
        elif rotamer_mover_id == MOVER_ID_ROTAMER_TRIALS_MOVER:
            rotamer_mover = RotamerTrialsMover(score_fxn, packer_task)

        # apply bb and rotamer
        print("Score " + str(mutations) + str(suffix) + " : " + str(score_fxn(pose)))
        if backbone_mover is not None:
            for x in range(number_of_backbone_moves):
                backbone_mover.apply(pose)
                print("Score " + str(mutations) + str(suffix) + " : " + str(score_fxn(pose)))
        print("Score " + str(mutations) + str(suffix) + " : " + str(score_fxn(pose)))
        if rotamer_mover is not None:
            for x in range(number_of_rotamer_moves):
                rotamer_mover.apply(pose)
                print("Score " + str(mutations) + str(suffix) + " : " + str(score_fxn(pose)))

        if extra_relax is not None:
            print("\nrelax structure...\n")
            if extra_relax == RELAX_FAST:
                f_relax = FastRelax()
                f_relax.set_scorefxn(score_fxn)
                if relax_max_iter > 0:
                    f_relax.max_iter(relax_max_iter)
                f_relax.apply(pose)
            if extra_relax == RELAX_CLASSIC:
                c_relax = ClassicRelax()
                c_relax.set_scorefxn(score_fxn)
                if relax_max_iter > 0:
                    c_relax.max_iter(relax_max_iter)
                c_relax.apply(pose)

        # save pdb
        if save_pdb:
            pose.dump_pdb(out_path + str(suffix) + ".pdb")
        mutagenesis_output.append(pose)

    return mutagenesis_output


def get_specific_output_from_pdb(mutant_out_path, use_existent_mutate_out_path):
    """
    Reads in paths to PDBs to recreate PyRosetta Pose objects, save these again as PDBs in new output folder
    (if desired) and
    :param mutant_out_path:
    :param use_existent_mutate_out_path:
    :return:
    """
    # extract mutant name (AAs combination) from out_path
    path_split = str(mutant_out_path).split('/')
    prefix = path_split[len(path_split) - 1]
    paths = []
    suffix = 0
    while True:
        suffix += 1
        try:
            path = use_existent_mutate_out_path + prefix + str(suffix) + ".pdb"
            # TODO more efficient way
            test_file = open(path, 'r')
            test_file.close()
            paths.append(path)
        except FileNotFoundError:
            try:
                path = use_existent_mutate_out_path + prefix + "_mut" + str(suffix) + ".pdb"
                # TODO more efficient way
                test_file = open(path, 'r')
                test_file.close()
                paths.append(path)
            except FileNotFoundError:
                try:
                    path = use_existent_mutate_out_path + prefix + "_apl" + str(suffix) + ".pdb"
                    # TODO more efficient way
                    test_file = open(path, 'r')
                    test_file.close()
                    paths.append(path)
                except FileNotFoundError:
                    break

    if not paths:
        # file loading failed
        return

    mutagenesis_output = []
    suffix = 0
    for path in paths:
        suffix += 1
        load_pose = pose_from_file(path)
        mutagenesis_output.append(load_pose)
        if save_pdb:
            load_pose.dump_pdb(mutant_out_path + str(suffix) + ".pdb")

    return mutagenesis_output


def get_compatibility_out():
    """
    Get object which is included in module specific output, in this case a pose.
    :return: A PyRosetta Pose Object.
    """
    pose = Pose()
    return pose


def save_pdb_file(pose, out_path):
    """
    Save a given pose as pdb. Used to save poses later after evaluation.
    :param pose: Input PyRosetta pose object.
    :param out_path: Outpath.
    :return: None.
    """
    pose.dump_pdb(out_path + ".pdb")
    return
