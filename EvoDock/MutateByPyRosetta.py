# try importing PyRosetta modules
try:
    from pyrosetta import init, toolbox, pose_from_file, dump_pdb, get_fa_scorefxn, standard_packer_task
    from pyrosetta.toolbox import mutate_residue, cleanATOM
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.core.pack.task import PackerTask
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.relax import ClassicRelax
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMinMover
    from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover
    from pyrosetta.rosetta.protocols.simple_moves import SmallMover
    from pyrosetta.rosetta.protocols.simple_moves import ShearMover

    from pyrosetta.rosetta.core.kinematics import MoveMap

except ImportError as e:
    exit("ImportError in the module \"pyrosetta\": " + str(e))

init()

NONE = "NONE"
TRUE = "TRUE"
FALSE = "FALSE"

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

RELAX = "-relax"
RELAX_FAST = "FAST-RELAX"
RELAX_CLASSIC = "CLASSIC-RELAX"
RELAX_BOTH = "RELAX-BOTH"
SAVE_PDB = "-save-pdb"
PACK_RADIUS = "-pack-radius"
ROTAMER_MOVES_NUMBER = "-rotamer-moves-number"
BACKBONE_MOVES_NUMBER = "-backbone-moves-number"
MAKE_POSES = "-make-poses"
ROTAMER_MOVER = "-rotamer-mover"
BACKBONE_MOVER = "-backbone-mover"

relax_mode = None
save_pdb = True
pack_radius = 8
number_of_rotamer_moves = 1
number_of_backbone_moves = 1
make_poses = 1
rotamer_mover_id = MOVER_ID_ROTAMER_TRIALS_MIN_MOVER
backbone_mover_id = NONE


def parameter_handling(params):
    if params[0] == RELAX:
        global relax_mode
        if params[1] == NONE:
            relax_mode = None
        elif params[1] == RELAX_FAST or params[1] == RELAX_CLASSIC or params[1] == RELAX_BOTH:
            relax_mode = params[1]
        else:
            exit("ERROR in MutateByPyRosetta: unknown argument(s): " + str(params))
    elif params[0] == SAVE_PDB:
        global save_pdb
        if params[1] == TRUE:
            save_pdb = True
        else:
            save_pdb = False
    elif params[0] == PACK_RADIUS:
        try:
            radius = float(params[1])
            if radius >= 0:
                global pack_radius
                pack_radius = radius
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
    elif params[0] == ROTAMER_MOVES_NUMBER:
        try:
            number = int(params[1])
            if number >= 0:
                global number_of_rotamer_moves
                number_of_rotamer_moves = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
    elif params[0] == BACKBONE_MOVES_NUMBER:
        try:
            number = int(params[1])
            if number >= 0:
                global number_of_backbone_moves
                number_of_backbone_moves = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
    elif params[0] == MAKE_POSES:
        try:
            number = int(params[1])
            if number >= 0:
                global make_poses
                make_poses = number
        except ValueError:
            exit("ERROR in MutateByPyRosetta: wrong value type of argument(s): " + str(params))
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


def is_mutation_module():
    """
    This function is essential to verify, that it is a mutation module and is required, so it cannot be used as
    a module of a different type.
    """
    return True


def preparation_result_path(protein_path, out_path):
    prot_split = protein_path.split('/')
    prot_name = prot_split[len(prot_split) - 1]
    if not relax_mode is None:
        return out_path + "/" + prot_name[:(len(prot_name) - 3)] + "relaxed.pdb"
    else:
        return prot_name


def prepare_files_for_tool(protein_path, out_path):
    """
    Prepares input files for the main task.
    :return: True, if preparation was successful.
    """

    pose = pose_from_file(protein_path)

    # set score function
    score_fxn = get_fa_scorefxn()

    # relax structure
    if not relax_mode is None:
        test_pose = Pose()
        test_pose.assign(pose)
        score1 = score_fxn(pose)
        content = "Scores:\n"
        content += "before: " + str(score1) + "\n"
        print("\nrelax structure...\n")
        if relax_mode == RELAX_FAST or relax_mode == RELAX_BOTH:
            f_relax = FastRelax()
            f_relax.set_scorefxn(score_fxn)
            f_relax.apply(pose)
            relax_out_path_name = preparation_result_path(protein_path, out_path)
            pose.dump_pdb(relax_out_path_name)
            print("\nrelax done! Saved relaxed strucutre as " + str(relax_out_path_name) + "\n")
            score2 = score_fxn(pose)
            content += "fast: " + str(score2) + "\n"
        if relax_mode == RELAX_CLASSIC or relax_mode == RELAX_BOTH:
            c_relax = ClassicRelax()
            c_relax.set_scorefxn(score_fxn)
            c_relax.apply(test_pose)
            score3 = score_fxn(test_pose)
            content += "classic: " + str(score3) + "\n"
        run_info_file = open(out_path + "/MUTATE_information.txt", 'w')
        run_info_file.write(content)
        run_info_file.close()

    return True


def generate_docking_input(protein_path, amino_acid_paths, mutations, out_path):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in a new pdb file. During this process a pml file is
    generated, containing the PyMOL script that performs the mutagenesis.
    :param protein_path: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param out_path: The path leading to the output files specific to the mutation.
    :param mutations: List of mutations relative to the original protein. An empty string represents no mutation while
    any substitution is represented by the given single letter code of the amino acid.
    :return: None. The generated files are of interest.
    """
    print(protein_path)

    # get score function
    score_fxn = get_fa_scorefxn()

    # load pose
    is_original = True
    load_pose = pose_from_file(protein_path)

    poses = []
    for x in range(make_poses):
        new_pose = Pose()
        new_pose.assign(load_pose)
        poses.append(new_pose)

    index = 0
    for mut in mutations:
        if mut != '':
            is_original = False
            for pose in poses:
                mutate_residue(pose, int(amino_acid_paths[index]), mut)
        index += 1

    if is_original:
        # skip time consuming calculations, just save copy to work with
        if save_pdb:
            load_pose.dump_pdb(out_path + ".pdb")
        docking_input = [load_pose]
        return docking_input

    # make residues of interest packable and calculate centers of change
    index = 0
    centers = []
    for mut in mutations:
        if mut != '':
            centers.append(load_pose.residue(int(amino_acid_paths[index])).nbr_atom_xyz())
        index += 1
    # make residues in range of these packable
    pack_list = []
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
    print(move_map)
    # for each pose
    docking_input = []
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
        kT = 1.0
        n_moves = 5
        if backbone_mover_id == MOVER_ID_SMALL_MOVER:
            backbone_mover = SmallMover(move_map, 1.0, 5)
        elif backbone_mover_id == MOVER_ID_SHEAR_MOVER:
            backbone_mover = ShearMover(move_map, 1.0, 5)

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
        print(score_fxn(pose))
        if backbone_mover is not None:
            print("BB!")
            for x in range(number_of_backbone_moves):
                backbone_mover.apply(pose)
        print(score_fxn(pose))
        if rotamer_mover is not None:
            print("ROTAMER!")
            for x in range(number_of_rotamer_moves):
                rotamer_mover.apply(pose)
                print(score_fxn(pose))

        # save pdb
        if save_pdb:
            if make_poses == 1:
                suffix = ""
            pose.dump_pdb(out_path + str(suffix) + ".pdb")
        docking_input.append(pose)


    return docking_input
