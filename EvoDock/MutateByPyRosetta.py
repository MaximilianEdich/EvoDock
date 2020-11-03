
# try importing PyRosetta modules
try:
    from pyrosetta import init, toolbox, pose_from_pdb, dump_pdb, get_fa_scorefxn, standard_packer_task
    from pyrosetta.toolbox import mutate_residue, cleanATOM
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.core.pack.task import PackerTask
    from pyrosetta.rosetta.core.pack.task import operation
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMinMover

except ImportError as e:
    exit("ImportError in the module \"pyrosetta\": " + str(e))

init()

#TODO use parameter to perform relaxation or not, if done before
do_relax = False
number_of_rotamer_moves = 3
pack_radius = 8


def is_mutation_module():
    """
    This function is essential to verify, that it is a mutation module and is required, so it cannot be used as
    a module of a different type.
    """
    return True


def preparation_result_path(protein_path, out_path):
    if do_relax:
        return out_path + "/" + protein_path[:(len(protein_path) - 3)] + "relaxed.pdb"
    else:
        return protein_path


def prepare_files_for_tool(protein_path, out_path):
    """
    Prepares input files for the main task.
    :return: None.
    """
    pose = pose_from_pdb(protein_path)

    if not verify_pose(pose):
        return False

    # set score function
    score_fxn = get_fa_scorefxn()

    # relax structure TODO ligand constraints
    if do_relax:
        print("\nrelax structure...\n")
        f_relax = FastRelax()
        f_relax.set_scorefxn(score_fxn)
        f_relax.apply(pose)
        name = preparation_result_path(protein_path, out_path)
        pose.dump_pdb(name)
        print("\nrelax done! Saved relaxed strucutre as " + str(name) + "\n")



    return True


def verify_pose(pose):
    """
   Check, if the input is valid. Do this before calling a mutation function.
   :return:
   """
    if not pose.is_fullatom():
        return False

    return True


def verify_input(protein_path):
    """
    Check, if the input is valid. Do this before calling a mutation function.
    :return:
    """

    pose = pose_from_pdb(protein_path)

    return verify_pose(pose)


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

    # get score function
    score_fxn = get_fa_scorefxn()

    # load pose
    is_original = True
    pose = pose_from_pdb(protein_path)
    test_pose = Pose()
    test_pose.assign(pose)

    index = 0
    for mut in mutations:
        if mut != '':
            is_original = False
            mutate_residue(pose, int(amino_acid_paths[index]), mut)
        index += 1

    if is_original:
        # skip time consuming calculations, just save copy to work with
        pose.dump_pdb(out_path + ".pdb")
        return

    # init PackerTask
    packer_task = standard_packer_task(pose)
    packer_task.restrict_to_repacking()
    packer_task.temporarily_fix_everything()

    # make residues of interest packable and calculate centers of change TODO do it only on first run and save list
    index = 0
    centers = []
    for mut in mutations:
        if mut != '':
            centers.append(pose.residue(int(amino_acid_paths[index])).nbr_atom_xyz())
        index += 1
    # make residues in range of these packable
    pack_list = []
    for i in range(1, pose.total_residue() + 1):
        for c in centers:
            if c.distance_squared(test_pose.residue(i).nbr_atom_xyz()) <= pack_radius * pack_radius:
                # residue is in range, pack residue
                if (i not in pack_list
                        and not pose.residue(i).name() == "CYS:disulfide" and not pose.residue(i).is_ligand()):
                    pack_list.append(i)
    print(pack_list)
    print(len(pack_list))
    for p in pack_list:
        packer_task.temporarily_set_pack_residue(p, True)

    # init RotamerTrialsMinMover and apply repacking
    print(packer_task)
    rtmm = RotamerTrialsMinMover(score_fxn, packer_task)
    print(score_fxn(pose))
    for x in range(number_of_rotamer_moves):
        rtmm.apply(pose)
        print(score_fxn(pose))

    # save pdb
    pose.dump_pdb(out_path + ".pdb")

    return
