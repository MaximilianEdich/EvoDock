"""
The Mutation, Application and Scoring Module for EvoDock Version 0.1.

This is responsible for performing all necessary steps to modify a pdb file which correspond to a mutation of the
protein. Next the protein will be used in a task specific application and get a score.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

VERSION = "0.21_01_10"

# region Imports and fixed vars
import importlib

mutate_mod = None
apply_mod = None
score_mod = None
fold_mod = None

MODULE_PARAM_MUTATE = "-module-param-mutate"
MODULE_PARAM_APPLY = "-module-param-apply"
MODULE_PARAM_SCORE = "-module-param-score"
MODULE_PARAM_FOLD = "-module-param-fold"

poses = []

# endregion


# region Initialization validation, and preparation of modules
def init(mutate, apply, score, fold):
    """
    Try to import desired modules. Modules are not validated in this step.
    :return: None.
    """
    mutate_mod_name = mutate
    apply_mod_name = apply
    score_mod_name = score
    fold_mod_name = fold
    global mutate_mod
    global apply_mod
    global score_mod
    global fold_mod

    try:
        mutate_mod = importlib.import_module(mutate_mod_name)
    except ImportError as e:
        mutate_mod = None
        exit(
            "Error: Import of \"" + mutate_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
                                                      "This was the import of your specified MutagenesisModule. Make sure it is in the correct Folder. "
                                                      "Check the documentation for more information.\nError Message: " + str(e))
    try:
        apply_mod = importlib.import_module(apply_mod_name)
    except ImportError as e:
        apply_mod = None
        exit(
            "Error: Import of \"" + apply_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
                                                     "This was the import of your specified ApplicationModule. Make sure it is in the correct Folder. "
                                                     "Check the documentation for more information.\nError Message: " + str(e))
    try:
        score_mod = importlib.import_module(score_mod_name)
    except ImportError as e:
        score_mod = None
        exit(
            "Error: Import of \"" + score_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
                                                     "This was the import of your specified EvaluationModule. Make sure it is in the correct Folder. "
                                                     "Check the documentation for more information.\nError Message: " + str(e))
    try:
        fold_mod = importlib.import_module(fold_mod_name)
    except ImportError as e:
        fold_mod = None
        exit(
            "Error: Import of \"" + fold_mod_name + "\" failed. "
                                                    "This was the import of your specified FoldingModule. Make sure it is in the correct Folder. "
                                                    "Check the documentation for more information.\nError Message: " + str(e))
    except ValueError as e:
        # no fold module specified, proceed with fold_mod set to None
        fold_mod = None
    return


def check_imported_modules():
    """
    Check if the imported modules are of the correct type (if mutagenesis module is indeed a mutagenesis module and
    so on).
    :return: True, if all specified modules are valid.
    """
    try:
        mutate_mod.is_mutation_module()
    except AttributeError:
        exit("Error, specified mutation module is not a mutation module!")
    try:
        apply_mod.is_application_module()
    except AttributeError:
        exit("Error, specified application module is not an application module!")
    try:
        score_mod.is_scoring_module()
    except AttributeError:
        exit("Error, specified evaluation module is not an evaluation module!")
    if fold_mod is not None:
        try:
            fold_mod.is_folding_module()
        except AttributeError:
            exit("Error, specified folding module is not a folding module!")
    return True


def print_documentations():
    """
    Print the documentation of each module.
    """
    try:
        mutate_mod.print_documentation()
    except AttributeError:
        print("Error, mutagenesis module is missing documentation print!")
    try:
        apply_mod.print_documentation()
    except AttributeError:
        print("Error, application module is missing documentation print!")
    try:
        score_mod.print_documentation()
    except AttributeError:
        print("Error, evaluation module is missing documentation print!")
    if fold_mod is not None:
        try:
            fold_mod.print_documentation()
        except AttributeError:
            print("Error, folding module is missing documentation print!")
    return


def validate_module_data(protein_path, out_path, skip_application):
    """
    Passes the validation task to the single modules and checks, if all settings and inputs for each module
    are validated.
    :param skip_application: If True, aplication module wont be valided.
    :param protein_path: Path to the input PDB file, which represents the wild type protein.
    :param out_path: Path to the output folder of this run.
    :return: True, if all modules have validated their settings and input.
    """
    mutate_mod.validate_data(protein_path, out_path)
    if not skip_application:
        apply_mod.validate_data(protein_path, out_path)
    score_mod.validate_data(protein_path, out_path)
    if fold_mod is not None:
        fold_mod.validate_data(protein_path, out_path)
    print("All modules validated!\n")
    return True


def check_compatibility_mutate_apply():
    mut_out = mutate_mod.get_compatibility_out()
    return apply_mod.get_compatibility_check(mut_out)


def preparation_result_path(protein_path, out_path):
    """

    :param protein_path:
    :param out_path:
    :return:
    """
    return mutate_mod.preparation_result_path(protein_path, out_path)


def prepare_tool(protein_path, out_path, skip_application):
    """

    :param out_path:
    :param protein_path:
    :return:
    """
    mutate_mod.prepare_files_for_tool(protein_path, out_path)
    if not skip_application:
        apply_mod.prepare_files_for_tool(protein_path, out_path)

    return


def handle_module_params(params):
    """
    Pass the module specific parameters from initial settings file to the modules.
    :param params: String with parameter.
    :return: Specific string values from the modules in a list.
    """
    if params[0] == MODULE_PARAM_MUTATE:
        return mutate_mod.parameter_handling(params[1:])
    if params[0] == MODULE_PARAM_APPLY:
        return apply_mod.parameter_handling(params[1:])
    if params[0] == MODULE_PARAM_SCORE:
        return score_mod.parameter_handling(params[1:])
    if params[0] == MODULE_PARAM_FOLD:
        if fold_mod is not None:
            return fold_mod.parameter_handling(params[1:])
        else:
            print("WARNING in MutateApplyScoreModule: Received Fold-parameters, but has not folding module.")
    exit("ERROR in MutateApplyScoreModule, handle_module_params: unknown parameter identifier: " + params)
    return


def get_ref_protein_amino_acids(protein_path, amino_acid_paths):
    """
    Get AAs from the reference protein corresponding to the specified res-ids.
    :param protein_path:
    :param amino_acid_paths:
    :return:
    """
    return mutate_mod.get_initial_amino_acids(protein_path, amino_acid_paths)


def get_reformatted_amino_acids(protein_path, amino_acid_paths):
    return mutate_mod.get_reformatted_amino_acids(protein_path, amino_acid_paths)


# endregion

# region Pipeline input-output flow
class MutateApplyScore:
    def __init__(self):
        self.original_individual = None
        self.out_path = None
        self.protein_path = None
        self.amino_acid_paths = None
        self.use_specific_mutate_out = True
        self.use_existent_mutate_out_path = None
        self.skip_application = False
        self.keep_improvements = True
        self.symmetry_paths = None
        return

    def set_values(self, original_individual, out_path, protein_path, amino_acid_paths, use_specific_mutate_out,
                   use_existent_mutate_out_path, skip_application, keep_improvements, symmetry_paths):
        self.original_individual = original_individual
        self.out_path = out_path
        self.protein_path = protein_path
        self.amino_acid_paths = amino_acid_paths
        self.use_specific_mutate_out = use_specific_mutate_out
        self.use_existent_mutate_out_path = use_existent_mutate_out_path
        self.skip_application = skip_application
        self.keep_improvements = keep_improvements
        self.symmetry_paths = symmetry_paths
        return


def generate_application_input(mutations, mutant_out_path, protein_path, amino_acid_paths, fold_instead_mutate,
                               use_existent_mutate_out_path):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in new generated pdb files.
    :param fold_instead_mutate:
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param mutant_out_path: The path leading to the output files.
    :return: None. The generated files are of interest.
    """
    if use_existent_mutate_out_path is not None:
        # if alternative file available, load it, otherwise, do regular mutagenesis
        results = mutate_mod.get_specific_output_from_pdb(mutant_out_path, use_existent_mutate_out_path)
        if results is not None:
            return results

    if fold_instead_mutate and fold_mod is not None:
        return fold_mod.generate_application_input(protein_path, mutant_out_path, amino_acid_paths, mutations)
    else:
        return mutate_mod.generate_application_input(protein_path, mutant_out_path, amino_acid_paths, mutations)


def run_application(mutant_out_path, application_input, use_specific_mutate_out, use_existent_mutate_out_path):
    """

    :return:
    """
    if use_specific_mutate_out:
        # simply use output from mutagenesis as application input
        return apply_mod.perform_application(application_input, mutant_out_path)
    else:
        if use_existent_mutate_out_path is not None:
            # generate new application input by loading pdb files from specified path.
            path_split = str(mutant_out_path).split('/')
            prefix = path_split[len(path_split) - 1]
            load_input = []
            suffix = 0
            while True:
                suffix += 1
                try:
                    path = use_existent_mutate_out_path + prefix + str(suffix) + ".pdb"
                    # TODO more efficient way
                    test_file = open(path, 'r')
                    test_file.close()
                    load_input.append(path)
                except FileNotFoundError:
                    try:
                        path = use_existent_mutate_out_path + prefix + "_mut" + str(suffix) + ".pdb"
                        # TODO more efficient way
                        test_file = open(path, 'r')
                        test_file.close()
                        load_input.append(path)
                    except FileNotFoundError:
                        try:
                            path = use_existent_mutate_out_path + prefix + "_apl" + str(suffix) + ".pdb"
                            # TODO more efficient way
                            test_file = open(path, 'r')
                            test_file.close()
                            load_input.append(path)
                        except FileNotFoundError:
                            break

            if load_input:
                # if there is input to load, use it, otherwise use regularly created PDBs from mutagenesis
                return apply_mod.perform_application_with_pdb(load_input, mutant_out_path)

        # generate new application input by loading pdb files of this run.
        new_input = []
        suffix = 0
        for _ in application_input:
            suffix += 1
            mutant_variant = mutant_out_path + str(suffix) + ".pdb"
            new_input.append(mutant_variant)

        return apply_mod.perform_application_with_pdb(new_input, mutant_out_path)


def calculate_fitness_score(specific_results):
    """
    Let the Evaluation module calculate the final fitness score, which in turn is used to determine the mutant's
    fitness in the context of the current optimization problem.
    :param specific_results: A list containing exactly two elements. The first one is a list, containing all results
    from the Mutagenesis module. The second one is a list, containing all results from the Application module.
    """
    # fetch used score functions from Mutagenesis and Application modules
    sfx_mut = mutate_mod.get_score_function()
    sfx_apply = apply_mod.get_score_function()
    score = score_mod.calculate_fitness_score(specific_results, sfx_mut, sfx_apply)

    return score


def get_fitness_score(target_individual, mas: MutateApplyScore, fold_instead_mutate):
    """
    Use the given input to generate a mutant and perform a task specific application, both via external software tools.
    Output information of both tools is then used to determine a score that is usable as a fitness score for evolution.
    :param fold_instead_mutate:
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.

    :return: The calculated score as the individuals fitness, that was written into the score value of the individual.
    """

    # catch error
    if fold_instead_mutate and fold_mod is None:
        fold_instead_mutate = False

    # define output-folder for this mutant
    mutant_out_path = mas.out_path + "/Mutant"
    for x in target_individual[0]:
        mutant_out_path += "_" + str(x)

    # identify real mutations
    mutations = []
    for k in range(len(mas.original_individual[0])):
        if mas.original_individual[0][k] == target_individual[0][k]:
            # no change in this position compared to target
            mutations.append("")
        else:
            # mutation detected
            mutations.append((target_individual[0][k]))

    # do error prevention
    if not mas.use_specific_mutate_out:
        # set PDB saving to true, since these are the input for the next step
        mutate_mod.parameter_handling([mutate_mod.SAVE_PDB, mutate_mod.TRUE])

    # set up amino acid paths
    aa_paths = mas.amino_acid_paths
    if mas.symmetry_paths is not None:
        mutations = mutations + mutations
        aa_paths = aa_paths + mas.symmetry_paths

    # generate mutagenesis output as input for the application and keep results for scoring
    # specific results stores in element 0 the mutagenesis results and in element 1 the application results
    # these results are specific to the module which generate them and may cannot be interpreted by the next pipeline
    # modules. Alternatively, each module generates PDBs in the form of:
    # Mutant_A_B_C1.pdb, where A, B, C, are placeholders for the amino acid on the respective residue and the number
    # at the end shows which pose of this mutant it is.
    specific_results = [[], []]

    specific_results[0] = generate_application_input(mutations, mutant_out_path, mas.protein_path,
                                                     aa_paths, fold_instead_mutate,
                                                     mas.use_existent_mutate_out_path)
    # perform application
    if not mas.skip_application:
        specific_results[1] = run_application(mutant_out_path, specific_results[0], mas.use_specific_mutate_out,
                                              mas.use_existent_mutate_out_path)

    # evaluate all scores for final fitness
    evaluation_results = calculate_fitness_score(specific_results)
    fitness_score = evaluation_results[0]
    mut_pose = evaluation_results[1]
    apply_pose = evaluation_results[2]

    # if pose is potentially saved later, keep it
    save_bests = True
    if save_bests:
        if mut_pose is not None:
            mutate_mod.save_pdb_file(mut_pose, mutant_out_path + "_mut")
        if apply_pose is not None:
            apply_mod.save_pdb_file(apply_pose, mutant_out_path + "_apl")
    elif mas.keep_improvements:
        pass
        # poses.append([target_individual, mut_pose, apply_pose])

    return fitness_score


def save_improvements(population, mas: MutateApplyScore):
    """
    Takes a population of all improvements and looks up their kept poses to save them via mutate/apply module.
    Afterwards, the pose list is cleared.
    :param mas:
    :param population: Already evaluated population, which contains only improvements.
    :return: None.
    """
    print("test1")
    print(poses)
    for individual in population:
        # for each individual look up the kept poses and save it, if found
        for result in poses:
            if individual[0] == result[0]:
                print(individual[0])
                print(result[0])
                # define output-folder for this mutant
                mutant_out_path = mas.out_path + "/Mutant"
                for x in individual[0]:
                    mutant_out_path += "_" + str(x)

                print(mutant_out_path)
                mutate_mod.save_pdb_file(result[1], mutant_out_path)
                apply_mod.save_pdb_file(result[2], mutant_out_path)
                break
    print("clear")
    poses.clear()
    return

# endregion
