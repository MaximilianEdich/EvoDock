"""
The Mutation, Docking and Scoring Module for EvoDock Version 0.1.

This is responsible for performing all necessary steps to modify a pdb file which correspond to a mutation of the
protein. Next the protein will be used in docking and get a score.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import importlib

global mutate_mod
global apply_mod
global score_mod
global fold_mod

MODULE_PARAM_MUTATE = "-module-param-mutate"
MODULE_PARAM_APPLY = "-module-param-apply"
MODULE_PARAM_SCORE = "-module-param-score"
MODULE_PARAM_FOLD = "-module-param-fold"


def init(mutate, apply, score, fold):
    """
    Import desired modules.
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
        exit("Error: Import of \"" + mutate_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
             "This was the import of your specified MutagenesisModule. Make sure it is in the correct Folder. "
             "Check the documentation for more information.\nError Message: " + str(e))
    try:
        apply_mod = importlib.import_module(apply_mod_name)
    except ImportError as e:
        apply_mod = None
        exit("Error: Import of \"" + apply_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
             "This was the import of your specified ApplicationModule. Make sure it is in the correct Folder. "
             "Check the documentation for more information.\nError Message: " + str(e))
    try:
        score_mod = importlib.import_module(score_mod_name)
    except ImportError as e:
        score_mod = None
        exit("Error: Import of \"" + score_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
             "This was the import of your specified EvaluationModule. Make sure it is in the correct Folder. "
             "Check the documentation for more information.\nError Message: " + str(e))
    try:
        fold_mod = importlib.import_module(fold_mod_name)
    except ImportError as e:
        fold_mod = None
        exit("Error: Import of \"" + fold_mod_name + "\" failed. Make sure to provide this module, since it is essential. "
             "This was the import of your specified FoldingModule. Make sure it is in the correct Folder. "
             "Check the documentation for more information.\nError Message: " + str(e))
    return


def check_imported_modules():
    try:
        mutate_mod.is_mutation_module()
    except AttributeError:
        exit("Error, specified mutation module is not a mutation module!")
    try:
        apply_mod.is_application_module()
    except AttributeError:
        exit("Error, specified docking module is not an application module!")
    try:
        score_mod.is_scoring_module()
    except AttributeError:
        exit("Error, specified scoring module is not a scoring module!")
    try:
        fold_mod.is_folding_module()
    except AttributeError:
        exit("Error, specified folding module is not a folding module!")

    return


class MutateDockScore:
    def __init__(self):
        self.original_individual = None
        self.out_path = None
        self.protein_path = None
        self.amino_acid_paths = None
        return

    def set_values(self, original_individual, out_path, protein_path, amino_acid_paths):
        self.original_individual = original_individual
        self.out_path = out_path
        self.protein_path = protein_path
        self.amino_acid_paths = amino_acid_paths
        return


def preparation_result_path(protein_path, out_path):
    """

    :param protein_path:
    :param out_path:
    :return:
    """
    return mutate_mod.preparation_result_path(protein_path, out_path)


def prepare_tool(protein_path, out_path):
    """

    :param out_path:
    :param protein_path:
    :return:
    """
    mutate_mod.prepare_files_for_tool(protein_path, out_path)

    return


def handle_module_params(params):
    if params[0] == MODULE_PARAM_MUTATE:
        return mutate_mod.parameter_handling(params[1:])

    return


def generate_docking_input(mutations, run_out_path, protein_path, amino_acid_paths, fold_instead_mutate):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in new generated pdb files.
    :param fold_instead_mutate:
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param run_out_path: The path leading to the output files.
    :return: None. The generated files are of interest.
    """
    results = []
    if fold_instead_mutate:
        results = fold_mod.generate_docking_input(protein_path, amino_acid_paths, mutations, run_out_path)
    else:
        results = mutate_mod.generate_docking_input(protein_path, amino_acid_paths, mutations, run_out_path)

    return results


def perform_docking(mutations, run_out_path, amino_acid_paths, protein_path, docking_input):
    """

    :return:
    """

    apply_mod.perform_docking(docking_input, run_out_path)
    #dock_mod.perform_docking(mutations, run_out_path, amino_acid_paths, protein_path)

    return


def calculate_fitness_score(mutations, run_out_path, amino_acid_paths, protein_path, docking_results):
    """
    Uses data from mutagenesis and the docking to calculate a fitness score.
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :return: the calculated score as the individuals fitness.
    """

    score = score_mod.calculate_fitness(docking_results)
    #score = score_mod.calculate_fitness(mutations, run_out_path, amino_acid_paths, protein_path)

    return score


def get_score(target_individual, mds: MutateDockScore, fold_instead_mutate):
    """
    Use the given input to generate a mutant and perform a ligand docking, both via external software tools.
    Output information of both tools is then used to determine a score that is usable as a fitness score for evolution.
    :param fold_instead_mutate:
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.

    :return: The calculated score as the individuals fitness, that was written into the score value of the individual.
    """

    # define output-folder for this mutant
    run_out_path = mds.out_path + "/Mutant"
    for x in target_individual[0]:
        run_out_path += "_" + str(x)

    # identify real mutations
    mutations = []
    for k in range(len(mds.original_individual[0])):
        if mds.original_individual[0][k] == target_individual[0][k]:
            # no change in this position compared to target
            mutations.append("")
        else:
            # mutation detected
            mutations.append((target_individual[0][k]))

    # generate new pdb as input for docking
    docking_input = generate_docking_input(mutations, run_out_path, mds.protein_path, mds.amino_acid_paths, fold_instead_mutate)

    # use fresh input for protein-ligand docking
    #docking_results = perform_docking(mutations, run_out_path, mds.amino_acid_paths, mds.protein_path, docking_input)

    # use results to calculate fitness score
    score = calculate_fitness_score(mutations, run_out_path, mds.amino_acid_paths, mds.protein_path, docking_input)

    return score
