"""
The Mutation, Docking and Scoring Module for EvoDock Version 0.1.

This is responsible for performing all necessary steps to modify a pdb file which correspond to a mutation of the
protein. Next the protein will be used in docking and get a score.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import MutateTemplate
import DockTemplate
import ScoreTemplate

import MutateByPyMOL


def generate_docking_input(mutations, run_out_path, protein_code, amino_acid_paths):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in new generated pdb files.
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param run_out_path: The path leading to the output files.
    :return: None. The generated files are of interest.
    """

    MutateTemplate.generate_docking_input(mutations, run_out_path, amino_acid_paths, protein_code)

    return


def perform_docking(mutations, run_out_path, amino_acid_paths, protein_code):
    """

    :return:
    """

    DockTemplate.perform_docking(mutations, run_out_path, amino_acid_paths, protein_code)

    return


def calculate_fitness_score(target_individual, mutations, run_out_path, amino_acid_paths, protein_code):
    """
    Uses data from mutagenesis and the docking to calculate a fitness score.
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :return: the calculated score as the individuals fitness.
    """

    ScoreTemplate.calculate_fitness(mutations, run_out_path, amino_acid_paths, protein_code)

    score = 0
    for aa in target_individual[0]:
        score += ord(aa)
    score = score / len(target_individual[0])

    return score


def get_score(target_individual, original_individual, out_path, protein_code, amino_acid_paths):
    """
    Use the given input to generate a mutant and perform a ligand docking, both via external software tools.
    Output information of both tools is then used to determine a score that is usable as a fitness score for evolution.
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :param original_individual: The original protein input representing the very first individual.
    :param out_path: The path leading to the output files.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.

    :return: The calculated score as the individuals fitness, that was written into the score value of the individual.
    """
    # define output-folder for this mutant
    run_out_path = out_path + "/Mutant"
    for x in target_individual[0]:
        run_out_path += "_" + str(x)

    # identify real mutations
    mutations = []
    for k in range(len(original_individual[0])):
        if original_individual[0][k] == target_individual[0][k]:
            # no change in this position compared to target
            mutations.append("")
        else:
            # mutation detected
            mutations.append((target_individual[0][k]))

    # generate new pdb as input for docking
    generate_docking_input(mutations, run_out_path, protein_code, amino_acid_paths)

    # use fresh input for protein-ligand docking
    perform_docking(mutations, run_out_path, amino_acid_paths, protein_code)

    # use results to calculate fitness score
    score = calculate_fitness_score(target_individual, mutations, run_out_path, amino_acid_paths, protein_code)

    return score
