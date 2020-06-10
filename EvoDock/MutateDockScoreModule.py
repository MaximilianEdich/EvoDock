"""
The Mutation, Docking and Scoring Module for EvoDock Version 0.1.

This is responsible for performing all necessary steps to modify a pdb file which correspond to a mutation of the
protein. Next the protein will be used in docking and get a score.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""


def generate_docking_input(targetIndividual):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change.
    :param targetIndividual:
    :return:
    """

    return


def perform_docking():
    """

    :return:
    """

    return


def calculate_fitness_score():
    """

    :return:
    """

    return 0


def get_score(targetIndividual):
    """

    :param targetIndividual:
    :return:
    """

    # generate new pdb as input for docking
    generate_docking_input(targetIndividual)

    # use fresh input for protein-ligand docking
    perform_docking()

    # use results to calculate fitness score
    score = calculate_fitness_score()

    return score
