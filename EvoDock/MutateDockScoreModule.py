"""
The Mutation, Docking and Scoring Module for EvoDock Version 0.1.

This is responsible for performing all necessary steps to modify a pdb file which correspond to a mutation of the
protein. Next the protein will be used in docking and get a score.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import MutateByPyMOL


def generate_docking_input(targetIndividual, runOutPath, aminoAcidPaths):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change.
    :param aminoAcidPaths:
    :param runOutPath:
    :param targetIndividual:
    :return:
    """
    # define outputfolder for this mutant
    outPath = runOutPath + "/Mutant"
    for x in targetIndividual[0]:
        outPath += "_" + str(x)

    MutateByPyMOL.generate_docking_input(targetIndividual, outPath, aminoAcidPaths)

    return


def perform_docking():
    """

    :return:
    """

    return


def calculate_fitness_score(targetIndividual):
    """

    :return:
    """
    score = 0
    for aa in targetIndividual[0]:
        score += ord(aa)
    score = score / len(targetIndividual[0])

    return score


def get_score(targetIndividual, outPath, aminoAcidPaths):
    """

    :param aminoAcidPaths:
    :param outPath:
    :param targetIndividual:
    :return:
    """

    # generate new pdb as input for docking
    generate_docking_input(targetIndividual, outPath, aminoAcidPaths)

    # use fresh input for protein-ligand docking
    perform_docking()

    # use results to calculate fitness score
    score = calculate_fitness_score(targetIndividual)

    return score
