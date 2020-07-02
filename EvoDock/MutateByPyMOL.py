"""


created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import os

FULL_TEMPLATE = ""


def get_3_letter_code(singleLetter):
    """
    Converts a single letter amino acid code to a three letter amino acid code.
    :param singleLetter: Single letter of the one letter amino acid code.
    :return: The corresponding three letter amino acid code, if entered correct letter, otherwise an empty string.
    """

    if singleLetter == 'A':
        return "ALA"
    if singleLetter == 'R':
        return "ARG"
    if singleLetter == 'N':
        return "ASN"
    if singleLetter == 'D':
        return "ASP"
    if singleLetter == 'C':
        return "CYS"
    if singleLetter == 'Q':
        return "GLN"
    if singleLetter == 'E':
        return "GLU"
    if singleLetter == 'G':
        return "GLY"
    if singleLetter == 'H':
        return "HIS"
    if singleLetter == 'I':
        return "ILE"
    if singleLetter == 'L':
        return "LEU"
    if singleLetter == 'K':
        return "LYS"
    if singleLetter == 'M':
        return "MET"
    if singleLetter == 'F':
        return "PHE"
    if singleLetter == 'P':
        return "PRO"
    if singleLetter == 'S':
        return "SER"
    if singleLetter == 'T':
        return "THR"
    if singleLetter == 'W':
        return "TRP"
    if singleLetter == 'Y':
        return "TYR"
    if singleLetter == 'V':
        return "VAL"

    return ""


def get_template(aminoAcidPaths):
    """
    :param aminoAcidPaths:
    :return:
    """

    # create template on the first run
    global FULL_TEMPLATE
    if FULL_TEMPLATE == "":
        # load head part
        template = open("mutagenesis_template_head.pml", 'r')
        templateContent = ""
        for line in template.readlines():
            templateContent += line
        template.close()
        # load loop part
        template = open("mutagenesis_template_loop.pml", 'r')
        loopContent = template.readlines()
        for path in aminoAcidPaths:
            for line in loopContent:
                templateContent += line
            templateContent = templateContent.replace("RESIDUE_N", str(path))
            templateContent = templateContent.replace("RESIDUE_SL", str(path)[1:].replace("/", "_"))
        template.close()
        # load tail part
        template = open("mutagenesis_template_tail.pml", 'r')
        for line in template.readlines():
            templateContent += line
        template.close()
        # insert general information
        templateContent = templateContent.replace("PROT_NAME", "1rgs")
        templateContent = templateContent.replace("CLEAN_RANGE", "5")
        # save run specific template
        FULL_TEMPLATE = templateContent

    return FULL_TEMPLATE


def generate_docking_input(targetIndividual, outPath, aminoAcidPaths):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change.
    :param aminoAcidPaths:
    :param outPath:
    :param targetIndividual:
    :return:
    """

    # load template for mutagenesis
    templateContent = get_template(aminoAcidPaths)

    # insert unique mutation information
    for aa in targetIndividual[0]:
        templateContent = templateContent.replace("MUT_AA", get_3_letter_code(aa), 1)
    templateContent = templateContent.replace("SAVE_PATH", outPath + ".pdb")

    # save pml file
    mutagenesisFile = open(outPath + "_mutagenesis.pml", 'w')
    mutagenesisFile.write(templateContent)
    mutagenesisFile.close()

    # run pymol and generate pdb as docking input
    cmd = "C:/ProgramData/PyMOL/PyMOLWin.exe -cq " + outPath + "_mutagenesis.pml"
    #os.system(cmd)

    return
