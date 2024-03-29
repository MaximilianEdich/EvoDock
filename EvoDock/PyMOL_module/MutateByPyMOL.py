"""
Performs a mutagenesis via the external Software PyMOL.

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import os

HEAD_TEMPLATE = ""
LOOP_TEMPLATE = ""
TAIL_TEMPLATE = ""


def is_mutation_module():
    """
    This function is essential to verify, that it is a mutation module and is required, so it cannot be used as
    a module of a different type.
    """
    return True


def get_3_letter_code(single_letter):
    """
    Converts a single letter amino acid code to a three letter amino acid code.
    :param single_letter: Single letter of the one letter amino acid code.
    :return: The corresponding three letter amino acid code, if entered correct letter, otherwise an empty string.
    """

    if single_letter == 'A':
        return "ALA"
    if single_letter == 'R':
        return "ARG"
    if single_letter == 'N':
        return "ASN"
    if single_letter == 'D':
        return "ASP"
    if single_letter == 'C':
        return "CYS"
    if single_letter == 'Q':
        return "GLN"
    if single_letter == 'E':
        return "GLU"
    if single_letter == 'G':
        return "GLY"
    if single_letter == 'H':
        return "HIS"
    if single_letter == 'I':
        return "ILE"
    if single_letter == 'L':
        return "LEU"
    if single_letter == 'K':
        return "LYS"
    if single_letter == 'M':
        return "MET"
    if single_letter == 'F':
        return "PHE"
    if single_letter == 'P':
        return "PRO"
    if single_letter == 'S':
        return "SER"
    if single_letter == 'T':
        return "THR"
    if single_letter == 'W':
        return "TRP"
    if single_letter == 'Y':
        return "TYR"
    if single_letter == 'V':
        return "VAL"

    return ""


def get_template(amino_acid_paths, mutations, protein_code):
    """
    Create the template for a specific protein.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param mutations: List of mutations relative to the original protein. An empty string represents no mutation while
    any substitution is represented by the given single letter code of the amino acid.
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :return: The specific template for the creation of the specific pml file to perform the mutagenesis.
    """

    # create template on the first run
    global HEAD_TEMPLATE
    global LOOP_TEMPLATE
    global TAIL_TEMPLATE

    if HEAD_TEMPLATE == "":
        # load head part
        template = open("mutagenesis_template_head.pml", 'r')
        head_content = ""
        for line in template.readlines():
            head_content += line
        template.close()
        head_content = head_content.replace("PROT_NAME", protein_code)
        HEAD_TEMPLATE = head_content

        # load loop part
        template = open("mutagenesis_template_loop.pml", 'r')
        loop_content = ""
        for line in template.readlines():
            loop_content += line
        template.close()
        loop_content = loop_content.replace("CLEAN_RANGE", "5")
        LOOP_TEMPLATE = loop_content

        # templateContent = templateContent.replace("RESIDUE_N", str(path))
        # templateContent = templateContent.replace("RESIDUE_SL", str(path)[1:].replace("/", "_"))

        # load tail part
        template = open("mutagenesis_template_tail.pml", 'r')
        tail_content = ""
        for line in template.readlines():
            tail_content += line
        template.close()
        tail_content = tail_content.replace("PROT_NAME", protein_code)
        TAIL_TEMPLATE = tail_content

    i = 0
    template = HEAD_TEMPLATE
    for aa in mutations:
        if not aa == "":
            template += LOOP_TEMPLATE
            template = template.replace("RESIDUE_N", str(amino_acid_paths[i]))
            template = template.replace("RESIDUE_SL", str(amino_acid_paths[i])[1:].replace("/", "_"))
        i += 1
    template += TAIL_TEMPLATE

    return template


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
    :return: Nothing.
    """

    # load template for mutagenesis
    template_content = get_template(amino_acid_paths, mutations, protein_path)

    # insert unique mutation information
    for aa in mutations:
        if not aa == "":
            # replace next "MUT_AA" with 3-letter code of mutated amino acid
            template_content = template_content.replace("MUT_AA", get_3_letter_code(aa), 1)
    template_content = template_content.replace("SAVE_PATH", out_path + ".pdb")

    # save pml file
    mutagenesis_file = open(out_path + "_mutagenesis.pml", 'w')
    mutagenesis_file.write(template_content)
    mutagenesis_file.close()

    # run pymol and generate pdb as docking input
    cmd = "C:/ProgramData/PyMOL/PyMOLWin.exe -cq " + out_path + "_mutagenesis.pml"
    os.system(cmd)

    return
