
def perform_docking(mutations, out_path, amino_acid_paths, protein_code):
    """
    Performs a mutagenesis on the original pdb file. Substitutes specific amino acids and optimizes rotamer and
    adapt the backbone to the change. Results are saved in a new pdb file. During this process a pml file is
    generated, containing the PyMOL script that performs the mutagenesis.
    :param protein_code: The protein accession code, by wich the protein structure can be fetched with.
    :param amino_acid_paths: Paths within the pdb file to the single amino acids of interest.
    :param out_path: The path leading to the output files specific to the mutation.
    :param mutations: List of mutations relative to the original protein. An empty string represents no mutation while
    any substitution is represented by the given single letter code of the amino acid.
    :return: None. The generated files are of interest.
    """

    return

