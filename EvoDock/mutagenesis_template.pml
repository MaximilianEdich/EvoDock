# PROT_NAME = name of the protein, it can be fetched by
# RESIDUE_N = residue chain-path in normal form with slashes, like "/1rgs/A/A/333"
# RESIDUE_SL = residue chain-path without slashes, like "1rgs_A_A_333"
# MUT_AA = Amino acid the chosen residue mutates into
# CLEAN_RANGE = range in angstrom around the mutated residue, that will be cleaned
# SAVE_PATH = path the resulting pdb will be saved to

# fetch protein
fetch PROT_NAME, async=0

# start mutagenesis wizard
cmd.wizard("mutagenesis")
cmd.do("refresh_wizard")

# mutate
cmd.get_wizard().set_mode("MUT_AA")
cmd.get_wizard().do_select("RESIDUE_N")

# Select the first rotamer, which is most probable
cmd.frame(1)

# Apply the mutation
cmd.get_wizard().apply()
# Close wizard
cmd.set_wizard("done")


# adapt reagion
select RESIDUE_N, RESIDUE_N expand CLEAN_RANGE
clean RESIDUE_SL

save SAVE_PATH, PROT_NAME, -1, pdb