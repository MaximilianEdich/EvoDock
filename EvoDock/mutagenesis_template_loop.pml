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

