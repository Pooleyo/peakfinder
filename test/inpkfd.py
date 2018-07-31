path = "path_1_stat" # Choose: "path_1_stat"
run_soh = True
num_cores = 20

# Input for "select_peak_positions"
gsqr_max = 3
negative_k = True
remove_000 = False

# Input for "use_soh_for_3DFT"
source_name = "uncompressed_300K.atom"
timestep = "10000" # Only used in moving soh output files.
mass = 63.546 # In amu
a_lattice = 3.615 # In Angstroms
k_steps = 3
N_atoms = [20,20,20]

# Input for use by "calc_debye_waller"
temperature = 300.0
