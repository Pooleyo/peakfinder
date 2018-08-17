path = "path_1_stat"  # Choose: "path_1_stat"
run_soh = True
make_plots = True
num_cores = 4

# Input for "select_peak_positions"
gsqr_max = 100
negative_k = True
remove_000 = False

# Input for "use_soh_for_3DFT"
source_name = "uncompressed_300K_10000.atom"
timestep = "10000"  # Only used in moving soh output files.
mass = 63.546  # In amu
a_lattice = 3.615  # In Angstroms
k_steps = 11
N_atoms = [5, 5, 5]

# Input for use by "calc_debye_waller"
uncompressed_debye_temperature = 320.0
temperature = 300.0

# These Cu model values are from Will Murphy et. al. PHYSICAL REVIEW B 78, 014109 (2008)
single_term_model_gamma_0_values = [1.98, 1.93, 2.008]
single_term_model_exponent_values = [1.0, 1.085, 1.33]

triple_term_model_gamma_0_values = [2.04]
triple_term_model_constants = [[-3.296, 10.493, -19.264]]
