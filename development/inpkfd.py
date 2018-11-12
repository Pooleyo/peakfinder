path = "path_2_dynamic_peakfinding"  # Choose from: "path_1_static_peakfinding", "path_2_dynamic_peakfinding"
run_soh = True
make_peak_plots = True
num_cores = 2

# Input for "select_peak_positions"
gsqr_max = 20
negative_k = False
remove_000 = False

# Input for "find_compression_ratio"
uncompressed_peak_positions = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
undershoot = 0.5
overshoot = 1.5
lineout_k_steps = 1e4 + 1

# Input for "use_soh_for_3DFT"
source_name = "uncompressed_cu_300K_5x5x60_10000.atom"
timestep = "10000"  # Only used in moving soh output files.
mass = 63.546  # In amu
a_lattice = 3.628  # In Angstroms
k_steps = 11
N_atoms = [5, 5, 60]

# Input for use by "calc_debye_waller"
uncompressed_debye_temperature = 320.0
temperature = 300.0

# These Cu model values are from Will Murphy et. al. PHYSICAL REVIEW B 78, 014109 (2008)
single_term_model_gamma_0_values = [1.98, 1.93, 2.008]
single_term_model_exponent_values = [1.0, 1.085, 1.33]

triple_term_model_gamma_0_values = [2.04]
triple_term_model_constants = [[-3.296, 10.493, -19.264]]
