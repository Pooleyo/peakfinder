path = "path_2_dynamic_peakfinding"  # Choose from: "path_1_static_peakfinding", "path_2_dynamic_peakfinding", "path_3_rotated_dynamic_peakfinding"
run_soh = True
make_final_peak_plots = False
#make_plots_peak_centre_fit = True # THIS DOESN


# 'T WORK YET
num_cores = 8
soh_command = "mpiexec -np " + str(num_cores) + " sonOfHoward "

N_atoms = [5, 5, 5]

# Input for "select_peak_positions"
crystal_type = 'FCC'  # Choose between 'FCC' and 'BCC'
gsqr_max = 11
negative_k = False
remove_000 = False

# Input for "find_compression_ratio"
uncompressed_peak_positions = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]] #[[1.4142135623730951, 0.81649658092772592, 1.1547005383792515], [-1.4142135623730949, 0.81649658092772603, 1.1547005383792517], [0.0, -1.6329931618554521, 1.1547005383792515]]
compression_ratio_undershoot = 0.9
compression_ratio_overshoot = 1.1
lineout_k_steps = 1e3 + 1

# Input for "fit_to_peak_centres"
k_steps_find_centre_1DFT = 21
k_steps_find_centre_3DFT = 21

# Input for "fit_to_peak_edges"
peak_edge_undershoot = [1.0/N_atoms[0], 1.0/N_atoms[1], 1.0/N_atoms[2]]
peak_edge_overshoot = [1.0/N_atoms[0], 1.0/N_atoms[1], 1.0/N_atoms[2]]
peak_edge_k_steps = 1e3 + 1

# Input for "use_soh_for_3DFT"
source_name = "uncompressed_450K_1000.atom"
timestep = "1000"  # Only used in moving soh output files.
mass = 63.546  # In amu
a_lattice = 3.6288 # In Angstroms
k_steps = 11

# Input for use by "calc_debye_waller"
uncompressed_debye_temperature = 311.0
temperature = 450.0

# Input for calc_md_temperature. Simulations must have "units metal" for this calculation to work properly.
calc_md_temperature_from_dump_file = True
calculated_temperature_dimensionality = 3  # Enter "3" for 3D temperature, and "2" for temperature calculated from vx
# and vy only.
velocity_columns = [5, 6, 7]  # The columns in the lammps dump file that contain vx, vy, vz, respectively. The first
# column is column number 0.
number_velocity_bins = 10  # The number of bins used to make the histogram of atom velocity.


# These Cu model values are from Will Murphy et. al. PHYSICAL REVIEW B 78, 014109 (2008)
single_term_model_gamma_0_values = [1.98, 1.93, 2.008]
single_term_model_exponent_values = [1.0, 1.085, 1.33]

triple_term_model_gamma_0_values = [2.04]
triple_term_model_constants = [[-3.296, 10.493, -19.264]]
