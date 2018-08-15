import test_units as tu
import test_bricks as tb
import numpy as np

########################
# Unit tests.

tu.run("build_all_k_values", [4, True], ([[-2, 0, 0], [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [0, -2, 0], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, 0, -2], [0, 0, -1], [0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, -1], [0, 1, 0], [0, 1, 1], [0, 2, 0], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 1, -1], [1, 1, 0], [1, 1, 1], [2, 0, 0]]))

tu.run("build_all_k_values", [4, False], ([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 0], [0, 1, 1], [0, 2, 0], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1], [2, 0, 0]]))

tu.run("remove_fcc_forbidden_reflections", [[[-2, 0, 0], [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [0, -2, 0], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, 0, -2], [0, 0, -1], [0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, -1], [0, 1, 0], [0, 1, 1], [0, 2, 0], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 1, -1], [1, 1, 0], [1, 1, 1], [2, 0, 0]]], ([[-2, 0, 0], [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [0, -2, 0], [0, 0, -2], [0, 0, 0], [0, 0, 2], [0, 2, 0], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1], [2, 0, 0]]))

tu.run("remove_000", [[[-2, 0, 0], [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [0, -2, 0], [0, 0, -2], [0, 0, 0], [0, 0, 2], [0, 2, 0], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1], [2, 0, 0]]], [[-2, 0, 0], [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [0, -2, 0], [0, 0, -2], [0, 0, 2], [0, 2, 0], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1], [2, 0, 0]])

tu.run("get_gsqr_values", [[[-2, 0, 0], [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [0, -2, 0], [0, 0, -2], [0, 0, 2], [0, 2, 0], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1], [2, 0, 0]]], [4, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4])

tu.run("calc_k_offset_with_N_atoms", [[10,10,10]], [0.1, 0.1, 0.1])

tu.run("find_k_start", [[1,1,1], [0.1,0.1,0.1]], [0.9, 0.9, 0.9])

tu.run("find_k_stop", [[1,1,1], [0.1,0.1,0.1]], [1.1, 1.1, 1.1])

tu.run("convert_to_per_angstrom", [[1, 1, 1], 3.0], [1.0 * (2*np.pi) / 3.0, 1.0 * (2 * np.pi) / (3.0), (2 * np.pi) / (3.0)])

tu.run("calc_dvol", [[[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1], [0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1,0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1,0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1], [0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1], [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]]], 0.001000000000000001 )

tu.run("calc_integrated_intensity", [[[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1], [0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1,0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1,0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1], [0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1,0.9,1.0,1.1], [0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0]], 0.001], 0.002)

tu.run("get_ln_intensity", [7], 1.9459101490553132)

tu.run("calc_line_slope_and_constant", [[0,2],[1,4]], (1.4999999999999996,0.99999999999999989) )

tu.run("calc_debye_waller_constant", [1.0], 145.5262505987811)

tu.run("calc_debye_temperature_from_xrd", [300,-2,150], 150)

tu.run("calc_debye_temperature_from_single_term_gruneisen_model", [200, 1.0, 0.9, 2.0, 2.0], 241.8499195314503)

tu.run("calc_debye_temperature_from_single_term_gruneisen_model", [320, 1.0, 0.9, 2.008, 1.33], 389.8377145278562)

tu.run("calc_volume_lattice_units", [2.0, [0.8, 0.99, 0.99]], 6.27264)

tu.run("calc_debye_temperature_from_triple_term_gruneisen_model", [320, 1.0, 0.9, 2.0, [-4.0, 10.0, -20.0]], 387.44039572033546)

tu.run("calc_temperature_xrd", [320.0, -1, 2.29], 44716.157205240175)

########################
# Brick tests.

tb.run( "select_peak_positions", [6, False, False], ([[0,0,0], [0,0,2], [0,2,0], [1,1,1], [2,0,0]], [0,4,4,3,4]))
