def run(run_soh, current_pos_est, raw_pos_est, source_name, timestep, undershoot, overshoot, source, mass, a_lattice, k_steps, num_cores):

    import units as un
    import copy
    import os

    print "Fitting to peak edges..."

    direction_str = ["kx", "ky", "kz"]

    k_start_accurate = copy.deepcopy(current_pos_est)
    k_stop_accurate = copy.deepcopy(current_pos_est)

    for i, pos in enumerate(raw_pos_est):

        peak_str = un.make_peak_str(pos)

        for j, direction in enumerate(direction_str):

            k_start, k_stop = un.calc_peak_edge_k_start_stop(current_pos_est[i], undershoot[j], overshoot[j])

            soh_location = un.determine_soh_edge_finding_input_file_location(direction, peak_str)

            un.write_soh_input_1DFT(source_name, soh_location, peak_str + "_find_edges_" + direction, mass, a_lattice, k_start, k_stop, k_steps)

    if run_soh is True:

        for i, pos in enumerate(raw_pos_est):

            peak_str = un.make_peak_str(pos)

            for j, direction in enumerate(direction_str):

                soh_location = un.determine_soh_edge_finding_input_file_location(direction, peak_str)

                un.run_soh(soh_location, num_cores)

                un.move_soh_rough_output_to_peak_folder(peak_str, peak_str + "_find_edges_" + direction, source_name, timestep)

    for i, pos in enumerate(raw_pos_est):

        peak_str = un.make_peak_str(pos)

        for j, direction in enumerate(direction_str):

            soh_location = un.determine_soh_edge_finding_output_file_location(peak_str, direction, source, timestep)

            soh_output = un.read_from_soh_output(soh_location)

            un.plot_pygnuplot(soh_output[j], soh_output[3], os.getcwd() + "/data/" + peak_str + "/" + direction + "_lineout.png", os.getcwd() + "/data/" + peak_str + "/" + direction + "_lineout.dat")

            k_start_accurate[i][j], k_stop_accurate[i][j] = un.find_k_start_stop_for_peak_from_first_minima(soh_output[j], soh_output[3])

    return k_start_accurate, k_stop_accurate