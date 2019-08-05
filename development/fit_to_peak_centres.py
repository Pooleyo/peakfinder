def run(run_soh, raw_pos_est, pos_est, gsqr_est, compression_ratio, source_name, N_atoms, mass, a_lattice, k_steps_find_centre, timestep, soh_command):

    import units as un
    import copy

    print "Fitting to peak centres..."

    fitted_pos_est = copy.deepcopy(pos_est)

    compressed_pos_est, compressed_gsqr_est = un.apply_compression_ratio_to_pos_est(pos_est, gsqr_est, compression_ratio)

    offset = un.calc_k_offset_with_N_atoms(N_atoms)

    for i, pos in enumerate(compressed_pos_est):

        k_start = un.find_simple_k_start(pos, offset)

        k_stop = un.find_k_stop(pos, offset)

        peak_str = un.make_peak_str(raw_pos_est[i])

        input_file_location = un.determine_rough_soh_input_file_location(peak_str)

        un.write_soh_input_3DFT(source_name, input_file_location, "find_centre_" + peak_str, mass, a_lattice, k_steps_find_centre, k_start, k_stop)

    if run_soh is True:

        for i, pos in enumerate(compressed_pos_est):

            peak_str = un.make_peak_str(raw_pos_est[i])

            input_file_location = un.determine_rough_soh_input_file_location(peak_str)

            un.run_soh(input_file_location, soh_command)

            un.move_soh_rough_output_to_peak_folder(peak_str, "find_centre_" + peak_str, source_name, timestep)

    for i, pos in enumerate(compressed_pos_est):

        peak_str = un.make_peak_str(raw_pos_est[i])

        soh_output_file_location = un.determine_rough_soh_output_file_location(peak_str, source_name, timestep)

        soh_output = un.read_from_soh_output(soh_output_file_location)

        fitted_pos_est[i] = un.find_point_of_max_height(soh_output)

    return fitted_pos_est
