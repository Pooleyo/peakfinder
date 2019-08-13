def run(run_soh, raw_pos_est, pos_est, gsqr_est, compression_ratio, source_name, N_atoms, mass, a_lattice,
        k_steps_find_centre_1D, k_steps_find_centre_3D, timestep, soh_command, plot):

    import units as un
    import copy

    print "Fitting to peak centres..."

    centre_guess_3DFT = copy.deepcopy(pos_est)  # This list will contain the first guesses of peak centres gained from
    # the three 1DFTs performed below. It is used as a first guess of the centre for the final 3DFT.

    fitted_pos_est = copy.deepcopy(pos_est)  # This list contains the final peak centres, determined from the final
    # 3DFT. These are used as the input for peak centres when calculating the intensity of a full peak.

    compressed_pos_est, compressed_gsqr_est = un.apply_compression_ratio_to_pos_est(pos_est, gsqr_est, compression_ratio)

    offset = un.calc_k_offset_with_N_atoms(N_atoms)

    for i, pos in enumerate(compressed_pos_est):

        k_start = un.find_simple_k_start(pos, offset)

        k_stop = un.find_simple_k_stop(pos, offset)

        kx_start, kx_stop, ky_start, ky_stop, kz_start, kz_stop = un.calc_lineout_k_start_stop_along_xyz(k_start,
                                                                                                         k_stop,
                                                                                                         pos)
        peak_str = un.make_peak_str(raw_pos_est[i])

        input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_kx_1DFT.in")

        un.write_soh_input_1DFT(source_name, input_file_location, "find_centre_kx_" + peak_str, mass, a_lattice, kx_start, kx_stop, k_steps_find_centre_1D)

        input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_ky_1DFT.in")

        un.write_soh_input_1DFT(source_name, input_file_location, "find_centre_ky_" + peak_str, mass, a_lattice, ky_start, ky_stop, k_steps_find_centre_1D)

        input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_kz_1DFT.in")

        un.write_soh_input_1DFT(source_name, input_file_location, "find_centre_kz_" + peak_str, mass, a_lattice, kz_start, kz_stop, k_steps_find_centre_1D)

    if run_soh is True:

        for i, pos in enumerate(compressed_pos_est):

            peak_str = un.make_peak_str(raw_pos_est[i])

            input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_kx_1DFT.in")

            un.run_soh(input_file_location, soh_command)

            un.move_soh_rough_output_to_peak_folder(peak_str, "find_centre_kx_" + peak_str, source_name, timestep)

            input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_ky_1DFT.in")

            un.run_soh(input_file_location, soh_command)

            un.move_soh_rough_output_to_peak_folder(peak_str, "find_centre_ky_" + peak_str, source_name, timestep)

            input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_kz_1DFT.in")

            un.run_soh(input_file_location, soh_command)

            un.move_soh_rough_output_to_peak_folder(peak_str, "find_centre_kz_" + peak_str, source_name, timestep)

    for i, pos in enumerate(compressed_pos_est):

        peak_str = un.make_peak_str(raw_pos_est[i])

        appended_string = "find_centre_kx_" + peak_str

        soh_output_file_location = un.determine_rough_soh_output_file_location(peak_str, source_name, timestep,
                                                                               appended_string)

        soh_output = un.read_from_soh_output(soh_output_file_location)

        k_max = un.find_point_of_max_height(soh_output)

        kx_centre = k_max[0]

        if plot is True:

            un.plot_pygnuplot(soh_output[0], soh_output[3], "./data/" + peak_str + "/find_centre_kx.png",
                              "./data/" + peak_str + "/find_centre_kx.dat")

        appended_string = "find_centre_ky_" + peak_str

        soh_output_file_location = un.determine_rough_soh_output_file_location(peak_str, source_name, timestep,
                                                                               appended_string)

        soh_output = un.read_from_soh_output(soh_output_file_location)

        k_max = un.find_point_of_max_height(soh_output)

        ky_centre = k_max[1]

        if plot is True:

            un.plot_pygnuplot(soh_output[1], soh_output[3], "./data/" + peak_str + "/find_centre_ky.png",
                              "./data/" + peak_str + "/find_centre_ky.dat")


        appended_string = "find_centre_kz_" + peak_str

        soh_output_file_location = un.determine_rough_soh_output_file_location(peak_str, source_name, timestep,
                                                                               appended_string)

        soh_output = un.read_from_soh_output(soh_output_file_location)

        k_max = un.find_point_of_max_height(soh_output)

        kz_centre = k_max[2]


        if plot is True:

            un.plot_pygnuplot(soh_output[2], soh_output[3], "./data/" + peak_str + "/find_centre_kz.png",
                              "./data/" + peak_str + "/find_centre_kz.dat")

        centre_guess_3DFT[i] = [kx_centre, ky_centre, kz_centre]

    #  Now that the initial guesses have been determined using the 1DFTs above, a 3DFT is performed to get as close as
    #  possible to the centre of the peak.

    dk = un.calc_dk_from_offset(offset, k_steps_find_centre_1D, k_steps_find_centre_1D, k_steps_find_centre_1D)

    for i, pos in enumerate(centre_guess_3DFT):

        k_start = un.find_simple_k_start(pos, dk)

        k_stop = un.find_simple_k_stop(pos, dk)

        peak_str = un.make_peak_str(raw_pos_est[i])

        input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_3DFT.in")

        un.write_soh_input_3DFT(source_name, input_file_location, "find_centre_3DFT_" + peak_str, mass, a_lattice,
                                k_steps_find_centre_3D, k_start, k_stop)

    if run_soh is True:

        for i, pos in enumerate(centre_guess_3DFT):

            peak_str = un.make_peak_str(raw_pos_est[i])

            input_file_location = un.determine_rough_soh_input_file_location(peak_str, "find_centre_3DFT.in")

            un.run_soh(input_file_location, soh_command)

            un.move_soh_rough_output_to_peak_folder(peak_str, "find_centre_3DFT_" + peak_str, source_name, timestep)

    for i, pos in enumerate(centre_guess_3DFT):

        peak_str = un.make_peak_str(raw_pos_est[i])

        appended_string = "find_centre_3DFT_" + peak_str

        soh_output_file_location = un.determine_rough_soh_output_file_location(peak_str, source_name, timestep,
                                                                               appended_string)

        soh_output = un.read_from_soh_output(soh_output_file_location)

        k_max = un.find_point_of_max_height(soh_output)

        fitted_pos_est[i] = k_max

    return fitted_pos_est
