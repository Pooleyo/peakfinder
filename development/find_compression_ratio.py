def run(run_soh, lineout_directions, undershoot, overshoot, source, mass, a_lattice, lineout_k_steps,
        timestep, soh_command):

    import units as un
    import os

    print "Finding compression ratios..."

    compression_ratio = [1.0, 1.0, 1.0]

    direction_str = ["kx", "ky", "kz"]

    if os.path.exists("./data/lineouts/"):

        pass

    else:

        un.make_lineout_directory()

    for i, direction in enumerate(lineout_directions):

        k_start, k_stop = un.calc_lineout_k_start_stop(direction, undershoot, overshoot)

        soh_location = un.determine_soh_compression_finding_input_file_location(direction_str[i])

        un.write_soh_input_1DFT(source, soh_location, "lineout_" + direction_str[i], mass, a_lattice, k_start, k_stop, lineout_k_steps)

    if run_soh is True:

        for i, direction in enumerate(lineout_directions):

            soh_location = un.determine_soh_compression_finding_input_file_location(direction_str[i])

            un.run_soh(soh_location, soh_command)

            un.move_soh_output_to_lineout_folder(direction_str[i], source, timestep)

    for i, direction in enumerate(lineout_directions):

        soh_location = un.determine_soh_1DFT_output_file_location(direction_str[i], source, timestep)

        soh_output = un.read_from_soh_output(soh_location)

        un.plot_pygnuplot(soh_output[i], soh_output[3], os.getcwd() + "/data/lineouts/" + direction_str[i] + "_lineout.png", os.getcwd() + "/data/lineouts/" + direction_str[i] + "_lineout.dat")

        k_of_max_height = un.find_point_of_max_height(soh_output)

        compression_ratio[i] = un.calc_compression_ratio(k_of_max_height[i], lineout_directions[i][i])

    return compression_ratio
