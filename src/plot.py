def run(peak_str, peak_centre, source_name, timestep):

    import units as un

    directions = ["kx", "ky", "kz"]

    constant_axes = [[1,2], [0,2], [0,1]]

    variable_axes = [0, 1, 2]

    for i, current_peak_str in enumerate(peak_str):

        data_filename = un.determine_soh_output_file_location(current_peak_str, source_name, timestep)

        soh_output = un.read_from_soh_output(data_filename)

        centre_point = peak_centre[i]

        for j, direction in enumerate(directions):

            k, intensity = un.find_line_data_from_3DFT(constant_axes[j], variable_axes[j], centre_point, soh_output)

            print k
            print intensity

            plot_filename = direction + ".png"

            un.plot(k, intensity, plot_filename)

        for direction in directions:

            un.move_plot_output_to_peak_folder(direction, current_peak_str)

    return
