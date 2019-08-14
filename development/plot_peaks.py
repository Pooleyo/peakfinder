def run(peak_str, peak_centre, source_name, timestep):

    import units as un
    import numpy as np

    print "Plotting each peak in kx, ky, kz, and |G^2|..."

    directions = ["kx", "ky", "kz"]

    constant_axes = [[1,2], [0,2], [0,1]]

    variable_axes = [0, 1, 2]

    for i, current_peak_str in enumerate(peak_str):

        data_filename = un.determine_accurate_soh_output_file_location(current_peak_str, source_name, timestep)

        soh_output = un.read_from_soh_output(data_filename)

        centre_point = peak_centre[i]

        # This for loop plots the peak through kx, ky, and kz.

        for j, direction in enumerate(directions):

            k, intensity = un.find_line_data_from_3DFT(constant_axes[j], variable_axes[j], centre_point, soh_output)

            plot_data_filename = "./data/" + current_peak_str + "/" + direction + ".dat"

            plot_filename = "./data/" + current_peak_str + "/" + direction + ".png"

            un.plot_pygnuplot(k, intensity, plot_filename, plot_data_filename)

        # This section plots every point in the 3DFT as intensity vs. G^2.

        gsqr = list((np.array(soh_output[0])**2) + (np.array(soh_output[1])**2) + (np.array(soh_output[2])**2))

        intensity = soh_output[3]

        plot_data_filename = "./data/" + current_peak_str + "/I_vs_gsqr.dat"

        plot_filename = "./data/" + current_peak_str + "/I_vs_gsqr_" + current_peak_str + ".png"

        un.plot_pygnuplot(gsqr, intensity, plot_filename, plot_data_filename)

    return

