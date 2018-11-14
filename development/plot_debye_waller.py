def run(gsqr, ln_intensity, raw_pos_est, temperature, mass):

    import units as un

    print "Plotting intensity vs. G^2 for all peaks, and in kx, ky, kz..."

    filename = "ln_intensity_vs_gsqr_per_angstrom.png"

    x_label = "$G^2$ (A$^-2$)"

    y_label = "ln(I) (arb.)"

    plot_title = "Intensity vs. G$^2$"

    un.plot_matplotlib(gsqr, ln_intensity, filename, x_label, y_label, plot_title)

    plotting_directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    direction_plot_filenames = ["parallel_x.png", "parallel_y.png", "parallel_z.png"]

    direction_result_filenames = ["parallel_x.pkfd", "parallel_y.pkfd", "parallel_z.pkfd"]

    for i, direction in enumerate(plotting_directions):

        current_peak_list = []
        current_intensity_list = []

        for j, pos in enumerate(raw_pos_est):

            if un.find_if_vectors_parallel(direction, pos) is True:

                current_peak_list.append(gsqr[j])
                current_intensity_list.append(ln_intensity[j])

            else:

                pass

        un.plot_matplotlib(current_peak_list, current_intensity_list, direction_plot_filenames[i], x_label, y_label, plot_title)

        slope, constant = un.calc_line_slope_and_constant(current_peak_list, current_intensity_list)

        debye_waller_constant = un.calc_debye_waller_constant(mass)

        debye_temperature_xrd = un.calc_debye_temperature_xrd(temperature, slope, debye_waller_constant)

        un.write_temperatures_to_file(debye_temperature_xrd, temperature, direction_result_filenames[i])

    return
