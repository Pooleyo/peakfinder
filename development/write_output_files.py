def run(debye_temperature, temperature, model_debye_temperatures, pos_est, peak_centre, gsqr, integrated_intensity, ln_intensity, slope):

    import units as un

    filename_temperatures = "data/results.pkfd"

    un.write_temperatures_to_file(slope, debye_temperature, temperature, model_debye_temperatures, filename_temperatures)

    filename_peaks = "data/integrated_intensity.dat"

    un.write_peak_intensities_to_file(pos_est, peak_centre, gsqr, integrated_intensity, ln_intensity, filename_peaks)

    return