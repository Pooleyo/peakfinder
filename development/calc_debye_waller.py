def run(peak_centre, intensity, a_lattice, mass, temperature, uncompressed_debye_temperature, single_term_model_gamma_0_values, single_term_model_q_values, triple_term_model_gamma_0_values, triple_term_constant_values):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)

    print "Calculating Debye-Waller effect..."
    
    ln_intensity = []

    for i in intensity:
        current_ln_intensity = un.get_ln_intensity(i)
        ln_intensity.append(current_ln_intensity)

    peak_centre_per_angstrom = []

    for i in peak_centre:
        g_per_angstrom = un.convert_to_per_angstrom(i, a_lattice)
        peak_centre_per_angstrom.append(g_per_angstrom)

    gsqr_per_angstrom = un.get_gsqr_values(peak_centre_per_angstrom)

    slope, constant = un.calc_line_slope_and_constant(gsqr_per_angstrom, ln_intensity)

    debye_waller_constant = un.calc_debye_waller_constant(mass)

    debye_temperature_xrd = un.calc_debye_temperature_xrd(temperature, slope, debye_waller_constant)

    initial_volume = un.calc_volume_lattice_units(a_lattice, [1.0, 1.0, 1.0])

    compression_factors = [1.0, 1.0, 1.0]

    final_volume = un.calc_volume_lattice_units(a_lattice, compression_factors)

    model_debye_temperatures = []

    for i in range(len(single_term_model_gamma_0_values)):

        debye_temperature = un.calc_debye_temperature_from_single_term_gruneisen_model( \
            uncompressed_debye_temperature, initial_volume, final_volume, single_term_model_gamma_0_values[i],
            single_term_model_q_values[i])

        model_debye_temperatures.append(debye_temperature)

    for i in range(len(triple_term_model_gamma_0_values)):

        debye_temperature = un.calc_debye_temperature_from_triple_term_gruneisen_model( \
            uncompressed_debye_temperature, initial_volume, final_volume, triple_term_model_gamma_0_values[i],
            triple_term_constant_values[i])

        model_debye_temperatures.append(debye_temperature)

    temperature_xrd = []

    for theta in model_debye_temperatures:

        current_temperature_xrd = un.calc_temperature_xrd(theta, slope, debye_waller_constant)

        temperature_xrd.append(current_temperature_xrd)

    log.debug("Brick %s finished.\n", __name__)

    return debye_temperature_xrd, temperature_xrd
