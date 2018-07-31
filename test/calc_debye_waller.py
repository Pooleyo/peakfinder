def run(peak_centre, intensity, mass, temperature):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)
    
    ln_intensity = []

    for i in intensity:
        current_ln_intensity = un.get_ln_intensity(i)
        ln_intensity.append(current_ln_intensity)

    
    gsqr_lattice_units = un.get_gsqr_values(peak_centre)
    gsqr_per_angstrom = un.convert_to_per_angstrom(gsqr_lattice_units)
    slope, constant = un.calc_line_slope_and_constant(gsqr_per_angstrom, ln_intensity)
    debye_waller_constant = un.calc_debye_waller_constant(mass)
    debye_temperature = un.calc_debye_temperature(temperature, slope, debye_waller_constant)
    #un.calc_temperature_XRD()
    
    log.debug("Brick %s finished.\n", __name__)

    return
