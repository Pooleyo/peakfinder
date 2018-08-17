def run():

    import inpkfd as ip

    import select_peak_positions
    import build_datafile_structure
    import use_soh_for_3DFT
    import calc_peak_intensities
    import calc_debye_waller
    import write_results
    import plot

    import logging as log
    
    log.info("Path %s started.\n", __name__)
    
    pos_est, gsqr_est = select_peak_positions.run(ip.gsqr_max, ip.negative_k, ip.remove_000)
    
    build_datafile_structure.run(pos_est)
    
    use_soh_for_3DFT.run(pos_est, ip.source_name, ip.timestep, ip.mass, ip.a_lattice, ip.N_atoms, ip.k_steps, ip.run_soh, ip.num_cores)
    
    peak_centre, integrated_intensity = calc_peak_intensities.run(pos_est, ip.source_name, ip.timestep)
    
    debye_temperature, temperature = calc_debye_waller.run(peak_centre, integrated_intensity, ip.a_lattice, ip.mass,
                                                           ip.temperature, ip.uncompressed_debye_temperature,
                                                           ip.single_term_model_gamma_0_values,
                                                           ip.single_term_model_exponent_values,
                                                           ip.triple_term_model_gamma_0_values,
                                                           ip.triple_term_model_constants)

    write_results.run(debye_temperature, temperature)

    if ip.make_plots is True:

        plot.run()
        
    log.info("Path %s finished.\n", __name__)
    
    return
