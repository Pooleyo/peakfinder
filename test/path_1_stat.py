def run():

    import select_peak_positions
    import build_datafile_structure
    import use_soh_for_3DFT
    import calc_peak_intensities
    import calc_debye_waller
    import inpkfd as ip
    import logging as log
    
    log.info("Path %s started.\n", __name__)
    
    pos_est, gsqr_est = select_peak_positions.run(ip.gsqr_max, ip.negative_k, ip.remove_000)
    
    build_datafile_structure.run(pos_est)
    
    use_soh_for_3DFT.run(pos_est, ip.source_name, ip.timestep, ip.mass, ip.a_lattice, ip.N_atoms, ip.k_steps, ip.run_soh, ip.num_cores)
    
    peak_centre, integrated_intensity = calc_peak_intensities.run(pos_est, ip.source_name, ip.timestep)
    
    calc_debye_waller.run(peak_centre, integrated_intensity, ip.mass, ip.temperature)
    
    #plot.run()
        
    log.info("Path %s finished.\n", __name__)
    
    return
