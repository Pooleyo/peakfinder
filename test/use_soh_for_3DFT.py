def run(pos_est, source_name, timestep, mass, a_lattice, N_atoms, k_steps, run_soh, num_cores):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)
    
    offset = un.calc_k_offset_with_N_atoms(N_atoms)
    
    for i in pos_est:
        k_start = un.find_k_start(i, offset)
        k_stop = un.find_k_stop(i, offset)
        peak_str = un.make_peak_str(i)
        input_file_location = un.determine_soh_input_file_location(peak_str)
        un.write_soh_input_3DFT(source_name, input_file_location, peak_str, mass, a_lattice, k_steps, k_start, k_stop)
    
    if run_soh == True:
    
        for i in pos_est:
            peak_str = un.make_peak_str(i)
            input_file_location = un.determine_soh_input_file_location(peak_str)
            un.run_soh(input_file_location, num_cores)
            un.move_soh_output_to_peak_folder(peak_str, source_name, timestep)
            
    log.debug("Brick %s finished.\n", __name__)
        
    return
