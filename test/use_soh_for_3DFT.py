def run(pos_est, source_location, mass, a_lattice, N_atoms, k_steps, run_soh):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)
    
    offset = un.calc_k_offset_with_N_atoms(N_atoms)
    
    for i in pos_est:
        peak_str = un.make_peak_str(i)
        input_file_location = un.determine_accurate_soh_input_file_location(peak_str)
        k_start = un.find_k_start(i, offset)
        k_stop = un.find_k_stop(i, offset)
        k_start = un.convert_to_per_angstrom(k_start, a_lattice)
        k_stop = un.convert_to_per_angstrom(k_stop, a_lattice)
        un.write_soh_input_3DFT(source_location, input_file_location, peak_str, mass, a_lattice, k_steps, k_start, k_stop)
    
    if run_soh == True:
    
        for i in pos_est:
            input_file_location = un.determine_accurate_soh_input_file_location(peak_str)
            un.run_soh(input_file_location)
    
    log.debug("Brick %s finished.\n", __name__)
        
    return
