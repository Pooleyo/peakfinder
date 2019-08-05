def run(current_pos_est, pos_est, source_name, timestep, mass, a_lattice, k_steps, run_soh, k_start, k_stop, soh_command):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)

    print "Performing accurate 3DFT of each peak..."

    for i, pos in enumerate(current_pos_est):

        peak_str = un.make_peak_str(pos_est[i])

        input_file_location = un.determine_accurate_soh_input_file_location(peak_str)

        un.write_soh_input_3DFT(source_name, input_file_location, peak_str, mass, a_lattice, k_steps, k_start[i], k_stop[i])

    if run_soh is True:
    
        for i in pos_est:

            peak_str = un.make_peak_str(i)

            input_file_location = un.determine_accurate_soh_input_file_location(peak_str)

            un.run_soh(input_file_location, soh_command)

            un.move_soh_accurate_output_to_peak_folder(peak_str, source_name, timestep)

    log.debug("Brick %s finished.\n", __name__)

    return
