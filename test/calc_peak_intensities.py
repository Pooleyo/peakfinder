def run(pos, source_name, timestep):

    import units as un
    import logging as log
    
    log.debug("Brick %s started.\n", __name__)
    
    peak_centre = []
    integrated_intensity = []

    for i in pos:
        peak_str = un.make_peak_str(i)
        soh_output_file_location = un.determine_soh_output_file_location(peak_str, source_name, timestep)
        soh_output = un.read_from_soh_output(soh_output_file_location)
        point_of_max_height = un.find_point_of_max_height(soh_output)
        dvol = un.calc_dvol(soh_output)
        current_integrated_intensity = un.calc_integrated_intensity(soh_output, dvol)
        
        peak_centre.append(point_of_max_height)
        integrated_intensity.append(current_integrated_intensity)

    log.debug("Brick %s finished.\n", __name__)

    return peak_centre, integrated_intensity
