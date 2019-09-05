def run():

    import inpkfd as ip

    import select_peak_positions
    import build_datafile_structure
    import use_static_peakfinding_for_3DFT
    import calc_peak_intensities
    import calc_debye_waller
    import write_output_files
    import plot_debye_waller
    import plot_peaks
    import find_compression_ratio
    import apply_compression_ratio
    import calc_md_temperature

    import logging as log

    log.info("Path %s started.\n", __name__)

    current_md_temperature = ip.temperature

    if ip.calc_md_temperature_from_dump_file is True:

        current_md_temperature = calc_md_temperature.run(ip.source_name, ip.temperature, ip.calculated_temperature_dimensionality, ip.mass, ip.velocity_columns, ip.number_velocity_bins)

    raw_pos_est, raw_gsqr_est = select_peak_positions.run(ip.gsqr_max, ip.negative_k, ip.remove_000, ip.crystal_type)

    current_pos_est = raw_pos_est

    current_gsqr_est = raw_gsqr_est

    peak_str = build_datafile_structure.run(raw_pos_est)

    compression_ratio = find_compression_ratio.run(ip.run_soh, ip.uncompressed_peak_positions,
                                                   ip.compression_ratio_undershoot,
                                                   ip.compression_ratio_overshoot, ip.source_name, ip.mass,
                                                   ip.a_lattice,
                                                   ip.lineout_k_steps, ip.timestep, ip.soh_command)

    compressed_pos_est, compressed_gsqr_est = apply_compression_ratio.run(current_pos_est, current_gsqr_est,
                                                                          compression_ratio)

    current_pos_est = compressed_pos_est

    current_gsqr_est = compressed_gsqr_est

    use_static_peakfinding_for_3DFT.run(current_pos_est, raw_pos_est, ip.source_name, ip.timestep, ip.mass, ip.a_lattice, ip.N_atoms, ip.k_steps,
                                        ip.run_soh, ip.num_cores)

    peak_centre, integrated_intensity = calc_peak_intensities.run(raw_pos_est, ip.source_name, ip.timestep)

    debye_temperature, temperature, model_debye_temperatures, gsqr_per_angstrom, ln_intensity = calc_debye_waller.run(
        peak_centre, integrated_intensity, ip.a_lattice, ip.mass,
        current_md_temperature, ip.uncompressed_debye_temperature,
        ip.single_term_model_gamma_0_values,
        ip.single_term_model_exponent_values,
        ip.triple_term_model_gamma_0_values,
        ip.triple_term_model_constants, compression_ratio, ip.polynomial_coeff)

    write_output_files.run(debye_temperature, temperature, model_debye_temperatures, raw_pos_est, peak_centre, gsqr_per_angstrom, integrated_intensity, ln_intensity)

    plot_debye_waller.run(gsqr_per_angstrom, ln_intensity, raw_pos_est, current_md_temperature, ip.mass, ip.uncompressed_peak_positions)

    if ip.make_final_peak_plots is True:

        plot_peaks.run(peak_str, peak_centre, ip.source_name, ip.timestep)

    log.info("Path %s finished.\n", __name__)

    return
