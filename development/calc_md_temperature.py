def run(lammps_file_location, user_input_temperature, temperature_dimensionality, atomic_mass, velocity_columns, number_velocity_bins):

    import units as un
    import numpy as np

    print "Calculating temperature from MD velocities..."

    md_temperature, velocity_squared = un.calc_MD_temperature(lammps_file_location, user_input_temperature, temperature_dimensionality, atomic_mass, velocity_columns)

    speed = np.sqrt(velocity_squared)

    histogram = un.bin_values(number_velocity_bins, speed)

    populations = list(histogram[0])
    total_population = float(sum(populations))

    normalised_populations = [0] * len(populations)

    for i, pop in enumerate(populations):

        normalised_populations[i] = pop / total_population

    bins = list(histogram[1])
    del(bins[-1])

    max_speed = np.max(speed)

    boltzmann_probability_list, boltzmann_speed_list = un.calc_maxwell_boltzmann_velocity_distribution(md_temperature, atomic_mass, max_speed, 10)

    un.plot_velocity_distribution(boltzmann_probability_list, boltzmann_speed_list, normalised_populations, bins)

    print "Temperature is: " + str(md_temperature) + " K"

    return md_temperature
