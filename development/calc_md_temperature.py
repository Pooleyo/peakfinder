def run(lammps_file_location, user_input_temperature, temperature_dimensionality, atomic_mass, velocity_columns, number_velocity_bins):

    import units as un
    import numpy as np

    print "Calculating temperature from MD velocities..."

    md_temperature, velocity_squared = un.calc_MD_temperature(lammps_file_location, user_input_temperature, temperature_dimensionality, atomic_mass, velocity_columns)

    histogram = un.bin_values(number_velocity_bins, np.sqrt(velocity_squared))

    populations = list(histogram[0])
    bins = list(histogram[1])
    del(bins[-1])

    un.plot_matplotlib(bins, populations, "histogram_of_atom_velocity.png", "Velocity Magnitude (Angstrom/ps)", "Bin Population", "Histogram of Atom Velocities")

    print "Temperature is: " + str(md_temperature) + " K"

    return md_temperature
