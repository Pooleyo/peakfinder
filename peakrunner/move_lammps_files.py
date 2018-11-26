def run(current_lammps_directory, lammps_prefix, lammps_suffix, lammps_variable, copy_locations):

    from subprocess import call

    for i, variable in enumerate(lammps_variable):

        lammps_filename = lammps_prefix + variable + lammps_suffix

        start_point = current_lammps_directory + "/" + lammps_filename
        print start_point
        end_point = copy_locations[i] + "/lammps/" + lammps_filename
        print end_point
        call("mv " + start_point + " " + end_point, shell=True)

    return
