import in_peakrunner as ip

def run():

    import create_copy_location_list
    import copy_peakfinder_to_locations
    import write_bash_script
    import move_lammps_files

    copy_location_list = create_copy_location_list.run(
        ip.root_directory, ip.subdirectory_prefix, ip.subdirectory_suffix, ip.subdirectory_variable_list)

    copy_peakfinder_to_locations.run(ip.peakfinder_src_location, copy_location_list)

    write_bash_script.run(ip.bash_script_filename, copy_location_list)

    return

run()
