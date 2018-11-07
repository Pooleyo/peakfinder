import in_peakrunner as ip

def run():

    import create_copy_location_list
    import copy_peakfinder_to_locations
    import rewrite_inpkfd

    copy_location_list = create_copy_location_list.run(
        ip.root_directory, ip.subdirectory_prefix, ip.subdirectory_suffix, ip.subdirectory_variable_list)

    print copy_location_list

    copy_peakfinder_to_locations.run(ip.peakfinder_src_location, copy_location_list)

    rewrite_inpkfd.run(copy_location_list)

    return

run()
