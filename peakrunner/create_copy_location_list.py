def run(root_directory, subdirectory_prefix, subdirectory_suffix, subdirectory_variable):

    copy_location_list = []

    for current_subdirectory_variable in subdirectory_variable:

        location = root_directory + "/" + subdirectory_prefix + current_subdirectory_variable + subdirectory_suffix + "/peakfinder"

        copy_location_list.append(location)

    return copy_location_list
