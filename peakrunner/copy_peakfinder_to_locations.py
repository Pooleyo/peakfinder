def run(src_location, copy_location_list):

    import shutil
    import os

    for copy_location in copy_location_list:

        if os.path.isdir(copy_location) is True:

            print "WARNING: Copy location already exists. The command to copy to this location has been aborted."
            pass

        else:

            shutil.copytree(src_location, copy_location)

    return
