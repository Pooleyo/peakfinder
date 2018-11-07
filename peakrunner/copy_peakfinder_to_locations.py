def run(src_location, copy_location_list):

    import shutil

    for copy_location in copy_location_list:

        shutil.copytree(src_location, copy_location)

    return
