def run(copy_location_list):

    import importlib

    for copy_location in copy_location_list:

        inpkfd_location = copy_location + "/inpkfd"

        inpkfd = importlib.import_module(inpkfd_location)

        print inpkfd

    return
