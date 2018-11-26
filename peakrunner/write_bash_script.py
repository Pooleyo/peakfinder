def run(bash_script_filename, peakfinder_directory_list):

    from subprocess import call

    f = open(bash_script_filename, 'w')

    f.write("#!/bin/bash\n")

    for peakfinder_directory in peakfinder_directory_list:

        cd_line = "cd " + peakfinder_directory + "\n"

        f.write(cd_line)

        python_line = "python peakfinder.py" + "\n"

        f.write(python_line)

    call("chmod 755 " + bash_script_filename, shell=True)

    return
