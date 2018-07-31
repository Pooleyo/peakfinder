def run():

    import units as un
    import inpkfd as ip
    import logging as log
    import shutil
    import os
    
    # Sets up the log system.
    log.basicConfig(filename = "log.pkfd", filemode = "w", level = log.DEBUG, format = "%(asctime)s\t\t%(filename)s\t\t%(funcName)s\n\t %(message)s")
    log.info("Peakfinder intialised.\n")
    
    # Removes previously created folders.
    if os.path.exists("data"):
        shutil.rmtree("data")

    return
