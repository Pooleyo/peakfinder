def run():

    import logging as log
    import shutil

    shutil.move("log.pkfd", "data/log.pkfd")
    shutil.move("speed_distribution_md_vs_boltzmann.png", "data/speed_distribution_md_vs_boltzmann.png")

    log.info("Peakfinder finalised.")

    print "Finished!"

    return
