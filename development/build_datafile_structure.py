def run(pos):
    
    import units as un
    import logging as log

    log.info("Brick %s started.\n", __name__)

    print "Building file structure..."

    peak_str = un.build_datafile_structure(pos)
    
    log.info("Brick %s finished.\n", __name__)
    
    return peak_str
