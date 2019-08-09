def run(gsqr_max, negative_k, remove_000, crystal_type):

    import units as un
    import logging as log
    
    log.info("Brick %s started.\n", __name__)

    print "Selecting peak positions..."
    
    pos_est = un.build_all_k_values(gsqr_max, negative_k)

    if crystal_type == 'FCC':

        pos_est = un.remove_fcc_forbidden_reflections(pos_est)

    elif crystal_type == 'BCC':

        pos_est = un.remove_bcc_forbidden_reflections(pos_est)

    if remove_000 == True:

        pos_est = un.remove_000(pos_est)
    
    gsqr_est = un.get_gsqr_values(pos_est)
    
    
    log.info("Brick %s finished.\n", __name__)

    return pos_est, gsqr_est
