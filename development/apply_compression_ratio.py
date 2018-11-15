def run(pos_est, gsqr_est, compression_ratio):

    import units as un

    compressed_pos_est, compressed_gsqr_est = un.apply_compression_ratio_to_pos_est(pos_est, gsqr_est,
                                                                                    compression_ratio)

    return compressed_pos_est, compressed_gsqr_est