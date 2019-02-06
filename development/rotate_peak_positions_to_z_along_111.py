def run(pos_est, gsqr_est):

    import units as un

    print "Rotating peak positions so that z lies along the 111 direction..."

    rot_x, rot_z = un.create_rotation_matrix_for_111_rotation()

    rotated_pos_est = un.rotate_pos_est_using_rotation_matrices(pos_est, rot_x, rot_z)

    rotated_gsqr_est = un.get_gsqr_values(rotated_pos_est)

    return rotated_pos_est, rotated_gsqr_est
