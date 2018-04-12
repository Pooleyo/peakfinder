import logging as log

def get_time():

    import time

    localtime = time.localtime()
    t0 = time.time()
    
    return t0, localtime
    
    
def build_all_k_values(gsqr_max, negative_k):

    # Builds all possible combinations of integers, up to the k-value given by sqrt(gsqr_max), rounded up. It only includes combinations that give gsqr < gsqr_max. Bool negative_k can be used to include/exclude negative k-values. 

    import numpy as np
    import math
    
    k_max = int(math.ceil(np.sqrt(gsqr_max)))
    
    if negative_k == True:
        k_values = range(-k_max, k_max + 1)
        
    elif negative_k == False:
        k_values = range(k_max + 1)
        
    else:
        print "Incorrect entry:\nnegative_k must be of bool type.\n"
        exit()
        
    pos = []
    
    for i in k_values:
        for j in k_values:
            for k in k_values:
                if i**2 + j**2 + k**2 <= gsqr_max:
                    pos.append([i, j, k])
                    
                else:

                    pass
    
    log.debug(pos)
                
    return pos
    

def remove_fcc_forbidden_reflections(old_pos):

    # Removes pos_est values that are forbidden in fcc crystals. Diffraction is allowed at positions where all the k-values are all-even or all-odd.

    new_pos = []

    for i in range(len(old_pos)):
        if old_pos[i][0] % 2 == 1 and old_pos[i][1] % 2 == 1 and old_pos[i][2] % 2 == 1: 
            new_pos.append(old_pos[i])
            
        elif old_pos[i][0] % 2 == 0 and old_pos[i][1] % 2 == 0 and old_pos[i][2] % 2 == 0: 
            new_pos.append(old_pos[i])
        
        else:
            pass
            
    log.debug(new_pos) 
    
    return new_pos
    

def remove_000(old_pos):
    
    # Removes [0, 0, 0] from pos.
    
    new_pos = []
    
    for i in old_pos:
        if i != [0, 0, 0]:
            new_pos.append(i)

    log.debug(new_pos)
    
    return new_pos
    
    
def get_gsqr_values(pos):

    # Calculates the value of G^2 for each position in pos.

    gsqr = []

    for i in pos:
        current_gsqr = (i[0] ** 2) + (i[1] ** 2) + (i[2] ** 2)
        
        gsqr.append(current_gsqr)
        
    log.debug(gsqr)
    
    return gsqr
    
    
def build_datafile_structure(pos):

    import os

    peak_str = []    
    
    for i in pos:
        current_peak_str = str(i[0]) + str(i[1]) + str(i[2])
        peak_str.append(current_peak_str)
        if not os.path.exists("data/" + current_peak_str):
            os.makedirs("data/" + current_peak_str)

    log.debug(peak_str)

    return peak_str


def calc_k_offset_with_N_atoms(N_atoms):

    offset = [1.0/N_atoms[0], 1.0/N_atoms[1], 1.0/N_atoms[2]]

    return offset
  
    
def convert_to_per_angstrom(element, a_lattice):
    
    import numpy as np

    element = np.asarray(element)
    converted_element = element * ( a_lattice / (2 *np.pi) )
    converted_element = list(converted_element)

    return converted_element
     
     
def make_peak_str(i):

    peak_str = str(i[0]) + str(i[1]) + str(i[2])
    
    return peak_str
        
        
def find_k_start(pos_element, offset):

    k_start = [pos_element[0] - offset[0], pos_element[1] - offset[1], pos_element[2] - offset[2]]

    return k_start
    

def find_k_stop(pos_element, offset):

    k_stop = [pos_element[0] + offset[0], pos_element[1] + offset[1], pos_element[2] + offset[2]]

    return k_stop


def determine_soh_input_file_location(peak_str):
    
    import os
    
    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/" + peak_str + ".in"
    
    return input_file_location
    
    
def write_soh_input_3DFT(source_location, file_destination, peak_str, mass, a_lattice, k_steps, k_start, k_stop):
    
    string_to_write = "VERBOSE 0\nFILE_TYPE lammps-multi\nDATA_FILE " + str(source_location) + "\nAPPEND_FILE_NAME " + str(peak_str) + "\nPLOT_OUTPUT pdf\nCOORDS_SCALED\nSET_MASS " + str(mass) + "\nSET_A_CELL " + str(a_lattice) + "\nCALC_3D_FT\nSET_KX " + str(k_start[0]) + " " + str(k_stop[0]) + " " + str(k_steps) + "\nSET_KY " + str(k_start[1]) + " " + str(k_stop[1]) + " " + str(k_steps) + "\nSET_KZ " + str(k_start[2]) + " " + str(k_stop[2]) + " " + str(k_steps)
    
    f= open(file_destination, "w") 
    f.write(string_to_write)
    f.close()
    
    log.debug("Unit finished, no internal output.")
    
    return string_to_write  

    
def run_soh(input_file_location):

    import subprocess

    shell_command = "sonOfHoward " + input_file_location
    
    subprocess.call(shell_command, shell=True)

    return  
