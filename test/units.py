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

    log.debug(offset)

    return offset
  
    
def convert_to_per_angstrom(element, a_lattice):
    
    import numpy as np

    element = np.asarray(element)
    converted_element = element * ( (2 * np.pi) / a_lattice )
    converted_element = list(converted_element)

    log.debug(converted_element)

    return converted_element
     
     
def make_peak_str(i):

    peak_str = str(i[0]) + str(i[1]) + str(i[2])
    
    log.debug(peak_str)
    
    return peak_str
        
        
def find_k_start(pos_element, offset):

    k_start = [pos_element[0] - offset[0], pos_element[1] - offset[1], pos_element[2] - offset[2]]

    log.debug(k_start)

    return k_start
    

def find_k_stop(pos_element, offset):

    k_stop = [pos_element[0] + offset[0], pos_element[1] + offset[1], pos_element[2] + offset[2]]

    log.debug(k_stop)

    return k_stop


def determine_soh_input_file_location(peak_str):
    
    import os
    
    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/" + peak_str + ".in"
    
    log.debug(input_file_location)
    
    return input_file_location
    
       
def write_soh_input_3DFT(source_name, file_destination, peak_str, mass, a_lattice, k_steps, k_start, k_stop):

    import os
	
    cwd = os.getcwd()
    source_location = cwd + "/lammps/" + source_name
    
    string_to_write = ("VERBOSE 0"
    + "\nFILE_TYPE lammps-multi"
    + "\nDATA_FILE " + str(source_location) 
   	+ "\nAPPEND_FILE_NAME " + str(peak_str) 
   	+ "\nPLOT_OUTPUT pdf"
   	+ "\nCOORDS_SCALED"
   	+ "\nSET_MASS " + str(mass) 
   	+ "\nSET_A_CELL " + str(a_lattice)
   	+ "\nCALC_3D_FT"
   	+ "\nSET_KX " + str(k_start[0]) + " " + str(k_stop[0]) + " " + str(k_steps)
   	+ "\nSET_KY " + str(k_start[1]) + " " + str(k_stop[1]) + " " + str(k_steps)
   	+ "\nSET_KZ " + str(k_start[2]) + " " + str(k_stop[2]) + " " + str(k_steps)
   	+ "\n"
    )
    
    f = open(file_destination, "w") 
    f.write(string_to_write)
    f.close()
    
    log.debug(string_to_write)
    
    return string_to_write  

    
def run_soh(input_file_location, num_cores):

    import subprocess
    
    shell_command = "mpiexec -np " + str(num_cores) + " sonOfHoward " + input_file_location
    
    subprocess.call(shell_command, shell=True)
    
    log.debug("sonOfHoward called using input file at " + input_file_location)
    
    return
    
    
def move_soh_output_to_peak_folder(peak_str, source_name, timestep):

    import shutil
    
    origin = "./lammps/" + source_name + "." + timestep + "." + peak_str + ".ft"
    destination = "./data/" + peak_str + "/"
    
    shutil.move(origin, destination)
    
    log.debug(origin + " moved to " + destination)
    
    return
    
    
def determine_soh_output_file_location(peak_str, source_name, timestep):
    
    import os
    
    cwd = os.getcwd()
    output_file_location = cwd + "/data/" + peak_str + "/" + source_name + "." + timestep + "." + peak_str + ".ft"
    
    log.debug(output_file_location)
    
    return output_file_location
    

def read_from_soh_output(filename):

    import numpy as np

    kx, ky, kz, intensity = np.loadtxt(filename, skiprows=1, usecols=(0,1,2,5), unpack=True)    
    soh_output = [kx, ky, kz, intensity]

    log.debug(soh_output)

    return soh_output
    
    
def find_point_of_max_height(soh_output):

    import numpy as np

    max_height_index = np.argmax(soh_output[3])

    point_of_max_height = [soh_output[0][max_height_index], soh_output[1][max_height_index], soh_output[2][max_height_index]]

    log.debug(point_of_max_height)
    
    return point_of_max_height
    

def calc_dvol(soh_output):
    
    k_step = len(soh_output[0]) ** (1.0/3.0)
    dkx = ( max(soh_output[0]) - min(soh_output[0]) ) / (k_step - 1)
    dky = ( max(soh_output[1]) - min(soh_output[1]) ) / (k_step - 1)
    dkz = ( max(soh_output[2]) - min(soh_output[2]) ) / (k_step - 1)
    dvol = dkx * dky * dkz

    log.debug(dvol)
    
    return dvol


def calc_integrated_intensity(soh_output, dvol):

    intensity_sum = sum(soh_output[3])
    integrated_intensity = dvol * intensity_sum

    log.debug(integrated_intensity)
    
    return integrated_intensity
    

def get_ln_intensity(intensity):

    import numpy as np

    ln_intensity = np.log(intensity)
    
    log.debug(ln_intensity)

    return ln_intensity
    
    
def calc_line_slope_and_constant(x, y):

    import numpy as np
    
    slope, constant = np.polyfit(x, y, 1, cov=False)

    log.debug("slope = " + str(slope) + "\nconstant = " + str(constant))

    return slope, constant
    

def calc_debye_waller_constant(m):

    from scipy.constants import h, N_A, k, pi

    debye_waller_constant = (10 ** 20) * 3 * (h ** 2) * N_A / (4 * (pi ** 2) * m * (10 ** -3) * k)
    
    return debye_waller_constant
            
        
def calc_debye_temperature(temperature, slope, debye_waller_constant):

    import numpy as np

    debye_temperature = np.sqrt( debye_waller_constant * temperature * abs( 1.0 / slope ) )
    
    print debye_temperature

    return debye_temperature
