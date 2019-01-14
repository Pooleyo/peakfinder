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


def make_lineout_directory():

    import os

    os.makedirs("data/lineouts")

    log.debug("Made directory for lineout data.")

    return

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
        
        
def find_simple_k_start(pos_element, offset):

    k_start = [pos_element[0] - offset[0], pos_element[1] - offset[1], pos_element[2] - offset[2]]

    log.debug(k_start)

    return k_start
    

def find_k_stop(pos_element, offset):

    k_stop = [pos_element[0] + offset[0], pos_element[1] + offset[1], pos_element[2] + offset[2]]

    log.debug(k_stop)

    return k_stop



def calc_lineout_k_start_stop(uncompressed_peak_centre, under_shoot, over_shoot):

    k_start = [uncompressed_peak_centre[0] * under_shoot, uncompressed_peak_centre[1] * under_shoot, uncompressed_peak_centre[2] * under_shoot]
    k_stop = [uncompressed_peak_centre[0] * over_shoot, uncompressed_peak_centre[1] * over_shoot, uncompressed_peak_centre[2] * over_shoot]

    log.debug(k_start)
    log.debug(k_stop)

    return k_start, k_stop


def calc_peak_edge_k_start_stop(predicted_peak_centre, under_shoot, over_shoot):

    k_start = [predicted_peak_centre[0] - under_shoot, predicted_peak_centre[1] - under_shoot, predicted_peak_centre[2] - under_shoot]
    k_stop = [predicted_peak_centre[0] + over_shoot, predicted_peak_centre[1] + over_shoot, predicted_peak_centre[2] + over_shoot]

    log.debug(k_start)
    log.debug(k_stop)

    return k_start, k_stop


def determine_soh_compression_finding_input_file_location(direction):

    import os

    cwd = os.getcwd()
    input_file_location = cwd + "/data/lineouts/" + direction + "_lineout.in"

    log.debug(input_file_location)

    return input_file_location


def determine_soh_edge_finding_input_file_location(direction, peak_str):

    import os

    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/" + direction + "_lineout.in"

    log.debug(input_file_location)

    return input_file_location


def determine_accurate_soh_input_file_location(peak_str):
    
    import os
    
    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/" + peak_str + ".in"
    
    log.debug(input_file_location)
    
    return input_file_location


def determine_rough_soh_input_file_location(peak_str):
    import os

    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/find_centre_" + peak_str + ".in"

    log.debug(input_file_location)

    return input_file_location


def determine_rough_lineout_soh_input_file_location(peak_str):
    import os

    cwd = os.getcwd()
    input_file_location = cwd + "/data/" + peak_str + "/find_edge_" + peak_str + ".in"

    log.debug(input_file_location)

    return input_file_location


def write_soh_input_1DFT(source_name, file_destination, appended_string, mass, a_lattice, k_start, k_stop, k_steps):

    import os
    cwd = os.getcwd()
    source_location = cwd + "/lammps/" + source_name
    string_to_write = ("VERBOSE 0"
                       + "\nFILE_TYPE lammps-multi"
                       + "\nDATA_FILE " + str(source_location)
                       + "\nAPPEND_FILE_NAME " + str(appended_string)
                       + "\nPLOT_OUTPUT pdf"
                       + "\nCOORDS_SCALED"
                       + "\nSET_MASS " + str(mass)
                       + "\nSET_A_CELL " + str(a_lattice)
                       + "\nCALC_1D_FT"
                       + "\nSET_K_START " + str(k_start[0]) + " " + str(k_start[1]) + " " + str(k_start[2]) + " "
                       + "\nSET_K_STOP " + str(k_stop[0]) + " " + str(k_stop[1]) + " " + str(k_stop[2]) + " "
                       + "\nSET_NK " + str(k_steps)
                       + "\n"
                       )

    f = open(file_destination, "w")
    f.write(string_to_write)
    f.close()

    log.debug(string_to_write)

    return string_to_write


def write_soh_input_3DFT(source_name, file_destination, appended_string, mass, a_lattice, k_steps, k_start, k_stop):

    import os
	
    cwd = os.getcwd()
    source_location = cwd + "/lammps/" + source_name
    
    string_to_write = ("VERBOSE 0"
    + "\nFILE_TYPE lammps-multi"
    + "\nDATA_FILE " + str(source_location) 
   	+ "\nAPPEND_FILE_NAME " + str(appended_string)
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
    
    shell_command = "mpiexec -np " + str(num_cores) + " sonOfHoward " + input_file_location + " >/dev/null"
    
    subprocess.call(shell_command, shell=True)
    
    log.debug("sonOfHoward called using input file at " + input_file_location)
    
    return
    
    
def move_soh_accurate_output_to_peak_folder(peak_str, source_name, timestep):

    import shutil
    
    origin = "./lammps/" + source_name + "." + timestep + "." + peak_str + ".ft"
    destination = "./data/" + peak_str + "/"
    
    shutil.move(origin, destination)
    
    log.debug(origin + " moved to " + destination)
    
    return


def move_soh_rough_output_to_peak_folder(peak_str, appended_string, source_name, timestep):

    import shutil

    origin = "./lammps/" + source_name + "." + timestep + "." + appended_string + ".ft"
    destination = "./data/" + peak_str + "/"

    shutil.move(origin, destination)

    log.debug(origin + " moved to " + destination)

    return


def move_soh_output_to_lineout_folder(lineout_str, source_name, timestep):
    import shutil

    origin = "./lammps/" + source_name + "." + timestep + ".lineout_" + lineout_str + ".ft"
    destination = "./data/lineouts/"

    shutil.move(origin, destination)

    log.debug(origin + " moved to " + destination)

    return


def move_plot_output_to_peak_folder(direction, peak_str):

    import shutil

    origin = direction + ".png"
    destination = "./data/" + peak_str + "/"

    shutil.move(origin, destination)

    log.debug(origin + " moved to " + destination)

    return
    
    
def determine_accurate_soh_output_file_location(peak_str, source_name, timestep):
    
    import os
    
    cwd = os.getcwd()
    output_file_location = cwd + "/data/" + peak_str + "/" + source_name + "." + timestep + "." + peak_str + ".ft"
    
    log.debug(output_file_location)
    
    return output_file_location


def determine_rough_soh_output_file_location(peak_str, source_name, timestep):
    import os

    cwd = os.getcwd()
    output_file_location = cwd + "/data/" + peak_str + "/" + source_name + "." + timestep + ".find_centre_" + peak_str + ".ft"

    log.debug(output_file_location)

    return output_file_location


def determine_soh_1DFT_output_file_location(direction_str, source_name, timestep):
    import os

    cwd = os.getcwd()
    output_file_location = cwd + "/data/lineouts/" + source_name + "." + timestep + ".lineout_" + direction_str + ".ft"

    log.debug(output_file_location)

    return output_file_location


def determine_soh_edge_finding_output_file_location(peak_str, direction_str, source_name, timestep):
    import os

    cwd = os.getcwd()
    output_file_location = cwd + "/data/" + peak_str + "/" + source_name + "." + timestep + "." + peak_str + "_find_edges_" + direction_str + ".ft"

    log.debug(output_file_location)

    return output_file_location


def determine_edge_finding_soh_output_file_location(direction_str, source_name, timestep):
    import os

    cwd = os.getcwd()
    output_file_location = cwd + "/data/lineouts/" + source_name + "." + timestep + ".lineout_" + direction_str + ".ft"

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

    log.debug(debye_waller_constant)

    return debye_waller_constant
            
        
def calc_debye_temperature_xrd(temperature, slope, debye_waller_constant):

    import numpy as np

    debye_temperature = np.sqrt(debye_waller_constant * temperature * abs( 1.0 / slope ))

    log.debug(debye_temperature)

    return debye_temperature


def calc_debye_temperature_from_single_term_gruneisen_model(debye_temperautre_300K_uncompressed, initial_volume, final_volume, gamma_uncompressed, exponent):

    import numpy as np

    #  See Will Murphy PHYSICAL REVIEW B 78, 014109 (2008) for the source of this equation.

    exponent_term = - (gamma_uncompressed / (exponent * initial_volume)) * ((final_volume ** exponent) - (initial_volume ** exponent))

    correction_factor = np.exp(exponent_term)

    model_debye_temperature = debye_temperautre_300K_uncompressed * correction_factor

    log.debug(model_debye_temperature)

    return model_debye_temperature


def calc_debye_temperature_from_triple_term_gruneisen_model(debye_temperature_300K_uncompressed, initial_volume, final_volume, gamma_uncompressed, constants):

    import numpy as np

    #  See Will Murphy PHYSICAL REVIEW B 78, 014109 (2008) for the source of this equation.

    constant_term_1 = gamma_uncompressed - constants[0] + constants[1] - constants[2]

    volume_term_1 = np.log(final_volume) - np.log(initial_volume)

    constant_term_2 = - constants[0] + 2 * constants[1] - 3 * constants[2]

    volume_term_2 = initial_volume * ((1 / final_volume) - (1 / initial_volume))

    constant_term_3 = - (constants[1] / 2.0) + (3 * constants[2] / 2.0)

    volume_term_3 = (initial_volume ** 2) * ((1 / (final_volume ** 2)) - (1 / (initial_volume ** 2)))

    constant_term_4 = - constants[2] / 3.0

    volume_term_4 = (initial_volume ** 3) * ((1 / (final_volume ** 3)) - (1 / (initial_volume ** 3)))

    exponent_term = (constant_term_1 * volume_term_1) + (constant_term_2 * volume_term_2) + (constant_term_3 * volume_term_3) + (constant_term_4 * volume_term_4)

    correction_factor = np.exp(- exponent_term)

    model_debye_temperature = debye_temperature_300K_uncompressed * correction_factor

    log.debug(model_debye_temperature)

    return model_debye_temperature


def calc_volume_lattice_units(a_lattice, compression_factors):

    volume = a_lattice ** 3 * (compression_factors[0] * compression_factors[1] * compression_factors[2])

    log.debug(volume)

    return volume


def calc_temperature_xrd(debye_temperature, slope, debye_waller_constant):

    temperature = (debye_temperature ** 2) * abs(slope) * (1.0 / debye_waller_constant)

    log.debug(temperature)

    return temperature


def plot_matplotlib(x, y, filename, x_label, y_label, plot_title):

    import matplotlib.pyplot as plt

    plt.scatter(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.savefig(filename)
    plt.close()

    return


def plot_pyqtgraph(x, y, filename):

    import pyqtgraph as pg
    import pyqtgraph.exporters

    class MyPlotClass():

        def __init__(self):
            self.windowplt = pg.plot()
            self.windowplt.win.hide()

        def savePlot(self, x, y, filename):
            self.windowplt.plot(x, y)
            exporter = pg.exporters.ImageExporter(self.windowplt.plotItem)
            exporter.params.param('width').setValue(256, blockSignal=exporter.widthChanged)
            exporter.params.param('height').setValue(256, blockSignal=exporter.heightChanged)
            exporter.export(filename)

    save_plot = MyPlotClass()
    save_plot.savePlot(x, y, filename)

    return


def plot_pygnuplot(x, y, filename, data_filename):

    import PyGnuplot as gnu

    gnu.s([x,y], data_filename)
    gnu.c('set terminal pngcairo size 350,262 enhanced font "Verdana,10"')
    gnu.c('set output "' + filename + '"')
    gnu.c('plot "' + data_filename + '" pt 1 ps 0.5')

    return


def find_line_data_from_3DFT(constant_axes, variable_axis, centre_point, soh_output):

    constant_value_0 = centre_point[constant_axes[0]]

    constant_value_1 = centre_point[constant_axes[1]]

    line_points = []

    line_intensity = []

    for i, intensity in enumerate(soh_output[3]):

        if soh_output[constant_axes[0]][i] == constant_value_0 and soh_output[constant_axes[1]][i] == constant_value_1:

            line_k = soh_output[variable_axis][i]

            line_points.append(line_k)

            line_intensity.append(intensity)

    return line_points, line_intensity


def write_temperatures_to_file(debye_temperature, temperature, filename_temperatures):

    f = open(filename_temperatures, "w")
    f.write(
        "Debye temperature\t\t\t\t" + str(debye_temperature) + "\n"
        "Temperature\t\t\t\t\t\t" + str(temperature)
    )
    f.close()

    log.debug(filename_temperatures)

    return


def write_peak_intensities_to_file(pos_est, peak_centre, gsqr, integrated_intensity, ln_intensity, filename):

    header_string = "peak_name peak_centre gsqr integrated_intensity ln_intensity\n"

    f = open(filename, "w")
    f.write(header_string)
    for i, pos in enumerate(pos_est):
        f.write("%s %s %s %s %s\n" % (pos, peak_centre[i], gsqr[i], integrated_intensity[i], ln_intensity[i]))
    f.close()

    log.debug(filename)

    return


def find_if_vectors_parallel(v_1, v_2):

    import numpy as np
    import math

    length_1 = np.linalg.norm(v_1)

    length_2 = np.linalg.norm(v_2)

    if length_1 < 0.001:

        result = False

    elif length_2 < 0.001:

        result = False

    else:

        normalised_1 = v_1 / length_1

        normalised_2 = v_2 / length_2

        dot_prod = np.dot(normalised_1, normalised_2)

        if math.isnan(dot_prod) is True:

            result = False

        elif int(dot_prod) is 1:

            result = True

        else:

            result = False

    log.debug(result)

    return result


def calc_compression_ratio(compressed_k, uncompressed_k):

    compression_ratio = compressed_k/uncompressed_k

    log.debug(compression_ratio)

    return compression_ratio


def apply_compression_ratio_to_pos_est(pos_est, gsqr_est, compression_ratio):

    import copy

    compressed_pos_est = copy.deepcopy(pos_est)
    compressed_gsqr_est = copy.deepcopy(gsqr_est)

    for i, pos in enumerate(compressed_pos_est):

        for j, compression in enumerate(compression_ratio):

            pos[j] = pos[j] * compression

        compressed_gsqr_est[i] = (pos[0] ** 2) + (pos[1] ** 2) + (pos[2] ** 2)

    return compressed_pos_est, compressed_gsqr_est


def find_k_start_stop_for_peak_from_first_minima(k_data, intensity):

    centre_index = len(k_data)/2

    for i, k in enumerate(k_data):

        if centre_index - i == 0:

            print "Couldn't find intensity minimum for k_start."

            exit()

        intensity_diff = intensity[centre_index - i] - intensity[centre_index - i - 1]

        if intensity_diff >= 0.0:

            continue

        elif intensity_diff < 0.0:

            k_start = k_data[centre_index - i]

            break

    for i, k in enumerate(k_data):

        if centre_index + i == (len(k_data) - 1):

            print "Couldn't find intensity minimum for k_stop."

            exit()

        intensity_diff = intensity[centre_index + i] - intensity[centre_index + i + 1]

        if intensity_diff >= 0.0:

            continue

        elif intensity_diff <= 0.0:

            k_stop = k_data[centre_index + i]

            break

    return k_start, k_stop


def calc_overstepped_k_start_k_stop(pos, undershoot, overshoot):

    k_start = [0.0, 0.0, 0.0]
    k_stop = [0.0, 0.0, 0.0]

    for i, k in enumerate(pos):

        k_start[i] = k - undershoot[i]
        k_stop[i] = k + overshoot[i]

    return k_start, k_stop


def calc_MD_temperature(lammps_file_location, user_input_temperature, temperature_dimensionality, atomic_mass, velocity_columns):

    import numpy as np
    from scipy.constants import k

    try:

        vx, vy, vz = np.loadtxt("lammps/" + lammps_file_location, skiprows=9, usecols=(velocity_columns[0], velocity_columns[1], velocity_columns[2]), unpack=True)

        number_of_atoms = len(vx)

        velocity_squared = [0] * number_of_atoms

        if temperature_dimensionality is 2:

            for i in range(len(vx)):

                velocity_squared[i] = (vx[i] ** 2) + (vy[i] ** 2)

        elif temperature_dimensionality is 3:

            for i in range(len(vx)):

                velocity_squared[i] = (vx[i] ** 2) + (vy[i] ** 2) + (vz[i] ** 2)

        velocity_squared_sum = sum(velocity_squared)

        MD_temperature = (1.660539040e-27 * 10000 * atomic_mass * velocity_squared_sum) / (temperature_dimensionality * number_of_atoms * k) # The number 1.66054e-27 is to convert the atomic mass from amu to kg. The factor of 100 is to conver the velocities from  Angstrom/ps to m/s.

    except:

        print "######### WARNING: calc_MD_temperature: Could not load values for velocity.\n\tThis could be due to " \
              + "the LAMMPS file not having enough columns. Check the LAMMPS file has the velocities set to column " \
              + str(velocity_columns[0]) + ", " + str(velocity_columns[1]) + ", and " \
              + str(velocity_columns[2]) + " (where 0 corresponds to the first column).\n\tThe user defined temperature" \
              + " will be used instead: " + str(user_input_temperature) + " K"

        MD_temperature = user_input_temperature

    log.debug(MD_temperature)

    return MD_temperature, velocity_squared

def bin_values(number_of_bins, list_to_bin):

    import numpy as np

    histogram = np.histogram(list_to_bin, number_of_bins)

    log.debug(histogram)

    return histogram

def plot_histogram(histogram, filename, x_label, y_label, plot_title):

    import matplotlib.pyplot as plt

    plt.hist(histogram)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.savefig(filename)
    plt.close()

    return
