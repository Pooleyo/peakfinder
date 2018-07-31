#Each function should have its own ability to write to the log file.
# I also want it to save plots of all the potentially interesting data.


########################################################################
# Opens a time log of the run.

def startwatch():

	import time

	start_time = time.localtime()
	tpy0 = time.clock()
	t0 = time.time()

	t = open('time.pkfd', 'w')
	t.write("Peakfinder started at " + str(start_time[3]) + ":" + str(start_time[4]) + ":" + str(start_time[5]) + " " + str(start_time[2]) + "/" + str(start_time[1]) + "/" + str(start_time[0]) + "\n")
	t.close()
	
	return t0, tpy0
	
	
########################################################################
# Closes the time log of the run.

def stopwatch(t0, tpy0):

	import time

	stop_time = time.localtime()
	tpyf = time.clock()
	tf = time.time()
	
	tpyt = tpyf - tpy0
	tt = tf - t0
	
	
	hours = int(tt/3600)
	minutes = int((tt - (hours * 3600))/60)
	seconds = int(tt - (hours * 3600) - (minutes * 60))
	
	time_elapsed = [hours, minutes, seconds]	

	t = open('time.pkfd', 'a')
	t.write("\n\nPeakfinder finished at " + str(stop_time[3]) + ":" + str(stop_time[4]) + ":" + str(stop_time[5]) + " " + str(stop_time[2]) + "/" + str(stop_time[1]) + "/" + str(stop_time[0]) + "\n")
	t.write("\nPeakfinder took " + str(tt) + " s (or " + str(time_elapsed[0]) + "h " + str(time_elapsed[1]) + "m " + str(time_elapsed[2]) + "s) to complete.")
	t.write("\n\nThe python portion took " + str(tpyt) + " s.")
	t.close()
	
	print "\nPeakfinder took " + str(tt) + " s (or " + str(time_elapsed[0]) + "h " + str(time_elapsed[1]) + "m " + str(time_elapsed[2]) + "s) to complete."

	return
	
	

########################################################################
#This function creates bcc positions. It takes as input range_num (int), negative_k (bool), and remove_000 (bool). As output it creates gsqr_est (list) and pos_est (list).

def make_bcc(range_num, negative_k, remove_000): 


	# Variables in this function:
	# negative_k -> this bool determines whether negative values of hkl are included.
	# x_est/y_est/z_est -> these are sets of integers that make up the estimates of each coordinate in reciprocal space.
	# range_num -> this integer determines how far into reciprocal space the peaks positions are estimated.
	# pos_est -> this list contains all the k coordinates of each reciprocal peak.
	# gsqr_est -> this list contains all the G^2 values for each peak contained in pos_est
	# h/k/l -> these are iteration variables used to cycle through each value of x_est/y_est/z_est.
	# remove_000 -> this bool determines whether the 000 peak is precluded from the list of pos_est.
	# gsqr_temp -> this temporary variable is used to hold the h^2 + k^ + l^2 calculation needed to create G^2, before being added to gsqr_est.
	# i -> used as a looping variable to write the log file and print to the console.



	if negative_k == True:
	
		x_est = range(-range_num+1, range_num)
		y_est = range(-range_num+1, range_num)
		z_est = range(-range_num+1, range_num)

		
	elif negative_k == False:
	
		x_est = range(0, range_num)
		y_est = range(0, range_num)
		z_est = range(0, range_num)



	pos_est = [] #This list will have our k coordinates for each peak.
	gsqr_est = [] #This list will have the G^2 values for each peak.

	
	for h in x_est:

		for k in y_est:

			for l in z_est:
			
				if remove_000 == False:
				
					#Here the positions are only accepted if h + k + l are even. If that is not the case, they are not entered into the pos_est list.

					if (h + k + l) % 2 == 0:

						gsqr_temp = (h * h) + (k * k) + (l * l)
						gsqr_est.append(gsqr_temp)
						pos_est.append([h, k, l]) 

					else:

						pass 
					
						
				elif remove_000 == True:
					
					#Here the positions are only accepted if h + k + l are even. If h = k = l = 0, the peak will not be written. If these conditions are not met the peak is not entered into the pos_est list.
					
					if (h + k + l) % 2 == 0 and abs(h) + abs(k) + abs(l) != 0:

						gsqr_temp = (h * h) + (k * k) + (l * l)
						gsqr_est.append(gsqr_temp)
						pos_est.append([h, k, l]) 

					else:

						pass 	

		
	
	for i in range(len(pos_est)):
		print "\nPeak " + str(i+1) + " of " + str(len(pos_est)) + " estimated: " + str(pos_est[i]) + " with G^2 = " + str(gsqr_est[i])
	
	
	#This part makes the log entry for the function.	
	f = open("log.pkfd", "w")
	
	f.write("\nFunction make_bcc called with input:\n"
	"range_num = " + str(range_num) + "\n"
	"negative_k = " + str(negative_k) + "\n"
	"remove_000 = " + str(remove_000) + "\n"
	"\nFunction make_bcc returned:\n")
	for i in range(len(pos_est)):
		f.write( "Peak " + str(i+1) + " of " + str(len(pos_est)) + " estimated: " + str(pos_est[i]) + " with G^2 = " + str(gsqr_est[i]) + "\n")
	
	f.close()
	
	
	return gsqr_est, pos_est
	
	

####################################################################
#This function creates bcc positions. It takes as input range_num (int), negative_k (bool), and remove_000 (bool). As output it creates gsqr_est (list) and pos_est (list).


def make_fcc(gsqr_max, negative_k, remove_000):

	
	# Variables in this function:
	# negative_k -> this bool determines whether negative values of hkl are included.
	# x_est/y_est/z_est -> these are sets of integers that make up the estimates of each coordinate in reciprocal space.
	# range_num -> this integer determines how far into reciprocal space the peaks positions are estimated.
	# pos_est -> this list contains all the k coordinates of each reciprocal peak (NOT in units of A^-1).
	# gsqr_est -> this list contains all the G^2 values for each peak contained in pos_est (note that this is NOT in units of A^-2)
	# h/k/l -> these are iteration variables used to cycle through each value of x_est/y_est/z_est.
	# remove_000 -> this bool determines whether the 000 peak is precluded from the list of pos_est.
	# gsqr_temp -> this temporary variable is used to hold the h^2 + k^ + l^2 calculation needed to create G^2, before being added to gsqr_est.
	# i -> used as a looping variable to write the log file and print to the console.



	import numpy as np
	import time
	
	t0 = time.time()
	
	
	range_num = int(np.sqrt(gsqr_max) + 1.0)

	
	
	if negative_k == True:
	
		x_est = range(-range_num+1, range_num)
		y_est = range(-range_num+1, range_num)
		z_est = range(-range_num+1, range_num)


	else:
	
		x_est = range(range_num)
		y_est = range(range_num)
		z_est = range(range_num)




	pos_est = [] # This list will have our k coordinates for each peak.
	gsqr_est = [] # This list will have the G^2 values for each peak.
	

	for i in x_est:

		for j in y_est:

			for k in z_est:
			
				#The values for i j k are only selected if they are all even or all odd. If remove_000 is true there is an extra condition that makes sure 000 is not included.
				if remove_000 == True:
				
					if i % 2 == 0 and j % 2 == 0 and k % 2 == 0 and abs(i) + abs(j) + abs(k) != 0:
				
						l = (i * i) + (j * j) + (k * k)
						gsqr_est.append(l)
						pos_est.append([i, j, k]) 

					elif i % 2 == 1 and j % 2 == 1 and k % 2 == 1:
				
						l = (i * i) + (j * j) + (k * k)
						gsqr_est.append(l)
						pos_est.append([i, j, k]) 

					else:
						pass   
				
				#This part is triggered if we want to keep the 000 peak.		
				elif remove_000 == False:
				
					if i % 2 == 0 and j % 2 == 0 and k % 2 == 0:
				
						l = (i * i) + (j * j) + (k * k)
						gsqr_est.append(l)
						pos_est.append([i, j, k]) 

					elif i % 2 == 1 and j % 2 == 1 and k % 2 == 1:
				
						l = (i * i) + (j * j) + (k * k)
						gsqr_est.append(l)
						pos_est.append([i, j, k]) 

					else:
					
						pass 

	
	i = 0
	
	print "Removing peaks with too large gsqr..."
	
	while i + 1 <= len(pos_est):
			
		if gsqr_est[i] <= gsqr_max:
		
			i += 1
			continue
		
			
		else:
		
			del gsqr_est[i]
			del pos_est[i]
			continue
  	

	
	# This section prints to the console.
	
	print "\nPeaks estimated for a fcc structure:"
	
	for i in range(len(pos_est)):
		print "Peak " + str(i+1) + " of " + str(len(pos_est)) + ": " + str(pos_est[i]) + " with G^2 = " + str(gsqr_est[i])
	
	
	#This part makes the log entry for the function.	
	f = open("log.pkfd", "w")
	
	f.write("\nFunction make_fcc called with input:\n"
	"range_num = " + str(range_num) + "\n"
	"negative_k = " + str(negative_k) + "\n"
	"remove_000 = " + str(remove_000) + "\n"
	"\nFunction make_fcc returned:\n")
	for i in range(len(pos_est)):
		f.write( "Peak " + str(i+1) + " of " + str(len(pos_est)) + " estimated: " + str(pos_est[i]) + " with G^2 = " + str(gsqr_est[i]) + "\n")
	
	f.close()
	
	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.make_fcc took \t\t\t\t" + str(tt) + " s to complete.")	
	
	return gsqr_est, pos_est
	
	
##################################################################

# This function rotates the estimated pos_est such that the 111 direction is parallel to the z-direction. This is necessary since, in LAMMPS, the way that compression along z is achieved is by rotating the crystal in this way and then simply compressing in the z-direction. The result is that all of our atom positions are rotated and thus the expected peak positions. It takes as input: pos_est (string).


def enforce_rotation_111(pos_est):



	# Variables:
	# pos_est -> estimated peak positions from, for example, the make_bcc function.
	# theta_x -> variable used to calculate the rotation matrix around x-axis.
	# theta_z -> variable used to calculate the rotation matrix around z-axis.
	# rot_x -> rotation matrix about x-axis.
	# rot_z -> rotation matrix about z-axis.
	# rot_pos_est -> contains the rotates position_estimates.
	# new, new_2 -> intermediate variable to hold the partially rotated pos_est.
	
	

	import numpy as np
	import time
	
	t0 = time.time()
	

	theta_x = np.pi/2.0 - np.arctan(1/np.sqrt(2)) # Calculated by looking at a 111 vector which has been rotated by 45 degrees around the z-axis.
	theta_z = np.pi/4.0 # Rotate around z-axis by 45 degrees.


	rot_x = np.array([[1,0,0], [0, np.cos(theta_x), -np.sin(theta_x)], [0, np.sin(theta_x), np.cos(theta_x)]]) # creates the array which rotates the positions around the x-axis i.e. rotation matrix
	rot_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0], [np.sin(theta_z), np.cos(theta_z), 0], [0, 0, 1]]) # same as above, this time around z. Note that we won't create a y-rotation matrix since it isn't needed in this instance of 111.



	rot_pos_est = [0] * len(pos_est)

	pos_est = np.asarray(pos_est) # Converts pos_est to an array so we can multiply it by our rotational matrices.


	print "\nPeak estimates rotated:"

	for i in range(len(pos_est)): # Loops over all peaks to populate the pos_est list with rotated peak positions.

		new = np.dot(rot_z, pos_est[i]) # First matrix multiply the z-rotation with the original position estimate to get "new", an intermediate variable.
		new_2 = np.dot(rot_x, new) # Then matrix multiply the x-rotation with "new" to get the array version of the rotated peak.
		rot_pos_est[i] = list(new_2) # Convert this to a list (for compatability with the rest of the code).
		print "Peak " + str(i + 1) + " of " + str(len(pos_est)) + " at " + str(pos_est[i]) + " rotated to " + str(rot_pos_est[i])





	#This part makes the log entry for the function.	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction enforce_rotation_111 called with input:\n"
	"pos_est = (listed below with rotated counterpart)\n"	
	"\nFunction enforce_rotation_111 returned:\n")
	
	for i in range(len(pos_est)):
		f.write("Peak " + str(i+1) + " of " + str(len(pos_est)) + ": " + str(pos_est[i]) + " rotated to " + str(rot_pos_est[i]) + "\n")

	
	f.close()
	
	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.enforce_rotation_111 took \t\t\t" + str(tt) + " s to complete.")	

	return rot_pos_est
	
	
	
################################################################
	
# This function cuts up a LAMMPS file to atoms contained within a certain volume, as defined by the user. It takes as input: source (string), xlo, xhi,ylo, yhi, zlo, zhi (all floats).


	

def cut_atoms(source, xlo, xhi, ylo, yhi, zlo, zhi):


	#Variables:
	# source -> the original LAMMPS file to be cut up.
	# intermediate_file -> this file is an intermediary to the final cut .atom file. It only differes from the cut .atom file by the number of atoms (which is reported incorrectly in the intermediate file).
	# cut_atom_filename -> the name of the file after it has been cut.
	# while_initialiser -> variable for looping through a while loop.
	# current_line_list -> list of all the words in a line.
	# xs_ind, ys_ind, zs_ind, vx_ind, vy_ind, vz_ind, fx_ind, fy_ind, fz_ind -> stores the indices of the columns which contain the atomic coordinates (xs, ys, zs), the velocity components of each atom (vx, vy, vz), and the force components on each atom (fx, fy, fz).
	# xlo, xhi, ylo, yhi, zlo, zhi -> floats that describe the boundary of the volume to be cut to. These values are between 0 and 1 (inclusive) and are fractions of the box dimensions.
	# atom_counter -> counts the number of atoms in the new cut .atom file.


	import time
	
	t0 = time.time()
	

	intermediate_file = "intermediate.file"

	cut_atom_filename = "cut_" + source



	with open(source, 'r') as f:


		


		g = open(intermediate_file, 'w')

			
		while_initialiser = True
		
		
		while while_initialiser == True:


			for line in f: 

			
				g.write(line)

			
				current_line_list = line.split()
				
								
				for i in range(len(current_line_list)):

		
					if current_line_list[i] == "xs":

					
						xs_ind = current_line_list.index("xs") - 2

						
						ys_ind = current_line_list.index("ys") - 2

						
						zs_ind = current_line_list.index("zs") - 2
						

						vx_ind = current_line_list.index("vx") - 2
					

						vy_ind = current_line_list.index("vy") - 2
						

						vz_ind = current_line_list.index("vz") - 2
						

						fx_ind = current_line_list.index("fx") - 2
					
	
						fy_ind = current_line_list.index("fy") - 2
						

						fz_ind = current_line_list.index("fz") - 2
						
					
						while_initialiser = False
						
						
						break

				
				
					else:


						continue


				break
			
		
		atom_counter = 0
			
		for line in f: 


			current_line_list = line.split()

			
			if float(current_line_list[fx_ind]) == 0 and float(current_line_list[fy_ind]) == 0 and float(current_line_list[fz_ind]) == 0: # This line gets rid of any of the piston and back wall atoms.
	
				continue
			
			# The following six conditions are triggered if the atoms are outside of the user-defined volume.
			
			if float(current_line_list[xs_ind]) < xlo:
			
				continue
				
			if float(current_line_list[ys_ind]) < ylo:
			
				continue
				
			if float(current_line_list[zs_ind]) < zlo:
			
				continue
			
			if float(current_line_list[xs_ind]) > xhi:
			
				continue
				
			if float(current_line_list[ys_ind]) > yhi:
			
				continue
				
			if float(current_line_list[zs_ind]) > zhi:
			
				continue	
			
			# If none of the above conditions are triggered, the current line in the lammps dump file is written to the intermediate file.
			
			else:
			
				g.write(line)
				
				atom_counter += 1
				
					
		g.close()


		
	with open(intermediate_file, 'r') as g:
	
	
		number_of_atoms = False
		
	
		h = open(cut_atom_filename, 'w')
		
	
		for line in g:
		
			current_line_list = line.split()
			
			
			
			if number_of_atoms == True:
			
				h.write(str(atom_counter) + "\n")
				
				number_of_atoms = False
				
				continue
			
		
		
		
			if len(current_line_list) == 4:
			
		
				if current_line_list[1] == "NUMBER" and current_line_list[2] == "OF" and current_line_list[3] == "ATOMS":
		
					h.write(line)
			
					number_of_atoms = True	
					

				
				else:
				
					h.write(line)
		
		
		
			else:
		
				h.write(line)



	h.close()

	g.close()




	lammps_file_name = cut_atom_filename
	
	
	
			
	print "\nThe lammps file has been cut to the user defined parameters. \nThe filename of this abridged list is: " + str(lammps_file_name)
	print "\nThe region to which the positions have been cut are bounded by: \n xlo = " + str(xlo) + "\n xhi = " + str(xhi) + "\n ylo = " + str(ylo) + "\n yhi = " + str(yhi) + "\n zlo = " + str(zlo) + "\n zhi = " + str(zhi)
	print "\nAll atoms outside of this volume have been removed."
	print "There are now " + str(atom_counter) + " atoms remaining after the cut."
	
	




	#This part makes the log entry for the function.	
	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction cut_atoms called with input:\n"
	"source = " + str(source) + "\n"
	"xlo = " + str(xlo) + "\n"	
	"xhi = " + str(xhi) + "\n"
	"ylo = " + str(ylo) + "\n"
	"yhi = " + str(yhi) + "\n"
	"zlo = " + str(zlo) + "\n"
	"zhi = " + str(zhi) + "\n"
	"\nFunction cut_atoms returned:\n"
	"lammps_file_name = " + str(lammps_file_name) + "\n"
	"atom_counter = " + str(atom_counter) + "\n")

	
	f.close()


	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.cut_atoms took \t\t\t\t" + str(tt) + " s to complete.")	


	return lammps_file_name, atom_counter


#################################################################

# This function calculates the exact temperature of the MD based on the velocities of the atoms. 



def get_md_temperature(source, mass, piston_velocity):


	import numpy as np
	import scipy.constants as codata
	import time
	
	t0 = time.time()
	
	atoms_vx, atoms_vy, atoms_vz = np.loadtxt(source, skiprows=9, usecols=(5, 6, 7), unpack = True) #Atom velocities are loaded from the .atom or .atomcut file.


	atoms_vz_minus_piston = atoms_vz - piston_velocity #The piston velocity is subtracted from the z-velocity. The piston velocity is given in angstroms per ps. 
	

	v_sqr_2d = ((atoms_vx ** 2) + (atoms_vy ** 2)) * (10e12 ** 2)/(10e10 ** 2) #Since the velocities are in units of angstroms/ps(assuming units = metals is used in lammps), the square of the velocities is adjusted here to be in m/s. 
		
	v_sqr_3d = ( (atoms_vx ** 2) + (atoms_vy ** 2) + (atoms_vz_minus_piston ** 2) ) * (10e12 ** 2)/(10e10 ** 2) 




	E_k_atoms_2d = v_sqr_2d * 0.5 * mass * (10 ** -3) / codata.value("Avogadro constant")

	E_k_atoms_3d = v_sqr_3d * 0.5 * mass * (10 ** -3) / codata.value("Avogadro constant") 



	E_k_average_2d = np.mean(E_k_atoms_2d)

	E_k_average_3d = np.mean(E_k_atoms_3d)

	

	md_temperature_2d = (2.0/2.0) * (E_k_average_2d / codata.value("Boltzmann constant"))
	
	md_temperature_3d = (2.0/3.0) * (E_k_average_3d / codata.value("Boltzmann constant"))



	print "\nThe temperature has been recalculated for the atom data in " + str(source)
	print "\nThe 2D temperature has been calculated to be: " + str(md_temperature_2d) + " K"
	print "\nThe 3D temperature has been calculated to be: " + str(md_temperature_3d) + " K"
	


	#This part makes the log entry for the function.	
	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction get_md_temperature called with input:\n"
	"source = " + str(source) + "\n"
	"mass = " + str(mass) + "\n"	
	"piston_velocity = " + str(piston_velocity) + "\n"

	"\nThe temperature has been recalculated for the atom data in " + str(source) + ""
	"\nThe 2D temperature has been calculated to be: " + str(md_temperature_2d) + " K"
	"\nThe 3D temperature has been calculated to be: " + str(md_temperature_3d) + " K\n"	

	"\nFunction get_md_temperature returned:\n"
	"md_temperature_2d = " + str(md_temperature_2d) + "\n"
	"md_temperature_3d = " + str(md_temperature_3d) + "\n"
	)
	

	
	f.close()

	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.get_md_temperature took \t\t\t" + str(tt) + " s to complete.")	


	return md_temperature_2d, md_temperature_3d



###############################################################

# This function compares the expected peak positions to the actual positions in each dimension (first peak along x, y, and z), and then adjusts the predicted positions of all the peaks based on this. As input it takes: source (string), rotated_to_111 (bool), run_soh (bool), k_steps (int), pos_est (list), a_lattice (float), mass (float), show_plot (bool), timestep (int).



def compensate_for_compression(source, initial_hkl_pos_est, rotated_to_111, run_soh, k_steps, pos_est, a_lattice, mass, show_plot, timestep):


	import numpy as np
	import subprocess
	import os
	from scipy.optimize import curve_fit
	import matplotlib.pyplot as plt
	import time
	
	t0 = time.time()


	# Variables:
	# rotated_to_111 -> this determines the points where soh will expect to find the uncompressed peak positions.
	# k_lineout_stop -> the point whre soh will stop looking for a peak in each direction.
	# k_stop -> list of k_lineout_stop coordinates (one coordinate for each direction).
	# k_lineout_start -> the point whre soh will start looking for a peak in each direction.
	# k_start -> list of k_lineout_start coordinates (one coordinate for each direction).
	# k_lineout_direction -> list of strings for x, y, z.
	# k_lineout_file -> string containing the name of the soh input file to be written.
	# command_lineout -> the bash command for running soh in each direction.
	# current_working_directory -> contains string of the cwd.
	# popt, pcov -> stores the results of the gaussian fit routine.
	# lineout_out -> string containing the name of the soh output.
	# x_lineout_k, y_lineout_intensity -> k-values and intensities for each dimensional lineout, loaded form the soh output.
	# A -> amplitude of peak (maximum value of intensity lineout values).
	# A_index -> the index of the highest value of intensity.
	# sigma -> initial guess for width of gaussian.
	# mu -> initial guess for k position of peak centre.
	# p0 -> list of sigma, mu, A to be input to the gaussian approximartion.
	# compression_peak_pos -> contains the estimated positions of peaks once compression has been taken into consideration.
	# compression_factor -> the number by which the pos_est will be multiplied in order to find the corrected peak position estimates.
	# compressed_gsqr_est -> G^2 values of the compressed peak positions.




	if rotated_to_111 == False:

		#First we find the k limits of the 1D ft we want to get soh to perform.
		k_lineout_stop = 1.5 * ((2*np.pi)/(a_lattice * 0.5))/(2*np.pi/a_lattice) # The result is k_lineout_stop = 3, which is equivalent to the lattice being compressed to a/1.5. This means we will only detect compressions up to this point.

		k_stop = [0] * 3

		for i in range(len(k_stop)):
			k_stop[i] = [0.0] * 3
			k_stop[i][i] = k_lineout_stop




		k_lineout_start = 0.5 #Defines where we start looking. If peaks are down here, then we have expansion, not compression. 

		k_start = [0] * 3

		for i in range(len(k_start)):
			k_start[i] = [0.0] * 3
			k_start[i][i] = k_lineout_start


		
		subprocess.call("mkdir soh_compression_lineouts", shell=True)



		k_lineout_direction = ["x", "y", "z"]

		#Then we write the soh input file to run the lineouts in each direction.
		for i in range(len(k_lineout_direction)):




			k_lineout_file = "k" + k_lineout_direction[i] + "_compression_lineout.soh"





			f = open(k_lineout_file, "w")


			f.write("VERBOSE\t\t\t\t\t0\n\n"
			"FILE_TYPE\t\t\t\tlammps-multi\n"
			"DATA_FILE\t\t\t\t" + str(source) + "\n"
			"APPEND_FILE_NAME\t\tk" + str(k_lineout_direction[i]) + "_compression_lineout\n\n"
			"PLOT_OUTPUT\t\t\t\tpdf\n\n"
			"COORDS_SCALED\n"
			"SET_MASS\t\t" + str(mass) + "\n\n"
			"SET_A_CELL\t\t\t\t" + str(a_lattice) + "\n\n"
			"CALC_1D_FT\n\n"
			"SET_K_START\t\t\t\t" + str(k_start[i][0]) + " " + str(k_start[i][1]) + " " + str(k_start[i][2]) + "\n"
			"SET_K_STOP\t\t\t\t" + str(k_stop[i][0]) + " " + str(k_stop[i][1]) + " " + str(k_stop[i][2]) + "\n"
			"SET_NK\t\t\t\t\t" + str(k_steps) + "\n")


			f.close() # Remember to close the file before you try to run it!







			command_lineout = 'mpiexec -np 24 sonOfHoward ' + k_lineout_file # Stores the bash command we will run. If we want to make it faster, we can increase the processor_number here, we just have to make sure we don't get in anyone else's way!
	
			if run_soh == True:
				subprocess.call(command_lineout, shell=True)
	
	
				
	
			
			current_working_directory = os.getcwd()
	
	

			
			
			
			subprocess.call("mv " + str(current_working_directory) + "/" + str(k_lineout_file) + " " + str(current_working_directory) + "/soh_compression_lineouts/" , shell=True)
			subprocess.call("mv " + str(current_working_directory) + "/" + source + "." + str(timestep) + ".k" + k_lineout_direction[i] + "_compression_lineout.ft " + str(current_working_directory) + "/soh_compression_lineouts/" , shell=True)
			
			
		

		# Next we fit a symmetric function to the peaks we just found in each direction in k (this should be one peak per direction).
		#This method fits a Gaussian to the peak. 


		popt = [0] * len(k_lineout_direction)
		pcov = [0] * len(k_lineout_direction)

		def gauss(x, A, mu, sigma):
			return A*np.exp(-(x-mu)**2/(2.*sigma**2)) 







		for i in range(len(k_lineout_direction)):
			lineout_out = str(current_working_directory) + "/soh_compression_lineouts/" + source + "." + str(timestep) + ".k" + k_lineout_direction[i] + "_compression_lineout.ft"
			x_lineout_k, y_lineout_intensity = np.loadtxt(lineout_out, skiprows = 1, usecols = (i, 5), unpack=True)
			A = max(y_lineout_intensity)
			A_index = np.argmax(y_lineout_intensity)
			mu = x_lineout_k[A_index]	
			sigma = 0.01
			p0 = [A, mu, sigma]
			popt[i], pcov[i] = curve_fit(gauss, x_lineout_k, y_lineout_intensity, p0)
			
			
			
			#The following section will show a plot of each of the lineout peaks with the fitted Gaussian.

			plt.plot(x_lineout_k, y_lineout_intensity)
			plt.plot(x_lineout_k, gauss(x_lineout_k, popt[i][0], popt[i][1], popt[i][2]))#	plt.plot(x_lineout_k, y_lineout_intensity)
			plt.xlabel("k" + k_lineout_direction[i])
			plt.ylabel("Intensity")
			plot_name = "k" + str(k_lineout_direction[i]) + "_compression_lineout.png"				
			plt.savefig(plot_name, bbox_inches='tight')
			
			if show_plot == True:
				plt.show()
				
			plt.close()
			
			
			print "Plot of compression compensation lineout in k" + str(k_lineout_direction[i]) + " created."



			
		#Now we put the positions of the peak in each direction into a list, then normalise the actual position with respect to the predicited position (i.e divide by 2). This method will only work for fcc since the predicted position is hardcoded in (divided by two).

		compression_peak_pos = [0] * len(k_lineout_direction)


		
		for i in range(len(k_lineout_direction)):
			compression_peak_pos[i] = popt[i][1]

		
		# This should work for both fcc and bcc.
		
		
		compression_factor = [1,1,1]
		
		for i in range(len(k_lineout_direction)):
			compression_factor[i] = compression_peak_pos[i]/2.0
			
			
			
			
			
			
			
			
			
			
		
		
	if rotated_to_111 == True:

		k_start  = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

		k_start[0] = [1.41421356237*2.0*0.8, 0, 0]
		k_start[1] = [0, 2.44948974278*2.0*0.8, 0]
		k_start[2] = [0, 0, 1.73205080757*0.8]
	
	
		k_stop = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
		
		k_stop[0] = [1.41421356237*2.0*1.5, 0, 0]
		k_stop[1] = [0, 2.44948974278*2.0*1.5, 0]
		k_stop[2] = [0, 0, 1.73205080757*1.5]
		
		subprocess.call("mkdir soh_compression_lineouts", shell=True)



		k_lineout_direction = ["x", "y", "z"]

		#Then we write the soh input file to run the lineouts in each direction.
		for i in range(len(k_lineout_direction)):




			k_lineout_file = "k" + k_lineout_direction[i] + "_compression_lineout.soh"





			f = open(k_lineout_file, "w")


			f.write("VERBOSE\t\t\t\t\t0\n\n"
			"FILE_TYPE\t\t\t\tlammps-multi\n"
			"DATA_FILE\t\t\t\t" + str(source) + "\n"
			"APPEND_FILE_NAME\t\tk" + str(k_lineout_direction[i]) + "_compression_lineout\n\n"
			"PLOT_OUTPUT\t\t\t\tpdf\n\n"
			"COORDS_SCALED\n"
			"SET_MASS\t\t" + str(mass) + "\n\n"
			"SET_A_CELL\t\t\t\t" + str(a_lattice) + "\n\n"
			"CALC_1D_FT\n\n"
			"SET_K_START\t\t\t\t" + str(k_start[i][0]) + " " + str(k_start[i][1]) + " " + str(k_start[i][2]) + "\n"
			"SET_K_STOP\t\t\t\t" + str(k_stop[i][0]) + " " + str(k_stop[i][1]) + " " + str(k_stop[i][2]) + "\n"
			"SET_NK\t\t\t\t\t" + str(k_steps) + "\n")


			f.close() # Remember to close the file before you try to run it!







			command_lineout = 'mpiexec -np 24 sonOfHoward ' + k_lineout_file # Stores the bash command we will run. If we want to make it faster, we can increase the processor_number here, we just have to make sure we don't get in anyone else's way!
	
			if run_soh == True:
				subprocess.call(command_lineout, shell=True)
	
	
				
	
			
			current_working_directory = os.getcwd()
	
	

			
			
			
			subprocess.call("mv " + str(current_working_directory) + "/" + str(k_lineout_file) + " " + str(current_working_directory) + "/soh_compression_lineouts/" , shell=True)
			subprocess.call("mv " + str(current_working_directory) + "/" + source + "." + str(timestep) + ".k" + k_lineout_direction[i] + "_compression_lineout.ft " + str(current_working_directory) + "/soh_compression_lineouts/" , shell=True)
			
			
		

		# Next we fit a symmetric function to the peaks we just found in each direction in k (this should be one peak per direction).
		#This method fits a Gaussian to the peak. 


		popt = [0] * len(k_lineout_direction)
		pcov = [0] * len(k_lineout_direction)

		def gauss(x, A, mu, sigma):
			return A*np.exp(-(x-mu)**2/(2.*sigma**2)) 







		for i in range(len(k_lineout_direction)):
			lineout_out = str(current_working_directory) + "/soh_compression_lineouts/" + source + "." + str(timestep) + ".k" + k_lineout_direction[i] + "_compression_lineout.ft"
			x_lineout_k, y_lineout_intensity = np.loadtxt(lineout_out, skiprows = 1, usecols = (i, 5), unpack=True)
			A = max(y_lineout_intensity)
			A_index = np.argmax(y_lineout_intensity)
			mu = x_lineout_k[A_index]	
			sigma = 0.01
			p0 = [A, mu, sigma]
			popt[i], pcov[i] = curve_fit(gauss, x_lineout_k, y_lineout_intensity, p0)
			
			
			
			#The following section will show a plot of each of the lineout peaks with the fitted Gaussian.

			plt.plot(x_lineout_k, y_lineout_intensity)
			plt.plot(x_lineout_k, gauss(x_lineout_k, popt[i][0], popt[i][1], popt[i][2]))#	plt.plot(x_lineout_k, y_lineout_intensity)
			plt.xlabel("k" + k_lineout_direction[i])
			plt.ylabel("Intensity")
			plot_name = "k" + str(k_lineout_direction[i]) + "_compression_lineout.png"				
			plt.savefig(plot_name, bbox_inches='tight')
			
			if show_plot == True:
				plt.show()
				
			plt.close()
			
			
			print "Plot of compression compensation lineout in k" + str(k_lineout_direction[i]) + " created."





			
		#Now we put the positions of the peak in each direction into a list, then normalise the actual position with respect to the predicited position (i.e divide by 2). This method will only work for fcc since the predicted position is hardcoded in (divided by two).

		compression_peak_pos = [0] * len(k_lineout_direction)


		
		for i in range(len(k_lineout_direction)):
			compression_peak_pos[i] = popt[i][1]

		
		# This should work for both fcc and bcc.
		
		
		compression_factor = [1,1,1]
		

		compression_factor[0] = compression_peak_pos[0]/(1.41421356237*2.0)
		compression_factor[1] = compression_peak_pos[1]/(2.44948974278*2.0)
		compression_factor[2] = compression_peak_pos[2]/1.73205080757

		
	
	
	
	

	
	
	
	
	
	
	
	
	
	compressed_pos_est = [0] * len(pos_est)
	compressed_gsqr_est = [0] * len(pos_est)

	
	

	for i in range(len(compressed_pos_est)):
	
		compressed_pos_est[i] = [ compression_factor[0] * pos_est[i][0], compression_factor[1] * pos_est[i][1], compression_factor[2] * pos_est[i][2] ]

		
		compressed_gsqr_est[i] = (compressed_pos_est[i][0] ** 2) + (compressed_pos_est[i][1] ** 2) + (compressed_pos_est[i][2] ** 2)
		


		print "\nPeak " + str(i+1) + " of " + str(len(compressed_pos_est)) + " " + str(initial_hkl_pos_est[i]) + ":\n" + str(pos_est[i]) +  " compensated for compression to \n" + str(compressed_pos_est[i]) + " with G^2 = " + str(compressed_gsqr_est[i]) + "."





	#This part makes the log entry for the function.	
	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction compensate_for_compression called with input:\n"
	"source = " + str(source) + "\n"
	"initial_hkl_pos_est = " + str(initial_hkl_pos_est) + "\n"	
	"rotated_to_111 = " + str(rotated_to_111) + "\n"
	"run_soh = " + str(run_soh) + "\n"	
	"a_lattice = " + str(a_lattice) + "\n"
	"k_steps = " + str(k_steps) + "\n"
	"mass = " + str(mass) + "\n"
	"show_plot = " + str(show_plot) + "\n"
	"timestep = " + str(timestep) + "\n"
	"pos_est = (given below as the uncompensated coordinates)\n")


	for i in range(len(pos_est)):
	
		f.write("\nPeak " + str(i+1) + " of " + str(len(compressed_pos_est)) + " " + str(initial_hkl_pos_est[i]) + ":\n" + str(pos_est[i]) +  " compensated for compression to \n" + str(compressed_pos_est[i]) + " with G^2 = " + str(compressed_gsqr_est[i]) + ".")

	f.write("\n\nFunction compensate_for_compression returned:\n"
	"compressed_pos_est = (shown above)\n"
	"compressed_gsqr_est = (shown above)\n"
	"compression_factor = " + str(compression_factor) + "\n"
	)
	
	f.close()

	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.compensate_for_compression took \t\t" + str(tt) + " s to complete.")	


	return compressed_pos_est, compressed_gsqr_est, compression_factor



##################################################################
# This function creates a box around each point in reciprocal space and performs a fourier transform of the atoms in a lammps .atom file with the reciprocal lattice vectors inside the box. This function creates input files for SoH, then runs SoH for each input file. It takes as input: source (str), pos_est (list), compression_factor (list), initial_hkl_pos_est (list), a_lattice (float), del_kx(float), del_ky(float), del_kz(float), k_steps (int), run_soh (bool), mass (float). It does not produce output variables. Instead it creates files which contain the SoH outputs, including intensities at each point.


def get_peak_intensities(source, pos_est, compression_factor, initial_hkl_pos_est, a_lattice, mass, del_kx, del_ky, del_kz, k_steps, k_steps_accurate, run_soh, timestep):



	# Variables:
	# source -> the name of the lamps file to be analysed.
	# pos_est -> the estimated positions of the peaks in reciprocal space.
	# compression_factor -> the ratio of expected peak position to actual peak position in each direction.
	# a_lattice -> the lattice constant, in Angstroms.
	# mass -> the mass of the material in g/mol.
	# del_kx -> determines the size of the reciprocal space box in x over which the FT will take place.
	# del_ky -> determines the size of the reciprocal space box in y over which the FT will take place.
	# del_kz -> determines the size of the reciprocal space box in z over which the FT will take place.
	# k_steps -> determines the resolution of the FT; this number sets the number of points in each dimension of the box such that k_steps^3 is the total number of reciprocal space points per box.
	# run_soh -> turns on/off SoH.
	# timestep -> for file location reasons.
	# current_working_directory -> holds a string of the current working directory.
	# source_location -> contains location of the lammps source file.
	# kx_start/ky_start/kz_start -> holds the start points of the FT in x/y/z.
	# kx_end/ky_end/kz_end -> holds the end points of the FT in x/y/z.
	# filenum -> the number to be appended to the soh in/out filenames.
	# soh_input -> stores the filename we want to write our soh input to.
	# t_start_peak/t_end_peak -> stores the start and time of each peak run.
	# time_peak -> the difference between t_end and t_start.
	# time_remaining -> approximates the time left for the FT based on the time it took for the last calculation.


	
	import time
	import subprocess
	import os
	import numpy
	import copy


	t0 = time.time()

	print "get_peak_intensities started..."

	current_working_directory = os.getcwd()

	
	source_location = str(current_working_directory) + "/" + str(source)


	
	subprocess.call("mkdir soh_input", shell=True)
	subprocess.call("mkdir soh_output", shell=True)
	
	
	def accurate_peak_centre_and_breadth(over_width, make_plots_accurate):
		
		print "Finding accurate centre and breadths for each peak..."
		
		acc_dir = "accurate_peak_lineouts"
		
		subprocess.call("rm -r " + acc_dir, shell = True)
		subprocess.call("mkdir " + acc_dir, shell=True)
		
		print "accurate_peak_lineouts directory made"
		
		for i in range(len(pos_est)):
		
			lineout_direction = ["kx", "ky", "kz"]
			del_k = [del_kx, del_ky, del_kz]

			peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])
			

			

			subprocess.call("mkdir " + current_working_directory + "/" + peak_dir, shell=True)
			subprocess.call("mv " + current_working_directory + "/" + peak_dir + "/ " + current_working_directory + "/" + acc_dir + "/", shell = True)	

		
			for j in range(len(lineout_direction)):
		

				width = del_k[j] * (1.0 + over_width)
				
				k_start = pos_est[i][j] - width
				k_end = pos_est[i][j] + width			
			
				print "\nPeak " + peak_dir

			


		
				filenum = lineout_direction[j]
		
				soh_input = current_working_directory + "/" + acc_dir + "/" + peak_dir + "/in_" + filenum + ".soh" 	
		
			
				if j == 0 :
	
					f = open(str(soh_input), "w")

					f.write("VERBOSE\t\t\t0\n\n"
					"FILE_TYPE\t\tlammps-multi\n"
					"DATA_FILE\t\t" + source_location + "\n"
					"APPEND_FILE_NAME\t\t" + filenum + "\n\n"
					"PLOT_OUTPUT\t\tpdf\n\n"
					"COORDS_SCALED\n"
					"SET_MASS\t\t" + str(mass) + "\n\n"
					"SET_A_CELL\t\t" + str(a_lattice) + "\n\n"
					"CALC_1D_FT\n\n"
					"SET_K_START\t\t\t" + str(k_start) + " " + str(pos_est[i][1]) + " " + str(pos_est[i][2]) + "\n"
					"SET_K_STOP\t\t\t" + str(k_end) + " " + str(pos_est[i][1]) + " " + str(pos_est[i][2]) + "\n"
					"SET_NK\t\t\t" + str(k_steps_accurate) + "\n")

					f.close()
	
	
				if j == 1 :
	
					f = open(str(soh_input), "w")

					f.write("VERBOSE\t\t\t0\n\n"
					"FILE_TYPE\t\tlammps-multi\n"
					"DATA_FILE\t\t" + source_location + "\n"
					"APPEND_FILE_NAME\t\t" + filenum + "\n\n"
					"PLOT_OUTPUT\t\tpdf\n\n"
					"COORDS_SCALED\n"
					"SET_MASS\t\t" + str(mass) + "\n\n"
					"SET_A_CELL\t\t" + str(a_lattice) + "\n\n"
					"CALC_1D_FT\n\n"
					"SET_K_START\t\t\t" + str(pos_est[i][0]) + " " + str(k_start) + " " + str(pos_est[i][2]) + "\n"
					"SET_K_STOP\t\t\t" + str(pos_est[i][0]) + " " + str(k_end) + " " + str(pos_est[i][2]) + "\n"
					"SET_NK\t\t\t" + str(k_steps_accurate) + "\n")

					f.close()


				if j == 2 :
	
					f = open(str(soh_input), "w")

					f.write("VERBOSE\t\t\t0\n\n"
					"FILE_TYPE\t\tlammps-multi\n"
					"DATA_FILE\t\t" + source_location + "\n"
					"APPEND_FILE_NAME\t\t" + filenum + "\n\n"
					"PLOT_OUTPUT\t\tpdf\n\n"
					"COORDS_SCALED\n"
					"SET_MASS\t\t" + str(mass) + "\n\n"
					"SET_A_CELL\t\t" + str(a_lattice) + "\n\n"
					"CALC_1D_FT\n\n"
					"SET_K_START\t\t\t" + str(pos_est[i][0]) + " " + str(pos_est[i][1]) + " " + str(k_start) + "\n"
					"SET_K_STOP\t\t\t" + str(pos_est[i][0]) + " " + str(pos_est[i][1]) + " " + str(k_end) + "\n"
					"SET_NK\t\t\t" + str(k_steps_accurate) + "\n")

					f.close()

	
	
				if run_soh == True:
					subprocess.call('cd soh_input ; mpiexec -np 24 sonOfHoward ' + soh_input, shell=True)	
	
			
				soh_output = str(source) + "." + str(timestep) + "." + str(filenum) + ".ft"
			
				subprocess.call("mv " + soh_output + " " + current_working_directory + "/" + acc_dir + "/" + peak_dir + "/", shell=True)
	

		if make_plots_accurate == True:
		
			lineout_direction = ["kx", "ky", "kz"]
			
			for i in range(len(pos_est)):
		

				peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])	
				
				
				for j in range(len(lineout_direction)):		
				
					soh_output = str(source) + "." + str(timestep) + "." + str(lineout_direction[j]) + ".ft"
					
					plot_datafile = current_working_directory + "/" + acc_dir + "/" + peak_dir + "/" + soh_output
					plot_name = lineout_direction[j] + ".png"
					gnuplot_input = "in_gnuplot_" + lineout_direction[j]
		
					g = open(gnuplot_input, "w")
					g.write(
					"set terminal png size 1600,1200 enhanced font 'Helvetica,20'"
					"\nset output '" + str(plot_name)  + "'"
					"\nplot '" + plot_datafile + "' using " + str(j+1) + ":6")
					g.close()
					
					print "Plotted " + peak_dir + " along " + str(lineout_direction[j]) 
					
					subprocess.call("gnuplot " + str(gnuplot_input), shell=True)
					
					subprocess.call("mv " + gnuplot_input + " " + current_working_directory + "/" + acc_dir + "/" + peak_dir + "/", shell=True)
					
					subprocess.call("mv " + plot_name + " " + current_working_directory + "/" + acc_dir + "/" + peak_dir + "/", shell=True)
					
	
		accurate_pos_est = [0] * len(pos_est)
		accurate_breadths = [0] * len(pos_est)
	
		for i in range(len(pos_est)):
		
			accurate_pos_est[i] = [0] * 3
			accurate_breadths[i] = [0] * 3
			
			
		print "\nFinding accurate peak centres and breadths for:"		
		
		for i in range(len(pos_est)):
		
			peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])
			
			lineout_direction = ["kx", "ky", "kz"]
	

			for j in range(len(lineout_direction)):
			
				print "\n" + str(peak_dir) + " along " + lineout_direction[j]
			
				datafile = current_working_directory + "/" + acc_dir + "/" + peak_dir + "/" + str(source) + "." + str(timestep) + "." + str(lineout_direction[j]) + ".ft"
		
				
		
				k_temp, intensity_temp = numpy.loadtxt(datafile, skiprows=1, usecols=(j,5), unpack=True)
		
				accurate_pos_est[i][j] = max(k_temp)
				
				ind = numpy.argmax(intensity_temp)
				accurate_pos_est[i][j] = k_temp[ind]

				
				
				for k in range(len(k_temp)):
					
					intensity_diff_left = intensity_temp[ind - k] - intensity_temp[ind - k - 1]


					if ind - k - 1 < 0:
						
						print "Lower bound for peak " + peak_dir + " could not be found."
						exit() 



					if intensity_diff_left <= 0.0:
					
						k_acc_start = k_temp[ind - k]
						print "\nIntensity diff left for " + peak_dir + " " + lineout_direction[j] + " = " + str(intensity_diff_left)
						break
						
					else:
						
						continue
						
					
					
				for k in range(len(k_temp)):
				
					if k + 1 + ind >= len(k_temp):
						
						print "Upper bound for peak " + peak_dir + " could not be found."
						exit() 
					
					intensity_diff_right = intensity_temp[ind + k] - intensity_temp[ind + k + 1]

					
					if intensity_diff_right <= 0.0:
					
						k_acc_end = k_temp[ind + k]
						print "Intensity diff right for " + peak_dir + " = " + str(intensity_diff_right)

						break
						
					else:
						
						continue
					
					
				accurate_breadths[i][j] = [k_acc_start, k_acc_end]
				

				
		
		return accurate_pos_est, accurate_breadths;
				
					

	accurate_pos_est, accurate_breadths = accurate_peak_centre_and_breadth(0.5, True)
	
	print "\nCreated accurate estimates of peak centres and breadths.\n"
	
	


	for i in range(len(pos_est)):


		peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])


		t_start_peak = time.time()



		kx_start = accurate_breadths[i][0][0]
		kx_end = accurate_breadths[i][0][1]
		
		ky_start = accurate_breadths[i][1][0]
		ky_end = accurate_breadths[i][1][1]
		
		kz_start = accurate_breadths[i][2][0]
		kz_end = accurate_breadths[i][2][1]





		filenum = peak_dir
		
		soh_input = "in_" + filenum + ".soh" 





		f = open(str(current_working_directory) + "/soh_input/" + str(soh_input), "w")

		f.write("VERBOSE\t\t\t0\n\n"
		"FILE_TYPE\t\tlammps-multi\n"
		"DATA_FILE\t\t" + source_location + "\n"
		"APPEND_FILE_NAME\t\t" + filenum + "\n\n"
		"PLOT_OUTPUT\t\tpdf\n\n"
		"COORDS_SCALED\n"
		"SET_MASS\t\t" + str(mass) + "\n\n"
		"SET_A_CELL\t\t" + str(a_lattice) + "\n\n"
		"CALC_3D_FT\n\n"
		"SET_KX\t\t\t" + str(kx_start) + " " + str(kx_end) + " " + str(k_steps) + "\n"
		"SET_KY\t\t\t" + str(ky_start) + " " + str(ky_end) + " " + str(k_steps) + "\n"
		"SET_KZ\t\t\t" + str(kz_start) + " " + str(kz_end) + " " + str(k_steps) + "\n")

		f.close()



		
		if run_soh == True:
			subprocess.call('cd soh_input ; mpiexec -np 24 sonOfHoward ' + soh_input, shell=True)

		

		
		soh_output = str(source) + "." + str(timestep) + "." + str(filenum) + ".ft"
			



		subprocess.call("mv " + soh_output + " " + str(current_working_directory) + "/soh_output/", shell=True)



		t_end_peak = time.time()



		time_peak = t_end_peak - t_start_peak
		print "\nTime for peak " + str(i + 1) + " of " + str(len(pos_est)) + " = " + str(time_peak) + " s"



		time_remaining = time_peak * ((int(len(pos_est))) - i)
		print "Approximate time remaining = " + str(time_remaining) + " s\n"
		
		



	#This part makes the log entry for the function.	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction get_peak_intensities called with input:\n"
	"source = " + str(source) + "\n"
	"initial_hkl_pos_est and pos_est =\n")
	
	for i in range(len(pos_est)):
		f.write(str(initial_hkl_pos_est[i]) + " sought at " + str(pos_est[i]) + "\n")
	
	f.write("a_lattice = " + str(a_lattice) + "\n"
#	"del_kx, del_ky, del_kz = " + str(del_kx) + ", " + str(del_ky) + " ," + str(del_kz) + "\n"
	"k_steps = " + str(k_steps) + "\n"
	"run_soh = " + str(run_soh) + "\n"
	"\nFunction get_peak_intensities returned:\n"
	"This function does not return any values. It produces fourier transforms of lammps .atom files.\n"
	"The SoH inputs are stored at " + str(current_working_directory) + "/soh_inputs/\n"
	"The SoH outputs are stored at " + str(current_working_directory) + "/soh_outputs/\n")
	
	f.close()
	
	t1 = time.time()
	tt = t1 - t0
	t = open('time.pkfd', 'a')
	t.write("\nmod.get_peak_intensities took \t\t\t" + str(tt) + " s to complete.")	


	return;


##################################################################

# This function


def get_ln_intensity(pos_est, initial_hkl_pos_est, miller_pos_est, source, show_plot, timestep, a_lattice, del_kx, del_ky, del_kz, k_steps, compression_factor, make_plots):

	import numpy as np
	import os
	import subprocess
	import matplotlib.pyplot as plt
	import time
	
	t0 = time.time()
	
	print "get_ln_intensity started..."
	
	cwd = os.getcwd()
	
	if make_plots == True:
			
		subprocess.call("mkdir " + str(cwd) + "/plots_of_data/", shell = True)


	simple_intensity_integrated = [0] * len(pos_est) # Stores a simple sum of intensities of each peak.
	complex_intensity_integrated = [0] * len(pos_est) # Stores the sum of intensity*volumes for each peak.

	
	gsqr_integrated = [0] * len(pos_est)
	
	pos_integrated = [[0, 0, 0]] * len(pos_est)

	kx = [0] * len(pos_est)
	ky = [0] * len(pos_est)
	kz = [0] * len(pos_est)
	

	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction get_ln_intensity called with input:\n"
	"pos_est = " + str(pos_est) + "\n"	
	"source = " + str(source) + "\n"
	"timestep = " + str(timestep) + "\n"
	"a_lattice = " + str(a_lattice) + "\n")

	first_peak_dir = str(initial_hkl_pos_est[0][0]) + str(initial_hkl_pos_est[0][1]) + str(initial_hkl_pos_est[0][2])
		

		
	first_soh_out = str(cwd) + "/soh_output/" + source + "." + str(timestep)+ "." + first_peak_dir+ ".ft" # Stores the name of the soh output file.

	kx_coord, ky_coord, kz_coord, first_intensity = np.loadtxt(first_soh_out, skiprows=1, usecols=(0, 1, 2, 5), unpack=True)
	#
	
	points_in_bulk = 0
	points_in_surface = 0 
	points_in_edge = 0
	points_in_corner = 0
	
	classification = ["bulk", "surface", "edge", "corner"]
	classification_ind = 0
	volume_fraction = list(first_intensity)

	#h = open(intensity_datafile, 'w')
	#h.write('#kx ky kz intensity_volume intensity classification')
	
	print "Sorting k-space points into corners, edges, surfaces, and bulk..."
	
	for j in range(len(first_intensity)):
		
		# This finds all of the corner, edge, and surface intensity points.
		if kx_coord[j] == min(kx_coord) or kx_coord[j] == max(kx_coord) or ky_coord[j] == min(ky_coord) or ky_coord[j] == max(ky_coord) or kz_coord[j] == min(kz_coord) or kz_coord[j] == max(kz_coord):
			# This finds all corner and edge intensity points.
			if kx_coord[j] == min(kx_coord) and ky_coord[j] == min(ky_coord) or kx_coord[j] == min(kx_coord) and kz_coord[j] == min(kz_coord) or ky_coord[j] == min(ky_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == max(kx_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == min(kx_coord) and kz_coord[j] == max(kz_coord) or kx_coord[j] == max(kx_coord) and kz_coord[j] == max(kz_coord) or ky_coord[j] == max(ky_coord) and kz_coord[j] == min(kz_coord) or ky_coord[j] == max(ky_coord) and kx_coord[j] == max(kx_coord) or kx_coord[j] == max(kx_coord) and ky_coord[j] == min(ky_coord) or ky_coord[j] == min(ky_coord) and kz_coord[j] == max(kz_coord) or ky_coord[j] == max(ky_coord) and kz_coord[j] == max(kz_coord) or kx_coord[j] == min(kx_coord) and ky_coord[j] == max(ky_coord):
						
				# This finds all the corner intensity points.
				if kx_coord[j] == min(kx_coord) and ky_coord[j] == min(ky_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == min(kx_coord) and ky_coord[j] == min(ky_coord) and kz_coord[j] == max(kz_coord) or kx_coord[j] == max(kx_coord) and ky_coord[j] == min(ky_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == min(kx_coord) and ky_coord[j] == max(ky_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == max(kx_coord) and ky_coord[j] == max(ky_coord) and kz_coord[j] == min(kz_coord) or kx_coord[j] == max(kx_coord) and ky_coord[j] == min(ky_coord) and kz_coord[j] == max(kz_coord) or kx_coord[j] == min(kx_coord) and ky_coord[j] == max(ky_coord) and kz_coord[j] == max(kz_coord) or kx_coord[j] == max(kx_coord) and ky_coord[j] == max(ky_coord) and kz_coord[j] == max(kz_coord):
				
					
					
					points_in_corner += 1
					classification_ind = 3
					volume_fraction[j] = 0.125

				
				# All the edge points must go here.
				else:
					
					points_in_edge += 1
					classification_ind = 2
					volume_fraction[j] = 0.25

			# All the surface points must go here.	
			else:

				points_in_surface += 1
				classification_ind = 1
				volume_fraction[j] = 0.5
				
		# All the bulk points must go here.
		else:

			points_in_bulk += 1
			classification_ind = 0
			volume_fraction[j] = 1.0

	print "Finished sorting k-space points."

	total_points = points_in_bulk + points_in_surface + points_in_edge + points_in_corner
	
	expected_bulk = (k_steps - 2) ** 3
	expected_surface = 6 * ((k_steps - 2) ** 2)
	expected_edge = 12 * (k_steps - 2)
	expected_corner = 8
	expected_total = expected_corner + expected_edge + expected_surface + expected_bulk

	print "\nPoints in bulk = " + str(points_in_bulk) + "			Expected " + str(expected_bulk)
	print "Points on a surface = " + str(points_in_surface) + "		Expected " + str(expected_surface)
	print "Points on an edge = " + str(points_in_edge) + "			Expected " + str(expected_edge)
	print "Point on a corner = " + str(points_in_corner) + "			Expected " + str(expected_corner)
	print "\nTotal points = "  + str(total_points) + "			Expected " + str(expected_total) + "\n"
	
	dk_vol_var = 1.0/(k_steps - 1.0) # Division is computationally expensive so best to do this outside the loop, then multiply it in.

	print "Integrating intensities..."

	for i in range(len(pos_est)):

		peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])
		
		print peak_dir
		
		soh_out = str(cwd) + "/soh_output/" + source + "." + str(timestep)+ "." + peak_dir+ ".ft" # Stores the name of the soh output file.

		kx_coord, ky_coord, kz_coord, tmp_intensity = np.loadtxt(soh_out, skiprows = 1, usecols = (0,1,2,5), unpack=True)   

		dk_vol = ( (max(kx_coord) - min(kx_coord)) * dk_vol_var) * ( (max(ky_coord) - min(ky_coord)) * dk_vol_var ) * ( (max(kz_coord) - min(kz_coord)) * dk_vol_var )

		peak_position_ind = np.argmax(tmp_intensity)

		tmp_gsqr_integrated = (kx_coord[peak_position_ind] * kx_coord[peak_position_ind]) + (ky_coord[peak_position_ind] * ky_coord[peak_position_ind]) + (kz_coord[peak_position_ind] * kz_coord[peak_position_ind])

		gsqr_integrated[i] = tmp_gsqr_integrated * (2 * np.pi / a_lattice ) * (2 * np.pi / a_lattice) # This is because of how soh handles the data. It is also the reason I was initially getting more peaks than Will.

		pos_integrated[i] = [ kx_coord[peak_position_ind], ky_coord[peak_position_ind],  kz_coord[peak_position_ind] ]

		subprocess.call("mkdir " + str(cwd) + "/plots_of_data/" + peak_dir, shell=True)
		
		intensity_datafile = str(cwd) + "/plots_of_data/" + peak_dir + "/intensity_vs_position.dat"		

		simple_intensity_integrated[i] = sum(tmp_intensity)
		
		intensity_volume = list(tmp_intensity)
	
		for j in range(len(tmp_intensity)):

			intensity_volume[j] = tmp_intensity[j] * dk_vol * volume_fraction[j]
		
		
		complex_intensity_integrated[i] = sum(intensity_volume)
	
		f.write("\nIntegrated intensity of " + str(miller_pos_est[i]) + " sought at " + str(pos_est[i]) + " = " + str(complex_intensity_integrated[i]) + "\n"
		"Simple sum of intensities = " + str(simple_intensity_integrated[i]) )
		
		
		if make_plots == True:

		
			plot_directory_name = str(miller_pos_est[i][0]) + str(miller_pos_est[i][1]) + str(miller_pos_est[i][2])
		
				
			subprocess.call("mkdir " + str(cwd) + "/plots_of_data/" + str(plot_directory_name), shell = True)
		
		
					
			kx_for_lineout_plot = []
				
			intensity_for_kx_lineout_plot = []
		
			
				
		
			for j in range(len(kx_coord)):
		
				
				if ky_coord[j] == ky_coord[peak_position_ind] and kz_coord[j] == kz_coord[peak_position_ind]:
				
				
					kx_for_lineout_plot.append(kx_coord[j])
		
				
					intensity_for_kx_lineout_plot.append(tmp_intensity[j])






			ky_for_lineout_plot = []
				
			intensity_for_ky_lineout_plot = []
		

			for j in range(len(kx_coord)):
		
				
				if kx_coord[j] == kx_coord[peak_position_ind] and kz_coord[j] == kz_coord[peak_position_ind]:
				
				
					ky_for_lineout_plot.append(ky_coord[j])
		
				
					intensity_for_ky_lineout_plot.append(tmp_intensity[j])





			kz_for_lineout_plot = []
				
			intensity_for_kz_lineout_plot = []
		

		
			for j in range(len(kx_coord)):
		
				
				if kx_coord[j] == kx_coord[peak_position_ind] and ky_coord[j] == ky_coord[peak_position_ind]:
				
				
					kz_for_lineout_plot.append(kz_coord[j])
		
				
					intensity_for_kz_lineout_plot.append(tmp_intensity[j])



			k_value = [kx_for_lineout_plot, ky_for_lineout_plot, kz_for_lineout_plot]
			
			intensity_value = [intensity_for_kx_lineout_plot, intensity_for_ky_lineout_plot, intensity_for_kz_lineout_plot]
			
			lineout_direction = ["kx", "ky", "kz"]
			
			for j in range(len(lineout_direction)):
			
				d = open(str(cwd) + "/plots_of_data/" + str(plot_directory_name) + "/I_vs_" + lineout_direction[j] + ".dat", "w")
				d.write("#" + lineout_direction[j] + "    intensity\n")
				
				for k in range(len(k_value[j])):
				
					d.write(str(k_value[j][k]) + " " + str(intensity_value[j][k]) + "\n")
					
				d.close()
	
	
	
	
	if make_plots == True:
	
		lineout_direction = ["kx", "ky", "kz"]
	
		for i in range(len(pos_est)):
		
			plot_directory_name = str(miller_pos_est[i][0]) + str(miller_pos_est[i][1]) + str(miller_pos_est[i][2])
			
			for j in range(len(lineout_direction)):
			
				datafile = cwd + "/plots_of_data/" + plot_directory_name + "/I_vs_" + lineout_direction[j] + ".dat"
			
				plot_name = cwd + "/plots_of_data/" + plot_directory_name + "/" + lineout_direction[j] + "_lineout.png"
			
				g = open(cwd + "/plots_of_data/" + plot_directory_name + "/" + lineout_direction[j] + "_gnuplot.in", "w")
				g.write(
				"set terminal png size 1600,1200 enhanced font 'Helvetica,20'"
				"\nset output '" + str(plot_name)  + "'"
				"\nplot '" + datafile + "' using 1:2"
				)
				g.close()
				
				
	if make_plots == True:
	
		lineout_direction = ["kx", "ky", "kz"]
	
		for i in range(len(pos_est)):
		
			plot_directory_name = str(miller_pos_est[i][0]) + str(miller_pos_est[i][1]) + str(miller_pos_est[i][2])
			
			for j in range(len(lineout_direction)):
				
				gnuplot_input = cwd + "/plots_of_data/" + plot_directory_name + "/" + lineout_direction[j] + "_gnuplot.in"
			
				subprocess.call("gnuplot " + gnuplot_input, shell=True)		

	# This section works on the complex_integrated_intensity.

	complex_intensity_integrated_max_ind = np.argmax(complex_intensity_integrated)


	ln_complex_intensity_integrated = np.log(complex_intensity_integrated)
	
	
	ln_norm_complex_intensity_integrated = np.log(complex_intensity_integrated/max(complex_intensity_integrated))
	
	
	ln_norm_complex_intensity_integrated.tolist()
	
	
	g = open("ln_complex_intensity_vs_g_squared.dat", "w")
	
	g.write("ln_complex_intensity g_squared h k l\n")
	
	for i in range(len(pos_est)):
		
		g.write(str(ln_complex_intensity_integrated[i]) + " " + str(gsqr_integrated[i]) + " " + str(miller_pos_est[i][0]) + " " + str(miller_pos_est[i][1]) + " " + str(miller_pos_est[i][2]) + "\n")
	
	g.close()
	
	
	# This section works on the simple_integrated_intensity.

	simple_intensity_integrated_max_ind = np.argmax(simple_intensity_integrated)


	ln_simple_intensity_integrated = np.log(simple_intensity_integrated)
	
	
	ln_norm_simple_intensity_integrated = np.log(simple_intensity_integrated/max(simple_intensity_integrated))
	
	
	ln_norm_simple_intensity_integrated.tolist()
	
	
	g = open("ln_simple_intensity_vs_g_squared.dat", "w")
	
	g.write("ln_simple_intensity g_squared h k l\n")
	
	for i in range(len(pos_est)):
		
		g.write(str(ln_simple_intensity_integrated[i]) + " " + str(gsqr_integrated[i]) + " " + str(miller_pos_est[i][0]) + " " + str(miller_pos_est[i][1]) + " " + str(miller_pos_est[i][2]) + "\n")
	
	g.close()
	
	#This part makes the final log entry for the function.	
	
	f.write("\nFunction get_ln_intensity returned:\n"
	"This function obtains the integrated intensity and then returns ln of the intensity."
	"\nIt also returns the gsqr values of the estimated peak centres.")
	
	f.close()
	
	t1 = time.time()
	tt = t1 - t0
	t = open('time.pkfd', 'a')
	t.write("\nmod.get_ln_intensity took \t\t\t" + str(tt) + " s to complete.")

	return pos_integrated, gsqr_integrated, ln_complex_intensity_integrated, ln_norm_complex_intensity_integrated, ln_simple_intensity_integrated, ln_norm_simple_intensity_integrated
	
	
################################################################

def get_slope_ln_intensity_vs_gsqr(gsqr, ln_intensity):

	import numpy as np
	import time
	
	t0 = time.time()

	print "\nget_slope_ln_intensity_vs_gsqr started..."

	slope_ln_intensity_vs_gsqr, constant_ln_intensity_vs_gsqr  = np.polyfit(gsqr, ln_intensity, 1)



	print "\n\nThe slope of ln(I) vs. G^2 = " + str(slope_ln_intensity_vs_gsqr)
	print "The line constant of ln(I) vs. G^2 = " + str(constant_ln_intensity_vs_gsqr)


	# The log entry for the function.

	f= open("log.pkfd", "a")
	
	f.write("\n\nFunction get_slope_ln_intensity called with input:\n"
	"ln_intensity, gsqr = the values written above for each peak.\n"

	"\nFunction get_slope_ln_intensity_vs_gsqr returned:\n"
	"slope_ln_intensity_vs_gsqr = " + str(slope_ln_intensity_vs_gsqr) + "\n"
	"constant_ln_intensity_vs_gsqr = " + str(constant_ln_intensity_vs_gsqr) + "\n"
	)
	
	f.close()

	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.get_slope_ln_intensity took \t\t" + str(tt) + " s to complete.")

	return slope_ln_intensity_vs_gsqr, constant_ln_intensity_vs_gsqr
	
	
################################################################



def calc_temperature_xrd(slope_ln_intensity_vs_gsqr, constant_ln_intensity_vs_gsqr, gruneisen_uncompressed, debye_temperature_uncompressed, a_lattice, compression_factor, mass, pos, gsqr, uncompressed_pos_est, uncompressed_gsqr_est, plot_name, show_plot, ln_intensity, md_temperature_3d, md_temperature_2d):



	import numpy as np
	import matplotlib.pyplot as plt
	import scipy.constants as codata
	from scipy.integrate import quad
	import time
	
	t0 = time.time()
		
	print "\ncalc_temperature_xrd started..."

	compressed_volume = (a_lattice ** 3)/( compression_factor[0] * compression_factor[1] * compression_factor[2])
	


	
	gruneisen_over_volume = gruneisen_uncompressed/((a_lattice ** 3) * (10 ** -30)) #This is the model that says gruneisen/V = constant.
	#Checked with google calculator.






	# This function is used in the Pandya and Ramakrishnan models.
	def gruneisen_over_volume_func(integrable_v, initial_gruneisen, q_power):

		return (1.0/(integrable_v)) * initial_gruneisen * ( ( (integrable_v)/((a_lattice * 1e-10) ** 3) ) ** q_power )



	# This function is used with Walsh's model. 
	def walsh_gruneisen_over_volume_func(integrable_v, initial_gruneisen):
		return (1.0/(integrable_v)) * ( initial_gruneisen +  ( -3.296 * ( (((a_lattice * 1e-10)** 3)/integrable_v) - 1.0) ) + ( 10.493 * ( ( (((a_lattice * 1e-10)** 3)/integrable_v) - 1.0) ** 2 ) ) + ( -19.264 * ( ( (((a_lattice * 1e-10)** 3)/integrable_v) - 1.0) ) ** 3)  ) 
	


	# The function is integrated between the uncompressed volume and the compressed volume. 
	integrated_gruneisen_1, err_integrated_gruneisen_1 = quad(gruneisen_over_volume_func, ((a_lattice * 1e-10)** 3), compressed_volume * 1e-30, args=(1.93, 1.085)) #Pandya
	#Checked with Wolfram Alpha
	
	integrated_gruneisen_2, err_integrated_gruneisen_2 = quad(gruneisen_over_volume_func, ((a_lattice * 1e-10)** 3), compressed_volume * 1e-30, args=(2.008, 1.33)) #Ramakrishnan
	#Checked with Wolfram Alpha



	#Again the function is integrated between the uncompressed volume and the compressed volume.
	integrated_gruneisen_3, err_integrated_gruneisen_3 = quad(walsh_gruneisen_over_volume_func, ((a_lattice * 1e-10)** 3), compressed_volume * 1e-30, args=(2.04)) #Walsh
	#Checked with integral-calculator.com





	# Here the Debye temperature at the compressed volume is calculated for each model.


	#This list will contain the Debye temperature at the given compression.
	estimated_debye_temperature = [0] * 4

	estimated_debye_temperature[0] = debye_temperature_uncompressed * np.exp(-gruneisen_over_volume * (10 ** -30) * (compressed_volume - (a_lattice ** 3))) #Uses the equation from Murphy et. al. 2008 which. Note that this is dependent on the material_debye_temperature, i.e the Debye temperature of the uncompressed crystal.
	# Checked with google calculator.

	estimated_debye_temperature[1] = debye_temperature_uncompressed * np.exp(-integrated_gruneisen_1)

	estimated_debye_temperature[2] = debye_temperature_uncompressed * np.exp(-integrated_gruneisen_2)

	estimated_debye_temperature[3] = debye_temperature_uncompressed * np.exp(-integrated_gruneisen_3)




	# Now we use the estimated Debye temperatures to calculate the temperature of the sample predicted by each model, using Debye-Waller.

	temperature_est = [0] * 4


	temperature_normalisation_factor =  ( mass * (10 ** -3) * codata.value("Boltzmann constant") * 4 * np.pi * np.pi) / ((10 ** 20) * 3 * codata.value("Planck constant") * codata.value("Planck constant") * codata.value("Avogadro constant"))  # Note that the factor of 10^20 is because our G^2 is in 1/Angstroms^2, so we convert it to meters here. This formulation is from Will Murphy's PhD thesis.
	#Checked with Google calculator


	for i in range(len(estimated_debye_temperature)):
		temperature_est[i] = - (estimated_debye_temperature[i] ** 2) * slope_ln_intensity_vs_gsqr * temperature_normalisation_factor
	
	



	
	
		
	



	# The following is all about plotting the ln(I) vs G^2 with a line fit, and an ideal line fit (as in, the line required to obtain the correct temperature).
	
	def line(x, m, c):
		return m*x + c


	line_point_x1 = 0
	line_point_x2 = max(gsqr)


	line_point_y1 = line(line_point_x1, slope_ln_intensity_vs_gsqr, constant_ln_intensity_vs_gsqr)
	line_point_y2 = line(line_point_x2, slope_ln_intensity_vs_gsqr, constant_ln_intensity_vs_gsqr)


	line_points_x = [line_point_x1, line_point_x2]
	line_points_y = [line_point_y1, line_point_y2]
	
	
	
	# This part calculates the slope required to perfectly calculate the temperature. It then creates the points necessary to plot this slope on the plot.
	
	
	
	ideal_slope_constant_model = - md_temperature_2d/( (estimated_debye_temperature[0] ** 2) * temperature_normalisation_factor)
	
	
	ideal_line_point_y1 = line(line_point_x1, ideal_slope_constant_model, constant_ln_intensity_vs_gsqr)
	ideal_line_point_y2 = line(line_point_x2, ideal_slope_constant_model, constant_ln_intensity_vs_gsqr)


	ideal_line_points_x = [line_point_x1, line_point_x2]
	ideal_line_points_y = [line_point_y1, line_point_y2]
	









	# This part calculates an approximation of the upper and lower bounds of the measured temperature from the peaks.
	
	
	
	gsqr_types = set(uncompressed_gsqr_est)
	


	
	minimum_ln_intensity = []
	
	minimum_gsqr = []
	
	minimum_actual_gsqr = [] # These values of actual_gsqr are in units of A^-2. The above are not.
	
	maximum_ln_intensity = []
	
	maximum_gsqr = []
	
	maximum_actual_gsqr = []
	
	
	
	

	
	for e in gsqr_types:
		
		
		
		intensities_for_each_gsqr_est = []
		
		
		actual_gsqrs_for_each_gsqr_est = []



				
		
		for i in range(len(uncompressed_gsqr_est)):
		
			if uncompressed_gsqr_est[i] == e:
			
				intensities_for_each_gsqr_est.append(ln_intensity[i])
				
				
				
				
				actual_gsqrs_for_each_gsqr_est.append(gsqr[i])
	
	
	
		
				
		minimum_ln_intensity.append(min(intensities_for_each_gsqr_est))
								
		maximum_ln_intensity.append(max(intensities_for_each_gsqr_est))
		

				
		index_of_minimum_intensity = intensities_for_each_gsqr_est.index(min(intensities_for_each_gsqr_est))
		
		index_of_maximum_intensity = intensities_for_each_gsqr_est.index(max(intensities_for_each_gsqr_est))
		

		
		minimum_actual_gsqr.append(actual_gsqrs_for_each_gsqr_est[index_of_minimum_intensity])
		
		maximum_actual_gsqr.append(actual_gsqrs_for_each_gsqr_est[index_of_maximum_intensity])
		

				
				
		
	

		
	slope_min_boundary, constant_min_boundary = np.polyfit(minimum_actual_gsqr, minimum_ln_intensity, 1)
	
	slope_max_boundary, constant_max_boundary = np.polyfit(maximum_actual_gsqr, maximum_ln_intensity, 1)
	
		
		
		


		
	# At this point we swap notation of the maximum to minimum and vice versa. This is because the upper boundary of the slope corresponds to the lower boundary of the temperature estimate and vice versa.
		
	maximum_temperature_est = [0] * len(estimated_debye_temperature)
		
	for i in range(len(estimated_debye_temperature)):
		maximum_temperature_est[i] = - (estimated_debye_temperature[i] ** 2) * slope_min_boundary * temperature_normalisation_factor
	
	
	
	minimum_temperature_est = [0] * len(estimated_debye_temperature)
		
	for i in range(len(estimated_debye_temperature)):
		minimum_temperature_est[i] = - (estimated_debye_temperature[i] ** 2) * slope_max_boundary * temperature_normalisation_factor
	








	# This part finds an average of the temperatures and standard deviation.



	central_temperature_mean = np.mean(temperature_est)
	
	
	central_temperature_stdev = np.std(temperature_est)
	
	
	minimum_temperature_mean = np.mean(minimum_temperature_est)
	
	
	maximum_temperature_mean = np.mean(maximum_temperature_est)
	
	temperature_error = np.mean([abs(central_temperature_mean - minimum_temperature_mean), abs(central_temperature_mean - maximum_temperature_mean)])



	# This part creates the plot.
	



	plt.plot(gsqr, ln_intensity, 'ko')	
	plt.plot(line_points_x, line_points_y, 'k')
	plt.plot(ideal_line_points_x, ideal_line_points_y, 'r')
	
	for i in range(len(gsqr)):
		label = "(" + str(uncompressed_pos_est[i][0]) + str(uncompressed_pos_est[i][1]) + str(uncompressed_pos_est[i][2]) + ")"
		plt.annotate(label, xy = (gsqr[i], ln_intensity[i]) )
	
	plt.xlabel('|$G^{2}$| / A$^{-2}$')
	plt.ylabel('$Ln(I/I_{0})$ / arb.')
	plt.title('$Ln(I/I_{0})$ vs. $G^{2}$')
	plt.xticks()
	plt.yticks()
	plt.savefig(plot_name, bbox_inches='tight')
	
	if show_plot == True:
		plt.show()
		
	plt.close()

	compression_ratio_x = 1/compression_factor[0]
	compression_ratio_y = 1/compression_factor[1]
	compression_ratio_z = 1/compression_factor[2]
	compression_ratio = compressed_volume/(a_lattice ** 3)


	# This part prints out our temperature calculations.


	print ("\nVolume compression ratio = " + str(compression_ratio) + "\n"
	"Compression ratio in x = " + str(compression_ratio_x) + "\n"
	"Compression ratio in y = " + str(compression_ratio_y) + "\n"
	"Compression ratio in z = " + str(compression_ratio_z) + "\n"	
	"\nThe best fit slope gives the following measurements of temperarature for each Gruneisen model:\n"
	"Gruneisen/vol constant -> T = " + str(temperature_est[0]) + "\n"
	"Pandya -> T = " + str(temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(temperature_est[2]) + "\n"
	"Walsh -> T = " + str(temperature_est[3]) + "\n"
	"Ideal slope for Gruneisen/vol constant model (in order to predict correct temperature) = " + str(ideal_slope_constant_model) + "\n"
	"\nA lower temperature boundary is given by the slope of the upper-most peaks :\n"
	"Gruneisen/vol constant -> T = " + str(minimum_temperature_est[0]) + "\n"
	"Pandya -> T = " + str(minimum_temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(minimum_temperature_est[2]) + "\n"
	"Walsh -> T = " + str(minimum_temperature_est[3]) + "\n"
	"\nAn upper temperature boundary is given by the slope of the lower-most peaks:\n"
	"Gruneisen/vol constant -> T = " + str(maximum_temperature_est[0]) + "\n"
	"Pandya -> T = " + str(maximum_temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(maximum_temperature_est[2]) + "\n"
	"Walsh -> T = " + str(maximum_temperature_est[3]) + "\n"
	"\nThe 3D MD temperature = " + str(md_temperature_3d) + "\n"
	"The 2D MD temperature = " + str(md_temperature_2d) + "\n"
	"From the best fit slope, the mean temperature is T = " + str(central_temperature_mean) + " +/- " + str(central_temperature_stdev) + " K."
	"\nIf we consider the boundary temperatures as well (something that might look like an experiment), we can get T = " + str(central_temperature_mean) + " +/- " + str(temperature_error) + " K.")






	#This part makes the log entry for the function.	
	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction calc_temperature_xrd called with input:\n"
	"slope_ln_intensity_vs_gsqr = " + str(slope_ln_intensity_vs_gsqr) + "\n"
	"constant_ln_intensity_vs_gsqr = " + str(constant_ln_intensity_vs_gsqr) + "\n"	
	"gruneisen_uncompressed = " + str(gruneisen_uncompressed) + "\n"
	"debye_temperature_uncompressed = " + str(debye_temperature_uncompressed) + "\n"
	"a_lattice = " + str(a_lattice) + "\n"
	"compressed_volume = " + str(compressed_volume) + "\n"
	"mass = " + str(mass) + "\n"
	"pos = " + str(pos) + "\n"
	"gsqr = " + str(gsqr) + "\n"
	"plot_name = " + str(plot_name) + "\n"
	"show_plot = " + str(show_plot) + "\n"
	"ln_intensity = " + str(ln_intensity) + "\n"
	"md_temperature_3d = " + str(md_temperature_3d) + "\n"
	"md_temperature_2d = " + str(md_temperature_2d) + "\n"
	
	"\nFunction calc_temperature_xrd returned:\n"
	"Compression ratio = " + str(compression_ratio) + "\n"
	"Ideal slope (required to exactly calculate T) = " + str(ideal_slope_constant_model) + "\n"
	"\nThe best fit slope gives the following measurements of temperarature for each Gruneisen model:\n"
	"Gruneisen/vol constant -> T = " + str(temperature_est[0]) + "\n"
	"Pandya -> T = " + str(temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(temperature_est[2]) + "\n"
	"Walsh -> T = " + str(temperature_est[3]) + "\n"
	"Ideal slope for Gruneisen/vol constant model (in order to predict correct temperature) = " + str(ideal_slope_constant_model) + "\n"
	"\nA lower temperature boundary is given by the slope of the upper-most peaks :\n"
	"Gruneisen/vol constant -> T = " + str(minimum_temperature_est[0]) + "\n"
	"Pandya -> T = " + str(minimum_temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(minimum_temperature_est[2]) + "\n"
	"Walsh -> T = " + str(minimum_temperature_est[3]) + "\n"
	"\nAn upper temperature boundary is given by the slope of the lower-most peaks:\n"
	"Gruneisen/vol constant -> T = " + str(maximum_temperature_est[0]) + "\n"
	"Pandya -> T = " + str(maximum_temperature_est[1]) + "\n"
	"Ramakrishnan -> T = " + str(maximum_temperature_est[2]) + "\n"
	"Walsh -> T = " + str(maximum_temperature_est[3]) + "\n"
	"\nFrom the best fit slope, the mean temperature is T = " + str(central_temperature_mean) + " +/- " + str(central_temperature_stdev) + " K."
	"\nIf we consider the boundary temperatures as well (something that might look like an experiment), we can get T = " + str(central_temperature_mean) + " +/- " + str(temperature_error) + " K.")
	
	f.close()

	t1 = time.time()
	tt = t1 - t0
	time_elapsed = time.localtime()
	t = open('time.pkfd', 'a')
	t.write("\nmod.calc_temperature_xrd took \t\t\t" + str(tt) + " s to complete.")

	return temperature_est, central_temperature_mean



#########################################################################



def calc_debye_temperature(slope_ln_intensity_vs_gsqr, mass, md_temperature):


	import scipy.constants as codata
	import numpy as np
	import time
	
	t0 = time.clock()


	debye_normalisation_factor = (10 ** 20) * md_temperature * 3 * codata.value("Planck constant") * codata.value("Planck constant") * codata.value("Avogadro constant") / (mass  * (10 ** -3) * codata.value("Boltzmann constant") * 4 * np.pi * np.pi)
	
	
	debye_temperature = np.sqrt(abs(debye_normalisation_factor/slope_ln_intensity_vs_gsqr))


	print "\n\nCalculated Debye temperature = " + str(debye_temperature) + " K\n"
	


	f = open("log.pkfd", "a")
	
	f.write("\n\nFunction calc_debye_temperature called with input:\n"
	"slope_ln_intensity_vs_gsqr = " + str(slope_ln_intensity_vs_gsqr) + "\n"
	"mass = " + str( mass) + "\n"
	"md_temperature = " + str(md_temperature) + "\n"
	"\nFunction calc_debye_temp returned:\n"
	"debye_temperature = " + str(debye_temperature) + "\n"	
	)
	
	f.close()

	t1 = time.clock()
	tt = t1 - t0
	t = open('time.pkfd', 'a')
	t.write("\nmod.calc_debye_temperature took \t\t" + str(tt) + " s to complete.")


	return debye_temperature
	
	
	
	
###########################################################################



def profile_peaks(source, timestep, initial_hkl_pos_est, make_plots):
		
	import time	
	import numpy as np
	import os
	import subprocess
	
	t0 = time.time()
	
	# This function profiles each peak by plotting each intensity point vs
	# distance from the central position.
	
	cwd = os.getcwd()
	
	intensity = [0] * len(initial_hkl_pos_est)
	k_diff_abs = [0] * len(initial_hkl_pos_est)
	
	for i in range(len(initial_hkl_pos_est)):

		peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])

		datafile = str(cwd) + "/soh_output/" + source + "." + str(timestep) + "." + peak_dir + ".ft"
	
		kx, ky, kz, intensity[i] = np.loadtxt(datafile, skiprows=1, usecols=(0,1,2,5), unpack=True)
	
		k_diff_abs[i] = [0] * len(intensity[i])
	
		peak_position = np.argmax(intensity[i])
		
		gsqr_centre = (kx[peak_position] ** 2) + (ky[peak_position] ** 2) + (kz[peak_position] ** 2) 
	
		for j in range(len(intensity[i])):
		
			gsqr = (kx[j] ** 2) + (ky[j] ** 2) + (kz[j] ** 2) 
			
			k_diff_abs[i][j] = np.sqrt(abs(gsqr - gsqr_centre))
			
		
	if make_plots == True:
	
		rm_command = "rm -r " + cwd + "/peak_histograms"
		subprocess.call(rm_command, shell=True)
		mkdir_command = "mkdir " + cwd + "/peak_histograms"
		subprocess.call(mkdir_command, shell=True)
		
		
		for i in range(len(initial_hkl_pos_est)):
		
			peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])
		
			location = str(cwd) + "/peak_histograms/" + peak_dir
		
			mkdir_command_2 = "mkdir " + location
		
			subprocess.call(mkdir_command_2, shell=True)
		
			dat_filename = location + "/histogram.dat"
		
			h = open(dat_filename, "w")
			h.write("#k_differential intensity")
			
			for j in range(len(intensity[i])):
			
				h.write("\n" + str(k_diff_abs[i][j]) + " " + str(intensity[i][j]) + "")
				
			h.close()
			
			in_filename = location + "/histogram_gnuplot.in"
			plot_name = "histogram_" + peak_dir + ".png"
			
			g = open(in_filename, "w")
			g.write(
			"set terminal png size 1600,1200 enhanced font 'Helvetica,20'"
			"\nset output '" + str(plot_name)  + "'"
			"\nplot '" + dat_filename + "' using 1:2"
			)
			
			g.close()
		
		
	if make_plots == True:
	
		for i in range(len(initial_hkl_pos_est)):
		
			peak_dir = str(initial_hkl_pos_est[i][0]) + str(initial_hkl_pos_est[i][1]) + str(initial_hkl_pos_est[i][2])
		
			location = str(cwd) + "/peak_histograms/" + peak_dir
			
			in_filename = location + "/histogram_gnuplot.in"
			
			plot_name = "histogram_" + peak_dir + ".png"
		
			subprocess.call("gnuplot < " + str(in_filename), shell=True)
			
			mv_command = "mv " + plot_name + " " + location
			
			subprocess.call(mv_command, shell=True)
		
					
	tf = time.time()
	tt = tf - t0			
	t = open('time.pkfd', 'a')
	t.write("\nmod.profile_peaks took \t\t" + str(tt) + " s to complete.")
	
	return;
	
				
############################################################################
# This function gives the final message at the end of a run.

def checkout(xrd_temperatures, xrd_temperature_labels, md_temperatures, md_temperature_labels):
	
	f = open("log.pkfd", "a")

	print "\n\n#########################\n\npeakfinder.py finished\n\n#########################"
	f.write("\n\n#########################\n\npeakfinder.py finished\n\n#########################")
	
	print "\n\nThe MD temperatures are:\n"
	f.write("\n\nThe MD temperatures are:\n")
	
	for i in range(len(md_temperatures)):
	
		print md_temperature_labels[i] + "\t\t= " + str(md_temperatures[i]) + " K"
		f.write(md_temperature_labels[i] + "\t\t= " + str(md_temperatures[i]) + " K\n")

	print "\nThe estimated x-ray diffraction temperatures are:\n"
	f.write("\nThe estimated x-ray diffraction temperatures are:\n")
	
	for i in range(len(xrd_temperatures)):
	
		percent_off = 100.0 * (abs(xrd_temperatures[i] - md_temperatures[0]))/md_temperatures[0]
		
		print xrd_temperature_labels[i] + "\t= " + str(xrd_temperatures[i]) + " K \t\t" + str(percent_off) + " % from the 2D MD temperature."
		f.write(xrd_temperature_labels[i] + "\t= " + str(xrd_temperatures[i]) + " K\n" + str(percent_off) + " % from the 2D MD temperature.")
		
		
		
	return;
		
		
