#####################################################
#
# Header last edited: 22/08/17
# This file contains the nuts and bolts of the new
# code. Each large process undertaken by peakfinder 
# is described by a "brick" function. These 
# bricks functions are made up of smaller functions called a 
# "unit"; units are accompanied by unit 
# tests, found in the testmodule.py file. 
#
# Unit functions are in the upper section of this
# file. Brick functions are in the lower section.
# 
#####################################################
# Units.

def LoadData(filename, columnnumbers):

    from numpy import loadtxt
    
    columnnumbers = tuple(columnnumbers)
 
    data = loadtxt(filename, usecols = columnnumbers, skiprows = 1, unpack = True)
    return data
   
   
   
def FindPeakCentre(data, intensityindex, kxkykzindex):
    
    from numpy import argmax
    
    intensity = data[intensityindex]
    peakintensityindex = argmax(intensity)
    
    kx = data[kxkykzindex[0]]
    ky = data[kxkykzindex[1]]
    kz = data[kxkykzindex[2]]
    
    peakcentre = [kx[peakintensityindex], ky[peakintensityindex], kz[peakintensityindex]]
    
    return peakcentre;



def FindOrthogonalLineout(data, intensityindex, kxkykzindex, directionindex, point):

    # directionindex = 0 for kx, = 1 for ky, and = 2 for kz.

    kindex = kxkykzindex[directionindex]
    kdirectiondata = data[kindex]
    otherdirectionsindex = list(kxkykzindex)
    otherdirectionsindex.remove(kindex)   

    orthogonallineout = []
    for i in range(len(kdirectiondata)):

        if data[otherdirectionsindex[0]][i] == point[otherdirectionsindex[0]] and data[otherdirectionsindex[1]][i] == point[otherdirectionsindex[1]]:
            lineoutpoint = [0,0,0,0]
            lineoutpoint[otherdirectionsindex[0]] = data[otherdirectionsindex[0]][i]
            lineoutpoint[otherdirectionsindex[1]] = data[otherdirectionsindex[1]][i]
            lineoutpoint[kindex] = kdirectiondata[i]
            lineoutpoint[intensityindex] = data[intensityindex][i]
            orthogonallineout.append(lineoutpoint)
     
    return orthogonallineout
    
    
    
def FindIntensityMinima1D(data, kxkykzindex, intensityindex):

    import numpy as np

    intensity = [0] * len(data)
    for i in range(len(data)):
        intensity[i] = data[i][intensityindex]
        
    maxintensityindex = np.argmax(intensity)
    maxpoint = data[maxintensityindex]
    minima = [0,0]

    imax = 1 + len(data)/2
    tempintensity = [0] * imax
    for i in range(imax):
        tempintensity[i] = intensity[i]

    minimumindex = np.argmin(tempintensity)
    minima[0] = data[minimumindex]
           
    tempintensity = [0] * imax
    for i in range(imax):
        tempintensity[i] = intensity[imax - 1 + i]

    minimumindex = np.argmin(tempintensity)
    minima[1] = data[imax - 1 + minimumindex]
       
    return minima
    
    
    
def BuildIntensityVolume(points, kxkykzindex, intensityindex):
    
    ############ Work in progress. #############

    for i in range(len(kxkykzindex)):
       coordinate1 = points[kxkykzindex[i]][0]
       coordinate2 = points[kxkykzindex[i]][1]
       intensity1 = points
       
    return
    
    

def GetCompressedGruneisenParameterModel1(uncompressedgruneisenparameter, volumetriccompressionratio):

    # Documentation edited: 22/08/17
    # Models the Gruneisen parameter under the
    # asumption that:
    #
    # Gruneisen_parameter/Volume = constant
    #
    # Since Gruneisen0/V0 = Gruneisencompressed / 
    # (V0 * volumetric_compression_ratio), the
    # calculation is simplified to:
    # Gruneisencompressed = Gruneisen0 * 
    # volumetric_compression_ratio.
    #
    # This function is only configured for cubic 
    # lattices. 
    #
    # Inputs:
    # 
    # uncompressedgruneisenparameter - The Gruneisen
    #   parameter of the material at ambient 
    #   conditions.
    # volumetriccompressionratio = V/V0. This value
    #   should be less than 1.
    # 
    #
    # Outputs:
    # compressedgruneisenparameter - the Gruneisen
    #   parameter at this compression, according to
    #   this model.

    compressedgruneisenparameter = uncompressedgruneisenparameter * volumetriccompressionratio
    
    return compressedgruneisenparameter
    
    
    
def GetDebyeTemperatureFromGruneisenParameter():

    ############ Work in progress. #############
    

    return
    
    
    
#######################################
# Bricks.

def RemoveTDS():
    
    #LoadData()
    #FindPeakCentre()
    #FindOrthogonalVector()
    #FindOrthogonalVector()
    #FindOrthogonalVector()
    #FindMinimumIntensity1D()
    #FindMinimumIntensity1D()
    #FindMinimumIntensity1D()
    #Build3DIntensityVolume()
    #SubtractIntensityVolumeFromIntensity()

    return;
