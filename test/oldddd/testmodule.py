########################################################
# These functions are multi-use functions used by the 
# tests.

def checktestresult(testresult, functionname):

    if testresult == 0:
            
        print "***FAIL***\tFunction '" + functionname + "' test complete."

    elif testresult == 1:

        print "PASS \t\tFunction '" + functionname + "' test complete."

    else:
        print "ERROR HANDLED:"
        print "Function '" + functionname + "' test returned non-binary value."
        print "Aborting test."
        exit() 
        
    return;

def testdatavariables():

    filename = "test.dat"
    datacolumnnumbers = [0, 1, 2, 3]
    intensityindex = 3
    kxkykzindex = [0,1,2]
    latticeparameter = 3.002
    #volumetriccompressionratio =
    
    return filename, datacolumnnumbers, intensityindex, kxkykzindex

##########################################################
# These test the smallest units of code.

def TestLoadData():

    from newmodule import LoadData

    def TEST():
    
        functionname = LoadData.__name__

        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()

        data = LoadData(filename, datacolumnnumbers)

        sumcolumn = [sum(data[0]), sum(data[1]), sum(data[2]), sum(data[3])]
        sumlastrow = data[0,-1] + data[1,-1] + data[2,-1] + data[3,-1]

        expectedsumcolumn = [54.0, 54.0, 54.0, 47.0]
        expectedsumlastrow = 10.0


        if sumcolumn == expectedsumcolumn and sumlastrow == expectedsumlastrow:   
            testresult = 1
        
        else:
            testresult = 0
            
        return testresult, functionname;

    testresult, functionname = TEST()

    checktestresult(testresult, functionname)

    return;



def TestFindPeakCentre():

    from newmodule import FindPeakCentre, LoadData
    
    def TEST():
  
        functionname = FindPeakCentre.__name__  

        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()       
        data = LoadData(filename, datacolumnnumbers)

        peakcentre = FindPeakCentre(data, intensityindex, kxkykzindex)
          
        expectedpeakcentre = [2.0,2.0,2.0]
        
        if peakcentre == expectedpeakcentre:
            testresult = 1
            
        else:
            testresult = 0
        
        return testresult, functionname;
        
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
        
    return;   
  
    

def TestFindOrthogonalLineout():

    from newmodule import FindOrthogonalLineout, LoadData
    
    def TEST():
  
        functionname = FindOrthogonalLineout.__name__  

        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()
        data = LoadData(filename, datacolumnnumbers)
        point = [2,2,2]

        directionindex = 0        
        orthogonallineout1 = FindOrthogonalLineout(data, intensityindex, kxkykzindex, directionindex, point)
        directionindex = 1
        orthogonallineout2 = FindOrthogonalLineout(data, intensityindex, kxkykzindex, directionindex, point)
        directionindex = 2
        orthogonallineout3 = FindOrthogonalLineout(data, intensityindex, kxkykzindex, directionindex, point)
        
        expectedorthogonallineout1 = [[1, 2, 2, 2],[2, 2, 2, 3],[3, 2, 2, 2]]
        expectedorthogonallineout2 = [[2, 1, 2, 2],[2, 2, 2, 3],[2, 3, 2, 2]]
        expectedorthogonallineout3 = [[2, 2, 1, 2],[2, 2, 2, 3],[2, 2, 3, 2]]
        
        if orthogonallineout1 == expectedorthogonallineout1 and orthogonallineout2 == expectedorthogonallineout2 and orthogonallineout3 == expectedorthogonallineout3:
            testresult = 1
        
        else:
            testresult = 0
        
        return testresult, functionname;
        
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
    
    return;     



def TestFindIntensityMinima1D():

    from newmodule import FindIntensityMinima1D

    def TEST():
  
        functionname = FindIntensityMinima1D.__name__
        
        lineout = [[2, 1, 2, 2],[2, 1.5, 2, 2.5],[2, 2, 2, 3],[2, 2.5, 2, 2.5],[2, 3, 2, 2]]
        kxkykzindex = [0,1,2]
        intensityindex = 3
        minima = FindIntensityMinima1D(lineout, kxkykzindex, intensityindex)
        expectedminima = [[2, 1, 2, 2],[2, 3, 2, 2]]
        if minima == expectedminima:
            testresult = 1
            
        else:    
            testresult = 0
            
        return testresult, functionname;
        
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
    
    return
    
    
    
def TestBuildIntensityVolume():
    
    ############ Work in progress. #############
    
    from newmodule import BuildIntensityVolume
    
    def TEST():
  
        functionname = BuildIntensityVolume.__name__
        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()
        
        intensityvolumecoefficients = 0
        
        points = [ [[1, 2, 2, 2],[3, 2, 2, 2]], [[2, 1, 2, 2],[2, 3, 2, 2]], [[2, 2, 1, 2],[2, 2, 3, 2]] ]
        
        BuildIntensityVolume(points, kxkykzindex, intensityindex)
        
        
        
        expectedintensityvolumecoefficients = [[0,2],[0,2],[0,2]]
        
        if intensityvolumecoefficients == expectedintensityvolumecoefficients:
            testresult = 1
            
        else:    
            testresult = 0
            
        return testresult, functionname;
    
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
        
    return
    


def TestGetCompressedGruneisenParameterModel1():

    from newmodule import GetCompressedGruneisenParameterModel1
    
    def TEST():
  
        functionname = GetCompressedGruneisenParameterModel1.__name__
        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()
        
        uncompresedgruneisenparameter = 1.75
        volumetriccompressionratio = 0.89
        
        compressedgruneisenparameter = GetCompressedGruneisenParameterModel1(uncompresedgruneisenparameter, volumetriccompressionratio)
 
        expectedcompressedgruneisenparameter = 1.5575
        
        if compressedgruneisenparameter == expectedcompressedgruneisenparameter:
            testresult = 1
            
        else:    
            testresult = 0
            
        return testresult, functionname;
    
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
        
    return    
    
    
    
def TestGetDebyeTemperatureFromGruneisenParameter():

    ############ Work in progress. #############

    from newmodule import GetDebyeTemperatureFromGruneisenParameter
    
    def TEST():
  
        functionname = GetDebyeTemperatureFromGruneisenParameter.__name__
        filename, datacolumnnumbers, intensityindex, kxkykzindex = testdatavariables()
        
        debyetemperature = GetDebyeTemperatureFromGruneisenParameter()
 
        expecteddebyetemperature = 300
        
        if debyetemperature == expecteddebyetemperature:
            testresult = 1
            
        else:    
            testresult = 0
            
        return testresult, functionname;
    
    testresult, functionname = TEST()
    
    checktestresult(testresult, functionname)
        
    return    
