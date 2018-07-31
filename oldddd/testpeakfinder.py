import testmodule

# Unit tests.

testmodule.TestLoadData()
testmodule.TestFindPeakCentre()         # Must come after TestLoadData.
testmodule.TestFindOrthogonalLineout()  # Must come after TestLoadData.
testmodule.TestFindIntensityMinima1D()
testmodule.TestBuildIntensityVolume()
testmodule.TestGetCompressedGruneisenParameterModel1()
testmodule.TestGetDebyeTemperatureFromGruneisenParameter()
