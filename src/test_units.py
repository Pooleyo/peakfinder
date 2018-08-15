def run(unit_name, test_input, expected_result):

    import numpy as np

    exec("from units import " + unit_name)
       
    string_test_input = str(test_input)
    string_test_input = string_test_input[1:-1]

    exec("actual_result = " + str(unit_name) + "(" + string_test_input + ")")

    if actual_result == expected_result:
    
        print "Pass\t\t"  + "UNIT\t\t" + unit_name
        
    elif actual_result != expected_result:

        print "\n####FAIL####\t\t" + "UNIT\t\t" + unit_name
        print "Expected result =\t%.12f" % expected_result
        print "Actual result \t=\t%.12f\n" % actual_result
        
    else:
    
        print "TEST FRAMEWORK BROKEN"

    return

