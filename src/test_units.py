def run(unit_name, test_input, expected_result):

    import importlib
    exec("from units import " + unit_name)
       
    string_test_input = str(test_input)
    string_test_input = string_test_input[1:-1]

    exec("actual_result = " + str(unit_name) + "(" + string_test_input + ")")
 
    if actual_result == expected_result:
    
        print "Pass\t\t"  + "UNIT\t\t" + unit_name
        
    else:
        
        print "####FAIL####\t\t" + "UNIT\t\t" + unit_name
        
    return

