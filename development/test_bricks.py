def run(brick_name, test_input, expected_result):

    exec("import " + brick_name)

    string_test_input = str(test_input)
    string_test_input = string_test_input[1:-1]

    exec("actual_result = " + str(brick_name) + ".run(" + string_test_input + ")")

    if actual_result == expected_result:
    
        print "Pass\t\t"  + "BRICK\t\t" + brick_name
        
    else:
        
        print "####FAIL####\t\t" + "BRICK\t\t" + brick_name
        
    return
