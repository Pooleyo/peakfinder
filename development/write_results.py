def run(debye_temperature, temperature):

    filename = "results.pkfd"

    f = open(filename, "w")
    f.write(
        "Debye temperature\t\t\t\t" + str(debye_temperature) + "\n"
        "Temperature\t\t\t\t\t\t" + str(temperature)
    )
    f.close()

    return