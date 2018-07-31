import inpkfd
import importlib
import initialise
import finalise

path = importlib.import_module(inpkfd.path) # Imports the path specified in inpkfd.

initialise.run()

path.run()

finalise.run()
