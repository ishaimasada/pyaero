### Example usage code for the Engine class ###
import json, os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("turbojet_parameters.json", "r") as file:
    parameters = json.load(file)

# Import all types from engine module
from pyaero.propulsion.engine import Engine

# Create an instance of the Engine class
engine = Engine(parameters)
off_design_results = engine.solve_offdesign_point({"altitude": 0, "Mach":0}, {"TET": 1300})
print(off_design_results)
