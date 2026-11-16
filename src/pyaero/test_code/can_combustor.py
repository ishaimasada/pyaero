### Example annular combustor ###
from pyaero.propulsion.engine import Burner
import json
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("can_burner_parameters.json", "r") as file:
    parameters = json.load(file)

burner = Burner(upstream=None, component_parameters=parameters)
results = burner.get_results()