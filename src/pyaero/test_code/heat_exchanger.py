### Example Heat Exchanger (Recuperator Class) ###
from pyaero.propulsion.engine import Recuperator
import json
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("heat_exchanger_parameters.json", "r") as file:
    parameters = json.load(file)

heat_exchanger = Recuperator(cold_inlet=None, component_parameters=parameters)
print(heat_exchanger.get_results())
