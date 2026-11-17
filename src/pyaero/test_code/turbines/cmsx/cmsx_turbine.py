### Free Power Turbine ###
from pyaero.propulsion.engine import Turbine
import json, os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("cmsx_turbine_parameters.json", "r") as file:
    parameters = json.load(file)

# Define the log path for the turbine
parameters["log path"] = os.path.join(directory) + r"\turbine_log.txt"

# Turbine Design
turbine = Turbine(upstream=None, component_parameters=parameters)
turbine.get_results()
print(turbine.stages[0].inner_radii*1e3)
print(turbine.stages[0].widths*1e3)