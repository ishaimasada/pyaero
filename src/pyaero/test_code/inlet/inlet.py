''' Example usage code for the Inlet component class '''
import json
import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("inlet_parameters.json", "r") as file:
    parameters = json.load(file)

sys.path.append(r"..\propulsion")

# Import all types from engine module
from engine import Inlet

inlet = Inlet(component_parameters=parameters)
print(inlet.radii)
print(inlet.z_coords)
print(inlet.N)
inlet.plot_contour()