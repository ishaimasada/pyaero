'''
Example usage code for component classes
'''
import json
import os
import sys

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("turbine_parameters.json", "r") as file:
    parameters = json.load(file)

sys.path.append(r"../propulsion")

# Import all types from engine module
from engine import *

turbine = Turbine(upstream=None, compressor=None, component_parameters=parameters)