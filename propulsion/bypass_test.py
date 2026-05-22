'''
Example usage code for the BypassEngine class
'''
import json
import os
import sys

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("bypass_parameters.json", "r") as file:
    parameters = json.load(file)["parameters"]

# Import all types from engine module
from engine import *

engine_parameters = parameters["engine"]
engine = BypassEngine(engine_parameters)
engine.set_components(parameters)
station_data = engine.get_station_data()
print(station_data)
engine.plot_thermo()
engine.display_performance()