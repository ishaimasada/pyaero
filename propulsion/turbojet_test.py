'''
Example usage code for the Engine class
'''
import json
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("parameters.json", "r") as file:
    parameters = json.load(file)["parameters"]

# Import all types from engine module
from engine import *

engine_parameters = parameters["engine"]
engine = Engine(engine_parameters)
engine.set_components(parameters)
station_data = engine.get_station_data()
print(station_data)
engine.plot_thermo()
engine.display_performance()