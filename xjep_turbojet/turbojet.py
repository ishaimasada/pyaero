''' 
Turbojet engine design for XJEP 

Replaces the role of the JetCats and KingTechs by allowing for easy instrumentation
'''

import json
import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("turbojet_parameters.json", "r") as file:
    parameters = json.load(file)["parameters"]

sys.path.append(r"..\propulsion")

# Import all types from engine module
from engine import *

# Create an instance of the Engine class
engine_parameters = parameters["engine"]
engine = Engine(engine_parameters)

# Pass the component parameters to the engine object for cycle analysis
engine.set_components(parameters)

# Retrieve the flow properties and full engine performance
station_data, performance = engine.get_performance()

# Change the current working directory to the file location
os.chdir(directory)

# Output the station data to an excel
station_data.to_excel("station_data.xlsx", index=False)

# Display performance & plot the temperatures and pressures
print(performance)