''' Example usage code for the Engine class '''
import json, os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("turbojet_parameters.json", "r") as file:
    parameters = json.load(file)["parameters"]

# Import all types from engine module
from pyaero.propulsion.engine import *

# Create an instance of the Engine class
engine_parameters = parameters["engine"]
engine = Engine(engine_parameters)

# Pass the component parameters to the engine object for cycle analysis
engine.set_components(parameters)

# Retrieve the flow properties and full engine performance
station_data, performance = engine.get_performance()

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Output the station data to an excel
engine.write_station_data("cycle.xlsx")

# Display performance & plot the temperatures and pressures
print(performance)
#engine.plot_thermo()