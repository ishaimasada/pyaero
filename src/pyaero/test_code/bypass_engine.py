''' Example usage code for the BypassEngine class '''
import json, os

# Import all types from engine module
from pyaero.propulsion.engine import * 

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("bypass_parameters.json", "r") as file:
    parameters = json.load(file)["parameters"]

engine_parameters = parameters["engine"]
engine = BypassEngine(engine_parameters) # type: ignore
engine.set_components(parameters)
station_data, performance = engine.get_performance()

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

station_data.to_excel("cycle.xlsx", index=False)
print(station_data)
print(performance)
engine.plot_thermo()