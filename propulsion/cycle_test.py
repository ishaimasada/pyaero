'''
Example driver code for testing the cycle analysis functionality of the engine module
'''

import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add to search locations
sys.path.append(r'..\aerodynamics')

# Import all types from engine module
from engine import *

# Initialize engine object
engine_parameters = {"PR": 5, "TET": 1510, "LHV": 43.15*10**6, "Minf": 0.85, "altitude": 12192, "front face diameter": 47.7, "bleed": Bleed(cooling=0, packing=0)}
engine = Engine(engine_parameters)

# Define custom architecture
intake_parameters = {"engine": engine, "total pressure recovery": 0.98, "M_inlet": 0.6, "M_exit": 0.55}
compressor_parameters = {"engine": engine, "e": 0.91, "M_exit": 0.3, "idx_inlet": 2, "idx_exit": 3}
burner_parameters = {"engine": engine, "total pressure loss": 0.05, "efficiency": 0.99, "M_exit": 0.4, "idx_exit": 4}
turbine_parameters = {"e": 0.95, "mechanical efficiency": 0.99, "M_exit": 0.5, "idx_exit": 5}
afterburner_parameters = {"engine": engine, "toggle": 1, "AET": 1800, "pi off": 0.98, "pi on": 0.94, "efficiency": 0.97}
nozzle_parameters = {"engine": engine, "C/CD": "C", "discharge coefficient": 0.97, "velocity coefficient": 0.98}
components = [
                inlet := Inlet(intake_parameters),
                compressor := Compressor(inlet.exit, compressor_parameters),
                burner := Burner(compressor.exit, burner_parameters),
                turbine :=  Turbine(burner.exit, compressor, turbine_parameters),
                afterburner := Afterburner(turbine.exit, afterburner_parameters),
                exhaust :=  Nozzle(afterburner.exit, nozzle_parameters)
]

# setting the engine architecture
engine.architecture = components

print(engine.get_station_data())
print(engine.get_performance())