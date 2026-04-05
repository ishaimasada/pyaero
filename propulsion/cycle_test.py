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

# Initialize Engine Object with Design Parameters
engine_parameters = {"PR": 5, "TET": 1510, "LHV": 43.15*10**6, "Minf": 0.85, "altitude": 12192, "front face diameter": 47.7, "bleed": Bleed(cooling=0, packing=0)}
engine = Engine(engine_parameters)

# Prepare cycle and component design parameters
intake_cycle_parameters = {"engine": engine, "total pressure recovery": 0.98, "M_inlet": 0.6, "M_exit": 0.55}
compressor_cycle_parameters = {"engine": engine, "e": 0.91, "M_exit": 0.3, "idx_inlet": 2, "idx_exit": 3}
burner_cycle_parameters = {"engine": engine, "total pressure loss": 0.05, "efficiency": 0.99, "M_exit": 0.4, "idx_exit": 4}
turbine_cycle_parameters = {"e": 0.8, "mechanical efficiency": 0.99, "M_exit": 0.5, "idx_exit": 5}
nozzle_cycle_parameters = {"engine": engine, "C/CD": "C", "discharge coefficient": 0.97, "velocity coefficient": 0.98}

# Build Engine Architecture
engine.components.append(inlet := Inlet(intake_cycle_parameters))
engine.components.append(compressor := Compressor(inlet.exit, compressor_cycle_parameters))
engine.components.append(burner := Burner(compressor.exit, burner_cycle_parameters))
engine.components.append(turbine :=  Turbine(burner.exit, compressor, turbine_cycle_parameters))
engine.components.append(exhaust :=  Nozzle(turbine.exit, nozzle_cycle_parameters))
print(engine.get_station_data())
print(engine.get_performance())