import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add to search locations
sys.path.append(r'..\aerodynamics')

from engine import *


engine_parameters = {
                     "OPR": 20, "TET": 1600, "BPR": 4, "FPR": 2.1, "LHV": 43.15*10**6, "Minf": 0.85, "altitude": 13106.4, "front face diameter": 43,
                     "fan tip total pressure recovery": 0.98, "fan hub total pressure recovery": 0.97, "OGV total pressure loss": 0.03,
                     "bypass duct total pressure loss": 0.01, "fan tip polytropic efficiency": 0.9, "fan hub polytropic efficiency": 0.9,
                     "swan neck loss": 0.02, "core comp. polytropic efficiency": 0.91, "burner efficiency": 0.99, "burner total pressure loss": 0.05,
                     "turbine polytropic efficiency": 0.8, "Nozzle Discharge Coefficienct": 0.97, "Nozzle Velocity Coefficienct": 0.99,
                     "shaft mechanical efficiency": 0.99, "M1": 0.6, "M1.2": 0.55, "M2": 0.55, "M3.1": 0.3, "M5": 0.5,
                     "high pressure cooling": 0, "high pressure packing": 0, "low pressure cooling": 0, "low pressure packing": 0
                     }

engine = Engine(engine_parameters)
print(engine.get_station_data())
