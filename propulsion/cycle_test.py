'''
Example driver code for testing the cycle analysis functionality of the engine module
'''
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Import all types from engine module
from engine import *

# Initialize engine object
engine_parameters = {"PR": 3.6, "mdotf": 0.0145, "LHV": 43.15*10**6, "Minf": 0, "altitude": 0, "mdot": 0.53, "bleed": Bleed(cooling=0, packing=0)}
engine = Engine(engine_parameters, burner="mdotf", size="mdot")

# Define custom architecture
intake_parameters = {"engine": engine, "total pressure recovery": 0.98, "M_inlet": 0.1, "M_exit": 0.3}
compressor_parameters = {"engine": engine, "e": 0.91, "M_exit": 0.4, "idx_inlet": 2, "idx_exit": 3}
recuperator_parameters = {"engine": engine, "pi_hot": 0.94, "pi_cold": 0.94, "delta_ht": 400, "Mexit_cold": 0.2, "Mexit_hot": 0.3}
burner_parameters = {"engine": engine, "total pressure loss": 0.05, "efficiency": 0.99, "M_exit": 0.4, "idx_exit": 4}
turbine_parameters = {"e": 0.95, "mechanical efficiency": 0.99, "M_exit": 0.5, "idx_exit": 5}
nozzle_parameters = {"engine": engine, "C/CD": "C", "discharge coefficient": 0.97, "velocity coefficient": 0.98}

parameters = {
    "inlet": intake_parameters,
    "compressor": compressor_parameters,
    "burner": burner_parameters,
    "turbine": turbine_parameters,
    "nozzle": nozzle_parameters
}

engine.set_components(parameters)
engine.toggle_recuperator(recuperator_parameters, toggle=True)
station_data = engine.get_station_data()
print(station_data)
engine.plot_thermo()
engine.display_performance()