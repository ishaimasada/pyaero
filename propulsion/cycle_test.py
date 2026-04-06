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
engine_parameters = {"PR": 5, "TET": 1510, "LHV": 43.15*10**6, "Minf": 0.85, "altitude": 12192, "front face diameter": 47.7, "bleed": Bleed(cooling=0, packing=0)}
engine = Engine(engine_parameters)

# Define custom architecture
intake_parameters = {"engine": engine, "total pressure recovery": 0.98, "M_inlet": 0.6, "M_exit": 0.55}
compressor_parameters = {"engine": engine, "e": 0.91, "M_exit": 0.4, "idx_inlet": 2, "idx_exit": 3}
recuperator_parameters = {"engine": engine, "pi_hot": 0.94, "pi_cold": 0.94, "delta_ht": 200, "Mexit_cold": 0.3, "Mexit_hot": 0.4}
burner_parameters = {"engine": engine, "total pressure loss": 0.05, "efficiency": 0.99, "M_exit": 0.4, "idx_exit": 4}
turbine_parameters = {"e": 0.95, "mechanical efficiency": 0.99, "M_exit": 0.5, "idx_exit": 5}
nozzle_parameters = {"engine": engine, "C/CD": "C", "discharge coefficient": 0.97, "velocity coefficient": 0.98}
components = [
                inlet := Inlet(intake_parameters),
                compressor := Compressor(inlet.exit, compressor_parameters),
                recuperator := Recuperator(compressor.exit, recuperator_parameters),
                burner := Burner(recuperator.cold_exit, burner_parameters),
                turbine :=  Turbine(burner.exit, compressor, turbine_parameters),
                hot_stream := recuperator.pass_hot_stream(turbine.exit),
                nozzle :=  Nozzle(hot_stream, nozzle_parameters),
]

# setting the engine architecture
engine.architecture = components

print(engine.get_station_data())
print(engine.get_performance())