import os, math

filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
parent_directory = "\\".join(directory.split("\\")[:-1])

from pyaero.aerodynamics.atmosphere import atmosphere 
os.chdir(directory)

# INPUTS
gamma = 1.4
R = 287
g = 9.8
[Tref, Pref, rhoref] = atmosphere(0)
AR = 10
e = 0.8
K1 = 1/(math.pi*AR*e)+0.02
K2 = -0.008
CDR = 0
acceleration = 0
wing_loading = 45 * 47.88
