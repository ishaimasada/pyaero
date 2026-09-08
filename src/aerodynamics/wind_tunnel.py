import numpy
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

from compressible import *
from atmosphere import *

"""
Wind Tunnel Calculations

ARCHITECTURE: Reservoir --> Converging Duct --> Throat --> Diverging Duct --> Constant Area (Test Section) --> Shock --> Atmosphere
"""

altitude = 0
gamma = 1.4
ambient = Ambient(altitude)
Mtest = 3
[_, _, Pt1_P1, _, Ae_Astar]= isentropic(Mtest, gamma, lookup_key="M")
[_, P2_P1, _, _, Ptexit_Pt1, _, Me] = normal_shock(Mtest, gamma, lookup_key="M")
[_, _, Ptexit_Pe, _, _]= isentropic(Me, gamma, lookup_key="M")
'''
Rt = 0.025 # throat radius (m)
At = numpy.pi * Rt**2
Ae = Ae_Astar * At
'''
Pt_reservoir = (1/Ptexit_Pt1) * Ptexit_Pe * ambient.P / ambient.P

print("Required Reservoir Pressure Ratio (kPa): ", Pt_reservoir * ambient.P / 10**3)