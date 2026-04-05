import math
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

from compressible import isentropic

# Only metric units
def atmosphere(height, unit = "metric"):
    if unit == "imperial": height = height * 0.3048

    # assuming height is in m
    troposphere = 11000 # m
    stratosphere = 25000 # m
    if height < troposphere:
        temperature = 15.04 - 0.00649 * height
        pressure = 101.29 * ((temperature + 273.1) / 288.08)**5.256
    elif troposphere < height < stratosphere:
        temperature = -56.46
        pressure = 22.65 * math.exp(1.73 - 0.000157 * height)
    elif height > stratosphere:
        temperature = -131.21 + 0.00299 * height
        pressure = 2.488 * ((temperature + 273.1) / 216.6)**-11.388
    
    pressure = pressure * 10**3 # Pa
    temperature = temperature + 273.15 # K
    R = 287
    rho = pressure / (R * temperature)

    return [temperature, pressure, rho]


class Ambient:
    R = 287
    gamma = 1.4

    def __init__(self, altitude, alpha=None, Minf=None, Vinf=None):
        self.T, self.P, self.rho = atmosphere(altitude)
        self.rhoinf = self.P / (self.R * self.T)
        self.a = math.sqrt(self.gamma * self.R * self.T)

        # Allow either velocity or Mach number as input
        if Minf != None:
            self.Minf = Minf
            self.Vinf = self.Minf * self.a
        elif Vinf != None:
            self.Vinf = Vinf
            self.Minf = self.Vinf / self.a
        else: self.Minf = self.Vinf = 0
        
        [_, Tt_T, Pt_P, rhot_rho, self.A_Astar] = isentropic(self.Minf, self.gamma, lookup_key="M")
        self.Tt = Tt_T * self.T
        self.Pt = Pt_P * self.P
        self.rhot = rhot_rho * self.rho
