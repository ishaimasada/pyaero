"""
Module for Cycle Analysis and preliminary component design 
"""

#NOTE: Axial Turbomachinery assumes constant radial work distribution (dht/dr = 0) and Free Vortex Design (dVax/dr = 0)

import matplotlib.pyplot as plt
import pandas
import numpy
import json
import copy
import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add to search locations
sys.path.append(r'..\aerodynamics')

from compressible import *
from atmosphere import *
from curves import *


class Station:
    """ Correlations for combustion products from Walsh and Fletcher """
    A = {"A0" : 992.313, "A1" : 236.688, "A2" : -1852.148, "A3" : 6083.152, "A4" : -8893.933, "A5" : 7097.112, 
         "A6" : -3234.725, "A7" : 794.571, "A8" : -81.873, "A9" : 422.178, "A10" : 1.053}

    B = {"B0" : -718.874, "B1": 8747.481, "B2": -15863.157, "B3": 17254.096, "B4": -10233.795, "B5": 3081.778, 
         "B6": -361.112, "B7": -3.919, "B8": 55.593, "B9": 1.6079}

    C = {"C0": 1.0001, "C1": 0.9248, "C2": -2.2078}
    
    REFH0 = 422.2202178
    Tref = 288.15
    Pref = 101325

    column_names = [
                    "Station",
                    "W [kg/sec]", 
                    "Tt [K]",
                    "Pt [kPa]",
                    "ht [kJ/kg]",
                    "Wc [kg/sec]", 
                    "Wf [kg/sec]", 
                    "FAR",
                    "M",
                    "Ts [K]",
                    "hs [kJ/kg]",
                    "Ps [kPa]",
                    "V [m/sec]",
                    "rho [kg/m^3]",
                    "Area [m^2]"
                   ]

    def __init__(self, W, Tt, Pt, FAR=None, M=None, idx=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.Wc = self.get_Wc(self.W, self.Tt, self.Pt)
        self.M = M
        self.idx = idx

        if FAR != None: 
            self.FAR = FAR
            self.Wa = self.W / (1 + FAR)
            self.Wf = FAR * self.Wa
        else:
            self.FAR = 0
            self.Wf = 0
        if self.M != None: 
            self.set_statics(self.M)
        self.ht = self.get_ht(self.Tt, self.FAR)

    def __str__(self):
        return f"Station {self.idx}\nW = {self.W}\nWc = {self.Wc}\nCorrected W = {self.Wc}\nTt = {self.Tt}\nPt = {self.Pt/1000}\nM = {self.M}"

    @property
    def Wc(self): return self.get_Wc(self.W, self.Tt, self.Pt)

    @Wc.setter
    def Wc(self, value): self._Wc = value

    def get_Wc(self, W, T, P): return W * numpy.sqrt(T / self.Tref) / (P / self.Pref)

    @property
    def R(self): return self.get_R(self.FAR)

    @R.setter
    def R(self, value): self._R = value

    def get_R(self, FAR): return 287.05 - 0.0099 * FAR + 0.0000001 * FAR**2 # Formula 3.22 for Kerosene

    @property
    def gamma(self): return self.get_gamma()

    @gamma.setter
    def gamma(self, value): self._gamma = value

    def get_gamma(self): return self.cp / (self.cp - self.R)

    @property
    def FAR(self): return self.get_FAR()

    @FAR.setter
    def FAR(self, value): self._FAR = value

    def get_FAR(self): return self.Wf / (self.W - self.Wf)

    @property
    def cp(self): return self.get_cp(self.Tt, self.FAR)

    @cp.setter
    def cp(self, value): self._cp = value

    def get_cp(self, T, FAR):
        TZ = T/1000
        cp = (
                self.A["A0"] + self.A["A1"] * TZ + (self.A["A2"] * TZ**2) + (self.A["A3"] * TZ**3) +
                (self.A["A4"] * TZ**4) + (self.A["A5"] * TZ**5) + (self.A["A6"] * TZ**6) +
                (self.A["A7"] * TZ**7) + (self.A["A8"] * TZ**8) + (self.B["B0"] + self.B["B1"] * TZ +
                (self.B["B2"] * TZ**2) + (self.B["B3"] * TZ**3) + (self.B["B4"] * TZ**4) +
                (self.B["B5"] * TZ**5) + (self.B["B6"] * TZ**6) + (self.B["B7"] * TZ**7)) * (FAR / (1 + FAR))
                ) # Formula 3.24
        return cp

    @property
    def ht(self): return self.get_ht(self.Tt, self.FAR)

    @ht.setter
    def ht(self, value): self._ht = value

    def get_ht(self, T, FAR):
        """ Formula 3.27 """
        TZ = T/1000
        ht = (
              self.A["A0"] * TZ + self.A["A1"] * TZ**2 / 2 + (self.A["A2"] * TZ**3) / 3 + 
              (self.A["A3"] * TZ**4) / 4 + (self.A["A4"] * TZ**5) / 5 + (self.A["A5"] * TZ**6) / 6 +
              (self.A["A6"] * TZ**7) / 7 + (self.A["A7"] * TZ**8) / 8 + (self.A["A8"] * TZ**9) / 9 +
              self.A["A9"] + (self.B["B0"] * TZ + self.B["B1"] * TZ**2 / 2 + (self.B["B2"] * TZ**3) / 3 +
              (self.B["B3"] * TZ**4) / 4 + (self.B["B4"] * TZ**5) / 5 + (self.B["B5"] * TZ**6) / 6 +
              (self.B["B6"] * TZ**7) / 7 + (self.B["B7"] * TZ**8) / 8 + self.B["B8"]) * (FAR / (1 + FAR))
            ) - self.REFH0
        return ht*1000

        
    def set_statics(self, M):
        [_, Tt_T, Pt_P, _, _] = isentropic(M, self.gamma, lookup_key="M")
        self.M = M
        self.T = (1 / Tt_T) * self.Tt
        self.h = self.T * self.cp
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        self.area = self.W / (self.rho * self.V)


    def T_from_H(self, h, FAR, Thi, Tlo):
        """ Finding temperature from enthalpy polynomial using the Bisection Method """
        h_lo = self.get_ht(Tlo, FAR)
        h_hi = self.get_ht(Thi, FAR)
        if not (h_lo <= h <= h_hi or h_hi <= h <= h_lo): raise ValueError("Bisection bounds do not bracket the solution")

        # Initial calculation of hmid to prevent uninitialized use
        Tmid = (Thi + Tlo) / 2
        hmid = self.get_ht(Tmid, FAR)
        error = abs(hmid - h) / h
        iterations = 0
        
        # Perform bisection
        while error > 0.001:
            Tmid = (Thi + Tlo) / 2
            hmid = self.get_ht(Tmid, FAR)
            if hmid < h: Tlo = Tmid
            elif hmid > h: Thi = Tmid
            iterations += 1
            T = Tmid
            error = abs(hmid - h) / h
        return T

    def get_properties(self):
        if self.M != None: station_data = [self.idx, self.W, self.Tt, self.Pt/1000, self.ht/1000, self.Wc, self.Wf, self.FAR, self.M, self.T, self.h/1000, self.P/1000, self.V, self.rho, self.area]
        else: station_data = [self.idx, self.W, self.Tt, self.Pt/1000, self.ht/1000, self.Wc, self.Wf, self.FAR, None, None, None, None, None, None, None]
        return station_data


class Inlet:
    def __init__(self, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        M_inlet = cycle_parameters["inlet M"]
        M_exit = cycle_parameters["exit M"]
        Pt_recovery = cycle_parameters["total pressure recovery"]
        engine = cycle_parameters["engine"]
        freestream = engine.ambient

        if hasattr(engine, "fan_diameter"):
            D1 = engine.fan_diameter * 2.54 / 100 # in to m
            [_, Tt_T, Pt_P, rhot_rho, _] = isentropic(M_inlet, freestream.gamma, lookup_key="M")
            T1 = (1/Tt_T) * freestream.Tt
            a1 = numpy.sqrt(freestream.gamma * freestream.R * T1)
            V1 = M_inlet * a1
            A1 = numpy.pi * (D1 / 2)**2
            rho1 = freestream.rhot * (1/rhot_rho)
            W1 = V1 * rho1 * A1
        elif hasattr(engine, "W"):
            W1 = engine.W

        FAR = 0
        self.freestream = Station(W1, freestream.Tt, freestream.Pt, FAR, M=freestream.Minf, idx=0) # Station 0
        self.inlet = Station(W1, freestream.Tt, freestream.Pt, FAR, M=M_inlet, idx=1) # Station 1
        
        # Check if it's a BypassEngine
        if isinstance(engine, BypassEngine):
            bypass_W = W1*engine.B / (engine.B+1)
            root_W = W1 - bypass_W
            self.bypass_exit = Station(bypass_W, freestream.Tt, freestream.Pt, FAR, M=M_exit, idx=1.2)
            self.root_exit = Station(root_W, freestream.Tt, freestream.Pt, FAR, M=M_exit, idx=2)
        else:
            self.exit = Station(W1, freestream.Tt, freestream.Pt*Pt_recovery, FAR, M=M_exit, idx=2) # Station 1

        # COMPONENT DESIGN
        if component_parameters != None: pass


class Compressor:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None, root_upstream=None):
        # CYCLE ANALYSIS
        engine = cycle_parameters["engine"]
        self.e_c = cycle_parameters["e"]

        # Handle fans and normal compressors
        if "fan" in cycle_parameters:
            machine = "axial"
            self.is_fan = True
            self.cooling = 0
            self.packing = 0
            self.B = engine.B
            self.rootPR = cycle_parameters["root PR"]
            self.tipPR = cycle_parameters["tip PR"]
            Mroot_exit = cycle_parameters["root exit M"]
            Mtip_exit = cycle_parameters["tip exit M"]
            self.root_inlet = copy.deepcopy(root_upstream)
            self.tip_inlet = copy.deepcopy(upstream)
            self.root_inlet.idx = "2"
            self.tip_inlet.idx = "1.2"
            self.root_exit = self.solve_exit(self.root_inlet, "2.05", self.rootPR, 0, 0, Mroot_exit)
            self.tip_exit = self.solve_exit(self.tip_inlet, "1.3", self.tipPR, 0, 0, Mtip_exit)
            root_delta_ht = self.root_exit.ht - self.root_inlet.ht
            tip_delta_ht = self.tip_exit.ht - self.tip_inlet.ht
            self.delta_ht = root_delta_ht + tip_delta_ht
            self.power = (self.root_inlet.W*root_delta_ht) + (self.tip_inlet.W*tip_delta_ht)
        else:
            machine = cycle_parameters["machine"]
            inlet_idx = cycle_parameters["inlet idx"]
            exit_idx = cycle_parameters["exit idx"]
            self.cooling = cycle_parameters["cooling"]
            self.packing = cycle_parameters["packing"]
            PR = cycle_parameters["PR"]
            M_exit = cycle_parameters["exit M"]
            self.inlet = upstream
            self.inlet.idx = inlet_idx
            self.exit = self.solve_exit(self.inlet, exit_idx, PR, self.cooling, self.packing, M_exit)
            self.delta_ht =  (self.exit.cp*self.exit.Tt) - (self.inlet.cp*self.inlet.Tt)
            self.power = self.inlet.W * self.delta_ht

        # COMPONENT DESIGN
        if component_parameters != None:
            match machine:
                case "axial":
                    pass
                case "radial":
                    pass


    def solve_exit(self, inlet, exit_idx, PR, cooling, packing, M_exit):
        exit = copy.deepcopy(inlet)
        exit.idx = exit_idx
        exit.M = M_exit
        exit.W = inlet.W * (1 - cooling - packing)
        exit.Pt = inlet.Pt * PR
        exit.Tt = inlet.Tt * (PR)**((inlet.gamma - 1) / (inlet.gamma*self.e_c))
        exit.set_statics(exit.M)
        return exit


class Burner:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        self.Ptloss_b = cycle_parameters["total pressure loss"]
        self.eta_b = cycle_parameters["efficiency"]
        self.M_exit = cycle_parameters["exit M"]
        self.engine = cycle_parameters["engine"]
        self.LHV = self.engine.LHV
        self.solve_exit(upstream)
        # COMPONENT DESIGN
        if component_parameters != None: pass
    
    def solve_exit(self, upstream):
        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = 4
        self.exit.M = self.M_exit
        if hasattr(self.engine, "TET"):
            TET = self.engine.TET
            FAR = self.get_FAR(TET, self.inlet.Tt, self.inlet.FAR, self.LHV, self.eta_b)
            self.exit.Wf = self.inlet.W * FAR
        elif hasattr(self.engine, "Wf"):
            self.exit.Wf = self.engine.Wf
            FAR = self.exit.Wf / upstream.W
            TET = bisection(self.get_FAR, FAR, 1800, 100, self.inlet.Tt, 0, self.LHV, self.eta_b)
        self.exit.W = self.inlet.W * (1 + self.exit.FAR)
        self.exit.Pt = self.inlet.Pt * (1 - self.Ptloss_b)
        self.exit.Tt = TET
        self.exit.set_statics(self.exit.M)

    def get_FAR(self, T2, T1, FAR1, LHV, eta):
        FARnew = 0.02
        FAR = -1 
        h1 = self.inlet.get_ht(T1, FAR1)
        tolerance = 0.00001
        error = (abs(FAR - FARnew) / FARnew) 
        while error > tolerance:
            FAR = FARnew
            h2 = self.exit.get_ht(T2, FAR)
            FARnew = (h2 - h1) / (LHV * eta)
            error = (abs(FAR - FARnew) / FARnew) 
        return FARnew
        

class Turbine:
    def __init__(self, upstream, compressor, cycle_parameters=None, component_parameters=None):
        if cycle_parameters != None:
            # CYCLE ANALYSIS
            machine = cycle_parameters["machine"]
            self.engine = cycle_parameters["engine"]
            self.upstream = upstream
            self.compressor = compressor
            self.e_t = cycle_parameters["e"]
            self.mdot_cool = self.compressor.cooling
            self.power = abs(self.compressor.power)
            self.eta_m = cycle_parameters["mechanical efficiency"]
            self.exit_idx = cycle_parameters["exit idx"]
            self.M_exit = cycle_parameters["exit M"]
            self.cycle_parameters = cycle_parameters
            self.solve_exit()
        elif component_parameters != None:
            # COMPONENT DESIGN
            machine = component_parameters["machine"]
            if machine != "turbine": raise ValueError("Parameters must be for a turbine!")
            flow = component_parameters["flow type"]
            self.component_parameters = component_parameters
            self.specification = self.component_parameters["specification"]
            self.ER = self.specification["ER"]
            self.power = self.specification["power"]
            self.rpm = self.specification["rpm"]
            self.delta_ht = self.power / self.specification["W"] * 1000
            self.stages = []
            match flow:
                case "axial":
                    for parameters in self.component_parameters["stages"]:
                        # Handle the upstream station of each stage
                        idx = parameters["idx"]
                        if idx == 1:
                            # First Stage
                            W = self.specification["W"]
                            Wc = self.specification["ER"]
                            Tt = self.specification["Tt"]
                            Pt = self.specification["Pt"]
                            R = self.specification["R"]
                            gamma = self.specification["gamma"]
                            Cp = self.specification["Cp"]
                            FAR = self.specification["FAR"]
                            upstream = Station(W, Tt, Pt, FAR=FAR)
                            parameters["upstream"] = upstream
                        else:
                            # Subsequent stages
                            parameters["upstream"] = self.stages[idx - 1].upstream
                        self.stages.append(TurbineStage(self, machine, flow, parameters))
                case "radial": pass

    # Turbine Cycle Analysis (design-point)
    def solve_exit(self):
        self.inlet = self.upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = self.exit_idx
        self.exit.M = self.M_exit
        self.exit.W = self.inlet.W + self.mdot_cool
        exit_ht = (self.inlet.ht * self.inlet.W - (self.power/self.eta_m)) / self.exit.W
        FAR = self.inlet.W * self.inlet.FAR / (self.exit.W - (self.inlet.W * self.inlet.FAR))
        self.exit.Wf = (self.inlet.W - self.inlet.Wf) * FAR
        self.exit.Tt = self.inlet.T_from_H(exit_ht, self.exit.FAR, self.engine.TET, 100)
        ER = (self.exit.Tt  / self.inlet.Tt)**(-self.inlet.gamma / ((self.inlet.gamma - 1)*self.e_t))
        self.exit.Pt = self.inlet.Pt / ER
        self.exit.set_statics(self.exit.M)


# Velocity Triangle for axial turbomachines
class VelocityTriangle:
    def __init__(self, radius, omega, Vu, Vax, alpha, station=None, Mabs=None, Mrel=None):
        self.station = station
        self.radius = radius
        self.omega = omega
        self.U = self.radius * self.omega

        # Absolute Velocities
        self.Vu = Vu
        self.Vax = Vax
        self.V = numpy.sqrt(self.Vu**2 + self.Vax**2)
        self.alpha = alpha

        # Relative Velocities
        self.Wax = self.Vax
        self.Wu = self.Vu - self.U
        if self.Wu == 0:
            self.W = 0
        else:
            self.W = numpy.sqrt(self.Wu**2 + self.Vax**2)
        self.beta = numpy.tan(self.Wu / self.Vax)

        # Mach Numbers
        self.Mabs = Mabs
        self.Mrel = Mrel

    @property
    def Mabs(self): return self.get_Mabs()

    @Mabs.setter
    def Mabs(self, value): self._Mabs = value

    def get_Mabs(self): 
        if self.station is None:
            return None
        else:
            return self.V / numpy.sqrt(self.station.gamma * self.station.R * self.station.T)

    @property
    def Mrel(self): return self.get_Mrel()

    @Mrel.setter
    def Mrel(self, value): self._Mrel = value

    def get_Mrel(self): 
        if self.station is None:
            return None
        else:
            return self.W / numpy.sqrt(self.station.gamma * self.station.R * self.station.T)


# Stations for axial turbomachines (adds radii and velocity triangles to the engine-flow stations)
class AxialStation(Station):
    def __init__(self, idx, W, Tt, Pt, omega, mid:VelocityTriangle, num_radii=None, FAR=None, M=None):
        super().__init__(W, Tt, Pt, FAR=FAR, M=M, idx=M)
        self.idx = idx
        self.mid = mid
        self.omega = omega
        self.num_radii = num_radii
        self.triangles = []
        
        # Solve for velocity triangles at every radius along the blade
        if self.num_radii is not None: self.define_velocity_triangles()
            

    def set_statics(self, M=None):
        if self.M is not None:
            self.M = M
        [_, Tt_T, Pt_P, _, _] = isentropic(self.M, self.gamma, lookup_key="M")
        self.T = (1 / Tt_T) * self.Tt
        self.h = self.T * self.cp
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        self.area = self.W / (self.rho * self.mid.Vax) # Calculate area based on the axial velocity since the velocity has a tangential component
        self.rhub = self.mid.radius - self.area/(4*numpy.pi*self.mid.radius)
        self.rtip = self.mid.radius + self.area/(4*numpy.pi*self.mid.radius)
        
    
    def define_velocity_triangles(self, num_radii=None):
        if num_radii is not None:
            self.num_radii = num_radii
        self.radii = list(numpy.linspace(self.rhub, self.rtip, self.num_radii, endpoint=True, dtype=float))
        for radius in self.radii:
            U = radius * self.omega
            Vu = self.mid.Vu # Free Vortex Equation
            alpha = numpy.atan(Vu / self.mid.Vax)
            self.triangles.append(VelocityTriangle(self, U, Vu, self.mid.Vax, alpha))

            
    def plot_triangles(self):
        #for triangle in self.triangles: triangle.plot_triangle
        pass


class BladeGeometry:
    def __init__(self, machine, flow, stage, blade:str, parameters:dict=None):
        self.machine = machine
        self.stage = stage
        self.blade = blade
        self.parameters = parameters

        if self.parameters is not None:
            match machine:
                case "turbine":
                    match flow:
                        case "axial":
                            self.AR = parameters["AR"]
                            self.zweiffel = parameters["zweiffel"]
                            match self.blade.lower():
                                case "stator":
                                    # AXIAL TURBINE STATOR 
                                    self.get_axial_turbine_stator(self.parameters)
                                case "rotor":
                                    # AXIAL TURBINE ROTOR 
                                    self.get_axial_turbine_rotor(self.parameters)
                        case "radial":
                            pass
                case "compressor":
                    match flow:
                        case "axial":
                            match self.blade.lower():
                                case "stator":
                                    # AXIAL COMPRESSOR STATOR 
                                    self.get_axial_compressor_stator(self.parameters)
                                case "rotor":
                                    # AXIAL COMPRESSOR ROTOR 
                                    self.get_axial_compressor_rotor(self.parameters)
                        case "radial":
                            pass

    def get_axial_turbine_stator(self, parameters):
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        S_rh_avg = (s1.rhub + s2.rhub) / 2
        S_rt_avg = (s1.rtip + s2.rtip) / 2
        self.HT_ratio = S_rh_avg / S_rt_avg
        self.h = S_rt_avg - S_rh_avg
        self.chord = self.h / self.AR
        
        # Stagger and Axial Chord
        self.stagger = list()
        self.cax = list()
        self.deflections = list()
        for radius_idx in range(self.stage.num_radii):
            self.stagger.append((s1.triangles[radius_idx].alpha + s2.triangles[radius_idx].alpha) / 2)
            self.cax.append(self.chord * numpy.cos(self.stagger[radius_idx]))
            self.deflections.append(s2.triangles[radius_idx].alpha - s1.triangles[radius_idx].alpha)
        
        # Stator specific parameters
        self.taper_ratio = self.cax[-1] / self.cax[0]
        self.axial_spacing = 0.25 * self.cax[(self.stage.num_radii + 1) / 2]
        
        # Solidity and Pitch
        self.solidity = (2/self.zweiffel) * numpy.cos(s2.mid.alpha)**2 * (numpy.tan(s2.mid.alpha) - numpy.tan(s1.mid.alpha))
        self.pitch = self.chord / self.solidity
        self.NOB = numpy.ceil((2 * numpy.pi * s1.mid.radius) / self.pitch)
        
        # Opening
        self.os = numpy.cos(s2.mid.alpha)
        self.opening = self.os * self.pitch

    def get_axial_turbine_rotor(self, parameters):
        # ROTOR 
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        R_rh_avg = (s2.rhub + s3.rhub) / 2
        R_rt_avg = (s2.rtip + s3.rtip) / 2
        self.HT_ratio = R_rh_avg / R_rt_avg
        self.h = R_rt_avg - R_rh_avg
        self.chord = self.h / self.AR
        
        # Stagger and Axial Chord
        self.stagger = list()
        self.cax = list()
        self.deflections = list()
        for radius_idx in range(self.stage.num_radii):
            self.stagger.append((s2.triangles[radius_idx].beta + s3.triangles[radius_idx].beta) / 2)
            self.cax.append(self.chord * numpy.cos(self.stagger[radius_idx]))
            self.deflections.append(s3.triangles[radius_idx].beta - s2.triangles[radius_idx].beta)
        
        # Stator specific parameters
        self.taper_ratio = self.cax[-1] / self.cax[0]
        self.axial_spacing = 0.25 * self.cax[(self.stage.num_radii + 1) / 2]
        
        # Solidity and Pitch
        self.solidity = (2/self.zweiffel) * numpy.cos(s2.mid.beta)**2 * (numpy.tan(s2.mid.beta) - numpy.tan(s1.mid.beta))
        self.pitch = self.chord / self.solidity
        self.NOB = numpy.ceil((2 * numpy.pi * s1.mid.radius) / self.pitch)
        
        # Opening
        self.os = numpy.cos(s2.mid.beta)
        self.opening = self.os * self.pitch


# General stage object to be used for turbines (handles stage calculations, not be used outside of component classes)
class TurbineStage:
    def __init__(self, turbine, flow, parameters):
        self.turbine = turbine
        self.parameters = parameters
        match flow: 
            case "axial":
                # AXIAL TURBINE
                self.solve_axial()
            case "radial": 
                # RADIAL TURBINE
                self.solve_radial()

    def solve_axial(self):
        # Load in stage parameters
        self.idx = self.parameters["idx"]
        self.aerodynamics = self.parameters["aerodynamics"]
        self.geometry = self.parameters["geometry"]
        self.cooling = self.parameters["cooling"]
        self.upstream = self.parameters["upstream"]
        self.num_radii = self.parameters["number of radii"]
        self.omega = self.turbine.rpm * (2*numpy.pi / 60)
        efficiency = self.parameters["specification"]["polytropic efficiency"]
        phi = self.aerodynamics["load coefficient"]
        psi = self.aerodynamics["work coefficient"]
        loss_coefficient = self.aerodynamics["loss coefficient"]
        M2m = self.aerodynamics["M2m"]
        U3m = self.aerodynamics["U3m"]
        Vu1m = self.aerodynamics["Vu1m"]
        Rm2_Rm1 = self.geometry["rm climb rate stator"]
        Rm3_Rm2 = self.geometry["rm climb rate rotor"]
        Vax2_Vax1 = self.geometry["Vax climb stator"]
        Vax3_Vax2 = self.geometry["Vax climb rotor"]
        AR_stator = self.geometry["AR stator"]
        AR_rotor = self.geometry["AR rotor"]
        zweiffel = self.geometry["Zweiffel"]

        # Ensure odd number of radii (must have a "mid-radius")
        if self.num_radii % 2 != 0 or self.num_radii < 3: raise ValueError("Number of radii per blade must be ")

        # Radii & Axial Velocities
        Rm3 = U3m * self.omega
        Rm2 = Rm3 / Rm3_Rm2
        Vax3 = U3m * phi
        Vax2 = Vax3 / Vax3_Vax2
        if self.idx == 1:
            # First Stage
            Rm1 = Rm2 / Rm2_Rm1
            Vax1 = Vax2 / Vax2_Vax1
        else:
            # Subsequent Stages
            Rm1 = self.upstream.radii[-1]
            Vax1 = self.upstream.mid.Vax
        
        # Stage Quantites
        self.delta_ht = psi * (Rm3*self.omega)**2
        self.work_split = self.delta_ht / self.turbine.delta_ht
        self.capacity = self.upstream.W * numpy.sqrt(self.upstream.Tt) / (self.upstream.Pt / 101.325)

        # Meanline and Radial Calculations (Velocity Triangles)
        # Station 1
        W1 = self.upstream.W
        Tt1 = self.upstream.Tt
        Pt1 = self.upstream.Pt
        alpha1 = numpy.atan(Vu1m / Vax1)
        mid1 = VelocityTriangle(Rm1, self.omega, Vu1m, Vax1, alpha1)
        s1 = AxialStation(W1, Tt1, Pt1, mid=mid1)
        s1.T = s1.Tt - s1.mid.V**2/(2*s1.cp)
        s1.M = numpy.sqrt((2/(s1.gamma - 1)) * (s1.Tt/s1.T) - 1)
        s1.set_statics()
        s1.define_velocity_triangles(self.num_radii)

        # Station 2
        s2 = copy.deepcopy(s1)
        s2.Pt = s1.Pt - self.loss_coefficient*(s1.Pt - s1.P)
        s2.set_statics(M2m)
        V2 = M2m * numpy.sqrt(s2.gamma * s2.R * s2.T)
        alpha2 = numpy.acos(Vax2 / V2)
        Vu2 = V2 * numpy.sin(alpha2)
        s2.mid = VelocityTriangle(Rm2, self.omega, Vu2, Vax2, alpha2)
        s2.define_velocity_triangles()

        # Station 3
        s3 = copy.deepcopy(s2)
        s3.Tt = s2.Tt - self.delta_ht/s2.cp
        s3.Pt = s2.Pt * (s3.Tt/s2.Tt)**(s3.gamma/(efficiency*(s3.gamma - 1)))
        Vu3 = (Rm2*self.omega*s2.mid.Vu - self.delta_ht) / (Rm3*self.omega) # Euler Turbine Equation
        alpha3 = numpy.atan(Vu3/Vax3)
        s3.mid = VelocityTriangle(Rm3*self.omega, Vu3, Vax3, alpha3)
        s3.T = s3.Tt - s3.mid.V**2/(2*s3.cp)
        s3.M = numpy.sqrt((2/(s3.gamma - 1)) * (s3.Tt/s3.T) - 1)
        s3.set_statics()
        s3.define_velocity_triangles()

        # Store stations in dictionary object
        self.stations = {1: s1, 2: s2, 3: s3}

        # Degrees of Reaction (DoR)
        self.get_DoR()

        # Deflections
        self.get_deflections()

        # Expansion Ratio
        self.ER = self.stations[1].Pt / self.stations[3].Pt

        # Blade Geometries
        self.stator = BladeGeometry(machine="turbine", flow="axial", stage=self, blade="stator", parameters={"AR": AR_stator, "zweiffel": zweiffel})
        self.rotor = BladeGeometry(machine="turbine", flow="axial", stage=self, blade="rotor", parameters={"AR": AR_rotor, "zweiffel": zweiffel})


    def solve_radial(self): pass

    
    def get_DoR(self):
        self.DoR = list()
        s1 = self.stations[1]
        s2 = self.stations[2]
        s3 = self.stations[3]
        for radius_idx in range(self.num_radii):
            h2 = s2.ht - s2.triangles[radius_idx].V**2/2
            h3 = s3.ht - s3.triangles[radius_idx].V**2/2
            ht1 = s1.ht
            ht3 = s3.ht
            self.DoR.append((h2 - h3) / (ht1 - ht3))


    def get_deflections(self):
        self.deflections = list()
        s1 = self.stations[1]
        s2 = self.stations[2]
        s3 = self.stations[3]
        for radius_idx in range(self.num_radii):
            stator_deflection = s2.triangles[radius_idx].alpha - s1.triangles[radius_idx].alpha
            rotor_deflection = s3.triangles[radius_idx].beta - s1.triangles[radius_idx].beta
            self.deflections.append({"stator": stator_deflection, "rotor": rotor_deflection})


    def axial_cooling(self, inlet:AxialStation, blade:str, parameters:dict):
        s1 = self.stations[1]
        s2 = self.stations[2]
        s3 = self.stations[3]
        OTDF = parameters["OTDF"]
        RTDF = parameters["RTDF"]
        CDT = parameters["CDT"]
        cp_coolant = parameters["Cp coolant"]
        metal_TtMax = parameters["Metal TtMax"]
        cooling_efficiency = parameters["cooling efficiency"]

        # STATOR
        # Peak gas temperature seen by stator (using Inlet Total Temp)
        T_peak_S = OTDF * (inlet.Tt - CDT) + inlet.Tt
        
        # Corrected temperature based on recovery factor (0.15 logic from original)
        T_peak_S_corr = T_peak_S - 0.15 * (s2.mid.V)^2 / (2 * cp_coolant)
        
        # Required cooling effectiveness
        eps_req_S = (T_peak_S_corr - metal_TtMax) / (T_peak_S_corr - CDT)
        
        # Non-dimensional coolant parameter
        K_cool_S = eps_req_S / (cooling_efficiency * (1 - eps_req_S))
        
        # Heat Transfer: Stator
        Area_S = 2 * self.stator.chord * self.stator.h * self.stator.NOB
        mu_S = 1.458e-6 * T_peak_S_corr^(3/2) / (T_peak_S_corr + 110.4)
        Re_S = (s1.rho * s1.mid.V * self.stator.chord) / mu_S
        k_S = 0.000053983 * T_peak_S_corr + 0.013568
        Nu_S = 0.488 * Re_S^0.592
        h_conv_S = Nu_S * k_S / self.stator.chord
        
        # Resulting Stator Coolant Flow
        self.stator.mdot_cool = K_cool_S * h_conv_S * Area_S / cp_coolant

        # ROTOR
        # Peak gas temperature seen by rotor (using RTDF)
        T_peak_R = RTDF * (inlet.Tt - CDT) + inlet.Tt
        
        # Corrected temperature based on relative velocity (W) for the rotor
        T_peak_R_corr = T_peak_R - 0.15 * (s2.mid.W)^2 / (2 * cp_coolant)
        
        # Required cooling effectiveness
        eps_req_R = (T_peak_R_corr - metal_TtMax) / (T_peak_R_corr - CDT)
        
        # Non-dimensional coolant parameter
        K_cool_R = eps_req_R / (cooling_efficiency * (1 - eps_req_R))
        
        # Heat Transfer: Rotor
        Area_R = 2 * self.rotor.chord * self.rotor.h * self.rotor.NOB
        mu_R = 1.458e-6 * T_peak_R_corr^(3/2) / (T_peak_R_corr + 110.4)
        Re_R = (s2.rho * s2.mid.W * self.rotor.Chord) / mu_R
        k_R = 0.000053983 * T_peak_R_corr + 0.013568
        Nu_R = 0.488 * Re_R^0.592
        h_conv_R = Nu_R * k_R / self.rotor.Chord
        
        # Resulting Rotor Coolant Flow
        self.rotor.mdot_cool = K_cool_R * h_conv_R * Area_R / cp_coolant
        
        # Total Stage Cooling
        self.total_mdot_cool = self.stator.mdot_cool + self.rotor.mdot_cool


    def display_stage_results(self, Flags):
        # Display Velocity Triangles
        if Flags["plot triangles"]:
            pass
            '''
            plot_triangles(stage.stations(1).Hub, stage.stations(1).Vax, 'Stator LE - Hub')
            plot_triangles(stage.stations(2).Hub, stage.stations(2).Vax, 'Stator TE / Rotor LE - Hub')
            plot_triangles(stage.stations(3).Hub, stage.stations(3).Vax, 'Rotor TE - Hub')
            plot_triangles(stage.stations(1).Mid, stage.stations(1).Vax, 'Stator LE - Mid')
            plot_triangles(stage.stations(2).Mid, stage.stations(2).Vax, 'Stator TE / Rotor LE - Mid')
            plot_triangles(stage.stations(3).Mid, stage.stations(3).Vax, 'Rotor TE - Mid')
            plot_triangles(stage.stations(1).Tip, stage.stations(1).Vax, 'Stator LE - Tip')
            plot_triangles(stage.stations(2).Tip, stage.stations(2).Vax, 'Stator TE / Rotor LE - Tip')
            plot_triangles(stage.stations(3).Tip, stage.stations(3).Vax, 'Rotor TE - Tip')
            '''
        # Display Meridional View
        if Flags["plot meridional"]:
            pass
            #plot_meridional_view(stage)
        # Display Blade Sections (and Export as CSV Files)
        if Flags["plot sections"]:
            pass
            '''
            export_section_coordinates(stage, 'Rotor', 'HPT_Rotor')
            export_section_coordinates(stage, 'Stator', 'HPT_Stator')
            '''
        # Display Velocity Tables
        if Flags["display velocities"]:
            pass
            '''
            display_velocity_table(stage, 'Hub')
            display_velocity_table(stage, 'Mid')
            display_velocity_table(stage, 'Tip')
            '''
        # Display Thermodynamics Table
        if Flags["display thermodynamics"]:
            pass
            #display_thermo_table(stage)
        # Display Geometry Table
        if Flags["display geometry"]:
            pass
            #display_geometry_table(stage)
        # Display Performance
        if Flags["display performance"]:
            pass
            #display_performance(stage, parameters.Specification)



class CompressorStage:
    def __init__(self, flow, parameters):
        match flow:
            case "axial":
                # AXIAL COMPRESSOR
                self.solve_axial()
            case "radial": 
                # RADIAL COMPRESSOR
                self.solve_axial()
    
    def solve_axial(): pass
    def solve_radial(): pass


class Mixer:
    def __init__(self, hot_inlet:Station, cold_inlet:Station, component_parameters=None): 
        # CYCLE ANALYSIS
        # NOTE The hot stream Mach Number is chosen, but the cold stream Mach Number can equivalently be chosen

        # Hot Stream (static properties)
        hot_inlet.T = hot_inlet.Tt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)**(-1)
        hot_inlet.P = hot_inlet.Pt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)**(-hot_inlet.gamma / (hot_inlet.gamma - 1))
        hot_inlet.rho = hot_inlet.P / (hot_inlet.R*hot_inlet.T)
        hot_inlet.V = hot_inlet.M * numpy.sqrt(hot_inlet.gamma * hot_inlet.R * hot_inlet.T)
        hot_inlet.area = hot_inlet.W / (hot_inlet.rho*hot_inlet.V)
        self.hot_inlet = hot_inlet
        momentum_hot = hot_inlet.W*hot_inlet.V + hot_inlet.area*hot_inlet.P

        # Cold Stream (static properties)
        cold_inlet.P = hot_inlet.P
        cold_inlet.M = numpy.sqrt(2 / (cold_inlet.gamma - 1) * ((cold_inlet.Pt/cold_inlet.P)**((cold_inlet.gamma - 1)/cold_inlet.gamma) - 1))
        cold_inlet.T = cold_inlet.Tt * (1 + (cold_inlet.gamma - 1)/2 * cold_inlet.M**2)**(-1)
        cold_inlet.rho = cold_inlet.P / (cold_inlet.R*cold_inlet.T)
        cold_inlet.V = cold_inlet.M * numpy.sqrt(cold_inlet.gamma * cold_inlet.R * cold_inlet.T)
        cold_inlet.area = cold_inlet.W / (cold_inlet.rho*cold_inlet.V)
        self.cold_inlet = cold_inlet
        momentum_cold = cold_inlet.W*cold_inlet.V + cold_inlet.area*cold_inlet.P

        # Design Parameter (hot stream Mach Number) Restriction
        if cold_inlet.P > cold_inlet.Pt: raise ValueError("Static pressure is larger than the total pressure (not physically possible). A larger Mhot must be chosen.")

        # Mixer Exit
        idx = "6"
        W = hot_inlet.W + cold_inlet.W
        Wf = hot_inlet.Wf
        A = hot_inlet.area + cold_inlet.area
        Tt = (hot_inlet.W*hot_inlet.Tt + cold_inlet.W*cold_inlet.Tt) / W
        FAR = Wf / (W - Wf)
        R = hot_inlet.get_R(FAR)
        Cp = hot_inlet.get_cp(Tt, FAR)
        gamma = Cp / (Cp - R)
        K = (momentum_hot + momentum_cold)**2 * (gamma / (W**2 * Tt * R))
        M = numpy.sqrt(((2*gamma - K) + numpy.sqrt((K - 2*gamma)**2 + 4*(((gamma-1)/2)*K-gamma**2))) / (2*((gamma-1)/2*K-gamma**2)))
        T = Tt * (1 + (gamma - 1)/2 * M**2)**(-1)
        P = (W*numpy.sqrt(R*T)/(A*M*numpy.sqrt(gamma)))
        Pt = P * (Tt/T)**(gamma / (gamma - 1))
        self.exit = Station(W, Tt, Pt, FAR=FAR, idx=idx, M=M)

        # COMPONENT DESIGN
        if component_parameters != None: pass


class Afterburner:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        engine = cycle_parameters["engine"]
        self.toggle = cycle_parameters["toggle"]
        self.Ttmax = cycle_parameters["AET"]
        self.pi_ab_on = cycle_parameters["pi_hot"]
        self.pi_ab_off = cycle_parameters["pi_cold"]
        self.eta = cycle_parameters["efficiency"]
        self.LHV = engine.LHV
        self.solve_exit(upstream)

        # COMPONENT DESIGN
        if component_parameters != None: pass


    def solve_exit(self, upstream):
        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = "7"

        # Handle different engine modes (afterburner on or off)
        match self.toggle: 
            case False:
                self.exit.Pt = self.inlet.Pt * self.pi_ab_off # Afterburner off
            case True: # Afterburner on
                self.exit.Pt = self.inlet.Pt * self.pi_ab_on
                Tt_Ttstar1 = get_Tt_Ttstar(self.inlet.M, self.inlet.gamma)
                Tt_Ttstar2 = self.Ttmax * (Tt_Ttstar1) / self.inlet.Tt
                Ttstar = 1 / Tt_Ttstar1 * self.inlet.Tt

                # Check if the afteburner chokes before reaching the given total
                # temperature (based on Rayleigh flow theory)
                if Tt_Ttstar2 > 1:
                    print(f"Given Afterburner Exit Temperature exceeds maximum allowed amount {Ttstar}.\n"
                        "This value has been set as the exit total temperatue and a lower value is recommended to be chosen.\n\n")
                    self.exit.Tt = Ttstar
                    self.exit.M = 1
                elif Tt_Ttstar2 <= 1:
                    self.exit.Tt = self.Ttmax
                    self.exit.M = bisection(get_Tt_Ttstar, Tt_Ttstar2, 1, 0, self.inlet.gamma)

                self.exit.FAR = Burner.get_FAR(self, self.exit.Tt, self.inlet.Tt, self.inlet.FAR, self.LHV, self.eta) 
                self.exit.Wf = (self.inlet.W - self.inlet.Wf) * self.exit.FAR
                self.exit.W = self.inlet.W + self.exit.Wf
                self.exit.set_statics(self.exit.M)


class Nozzle:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        self.geometry = cycle_parameters["C/CD"]
        self.engine = cycle_parameters["engine"]
        self.Pinf = self.engine.ambient.P
        self.solve_exit(upstream)
        
        # COMPONENT DESIGN
        if component_parameters != None: pass

    def solve_exit(self, upstream):
        self.inlet = upstream

        # Handle Converging or Converging-Diverging nozzle geometries
        match self.geometry:
            case "C": self.converging(self.Pinf)
            case "CD": self.CD(self.Pinf)

    def converging(self,  Pinf):
        gamma = self.inlet.gamma
        cp = self.inlet.cp
        R = self.inlet.R
        critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
        NPR = self.inlet.Pt / Pinf
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = 8 # Throat station

        if NPR >= critical_NPR:
            # Choked
            self.exit.M = 1
            self.exit.T = (2 / (gamma + 1)) * self.exit.Tt
            self.exit.V = numpy.sqrt(gamma * R * self.exit.T)
            self.exit.P = self.exit.Pt * (1 / critical_NPR)
        else:
            # Unchoked
            self.exit.P = Pinf
            self.exit.T = self.exit.Tt * (self.exit.P/self.exit.Pt)**((gamma - 1) / gamma)
            self.exit.V = numpy.sqrt(2 * cp * (self.exit.Tt - self.exit.T))
            self.exit.M = self.exit.V / numpy.sqrt(gamma * R * self.exit.T)

        self.exit.set_statics(self.exit.M)
    

    def CD(self, Pinf):
        # This function forces the throat to always be choked (CHECK THIS)
        gamma = self.inlet.gamma
        cp = self.inlet.cp
        R = self.inlet.R
        self.exit = self.throat = copy.deepcopy(self.inlet)
        self.throat.idx = 8
        self.exit.idx = 9

        # Throat Calculations
        if self.inlet.M < 1:
            self.throat.M = 1
            self.throat.Tt = self.inlet.Tt
            self.throat.ht = self.throat.Tt * cp
            self.throat.Pt = self.inlet.Pt
            self.throat.mdot = self.inlet.mdot
            self.throat.set_statics(self.throat.M)

        # Exit Plane Calculations
        self.exit.P = Pinf; # Assuming perfectly expanded flow
        self.exit.M = numpy.sqrt((2/(gamma - 1)) * ((self.exit.Pt/self.exit.P)^((gamma-1)/gamma) - 1)); # Isentropic Relation
        self.exit.set_statics(self.exit.M)


class Recuperator:
    def __init__(self, cold_inlet:Station, cycle_parameters, component_parameters=None):
        self.cyce_parameters = cycle_parameters
        self.engine = cycle_parameters["engine"]
        self.Mexit_cold = cycle_parameters["cold exit M"]
        self.Mexit_hot = cycle_parameters["hot exit M"]
        self.pi_cold = cycle_parameters["pi_cold"]
        self.pi_hot = cycle_parameters["pi_hot"]
        self.delta_ht = cycle_parameters["delta_ht"] * 10**3 # kJ/kg to J/kg

        self.cold_inlet = copy.deepcopy(cold_inlet)
        self.cold_inlet.idx = "3.06"
        self.cold_exit = copy.deepcopy(self.cold_inlet)
        self.cold_exit.idx = "3.07"
        self.cold_exit.Tt = self.cold_exit.T_from_H(self.cold_inlet.ht+self.delta_ht, self.cold_exit.FAR, 1800, self.cold_inlet.Tt)
        self.cold_exit.Pt *= self.pi_cold
        self.cold_exit.M = self.Mexit_cold
        self.cold_exit.set_statics(self.cold_exit.M)
    
        # COMPONENT DESIGN
        if component_parameters != None: pass

    def pass_hot_stream(self, hot_upstream):
        self.hot_inlet = hot_upstream
        self.hot_inlet.idx = "6.07"
        self.hot_exit = copy.deepcopy(self.hot_inlet)
        self.hot_exit.idx = "6.08"
        self.hot_exit.Tt = self.cold_exit.T_from_H(self.hot_inlet.ht-self.delta_ht, self.hot_exit.FAR, self.hot_inlet.Tt, 0)
        self.hot_exit.Pt *= self.pi_hot
        self.hot_exit.M = self.Mexit_hot
        self.hot_exit.set_statics(self.hot_exit.M)

        return self.hot_exit


class Bleed:
    def __init__(self, cooling, packing):
        self.cooling = cooling
        self.packing = packing


def format_axes(ax, title, ylabel):
    ax.set_title(title, fontsize=14, weight='bold')
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(frameon=False, loc="best")
    # Grid control
    ax.minorticks_on()
    ax.grid(True, which='major', linestyle='--', linewidth=0.7, alpha=0.7)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
    # Clean spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


'''
Essentially a turbojet core engine that will serve as a parent class to other architectures like turbofan, turboshaft, ramjet, etc.
Handles a recuperator and afterburner as well, but the default is just a core
'''
class Engine:
    def __init__(self, engine_parameters):
        self.spools = engine_parameters["spools"]
        self.LHV = engine_parameters["LHV"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.ambient = Ambient(self.altitude, Minf=self.Minf)

        # Handle user errors in engine parameters
        if "mdot" in engine_parameters and "fan diameter" in engine_parameters:
            raise ValueError("Error. Cannot have two sizing parameters. Must choose either inlet mass flow or fan face diameter.")
        else:
            if "mdot" in engine_parameters: self.W = engine_parameters["mdot"] # kg/sec
            elif "fan diameter" in engine_parameters: self.fan_diameter = engine_parameters["front diameter"] # inches
        if "TET" in engine_parameters and "mdotf" in engine_parameters:
            raise ValueError("Error. Cannot have two burner parameters. Must choose either TET or fuel mass flow rate.")
        else:
            if "TET" in engine_parameters: self.TET = engine_parameters["TET"] # K
            if "mdotf" in engine_parameters: self.Wf = engine_parameters["mdotf"] # kg/sec

    @property
    def components(self): return self._components

    @components.setter
    def components(self, value): self._components = value

    def set_components(self, parameters):
        self.parameters = parameters

        # Include the engine object in each set of component parameters
        for component_name in parameters:
            if component_name != "engine":
                parameters[component_name]["engine"] = self

        # Handle spool counts
        match self.spools:
            case 1:
                self.PR = parameters["engine"]["PR"]
                parameters["compressor"]["PR"] = self.PR
                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    compressor := Compressor(inlet.exit, parameters["compressor"]),
                    burner := Burner(compressor.exit, parameters["burner"]),
                    turbine := Turbine(burner.exit, compressor, parameters["turbine"]),
                    exhaust := Nozzle(turbine.exit, parameters["nozzle"])
                ]
                '''
                self.components = [inlet := Inlet(parameters["intake"])]
                parameters["compressor"]["upstream"] = inlet.exit
                parameters["burner"]["upstream"] = compressor.exit
                parameters["turbine"]["upstream"] = burner.exit
                parameters["turbine"]["compressor"] = compressor
                parameters["nozzle"]["upstream"] = turbine.exit
                self.components.append(compressor := Compressor(parameters["compressor"]))
                self.components.append(burner := Burner(parameters["burner"]))
                self.components.append(turbine := Turbine(parameters["turbine"]))
                self.components.append(exhaust := Nozzle(turbine.exit, parameters["nozzle"]))
                self.compressor = compressor
                self.turbine = turbine
                '''
                self.compressor = compressor
                self.turbine = turbine
            case 2:
                self.OPR = parameters["engine"]["OPR"]
                self.HPR = parameters["engine"]["HPR"]
                self.LPR = self.OPR / self.HPR
                parameters["hpc"]["PR"] = self.HPR
                parameters["lpc"]["PR"] = self.LPR

                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    lpc := Compressor(inlet.exit, parameters["lpc"]),
                    hpc := Compressor(lpc.exit, parameters["hpc"]),
                    burner := Burner(hpc.exit, parameters["burner"]),
                    hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                    lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                    exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                ]
                self.lpc = lpc
                self.hpc = hpc
                self.hpt = hpt
                self.lpt = lpt
            case 3:
                self.OPR = parameters["engine"]["OPR"]
                self.HPR = parameters["engine"]["HPR"]
                self.IPR = parameters["engine"]["IPR"]
                self.LPR = self.OPR / self.HPR / self.IPR
                parameters["hpc"]["PR"] = self.HPR
                parameters["ipc"]["PR"] = self.IPR
                parameters["lpc"]["PR"] = self.LPR

                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    lpc := Compressor(inlet.exit, parameters["lpc"]),
                    ipc := Compressor(lpc.exit, parameters["hpc"]),
                    hpc := Compressor(ipc.exit, parameters["hpc"]),
                    burner := Burner(hpc.exit, parameters["burner"]),
                    hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                    ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                    lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                    exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                ]
                self.lpc = lpc
                self.ipc = ipc
                self.hpc = hpc
                self.hpt = hpt
                self.ipt = ipt
                self.lpt = lpt
            case _:
                raise ValueError("Number of spools must be either 1, 2, or 3.")

        self.inlet = inlet
        self.burner = burner
        self.exhaust = exhaust
        # Check for a recuperator in the user input parameters
        if "recuperator" in parameters: self.add_recuperator(parameters["recuperator"])
        # Check for an afterburner in the user input parameters
        if "afterburner" in parameters: self.add_afterburner(parameters["afterburner"])


    # Add an afterburner to the engine
    def add_afterburner(self, parameters):
        insert_idx = len(self.components) - 1
        inlet_station = self.components[insert_idx - 1].exit
        afterburner = Afterburner(inlet_station, parameters)
        self.afterburner = afterburner
        self.components[insert_idx].solve_exit(afterburner.exit)
        self.components.insert(insert_idx, afterburner) 


    # Add a recuperator to the engine
    def add_recuperator(self, parameters):
        for idx, component in enumerate(self.components):
            if isinstance(component, Burner): 
                insert_idx = idx
                recuperator = Recuperator(self.components[idx - 1].exit, parameters)
                component.solve_exit(recuperator.cold_exit)
            elif isinstance(component, Turbine):
                component.solve_exit(self.components[idx-1].exit)
            elif not isinstance(component, Turbine) and isinstance(self.components[idx - 1], Turbine):
                component.solve_exit(recuperator.pass_hot_stream(self.components[idx-1].exit))
            elif idx > 2:
                component.solve_exit(self.components[idx - 1].exit)
        self.components.insert(insert_idx, recuperator)
        self.recuperator = recuperator

    
    # Add a mixer to the engine
    def add_mixer(self, parameters):
        pass


    # Retrieve flow properties at every station
    def get_station_data(self): 
        raw_data = list()
        # Handle Recuperator station data (afterburner is treated like the other components)
        if "recuperator" in self.parameters:
            for component in self.components:
                if isinstance(component, Inlet):
                    raw_data.append(component.freestream.get_properties())
                    raw_data.append(component.inlet.get_properties())
                    raw_data.append(component.exit.get_properties())
                elif isinstance(component, Recuperator):
                    raw_data.append(component.cold_exit.get_properties())
                elif isinstance(component, Turbine):
                    raw_data.append(component.exit.get_properties())
                    raw_data.append(self.recuperator.pass_hot_stream(component.exit).get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
        else:
            for component in self.components:
                if isinstance(component, Inlet):
                    raw_data.append(component.freestream.get_properties())
                    raw_data.append(component.inlet.get_properties())
                    raw_data.append(component.exit.get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
                    
        rounded_data = numpy.round(numpy.array(raw_data, dtype=float), 3).tolist()
        station_data = pandas.DataFrame(rounded_data, columns=Station.column_names)
        
        return station_data


    # Display the full engine performance
    def get_performance(self):
        """ Performance Parameters """
        station_data = self.get_station_data()
        W = station_data["W [kg/sec]"].values
        V = station_data["V [m/sec]"].values
        P = station_data["Ps [kPa]"].values
        A = station_data["Area [m^2]"].values
        FAR = station_data["FAR"].values
        Wf_in = self.burner.exit.Wf
        FARexit = self.exhaust.exit.FAR
        pressure_thrust = A[-1] * ((P[-1]*1000) - self.ambient.P)
        pressure_power = pressure_thrust * V[-1]
        # Specific Thrust
        T_ma = V[-1]*(1 + FAR[-1]) - V[0] + pressure_thrust/W[0] 
        # Thrust Specific Fuel Consumption (TSFC)
        TSFC = (FAR[-1] / T_ma) * 10**6 
        # Thrust
        thrust = T_ma * W[0] 
        # Propulsive Efficiency
        eta_p = thrust*V[0] / (0.5*W[-1]*V[-1]**2 - 0.5*W[0]*V[0]**2 + pressure_power) 
        # Thermal Efficiency
        eta_th = (0.5*W[-1]*V[-1]**2 - 0.5*W[0]*V[0]**2 + pressure_power) / (Wf_in*self.LHV) 
        # Overall Efficiency
        eta_o = eta_p * eta_th 
        performance = pandas.DataFrame({
                                         "Specific Thrust [m/sec]": T_ma,
                                         "TSFC [g/(sec*kN)]": TSFC,
                                         "Thrust [kN]": thrust/1000,
                                         "Propulsive Efficiency": eta_p,
                                         "Thermal Efficiency": eta_th,
                                         "Overall Efficiency": eta_o
                                        }, index=[0])
        return station_data, performance


    # Plot the temperatures and pressures throughout the whole engine
    def plot_thermo(self):
        Tt = self.get_station_data()["Tt [K]"]
        Ts = self.get_station_data()["Ts [K]"]
        Pt = self.get_station_data()["Pt [kPa]"]
        Ps = self.get_station_data()["Ps [kPa]"]
        stations = self.get_station_data()["Station"]
        marker = "o"
        size = 4
        color = "white"
        edge = 1.5
        width = 1
        plt.close("all")
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, axes = plt.subplots(2, 1, figsize=(7, 7), dpi=120, sharex=True)
        # Temperatures
        axes[0].plot(stations, Tt, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Total Temperature")
        axes[0].plot(stations, Ts, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Static Temperature")
        format_axes(axes[0], "Temperature", "Temperature [K]")
        # Pressures
        axes[1].plot(stations, Pt, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Total Pressure")
        axes[1].plot(stations, Ps, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Static Pressure")
        format_axes(axes[1], "Pressure", "Pressure [kPa]")
        axes[1].set_xlabel("Station", fontsize=11)
        plt.tight_layout()
        plt.show()


    def optimize(self, performance_parameter): pass
    def sensitivity_study(self): pass
    def off_design(self): pass


'''
handles turbofan (mixed/unmixed) and variable bypass ramjet
'''
class BypassEngine(Engine):
    def __init__(self, engine_parameters):
        self.spools = engine_parameters["spools"]
        self.B = engine_parameters["bypass ratio"]
        self.mixed = engine_parameters["mixed"]
        self.is_fan = engine_parameters["fan"]
        self.TET = engine_parameters["TET"]
        self.LHV = engine_parameters["LHV"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.ambient = Ambient(self.altitude, Minf=self.Minf)

        # Handle user errors in engine parameters
        if "mdot" in engine_parameters and "fan diameter" in engine_parameters:
            raise ValueError("Error. Cannot have two sizing parameters. Must choose either inlet mass flow or fan face diameter.")
        else:
            if "mdot" in engine_parameters: self.W = engine_parameters["mdot"] # kg/sec
            elif "fan diameter" in engine_parameters: self.fan_diameter = engine_parameters["fan diameter"] # inches


    def set_components(self, parameters):
        self.parameters = parameters

        # Include the engine object in each set of component parameters
        for component_name in parameters:
            if component_name != "engine":
                parameters[component_name]["engine"] = self

        # Handle spool counts
        match self.spools:
            case 1:
                self.PR = parameters["engine"]["PR"]
                parameters["compressor"]["PR"] = self.PR

                # Check for user errors
                if self.is_fan: raise ValueError("Cannot have a fan with only 1 spool. Must be at least two spools.")
                if self.mixed == False: raise ValueError("Must have a mixed exhaust with only 1 spool.")
                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    compressor := Compressor(inlet.root_exit, parameters["compressor"]),
                    burner := Burner(compressor.exit, parameters["burner"]),
                    turbine := Turbine(burner.exit, compressor, parameters["turbine"]),
                    mixer := Mixer(turbine.exit, inlet.bypass_exit),
                    exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                ]
                self.compressor = compressor
                self.turbine = turbine
            case 2:
                # Turbofan or Turbojet
                if self.is_fan: 
                    parameters["lpc"]["fan"] = True
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["lpc"]["root PR"]
                    self.HPR = self.OPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR

                    # Mixed Exhaust
                    match self.mixed:
                        case True: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.bypass_exit, parameters["lpc"], root_upstream=inlet.root_exit),
                                hpc := Compressor(lpc.root_exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                                mixer := Mixer(lpt.exit, lpc.tip_exit),
                                exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                            ]
                            self.mixer = mixer
                        case False: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.bypass_exit, parameters["lpc"], root_upstream=inlet.root_exit),
                                cold_nozzle := Nozzle(lpc.tip_exit, parameters["cold nozzle"]),
                                hpc := Compressor(lpc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                                exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                            ]
                            self.cold_nozzle = cold_nozzle
                else: 
                    parameters["lpc"]["fan"] = False
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["engine"]["LPR"]
                    self.HPR = self.OPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["lpc"]["PR"] = self.LPR

                    # Mixed Exhaust
                    if self.mixed == False: raise ValueError("Must have a mixed exhaust if there's no fan.")

                    self.components = [
                        inlet := Inlet(parameters["intake"]),
                        lpc := Compressor(inlet.root_exit, parameters["lpc"]),
                        hpc := Compressor(lpc.exit, parameters["hpc"]),
                        burner := Burner(hpc.exit, parameters["burner"]),
                        hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                        lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                        mixer := Mixer(lpt.exit, inlet.bypass_exit),
                        exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                    ]
                    self.mixer = mixer
                self.lpc = lpc
                self.hpc = hpc
                self.hpt = hpt
                self.lpt = lpt
            case 3:
                # Check if the engine is a turbofan
                if self.is_fan:
                    parameters["lpc"]["fan"] = True
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["lpc"]["root PR"]
                    self.HPR = parameters["engine"]["HPR"]
                    self.IPR = self.OPR / self.HPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["ipc"]["PR"] = self.IPR
                    # Mixed Exhaust
                    match self.mixed:
                        case True: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.exit, parameters["lpc"]),
                                ipc := Compressor(lpc.exit, parameters["hpc"]),
                                hpc := Compressor(ipc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                                lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                                mixer := Mixer(lpt.exit, lpc.tip_exit),
                                exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                            ]
                            self.mixer = mixer
                        case False:
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.exit, parameters["lpc"]),
                                cold_nozzle := Nozzle(lpc.tip_exit, parameters["cold nozzle"]),
                                ipc := Compressor(lpc.exit, parameters["hpc"]),
                                hpc := Compressor(ipc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                                lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                                exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                            ]
                            self.cold_nozzle = cold_nozzle
                else:
                    parameters["lpc"]["fan"] = False
                    self.OPR = parameters["engine"]["OPR"]
                    self.HPR = parameters["engine"]["HPR"]
                    self.IPR = parameters["engine"]["IPR"]
                    self.LPR = self.OPR / self.HPR / self.IPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["ipc"]["PR"] = self.IPR
                    parameters["lpc"]["PR"] = self.LPR
                    # Mixed Exhaust
                    if self.mixed == False: raise ValueError("Must have a mixed exhaust if there's no fan.")
                    self.components = [
                        inlet := Inlet(parameters["intake"]),
                        lpc := Compressor(inlet.root_exit, parameters["lpc"]),
                        ipc := Compressor(lpc.exit, parameters["hpc"]),
                        hpc := Compressor(ipc.exit, parameters["hpc"]),
                        burner := Burner(hpc.exit, parameters["burner"]),
                        hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                        ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                        lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                        mixer := Mixer(lpt.exit, inlet.bypass_exit),
                        exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                    ]
                    self.mixer = mixer
                self.lpc = lpc
                self.ipc = ipc
                self.hpc = hpc
                self.hpt = hpt
                self.ipt = ipt
                self.lpt = lpt
            case _:
                raise ValueError("Number of spools must be either 1, 2, or 3.")

        self.inlet = inlet
        self.burner = burner
        self.exhaust = exhaust
        # Check for an afterburner in the user input parameters
        if "afterburner" in parameters: self.add_afterburner(parameters["afterburner"])

    # Retrieve flow properties at every station
    def get_station_data(self): 
        raw_data = list()
        
        for component in self.components:
            if isinstance(component, Inlet):
                raw_data.append(component.freestream.get_properties())
                raw_data.append(component.inlet.get_properties())
                raw_data.append(component.bypass_exit.get_properties())
                raw_data.append(component.root_exit.get_properties())
            elif isinstance(component, Compressor):
                if hasattr(component, "is_fan"):
                    raw_data.append(component.root_exit.get_properties())
                    raw_data.insert(-2, component.tip_exit.get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
            else:
                raw_data.append(component.exit.get_properties())
                    
        rounded_data = numpy.round(numpy.array(raw_data, dtype=float), 3).tolist()
        station_data = pandas.DataFrame(rounded_data, columns=Station.column_names)
        return station_data

