"""
Module for Cycle Analysis and preliminary component design 
"""

import pandas
import numpy
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
                    "Pt [Pa]",
                    "ht [kJ/kg]",
                    "Wc [kg/sec]", 
                    "Wf [kg/sec]", 
                    "FAR",
                    "M",
                    "Ts [K]",
                    "hs [kJ/kg]",
                    "Ps [Pa]",
                    "V [m/sec]",
                    "rho [kg/m^3]",
                    "Area [m^2]"
                   ]

    def __init__(self, W, Wf, Tt, Pt, FAR=None, M=None, idx=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.Wc = self.get_Wc(self.W, self.Tt, self.Pt)
        self.Wf = 0
        self.Wa = self.W - self.Wf
        self.FAR = FAR
        self.ht = self.get_ht(self.Tt, self.FAR)
        self.M = M
        self.idx = idx

        if self.M != None: self.set_statics(self.M)
        if self.FAR != None: self.FAR = self.get_FAR()

    def __str__(self):
        return f"Station {self.idx}\nW = {self.W}\nWc = {self.Wc}\nCorrected W = {self.Wc}\nTt = {self.Tt}\nPt = {self.Pt}\nM = {self.M}"

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
        if self.M != None: station_data = [self.idx, self.W, self.Tt, self.Pt, self.ht/1000, self.Wc, self.Wf, self.FAR, self.M, self.T, self.h/1000, self.P, self.V, self.rho, self.area]
        else: station_data = [self.idx, self.W, self.Tt, self.Pt, self.ht/1000, self.Wc, self.Wf, self.FAR, None, None, None, None, None, None, None]
        return station_data


class Component:
    def __init__(self, *kwargs):
        self.attributes = kwargs

class Inlet:
    def __init__(self, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        M_inlet = cycle_parameters["M_inlet"]
        M_exit = cycle_parameters["M_exit"]
        Pt_recovery = cycle_parameters["total pressure recovery"]
        engine = cycle_parameters["engine"]
        freestream = engine.ambient
        D1 = engine.front_diameter * 2.54 / 100 # in to m
        [_, _, _, rhot_rho, _] = isentropic(freestream.Minf, freestream.gamma, lookup_key="M")
        A1 = numpy.pi * (D1 / 2)**2
        rho1 = freestream.rhot * (1/rhot_rho)
        W1 = freestream.Vinf * rho1 * A1
        freestream_FAR = freestream_Wf = inlet_FAR = inlet_Wf = exit_FAR = exit_Wf = 0

        self.freestream = Station(W1, freestream_Wf, freestream.Tt, freestream.Pt, freestream_FAR, M=freestream.Minf, idx=0) # Station 0
        self.inlet = Station(W1, inlet_Wf, freestream.Tt, freestream.Pt, inlet_FAR, M=M_inlet, idx=1) # Station 1
        self.exit = copy.deepcopy(self.inlet)
        self.exit.M = M_exit
        self.exit.Pt = self.inlet.Pt * Pt_recovery

        # COMPONENT ANALYSIS
        if component_parameters != None: pass


class Compressor:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        e_c = cycle_parameters["e"]
        engine = cycle_parameters["engine"]
        inlet_idx = cycle_parameters["idx_inlet"]
        exit_idx = cycle_parameters["idx_exit"]
        M_exit = cycle_parameters["M_exit"]
        PR = engine.PR
        bleed = engine.bleed

        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = exit_idx
        self.exit.M = M_exit
        self.inlet.idx = inlet_idx
        self.coolant = bleed.cooling * self.inlet.W
        self.exit.W = self.inlet.W - self.coolant
        self.exit.Pt = self.inlet.Pt * PR
        self.exit.Tt = self.inlet.Tt * (PR)**((self.inlet.gamma - 1)*e_c / self.inlet.gamma)
        self.delta_h =  self.inlet.cp * self.inlet.Tt - (self.exit.cp * self.exit.Tt)
        self.power = self.inlet.W * self.delta_h


class Burner:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        Ptloss_b = cycle_parameters["total pressure loss"]
        eta_b = cycle_parameters["efficiency"]
        exit_idx = cycle_parameters["idx_exit"]
        M_exit = cycle_parameters["M_exit"]
        engine = cycle_parameters["engine"]
        LHV = engine.LHV
        TET = engine.TET

        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = exit_idx
        self.exit.M = M_exit
        FAR = self.get_FAR(TET, self.inlet.Tt, 0, LHV, eta_b)
        self.exit.Wf = self.inlet.W * FAR
        self.exit.W = self.inlet.W * (1 + self.exit.FAR)
        self.exit.Pt = self.inlet.Pt * (1 - Ptloss_b)
        self.exit.Tt = TET

        # COMPONENT DESIGN
        if component_parameters != None: pass

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
    def __init__(self, upstream, compressor, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        e_t = cycle_parameters["e"]
        mdot_cool = compressor.coolant
        power = abs(compressor.power)
        eta_m = cycle_parameters["mechanical efficiency"]
        exit_idx = cycle_parameters["idx_exit"]
        M_exit = cycle_parameters["M_exit"]

        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = exit_idx
        self.exit.M = M_exit
        self.exit.W = self.inlet.W + mdot_cool
        exit_ht = (self.inlet.ht * self.inlet.W - (power/eta_m)) / self.exit.W
        FAR = self.inlet.W * self.inlet.FAR / (self.exit.W - (self.inlet.W * self.inlet.FAR))
        self.exit.Wf = (self.inlet.W - self.inlet.Wf) * FAR
        self.exit.Tt = self.inlet.T_from_H(exit_ht, self.exit.FAR, self.inlet.Tt, 100)
        ER = (self.exit.Tt  / self.inlet.Tt)**(-self.inlet.gamma / ((self.inlet.gamma - 1)*e_t))
        self.exit.Pt = self.inlet.Pt / ER

        # COMPONENT DESIGN
        if component_parameters != None:
            self.phi = component_parameters["load coefficient"]
            self.psi = component_parameters["work coefficient"]
            self.R = component_parameters["radius"]
            self.rpm = component_parameters["rpm"]
            self.Vum3 = component_parameters["tangential velocity at mid radius station 3"]
            self.Rm3 = component_parameters["mid radius at station 3"]
            self.AR = component_parameters[" aspect ratio"]
            self.Vax3 = component_parameters["axial velocity"]

            omega = self.rpm * (2*numpy.pi / 60)

    def solve_stage(self, R: list, omega, Vax, Vu_mid, Rm3):
        # Aerodynamics
        U = omega * numpy.array(R)
        Vu = Vu_mid * (R / Rm3) # Free vortex assumption
        Wu = Vu - U
        V = numpy.sqrt(Vax**2 + Vu**2)
        W = numpy.sqrt(Vax**2 + Wu**2)
        alpha = numpy.atan(Vu / Vax)
        beta = numpy.atan(Wu / Vax)
        reaction = (W[2]**2 - W[1]**2) / (V[2]**2 - V[1]**2 + W[2]**2 - W[1]**2)

        # Thermodynamic Properties
        delta_ht = self.psi * U[2]**2
        Tt1 = Tt2 = self.upstream.Tt
        Tt3 = self.upstream.T_from_H(ht[2], self.upstream.FAR, self.upstream.Tt, 100)
        ER = (self.exit_Tt  / self.upstream.Tt)**(-self.upstream.gamma / ((self.upstream.gamma - 1)*e_tt))
        Pt1 = Pt2 = self.upstream.Pt
        Pt3 = Pt2 / ER
        mdot = numpy.array([self.upstream.W, self.upstream.W, self.upstream.W + self.mdot_cool])
        Tt = numpy.array([Tt1, Tt2, Tt3])
        Pt = numpy.array([Pt1, Pt2, Pt3])
        mdot_corrected = numpy.vectorize(self.upstream.get_Wc)(mdot, Tt, Pt)
        ht = numpy.array([self.upstream.ht, self.upstream.ht, self.upstream.ht - delta_ht])
        T = ht - V**2 / self.upstream.cp
        P = Pt * (T / Tt)**(self.upstream.gamma / ((self.upstream.gamma - 1)*self.e_tt))
        rho = P / (self.upstream.R * T)
        a = numpy.sqrt(self.upstream.R * self.upstream.gamma * T)
        M_absolute = V / a
        M_relative = W / a

    def blade_geometry(self):
        pass
        '''
        sections = solve_section()
        # Geometry
        hub = sections[0]
        tip = sections[-1]
        h = (sections.R[-1][1] + Rt[0])/2 - (Rh[2] + Rh[1])/2
        c = h / AR
        stagger_angle = (alpha + alpha) / 2
        cax_t = c * numpy.cos(stagger_angle)
        cax_t = c * numpy.cos(stagger_angle)
        taper_ratio = 
        '''


class Mixer:
    def __init__(self, hot_inlet:Station, cold_inlet:Station, B): 
        # Mhot is a design choice that can alternatively be cold_inlet.M

        if B == 0: self.exit = copy.deepcopy(hot_inlet)
        
        # Hot Stream (static properties)
        hot_inlet.T = hot_inlet.Tt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)^(-1)
        hot_inlet.P = hot_inlet.Pt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)^(-hot_inlet.gamma / (hot_inlet.gamma - 1))
        hot_inlet.rho = hot_inlet.P / (hot_inlet.R*hot_inlet.T)
        hot_inlet.V = hot_inlet.M * numpy.sqrt(hot_inlet.gamma * hot_inlet.R * hot_inlet.T)
        hot_inlet.A = hot_inlet.W / (hot_inlet.rho*hot_inlet.V)
        self.hot_inlet = hot_inlet
        momentum_hot = hot_inlet.W*hot_inlet.V + hot_inlet.A*hot_inlet.P

        # Cold Stream (static properties)
        cold_inlet.P = hot_inlet.P
        cold_inlet.M = numpy.sqrt(2 / (cold_inlet.gamma - 1) * ((cold_inlet.Pt/cold_inlet.P)**((cold_inlet.gamma - 1)/cold_inlet.gamma) - 1))
        cold_inlet.T = cold_inlet.Tt * (1 + (cold_inlet.gamma - 1)/2 * cold_inlet.M**2)**(-1)
        cold_inlet.rho = cold_inlet.P / (cold_inlet.R*cold_inlet.T)
        cold_inlet.V = cold_inlet.M * numpy.sqrt(cold_inlet.gamma * cold_inlet.R * cold_inlet.T)
        cold_inlet.A = cold_inlet.W / (cold_inlet.rho*cold_inlet.V)
        cold_inlet.set_statics(cold_inlet.M)
        self.cold_inlet = cold_inlet
        momentum_cold = cold_inlet.W*cold_inlet.V + cold_inlet.A*cold_inlet.P

        if cold_inlet.P > cold_inlet.Pt:
            raise ValueError("Static pressure is larger than the total pressure (not physically possible). A larger Mhot must be chosen.")

        # Mixer Exit Stream
        exit_idx = "6"
        exit_W = hot_inlet.W + cold_inlet.W
        exit_Wf = hot_inlet.Wf
        exit_A = hot_inlet.A + cold_inlet.A
        exit_Tt = (hot_inlet.W*hot_inlet.Tt + cold_inlet.W*cold_inlet.Tt) / (exit_station.W)
        exit_FAR = exit_Wf / (exit_station.W - exit_Wf)
        exit_R = Station.get_R(exit_FAR)
        exit_Cp = Station.get_cp(exit_Tt, exit_FAR)
        exit_gamma = exit_Cp / (exit_Cp - exit_R)
        K = (momentum_hot + momentum_cold)^2 * (exit_gamma / (exit_W^2 * exit_Tt * exit_R))
        exit_M = numpy.sqrt(((2*exit_gamma - K) + numpy.sqrt((K - 2*exit_gamma)^2 + 4*(((exit_gamma-1)/2)*K-exit_gamma^2))) / (2*((exit_gamma-1)/2*K-exit_gamma^2)))
        exit_T = exit_Tt * (1 + (exit_gamma - 1)/2 * exit_M^2)^(-1)
        exit_P = (exit_W*numpy.sqrt(exit_R*exit_T)/(exit_A*exit_M*numpy.sqrt(exit_gamma)))
        exit_Pt = exit_P * (exit_Tt/exit_T)^(exit_gamma / (exit_gamma - 1))
        exit_station = Station(exit_W, exit_Wf, exit_Tt, exit_Pt, exit_FAR, idx=exit_idx, M=exit_M)
        exit_station.set_statics(exit_station.M)


class Afterburner:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        engine = cycle_parameters["engine"]
        toggle = cycle_parameters["toggle"]
        Ttmax = cycle_parameters["AET"]
        pi_ab_on = cycle_parameters["pi off"]
        pi_ab_off = cycle_parameters["pi on"]
        eta = cycle_parameters["efficiency"]
        LHV = engine.LHV

        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = "7"

        # Handle different engine modes (afterburner on or off)
        if toggle == 0: self.exit.Pt = self.inlet.Pt * pi_ab_off # Afterburner off
        elif toggle == 1: # Afterburner on
            self.exit.Pt = self.inlet.Pt * pi_ab_on
            Tt_Ttstar1 = get_Tt_Ttstar(self.inlet.M, self.inlet.gamma)
            Tt_Ttstar2 = Ttmax * (Tt_Ttstar1) / self.inlet.Tt
            Ttstar = 1 / Tt_Ttstar1 * self.inlet.Tt

            # Check if the afteburner chokes before reaching the given total
            # temperature (based on Rayleigh flow theory)
            if Tt_Ttstar2 > 1:
                print(f"Given Afterburner Exit Temperature exceeds maximum allowed amount {Ttstar}.\n"
                       "This value has been set as the exit total temperatue and a lower value is recommended to be chosen.\n\n")
                self.exit.Tt = Ttstar
                self.exit.M = 1
            elif Tt_Ttstar2 <= 1:
                self.exit.Tt = Ttmax
                self.exit.M = bisection(get_Tt_Ttstar, Tt_Ttstar2, 1, 0, self.inlet.gamma)

            self.exit.FAR = Burner.get_FAR(self, self.exit.Tt, self.inlet.Tt, self.inlet.FAR, LHV, eta) 
            self.exit.Wf = (self.inlet.W - self.inlet.Wf) * self.exit.FAR
            self.exit.W = self.inlet.W + self.exit.Wf
            self.exit.set_statics(self.exit.M)
        elif toggle > 1: raise ValueError("Error: bypass factor cannot be greater than 1 (0<B<1).\n")
        elif toggle < 0: raise ValueError("Error: bypass factor cannot be negative (0<B<1).\n")


class Nozzle:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        geometry = cycle_parameters["C/CD"]
        engine = cycle_parameters["engine"]
        Pinf = engine.ambient.P
        self.inlet = upstream

        # Handles Converging or Converging-Diverging Nozzles
        match geometry:
            case "C":
                self.exit = self.converging(Pinf)
            case "CD":
                self.exit, self.throat = self.CD(Pinf)
        
        # COMPONENT DESIGN
        if component_parameters != None: pass

    def converging(self,  Pinf):
        gamma = self.inlet.gamma
        cp = self.inlet.cp
        R = self.inlet.R
        critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
        NPR = self.inlet.Pt / Pinf
        exit_station = copy.deepcopy(self.inlet)
        exit_station.idx = 8 # Throat station

        if NPR >= critical_NPR:
            # Choked
            exit_station.M = 1
            exit_station.T = (2 / (gamma + 1)) * exit_station.Tt
            exit_station.V = numpy.sqrt(gamma * R * exit_station.T)
            exit_station.P = exit_station.Pt * (1 / critical_NPR)
        else:
            # Unchoked
            exit_station.P = Pinf
            exit_station.T = exit_station.Tt * (exit_station.P/exit_station.Pt)**((gamma - 1) / gamma)
            exit_station.V = numpy.sqrt(2 * cp * (exit_station.Tt - exit_station.T))
            exit_station.M = exit_station.V / numpy.sqrt(gamma * R * exit_station.T)

        exit_station.set_statics(exit_station.M)

        return exit_station
    
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
        #self.downstream.mdot_fuel = self.throat.mdot_fuel
        self.exit.P = Pinf; # Assuming perfectly expanded flow
        self.exit.M = numpy.sqrt((2/(gamma - 1)) * ((self.exit.Pt/self.exit.P)^((gamma-1)/gamma) - 1)); # Isentropic Relation
        self.exit.set_statics(self.exit.M)


class Recuperator:
    def __init__(self, cold_inlet:Station, cycle_parameters, component_parameters=None):
        self.cyce_parameters = cycle_parameters
        self.engine = cycle_parameters["engine"]
        self.Mexit_cold = cycle_parameters["Mexit_cold"]
        self.Mexit_hot = cycle_parameters["Mexit_hot"]
        self.pi_cold = cycle_parameters["pi_cold"]
        self.pi_hot = cycle_parameters["pi_hot"]
        self.delta_ht = cycle_parameters["delta_ht"] * 10**3 # kJ/kg to J/kg

        self.cold_inlet = copy.deepcopy(cold_inlet)
        self.cold_inlet.idx = "3.06"
        self.cold_exit = copy.deepcopy(self.cold_inlet)
        self.cold_exit.idx = "3.07"
        self.cold_exit.Tt = self.cold_exit.T_from_H(self.cold_inlet.ht+self.delta_ht, self.cold_exit.FAR, self.engine.TET, self.cold_inlet.Tt)
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

class Engine:
    '''
    Essentially a turbojet core engine will serve as a parent class to other architectures like turbofan, turboshaft, ramjet, etc.
    Handles a recuperator and afterburner as well, but the default is just a core
    '''
    def __init__(self, engine_parameters, parameters=None):
        # Design Parameters (metric units)
        self.PR = engine_parameters["PR"]
        self.TET = engine_parameters["TET"]
        self.LHV = engine_parameters["LHV"]
        self.bleed = engine_parameters["bleed"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.front_diameter = engine_parameters["front face diameter"] # inches

        # Ambient/Freestream
        self.ambient = Ambient(self.altitude, Minf=self.Minf)

        # Architecture (Turbojet by Default)
        self.set_components(parameters)


    @property
    def components(self): return self._components

    @components.setter
    def components(self, value): self._components = value

    def set_components(self, parameters):
        # Defualt turbojet design parameters
        if parameters == None:
            intake_parameters = {"engine": self, "total pressure recovery": 0.98, "M_inlet": 0.6, "M_exit": 0.55}
            compressor_parameters = {"engine": self, "e": 0.91, "M_exit": 0.3, "idx_inlet": 2, "idx_exit": 3}
            burner_parameters = {"engine": self, "total pressure loss": 0.05, "efficiency": 0.99, "M_exit": 0.4, "idx_exit": 4}
            turbine_parameters = {"e": 0.8, "mechanical efficiency": 0.99, "M_exit": 0.5, "idx_exit": 5}
            nozzle_parameters = {"engine": self, "C/CD": "C", "discharge coefficient": 0.97, "velocity coefficient": 0.98}
        else:
            intake_parameters = parameters["inlet"]
            compressor_parameters = parameters["compressor"]
            burner_parameters = parameters["burner"]
            turbine_parameters = parameters["turbine"]
            nozzle_parameters = parameters["nozzle"]
        
        self.components = [
                        inlet := Inlet(intake_parameters),
                        compressor := Compressor(inlet.exit, compressor_parameters),
                        burner := Burner(compressor.exit, burner_parameters),
                        turbine := Turbine(burner.exit, compressor, turbine_parameters),
                        exhaust := Nozzle(turbine.exit, nozzle_parameters)
        ]

        self.inlet = inlet
        self.compressor = compressor
        self.burner = burner
        self.turbine = turbine
        self.exhaust = exhaust


    def toggle_afterburner(self, parameters=None, toggle=False):
        match toggle:
            case False:
                for idx, component in enumerate(self.components):
                    if isinstance(component, Afterburner):
                        delattr(self, "afterburner")
                        self.components.remove(idx)
                self.afterburner = False
            case True:
                if parameters == None: ValueError("No afterburner parameters given. Must provide component parameters.\n")
                if any(isinstance(component, Afterburner) for component in self.components): 
                    return
                else:
                    for idx, component in enumerate(self.components):
                        if isinstance(component, Nozzle):
                            insert_idx = idx
                    self.components.insert(insert_idx, afterburner := Afterburner(self.components[insert_idx-1].exit, parameters)) 
                self.afterburner = afterburner


    def toggle_recuperator(self, parameters=None, toggle=False):
        match toggle:
            case False:
                for idx, component in enumerate(self.components):
                    if isinstance(component, Recuperator):
                        delattr(self, "recuperator")
                        self.components.remove(idx)
                self.recuperator = False
            case True:
                if parameters == None: ValueError("No recuperator parameters given. Must provide component parameters.\n")
                recuperator = Recuperator(self.compressor.exit, parameters)
                if any(isinstance(component, Recuperator) for component in self.components): 
                    return
                else:
                    for idx, component in enumerate(self.components):
                        if isinstance(component, Burner): 
                            insert_idx = idx
                    self.components.insert(insert_idx, recuperator)
                self.recuperator = recuperator


    def get_station_data(self): 
        raw_data = list()

        # Handle Recuperator station data (afterburner is treated like the other components)
        if hasattr(self, "recuperator"):
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


    def get_performance(self):
        """ Performance Parameters """
        station_data = self.get_station_data()

        W = station_data["W [kg/sec]"].values
        V = station_data["V [m/sec]"].values
        FAR = station_data["FAR"].values

        T_ma = V[-1]*(1 + FAR[-1]) - V[0] # Specific Thrust
        TSFC = (FAR[-1] / T_ma) * 10**3 # Thrust Specific Fuel Consumption (TSFC)
        thrust = T_ma * W[0] # Thrust
        eta_p = 2 / (1 + (V[-1]/V[0])) # Propulsive Efficiency
        eta_th = (((1 + FAR[-1])*V[-1]**2) - V[0]**2) / (FAR[-1]*self.LHV) # Thermal Efficiency
        eta_o = eta_p * eta_th # Overall Efficiency

        performance = pandas.DataFrame({
                                         "Specific Thrust [m/sec]": T_ma,
                                         "TSFC [g/sec/N]": TSFC,
                                         "Thrust [N]": thrust,
                                         "Propulsive Efficiency": eta_p,
                                         "Thermal Efficiency": eta_th,
                                         "Overall Efficiency": eta_o
                                        }, index=[0])

        return performance


    def optimize(self, performance_parameter):
        pass
    def sensitivity_study(self):
        pass
    def off_design(self):
        pass


class BypassEngine(Engine): pass
# handles turbofan (mixed/unmixed) and variable bypass ramjet