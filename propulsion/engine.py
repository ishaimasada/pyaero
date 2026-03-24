"""
Module for Cycle Analysis and preliminary component design 
"""

import pandas
import numpy
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
                    "W [kg/sec]", 
                    "Wc [kg/sec]", 
                    "Tt [K]",
                    "Pt [Pa]",
                    "ht [kJ/kg]",
                    "M",
                    "Ts [K]",
                    "hs [kJ/kg]",
                    "Ps [Pa]",
                    "V [m/sec]",
                    "rho [kg/m^3]",
                    "Area [m^2]"
                   ]

    def __init__(self, W, Tt, Pt, M=None, idx=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.Wc = self.get_Wc(self.W, self.Tt, self.Pt)
        self.FAR = self.Wf = 0
        self.ht = self.get_ht(self.Tt, self.FAR)
        self.M = M
        self.idx = idx

        if self.M != None: self.set_statics(self.M)

    def __str__(self):
        return f"Station {self.idx}\nW = {self.W}\nWc = {self.Wc}\nCorrected W = {self.Wc}\nTt = {self.Tt}\nPt = {self.Pt}"

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
        return ht

        
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

    def get_station_data(self):
        if self.M != None: station_data = [self.W, self.Wc, self.Tt, self.Pt, self.ht, self.M, self.T, self.h, self.P, self.V, self.rho, self.area]
        else: station_data = [self.W, self.Wc, self.Tt, self.Pt, self.ht, None, None, None, None, None, None, None]

        return station_data


class Inlet:
    def __init__(self, freestream:Ambient, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        D1 = cycle_parameters["front face diameter"]
        M1 = cycle_parameters["M1"]
        [_, _, _, rhot_rho, _] = isentropic(freestream.Minf, freestream.gamma, lookup_key="M")
        A1 = numpy.pi * (D1 / 2)**2
        rho1 = freestream.rhot * (1/rhot_rho)
        W1 = freestream.Vinf * rho1 * A1
        self.upstream = Station(W1, freestream.Tt, freestream.Pt, idx=1) # Station 1
        self.upstream.set_statics(M1)
        self.freestream = Station(W1, freestream.Tt, freestream.Pt, idx=0) # Station 0
        self.freestream.set_statics(freestream.Minf)
        
        self.exit_W = self.upstream.W
        self.exit_Pt = self.upstream.Pt
        self.exit_Tt = self.upstream.Tt
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt, 2) # Must be split into tip and root if going into a fan

        # COMPONENT ANALYSIS
        if component_parameters != None: pass


class Compressor:
    def __init__(self, upstream:Station, cycle_parameters):
        # CYCLE ANALYSIS
        PR = cycle_parameters["PR"]
        hp_bleed = cycle_parameters["bleed"]
        e_t = cycle_parameters["efficiency"]

        self.upstream = upstream
        self.coolant = hp_bleed.cooling * self.upstream.W
        self.exit_W = self.upstream.W - self.coolant
        self.exit_Pt = self.upstream.Pt * PR
        self.exit_Tt = self.upstream.Tt * (PR)**((self.upstream.gamma - 1)*e_t / self.upstream.gamma)
        self.downstream = Station(self.upstream.W, self.exit_Tt, self.exit_Pt)
        self.delta_h =  self.upstream.cp * self.upstream.Tt - (self.downstream.cp * self.downstream.Tt)
        self.power = self.upstream.W * self.delta_h


class Burner:
    def __init__(self, upstream:Station, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        Ptloss_b = cycle_parameters["burner total pressure loss"]
        TET = cycle_parameters["TET"]
        LHV = cycle_parameters["LHV"]
        self.eta_b = cycle_parameters["burner efficiency"]

        self.upstream = upstream
        self.exit_FAR = self.burn(TET, self.upstream.Tt, LHV)
        self.exit_W = self.upstream.W * (1 + self.exit_FAR)
        self.exit_Pt = upstream.Pt * (1 - Ptloss_b)
        self.exit_Tt = TET
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)

        # COMPONENT DESIGN
        if component_parameters != None: pass

    def burn(self, T2, T1, LHV):
        FARnew = 0.02
        FAR = -1 
        h1 = self.upstream.get_ht(T1, 0)
        tolerance = 0.00001
        error = (abs(FAR - FARnew) / FARnew) 
        
        while error > tolerance:
            FAR = FARnew
            h2 = self.upstream.get_ht(T2, FAR)
            FARnew = (h2 - h1) / (LHV * self.eta_b)
            error = (abs(FAR - FARnew) / FARnew) 

        return FARnew
        

class Turbine:
    def __init__(self, upstream, cycle_parameters, component_parameters=None):
        # CYCLE ANALYSIS
        TET = cycle_parameters["TET"]
        e_t = cycle_parameters["polytropic efficiency"]
        mdot_cool = cycle_parameters["compressor"].coolant
        power = cycle_parameters["compressor"].power 

        self.upstream = upstream
        self.exit_W = self.upstream.W + mdot_cool
        self.exit_ht = (self.upstream.ht * self.upstream.W - power) / self.upstream.W
        FAR = self.upstream.W * self.upstream.FAR / (self.exit_W - (self.upstream.W * self.upstream.FAR))
        self.exit_Tt = self.upstream.T_from_H(self.exit_ht/1000, FAR, TET, 0)
        ER = (self.exit_Tt  / self.upstream.Tt)**(-self.upstream.gamma / ((self.upstream.gamma - 1)*e_t))
        self.exit_Pt = self.upstream.Pt / ER
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)

        # COMPONENT DESIGN
        if component_parameters != None:
            self.phi = component_parameters["load coefficient"]
            self.psi = component_parameters["work coefficient"]
            self.R = component_parameters["radius"]
            self.rpm = component_parameters["rpm"]
            self.Vum3 = component_parameters["tangential velocity at mid radius station 3"]
            self.Rm3 = component_parameters["mid radius at station 3"]
            self.AR = component_parameters[" aspect ratio"]
    
            omega = self.rpm * (2*numpy.pi / 60)
            self.Vax3 = component_parameters["axial velocity"]

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


class Nozzle:
    def __init__(self, upstream:Station):
        self.upstream = upstream
        self.downstream = Station(self.upstream.W, self.upstream.Pt, self.upstream.Tt, self.upstream.ht)

        def statics(self,  Pinf):
            gamma = self.downstream.gamma
            cp = self.downstream.cp
            R = self.downstream.R
            critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
            NPR = self.upstream.Pt / Pinf
            if NPR > critical_NPR:
                self.downstream.M = 1
                self.downstream.V = numpy.sqrt(gamma * R * self.downstream.T)
                self.downstream.P = self.downstream.Pt * critical_NPR
            else:
                self.downstream.P = Pinf
                self.downstream.T = self.downstream * (self.downstream.P/self.downstream.Pt)**((gamma - 1) / gamma)
                self.downstream.V = numpy.sqrt(2 * cp * (self.downstream.Tt - self.downstream.T))
                self.downstream.M = self.downstream.V / numpy.sqrt(gamma * R * self.downstream.T)

class Diffuser:
    def __init__(self, upstream:Station):
        self.upstream = upstream
        self.downstream = Station(self.upstream.W, self.upstream.Pt, self.upstream.Tt, self.upstream.ht)

        def set_statics(self,  Pa):
            pass

class Bleed:
    def __init__(self, cooling, packing):
        self.cooling = cooling
        self.packing = packing

class Duct:
    def init(self):
        pass

class Engine:
    def __init__(self, engine_parameters):
        # Design Parameters (metric units)
        self.OPR = engine_parameters["OPR"]
        self.TET = engine_parameters["TET"]
        self.BPR = engine_parameters["BPR"]
        self.FPR = engine_parameters["FPR"]
        self.LHV = engine_parameters["LHV"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.front_diameter = engine_parameters["front face diameter"] # inches

        # Component Efficiencies
        self.Ptrec_ft = engine_parameters["fan tip total pressure recovery"] 
        self.Ptrec_fh = engine_parameters["fan hub total pressure recovery"]
        self.Ptloss_OGV = engine_parameters["OGV total pressure loss"]
        self.Ptloss_bypass = engine_parameters["bypass duct total pressure loss"]
        self.e_ft = engine_parameters["fan tip polytropic efficiency"]
        self.e_fh = engine_parameters["fan hub polytropic efficiency"]
        self.swan_loss = engine_parameters["swan neck loss"]
        self.e_core = engine_parameters["core comp. polytropic efficiency"]
        self.eta_b = engine_parameters["burner efficiency"]
        self.Ptloss_b = engine_parameters["burner total pressure loss"]
        self.e_t = engine_parameters["turbine polytropic efficiency"]
        self.CD = engine_parameters["Nozzle Discharge Coefficienct"]
        self.CV =engine_parameters["Nozzle Velocity Coefficienct"]
        self.eta_m = engine_parameters["shaft mechanical efficiency"]

        # Secondary Air System
        self.hp_bleed = Bleed(engine_parameters["high pressure cooling"], engine_parameters["high pressure packing"])
        self.lp_bleed = Bleed(engine_parameters["low pressure cooling"], engine_parameters["low pressure packing"])

        # Mach Numbers
        self.M1 =engine_parameters["M1"]
        self.M12 = engine_parameters["M1.2"]
        self.M2 = engine_parameters["M2"]
        self.M31 = engine_parameters["M3.1"]
        self.M5 = engine_parameters["M5"]

        # Prepare cycle analysis and component design inputs
        intake_cycle_parameters = {"front face diameter": self.front_diameter, "M1": self.M1}
        burner_cycle_parameters = {"TET": self.TET, "LHV": self.LHV, "burner total pressure loss": self.Ptloss_b, "burner efficiency": self.eta_b}
        compressor_cycle_parameters = {"PR": self.OPR, "bleed": self.hp_bleed, "efficiency": self.e_core}

        # Architecture
        self.ambient = Ambient(self.altitude, Minf=self.Minf)
        self.inlet = Inlet(self.ambient, intake_cycle_parameters)
        self.compressor = Compressor(self.inlet.downstream, compressor_cycle_parameters)
        self.burner = Burner(self.compressor.downstream, burner_cycle_parameters)

        turbine_cycle_parameters = {"TET": self.TET, "polytropic efficiency": self.e_t, "compressor": self.compressor}

        self.turbine =  Turbine(self.burner.downstream, turbine_cycle_parameters)
        self.exhaust =  Nozzle(self.turbine.downstream)
        
        self.components = [self.inlet, self.compressor, self.burner, self.turbine, self.exhaust]

    '''
    @property
    def components(self): return self.validate_cycle(engine_parameters)

    @components.setter
    def components(self, components): self._components = components

    def validate_cycle(self, engine_parameters): return self.cp / (self.cp - self.R)
    self.components = [self.ambient, self.inlet, self.compressor, self.burner, self.turbine, self.exhaust]
    '''


    def get_station_data(self): 
        raw_data = list()

        for component in self.components:
            if isinstance(component, Inlet): raw_data.append(component.upstream.get_station_data())
            else: raw_data.append(component.downstream.get_station_data())

        station_data = pandas.DataFrame(raw_data, columns=Station.column_names)

        return station_data


    def display_performance(self):
        """ Performance Parameters """
        '''
        # Specific Thrust
        T_ma = u_e * (1 + f44) - u_i
        # Thrust Specific Fuel Consumption (TSFC)
        TSFC = (f44 / T_ma) * 10**3
        # Thrust
        thrust = T_ma * mdota

        # Propulsive Efficiency
        eta_p = 2 / (1 + (u_e / u_i))
        # Thermal Efficiency
        eta_th = (((1 + f44) * u_e ** 2) - u_i**2) / (f44 * QR)
        # Overall Efficiency
        eta_o = eta_p * eta_th
        '''
        pass

    def optimize(self, performance_parameter):
        pass
    def sensitivity_study(self):
        pass
    def off_design(self):
        pass