### Module for shape design ###
import math
import matplotlib.pyplot as plt
import copy
import numpy


# This is a point, but the truer representation of the object is a position vector. Will manipulate it as such.
class Point:
    def __repr__(self):
            return f"Point: {round(float(self.x_coord), 5)}, {round(float(self.y_coord), 5)}, {round(float(self.z_coord), 5)}"

    def __init__(self, x, y, z):
        self.x_coord = x
        self.y_coord = y
        self.z_coord = z

    def __add__(self, other_point):
        return Point(self.x_coord + other_point.x_coord, self.y_coord + other_point.y_coord, self.z_coord + other_point.z_coord)

    def __sub__(self, other_point):
        return Point(self.x_coord - other_point.x_coord, self.y_coord - other_point.y_coord, self.z_coord - other_point.z_coord)

    def __mul__(self, num):
        return Point(self.x_coord * num, self.y_coord * num, self.z_coord * num)

    def __rmul__(self, num):
        return self.__mul__(num)
    
    def __truediv__(self, scalar):
        return Point(self.x_coord / scalar, self.y_coord / scalar, self.z_coord / scalar)
    
    def make_copy(self):
        return Point(copy.deepcopy(self.x_coord), copy.deepcopy(self.y_coord), copy.deepcopy(self.z_coord))

    def set_x(self, new_x):
        self.x_coord = new_x

    def set_y(self, new_y):
        self.y_coord = new_y

    def set_z(self, new_z):
        self.z_coord = new_z
    
    @property
    def magnitude(self):
        return numpy.sqrt(self.x_coord**2 + self.y_coord**2 + self.z_coord**2)


class Panel:
    def __init__(self, start_point:Point, end_point:Point):
        self.start_point = copy.deepcopy(start_point)
        self.end_point = copy.deepcopy(end_point)

    def get_position(self, parameter, data="point"):
        v = self.end_point - self.start_point
        position = self.start_point + parameter * v
        match data:
            case "point":
                return position
            case "coordinates":
                return position.x_coord, position.y_coord
    
    @property
    def length(self):
        return numpy.sqrt((self.end_point.x_coord - self.start_point.x_coord)**2 + (self.end_point.y_coord - self.start_point.y_coord)**2 + (self.end_point.z_coord - self.start_point.z_coord)**2)
    
    @property
    def tangent_vector(self):
        return (self.end_point - self.start_point) / self.length

    @property
    def normal_vector(self):
        pass

    @property
    def phi(self):
        if self.end_point.x_coord - self.start_point.x_coord != 0: 
            phi = numpy.atan2((self.end_point.y_coord - self.start_point.y_coord), (self.end_point.x_coord - self.start_point.x_coord)) 
        else: phi = 0
        return phi
    
    @property
    def beta(self):
        return self.phi + (numpy.pi / 2)


# Bezier Curve class for detailed curve design
class Bezier:
    def __init__(self, control_points, resolution):
        self.control_points = control_points
        self.resolution = resolution
        self.degree = len(self.control_points) - 1

    @property
    def positions(self): return self.get_positions()

    def get_positions(self, control_points=None):
        # Returns positions based on the parameters given 
        def basis_polynomial(parameter, iterator):
            # Returns the bernstein polynomial coefficient
            binomial_coefficient = math.factorial(self.degree) / (math.factorial(iterator) * math.factorial(self.degree - iterator))
            return binomial_coefficient * (parameter**iterator) * ((1 - parameter)**(self.degree - iterator))
        if control_points is not None:
            self.control_points = control_points
            self.degree = len(self.control_points) - 1
        positions = []

        # Iterate Through Each Parameter Step
        for t in numpy.linspace(0, 1, self.resolution):
            position = Point(0, 0, 0)
            # Apply the effects of each control point to the parameter
            for idx, point in enumerate(self.control_points):
                position += point * basis_polynomial(t, idx)
            # Store the position
            positions.append(position)
        return positions

    # Plots the positions on the curve 
    def plot_points(self, cp=False):
        if len(self.positions) == 0:
            print("No positions have been generated. Parameter(s) must be provided")
            return
        x_positions = numpy.array([point.x_coord for point in self.positions])
        y_positions = numpy.array([point.y_coord for point in self.positions])
        x_control = numpy.array([point.x_coord for point in self.control_points])
        y_control = numpy.array([point.y_coord for point in self.control_points])
        match cp:
            case True:
                plt.scatter(x_control, y_control, marker="o")
                plt.plot(x_control, y_control, color="red")
        plt.plot(x_positions, y_positions)


# Create a BSpline curve using the Cox-de Boor recursion formula based on basis coefficients
# NOTE: This program currently assumes a uniform knot vector for periodic curves.
class BSpline:
    def __init__(self, control_points, degree, num_points):
        self.control_points = control_points
        self.degree = degree
        self.positions = []

        # Uniform knot vector (Periodic Curve)
        num_knots = degree + len(self.control_points) + 1
        self.knots = numpy.linspace(0, num_knots, num=num_knots)

        # Parameter values
        min_parameter = self.knots[0]
        max_parameter = self.knots[-1]
        self.parameters = numpy.linspace(min_parameter, max_parameter, num_points)
        self.get_positions()

    # Returns the positions of each point based on the control points
    def get_positions(self):
        #### Basis Function Definition ####
        def basis_function(parameter, degree, i):
            # Check for curve simplicity
            if degree == 0:
                if self.knots[i] <= parameter < self.knots[i + 1]: return 1
                return 0
            # First Term
            if self.knots[i + degree] == self.knots[i]: 
                term1 = 0
            else: 
                term1 = ((parameter - self.knots[i]) / (self.knots[i + degree] - self.knots[i])) * basis_function(parameter, degree - 1, i)
            # Second Term
            if self.knots[i + degree + 1] == self.knots[i + 1]: 
                term2 = 0
            else: 
                term2 = ((self.knots[i + degree + 1] - parameter) / (self.knots[i + degree + 1] - self.knots[i + 1])) * basis_function(parameter, degree - 1, i + 1)
            return term1 + term2
        #### Basis Function Definition End ####
        for parameter in self.parameters:
            point = Point(0,0,0)
            for i in range(len(self.control_points)):
                basis_coefficient = basis_function(parameter, self.degree, i)
                point += self.control_points[i] * basis_coefficient
            self.positions.append(point)

    # Plots the spline points.
    def plot_points(self):
        if len(self.positions) == 0:
            print("No positions have been generated. Parameter(s) must be provided")
            return

        # Retrieve the x-components and y-components of the spline points and control points
        spline_x = [point.x_coord for point in self.positions]
        spline_y = [point.y_coord for point in self.positions]
        control_x = [point.x_coord for point in self.control_points]
        control_y = [point.y_coord for point in self.control_points]

        # Plot all points
        plt.scatter(spline_x, spline_y, label = "Resultant BSpline Curve")
        plt.plot(spline_x, spline_y, label = "Resultant BSpline Curve")
        plt.scatter(control_x, control_y)
        plt.axis("equal")
        plt.show()

class BezierSection:
    def __init__(self, parameters):
        self.parameters = parameters

    def get_positions(self):
        # Load parameters
        upper_control_x = numpy.array(self.parameters["upper"]["x"]).flatten()
        upper_control_y = numpy.array(self.parameters["upper"]["y"]).flatten()
        upper_end_control_x = numpy.array(self.parameters["upper end"]["x"]).flatten()
        upper_end_control_y = numpy.array(self.parameters["upper end"]["y"]).flatten()
        lower_control_x = numpy.array(self.parameters["lower"]["x"]).flatten()
        lower_control_y = numpy.array(self.parameters["lower"]["y"]).flatten()
        self.alpha = numpy.deg2rad(self.parameters["inlet angle"])
        self.beta = numpy.deg2rad(self.parameters["exit angle"])

        # Create Bezier curves for upper and lower surfaces
        self.upper = Bezier(points = [Point(x,y,0) for (x,y) in list(zip(upper_control_x, upper_control_y))], resolution=50) 
        self.upper_end = Bezier(points=[Point(x,y,0) for (x,y) in list(zip(upper_end_control_x, upper_end_control_y))], resolution=2)
        self.lower = Bezier(points=[Point(x,y,0) for (x,y) in list(zip(lower_control_x, lower_control_y))], resolution=52)
        # LE 
        AR = section["LE AR"]
        half_wedge = numpy.deg2rad(self.parameters["LE angle"])
        # Check for horizontal or vertical ellipse
        if AR < 1:
            e = numpy.sqrt(1 - AR**2)
        else:
            e = numpy.sqrt(1 - (1/AR)**2)
        if e == 0:
            # Circle
            theta = numpy.linspace(numpy.pi/2 + (half_wedge+self.alpha), 3*numpy.pi/2 - (half_wedge+self.alpha), 25)
            x1 = upper_x[0]
            y1 = upper_y[0]
            x0 = 0
            y0 = 0
            m1 = - 1 / slope_upper
            m0 = numpy.tan(alpha)
            center_x = (m1*x1 - m0*x0 + y0 - y1) / (m1 - m0)
            center_y = center_x * m0
            R = numpy.sqrt((x1 - center_x)**2 + (y1 - center_y)**2)
            r = R / (1 - e*numpy.cos(theta - alpha))
            LE_x = r*numpy.cos(theta) + center_x
            LE_y = r*numpy.sin(theta) + center_y
        elif e > 0 and e < 1:
            # Ellipse
            x0 = upper_x[0]
            y0 = upper_y[0]
            m0 = slope_upper
            x1 = lower_x[0]
            y1 = lower_y[0]
            m1 = slope_lower
            k = ((x0 - x1) + AR**2*(m0*y0 - m1*y1)) / (AR**2 * (m0 - m1))
            h = x0 + m0*AR**2*(y0 - k)
            angle = numpy.pi - abs(numpy.atan((h - x0)/(k - y0)))
            theta = numpy.linspace(numpy.pi/2 + (angle+alpha), 3*numpy.pi/2 - (angle+alpha), 25)
            a = numpy.sqrt(x0**2 - 2*x0*h + h**2 + y0**2*AR**2 - 2*y0*k*AR**2 + k**2*AR**2)
            b = a / AR
            r = a*b / numpy.sqrt((b*numpy.cos(theta-alpha))**2 + (a*numpy.sin(theta-alpha))**2)
            LE_x = r*numpy.cos(theta) + h
            LE_y = r*numpy.sin(theta) + k
        # TE
        mid_x = (shift[-1].x_coord + shift[0].x_coord) / 2
        mid_y = (shift[-1].y_coord + shift[0].y_coord) / 2
        AR = section["TE AR"] 
        # Check for horizontal or vertical ellipse
        if AR < 1:
            e = numpy.sqrt(1 - AR**2)
        else:
            e = numpy.sqrt(1 - (1/AR)**2)
        if e == 0:
            # Circle
            x1 = upper_end_x[-1]
            y1 = upper_end_y[-1]
            x0 = mid_x
            y0 = mid_y
            m1 = - 1 / ((upper_end_y[-1] - upper_end_y[-2]) / (upper_end_x[-1] - upper_end_x[-2]))
            m0 = numpy.tan(beta)
            center_x = (m1*x1 - m0*x0 + y0 - y1) / (m1 - m0)
            center_y = y0 + m0 * (center_x - x0)
            R = numpy.sqrt((x1 - center_x)**2 + (y1 - center_y)**2)
            theta = numpy.linspace(1.5*numpy.pi - beta, 2.5*numpy.pi - beta, 25)
            r = R / (1 - e*numpy.cos(theta - beta))
            TE_x = r*numpy.cos(theta) + center_x
            TE_y = r*numpy.sin(theta) + center_y
        elif e > 0 and e < 1:
            # Ellipse
            x0 = upper_end_x[-1]
            y0 = upper_end_y[-1]
            m0 = (upper_end_y[-1] - upper_end_y[-2]) / (upper_end_x[-1] - upper_end_x[-2])
            x1 = lower_x[-1]
            y1 = lower_y[-1]
            m1 = (lower_y[-1] - lower_y[-2]) / (lower_x[-1] - lower_x[-2])
            theta = numpy.linspace(-beta - numpy.pi/2, - beta + numpy.pi/2, 25)
            l = numpy.sqrt((x1 - x0)**2 + (y1 - y0)**2)
            r = (l/2) / (1 + e*numpy.cos(theta + beta))
            TE_x = r * numpy.cos(theta) + mid_x
            TE_y = r * numpy.sin(theta) + mid_y
        x = numpy.concatenate((numpy.flip(LE_x), upper_x, upper_end_x, numpy.flip(TE_x[2:-2]), numpy.flip(lower_x)))
        y = numpy.concatenate((numpy.flip(LE_y), upper_y, upper_end_y, numpy.flip(TE_y[2:-2]), numpy.flip(lower_y)))
        positions = [Point(x,y,0) for (x,y) in list(zip(x, y))]
        return positions

    def plot_section(self):
        x = [point.x_coord for point in self.positions]
        y = [point.y_coord for point in self.positions]
        plt.plot(x, y)
        plt.axis("equal")
        plt.show()

# Classes for thin airfoil theory
# Numerical Thin Airfoil Theory using Weissinger's Approximation
class NACA:
    def __init__(self, designation, chord=None, num_points=None):
        self.chord = chord
        self.num_points = num_points
        self.designation = designation
        if len(self.designation) == 4:
            self.digit = 4
            self.M = int(self.designation[0]) / 100
            self.P = int(self.designation[1]) / 10
            self.XX = int(self.designation[2:4]) / 100
        elif len(self.designation) == 5:
            self.digit = 5
            self.L = int(self.designation[0]) / 100
            self.P = int(self.designation[1]) / 100
            self.Q = int(self.designation[2]) / 100
            self.XX = int(self.designation[3:5]) / 100

    def get_upper(self, x_c):
        match self.digit:
            case 4:
                z_c, dz_dx = self.get_camber(x_c)
                z_t = self.XX * 5 * (0.2969 * numpy.sqrt(x_c) - 0.126 * x_c - 0.3516*x_c**2 + 0.2843*x_c**3 - 0.1015*x_c**4)
                theta = numpy.atan(dz_dx)
                x_u = x_c - z_t/2 * numpy.sin(theta)
                z_u = z_c + z_t/2 * numpy.cos(theta)

        return x_u, z_u

    def get_lower(self, x_c):
        match self.digit:
            case 4:
                z_c, dz_dx = self.get_camber(x_c)
                z_t = self.XX * 5 * (0.2969 * numpy.sqrt(x_c) - 0.126 * x_c - 0.3516*x_c**2 + 0.2843*x_c**3 - 0.1015*x_c**4)
                theta = numpy.atan(dz_dx)
                x_u = x_c + z_t/2 * numpy.sin(theta)
                z_u = z_c - z_t/2 * numpy.cos(theta)

        return x_u, z_u

    def get_camber(self, x_c):
        # input x/c instead of x-coordinates to get the chord-normalized z-coordinates
        match self.digit:
            case 4:
                ''' 4-digit Airfoil equations referenced from Airfoil Tools '''
                z_c = (self.M/self.P**2) * (2*self.P*x_c - x_c**2)
                dz_dx = (2*self.M/self.P**2) * (self.P - x_c)
            case 5:
                pass
        
        return z_c, dz_dx

    def plot_airfoil(self):
        if self.num_points is not None:
            x_values = numpy.linspace(0, 1, self.num_points)
            upper_z = numpy.vectorize(self.get_upper)(x_values)[1]
            upper_x = numpy.vectorize(self.get_upper)(x_values)[0]
            lower_x = numpy.vectorize(self.get_lower)(x_values)[0]
            lower_z = numpy.vectorize(self.get_lower)(x_values)[1]
            camber_z = numpy.vectorize(self.get_camber)(x_values)[0]
            fig, ax = plt.subplots()
            ax.plot(upper_x, upper_z)
            ax.plot(lower_x, lower_z)
            ax.plot(x_values, camber_z)
            ax.set_aspect("equal")
            plt.show()
        else:
            print("No point resolution provided (""num_points"" attribute).")
            return

    
class Thin_Airfoil(NACA):
    def __init__(self, designation:str, chord, num_panels, flap_angle=0, flap_start=None, ambient=None):
        super().__init__(designation)
        self.num_panels = num_panels
        self.chord = chord
        self.flap_angle = flap_angle
        self.flap_start = flap_start
        self.ambient = ambient
        self.panels = self.generate_panels()
        
    def generate_panels(self):
        x_positions = numpy.linspace(0, 1, self.num_panels + 1)
        z_positions = numpy.array([self.get_camber(x)[0] for x in x_positions])
        panels = [Panel(Point(x_positions[i], z_positions[i], 0), Point(x_positions[i + 1], z_positions[i + 1], 0)) for i in range(self.num_panels)]

        common_length = 1 / (self.num_panels + 1)

        precision = 5
        error_threshold = 10**-precision
        panel_lengths = set(round(panel.length, precision) for panel in panels)
        while len(panel_lengths) > 1:
            for idx, panel in enumerate(panels):
                if idx > 0:
                    panel.start_point = copy.deepcopy(panels[idx - 1].end_point)
                if idx == len(panels) - 1:
                    continue
                length_error = abs(common_length - panel.length)
                while length_error > error_threshold:
                    #print(panel.get_length(), common_length, panel.start_point.x_coord, panel.end_point.x_coord)
                    if panel.length < common_length:
                        panel.end_point.set_x(panel.end_point.x_coord + error_threshold/2)
                    elif panel.length > common_length and panel.end_point.x_coord > panel.start_point.x_coord:
                        panel.end_point.set_x(panel.end_point.x_coord - error_threshold/2)
                    if panel.end_point.x_coord > 1: panel.end_point.set(1)
                    else: break
                    panel.end_point.set_y(self.get_camber(panel.end_point.x_coord)[0])
                    length_error = abs(common_length - panel.length)
            if panels[-1].length > common_length:
                common_length += abs(common_length - panels[-1].length) / self.num_panels
            elif panels[-1].length < common_length:
                common_length -= abs(common_length - panels[-1].length) / self.num_panels
            panel_lengths = set(round(panel.length, precision - 2) for panel in panels)

        # Each panel's vortex influences every panel's 3/4 flow tangency point
        self.quarter_points = numpy.tile(numpy.array([panel.get_position(0.25) for panel in panels]).reshape(self.num_panels, 1), (1, self.num_panels)) # vortex position vectors
        self.three_quarter_points = numpy.tile(numpy.array([panel.get_position(0.75) for panel in panels]).reshape(1, self.num_panels), (self.num_panels, 1)) # flow tangency position vectors

        return panels

    def get_camber(self, x_c):
        # input x/c instead of x-coordinates to get the chord-normalized z-coordinates
        match self.digit:
            case 4:
                ''' 4-digit Airfoil equations referenced from Airfoil Tools '''
                if 0 <= x_c < self.P:
                    z_c = (self.M/self.P**2) * (2*self.P*x_c - x_c**2)
                    dz_dx = (2*self.M/self.P**2) * (self.P - x_c)
                elif self.flap_start == None and self.P <= x_c <= 1:
                    z_c = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*x_c - x_c**2)
                    dz_dx = (2*self.M/(1 - self.P)**2) * (self.P - x_c)
                elif self.flap_start != None and self.P <= x_c <= self.flap_start:
                    z_c = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*x_c - x_c**2)
                    dz_dx = (2*self.M/(1 - self.P)**2) * (self.P - x_c)
                elif self.flap_start != None and self.flap_start < x_c <= 1:
                    flap_start_height = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*self.flap_start - self.flap_start**2) # b
                    dz_dx = numpy.tan(-self.flap_angle) # m
                    z_c = dz_dx*(x_c - self.flap_start) + flap_start_height # mx + b
                else: return None
            case 5:
                pass
        
        return z_c, dz_dx


    @property
    def lift(self):

        # Camber slopes at the points of tangency (three quarter point on panel)
        slopes = numpy.array([self.get_camber(position.x_coord)[1] for position in self.three_quarter_points[0]]) 
        
        # Connect every quarter-panel point to every three-quarter-panel point of every panel
        panel_connectors = numpy.vectorize(Panel)(self.quarter_points, self.three_quarter_points)
        r = numpy.vectorize(lambda panel: panel.length)(panel_connectors) # distances between control points and flow tangency points

        # Induced velocity with an absolute angle determined by a clockwise circulation
        theta = numpy.vectorize(lambda panel: panel.phi + 3*numpy.pi/2)(panel_connectors) # direction of induced velocity
        panel_lengths = numpy.tile(numpy.array([panel.length for panel in self.panels]), (self.num_panels, 1))

        # Apply Weissinger's Approximation in matrix form
        A = panel_lengths / (2*numpy.pi*r) * (numpy.sin(theta) - numpy.cos(theta)*slopes) # vertical velocities minus horizontal velocities
        b = numpy.transpose(self.ambient.Vinf*(numpy.cos(self.ambient.alpha)*slopes - numpy.sin(self.ambient.alpha)))

        # Lift
        gamma_distribution = numpy.linalg.solve(A, b) # distribution of length-specific vortex strengths
        Vz = self.ambient.Vinf*numpy.sin(self.ambient.alpha) + numpy.sum(((panel_lengths*gamma_distribution) / (2*numpy.pi*r)) * numpy.sin(slopes), axis=0)
        Vx = self.ambient.Vinf*numpy.cos(self.ambient.alpha) + numpy.sum(((panel_lengths*gamma_distribution) / (2*numpy.pi*r)) * numpy.cos(slopes), axis=0)
        V = numpy.vectorize(numpy.sqrt)(Vx**2 + Vz**2)
        Cp = 1 - (V/self.ambient.Vinf)**2
        L = float(self.ambient.rhoinf * self.ambient.Vinf * sum(gamma_distribution * panel_lengths[0]))
        L_prime = self.ambient.rhoinf * self.ambient.Vinf * gamma_distribution * panel_lengths[0]
        Cl = float(L / ((1/2)*self.ambient.rhoinf*self.ambient.Vinf**2))

        return L, Cl, L_prime, Cp
    
    @property
    def moment(self):
        if self.flap_start != None:
            theta_f = numpy.acos(1 - (2*self.flap_start))
            A1 = 4*self.M / (numpy.pi * (1 - self.P)**2) * ((1 - self.chord/2)*numpy.sin(theta_f) + (theta_f*self.chord/2) + (numpy.sin(2*theta_f)/4)) - ((2/numpy.pi)*numpy.tan(-self.flap_angle)*numpy.sin(theta_f))
            A2 = 4*self.M / (numpy.pi * (1 - self.P)**2) * ((1 - self.chord/2)*(theta_f/2 + numpy.sin(2*theta_f)/4) + (self.chord/2)*(numpy.sin(theta_f) - (1/3)*(numpy.sin(theta_f)**3))) + (2/numpy.pi*numpy.tan(-self.flap_angle)*((numpy.pi - theta_f)/2 - numpy.sin(2*theta_f)/4))
        else:
            A1 = self.chord * self.M / (1 - self.P)**2
            A2 = - self.chord * self.M / (2 * (1 - self.P)**2)
        cmc4 = (numpy.pi/4) * (A2 - A1)
        return cmc4

    def plot_results(self):
        fig, ax = plt.subplots(4, 1, figsize=(12, 5))

        # Camberline
        x_positions = list()
        y_positions = list()
        for panel in self.panels:
            x_positions.append(float(panel.start_point.x_coord))
            x_positions.append(float(panel.end_point.x_coord))
            y_positions.append(float(panel.start_point.y_coord))
            y_positions.append(float(panel.end_point.y_coord))
        
        control_x = [point.x_coord for point in self.quarter_points.flatten()]
        control_y = [point.y_coord for point in self.quarter_points.flatten()]
        vortex_x = [point.x_coord for point in self.three_quarter_points.flatten()]
        vortex_y = [point.y_coord for point in self.three_quarter_points.flatten()]
        ax[0].scatter(control_x, control_y, color="blue", label="Control Points")
        ax[0].scatter(vortex_x, vortex_y, color="red", label="Flow Tangency Points")
        ax[0].scatter(x_positions, y_positions, color='green', label="Endpoints")
        ax[0].plot(x_positions, y_positions)
        ax[0].set_xlim(0, 1)
        ax[0].set_ylim(-0.2, 0.2)
        ax[0].legend()

        _, _, L_prime, Cp = self.lift
        x_c = numpy.transpose(numpy.array([point.x_coord for point in self.three_quarter_points[0]]))
        
        # Lift Distribution
        ax[1].plot(x_c, L_prime, label="L_prime")
        if self.flap_start != None: ax[1].title("Flapped " + self.designation + " Lift Distribution")
        else: ax[1].set_title(self.designation + " Lift Distribution")

        # Pressure Distribution
        ax[2].plot(x_c, Cp, label="Cp")
        if self.flap_start != None: ax[2].title("Flapped " + self.designation + " Pressure Distribution")
        else: ax[2].set_title(self.designation + " Pressure Distribution")

        plt.show()
