from pyaero.aerodynamics.curves import *
from pyaero.aerodynamics.potential import *
import os

# NACA 
airfoil = NACA(designation="2412", num_points=100)
airfoil.plot_airfoil()

# Bezier 
control_points = [Point(0,0,0), Point(5,0.25,0), Point(7,0,0), Point(5, -0.25,0)]
bezier = Bezier(control_points, resolution=25)
bezier.plot_points()

# Spline 
control_points = [Point(0,0,0), Point(5,0.25,0), Point(7,0,0), Point(5, -0.25,0)]
spline = BSpline(control_points, degree=3, num_points=100)
spline.plot_points()
