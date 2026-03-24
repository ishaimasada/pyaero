import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Import the Position Vector class
from curves import *
from potential import *
from atmosphere import Ambient

# NACA class
designation = "2412"
num_points = 100
airfoil = NACA(designation, num_points=num_points)
#airfoil.plot_airfoil()

# Spline Class
control_points = [Point(0,0,0), Point(5,0.25,0), Point(7,0,0), Point(5, -0.25,0)]
degree = 3
num_points = 100
spline = BSpline(control_points, degree, num_points)
#spline.plot_points()

# composite airfoil (Bezier or Cubic Spline)