### ROSS Sample Code (Requires Python versions later than 3.13) ###

import os, numpy, subprocess
import ross 
from ross.units import Q_
from ross.probe import Probe

# clears the terminal (package warnings from ROSS dependencies)
subprocess.run(["cls"], shell=True)

# Inputs
L1 = 0.044
L2 = 0.08
L3 = 0.044
od1 = 0.0085374
od2 = 0.0227236
od3 = 0.0085374
num_elements1 = 3
num_elements2 = 5
num_nodes = 2*num_elements1 + num_elements2
disk_mass = 0.63220218 # kg from inconel 718 onshape custom material (density 8190kg/m^3)
disk_ID = 782.79312 * 1e-6 # kg*m^2
disk_IP = 1561.88144 * 1e-6 # kg*m^2
bearing1_node = round(44/168 * num_nodes)
bearing2_node = round((168-44)/168 * num_nodes)
num_balls = 10
ball_diameter = 0.008
static_loading = 500
shaft_speed = 48000 * numpy.pi/30
speed_range = numpy.linspace(0, shaft_speed, 50)
steel = ross.steel

# Assembling Shaft Parameters
lengths = numpy.concatenate((numpy.full(num_elements1,L1/num_elements1), numpy.full(num_elements2,L2/num_elements2), numpy.full(num_elements1,L3/num_elements1)))
od = numpy.concatenate((numpy.full(num_elements1,od1), numpy.full(num_elements2,od2), numpy.full(num_elements1,od3)))

# Objects (shaft elements, disk, bearings)
elements = [ross.ShaftElement(L=lengths[i], idl=0, odl=od[i], material=steel, shear_effects=True, rotary_inertia=True, gyroscopic=True) for i in range(num_nodes)]
disk = ross.DiskElement(n=0, m=disk_mass, Ip=disk_IP, Id=disk_ID)
bearing1 = ross.BallBearingElement(n=bearing1_node, n_balls=num_balls, d_balls=ball_diameter, fs=static_loading, alpha=0)
bearing2 = ross.BallBearingElement(n=bearing2_node, n_balls=num_balls, d_balls=ball_diameter, fs=static_loading, alpha=0)

# Rotor Object & Analyses
rotor = ross.Rotor(shaft_elements=elements, disk_elements=[disk], bearing_elements = [bearing1, bearing2])
static_analysis = rotor.run_static()
modal_analysis = rotor.run_modal(speed=shaft_speed, num_modes=16)
torsional_modes = [mode for mode, shape in enumerate(modal_analysis.shapes) if shape.mode_type == "Torsional"]

### RESULTS ###

# Geometry
# rotor.plot_rotor().show()
# print(round(rotor.m, 2))
# print(rotor.CG)

# Static Analysis
# static_analysis.plot_free_body_diagram().show()
# static_analysis.plot_deformation().show()
# static_analysis.plot_shearing_force().show()
# static_analysis.plot_bending_moment().show()

# Modal Analysis
# print(len(modal_analysis.wn))
# modal_analysis.plot_mode_2d(1).show()
# modal_analysis.plot_mode_2d(4).show()
# modal_analysis.plot_mode_2d(6).show()
# modal_analysis.plot_mode_2d(7).show()
# modal_analysis.plot_mode_3d(1).show()
# modal_analysis.plot_mode_3d(4).show()
# modal_analysis.plot_mode_3d(6).show()
# modal_analysis.plot_mode_3d(7).show()

# Torsional Analysis
# modal_analysis.plot_mode_3d(torsional_modes[0], animation=True).show()

# Campbell Diagram
campbell = rotor.run_campbell(speed_range)
# print(campbell.log_dec.shape[1])
# print(campbell.log_dec[:,1])
# print(campbell.log_dec[:,3])
# print(campbell.log_dec[:,5])
# campbell.plot(harmonics=[0.5, 1]).show()

# Unbalance Response (Frequency Domain)
frequency_response = rotor.run_freq_response(speed_range)
global_dof = 0*rotor.number_dof + 1 # node 0, dof: y
# frequency_response.plot(inp=global_dof, out=global_dof).show()
unbalance_response = rotor.run_unbalance_response([0], [0.001], [0], speed_range)
probe = Probe(0, Q_(45, "deg"))  # node 0, orientation 45°
# unbalance_response.plot(probe=[probe]).show()
# unbalance_response.plot_deflected_shape(speed=shaft_speed).show()

# Unbalance Response (Time Domain)
node = 0
num_samples = 100
num_revolutions = 50
num_times = num_revolutions * num_samples
times = numpy.linspace(0, num_revolutions*2*numpy.pi/shaft_speed, num_times)
F = numpy.zeros((num_times, rotor.ndof))
F[:, node * rotor.number_dof + 0] = 10 * numpy.cos(shaft_speed * times)
F[:, node * rotor.number_dof + 1] = 10 * numpy.sin(shaft_speed * times)
time_response = rotor.run_time_response(shaft_speed, F, times)
# time_response.plot_1d(probe=[probe], one_revolution=True).show()
# time_response.plot_2d(node=0, one_revolution=True).show()
# time_response.plot_3d(one_revolution=True).show()

# Undamped Critical Speed (UCS) Map
bearing_k_range = (6, 11)
ucs_results = rotor.run_ucs(stiffness_range=bearing_k_range, num=20, num_modes=16)
# ucs_fig = ucs_results.plot().show()