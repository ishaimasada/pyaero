Programs for aerospace applications.

Split into aerodynamics, stability, and propulsion categories (control code may be added, but probably not cause it's harder than using simulink).

- Aerodynamics contains incompressible (potential flow theory) and compressible applications.
  - curve design (useful for wing sections and similar geometries)
  - potential flow solvers
  - 1D compressible flow lookup tables
  - functions for hypersonic vehicle design (Taylor-Maccoll integrator and Busemann Inlet) 
- Propulsion contains programs that automate the jet engine design process.
  - the engine class is intended to wrap the preliminary design steps into a simple-to-manipulate object
    - flight performance (constraint & mission analysis) --> Not Implemented 
    - design point cycle analysis --> Functional
    - preliminary component design functionality --> Not Implemented
- Stability contains a calculator for stability derivatives based on aircraft geometry
  - the geometry must be provided first and each value input into the program
  - the process is laborious because it acts as just a calculator right now
  - lacks a control analysis
