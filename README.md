Programs for aerospace applications.

Split into aerodynamics, stability, and propulsion categories (control code may be added, but probably not cause it's harder than using simulink).

- "Aerodynamics" contains incompressible (potential flow theory) and compressible tools.
  - curve design (useful for wing sections and turbomachinery blade design)
  - potential flow solvers
  - compressible flow lookup tables
  - functions for hypersonic vehicle design (Taylor-Maccoll integrator and Busemann Inlet) 
- "Propulsion" contains programs that make the gas turbine design process a little easier by integrating different parts of the design process.
  - the engine class is intended to wrap the preliminary design steps into a simple-to-manipulate object
    - flight performance (constraint & mission analysis) 
    - cycle analysis (design-point and off-design)
    - preliminary component design
- "Stability" contains a calculator for stability derivatives based on aircraft geometry
  - the geometry must be provided by the user
  - the process is laborious because it acts as just a calculator right now
  - lacks a control analysis
