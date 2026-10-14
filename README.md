# 3D Wing Aerodynamic Solver

A MATLAB implementation of **Prandtl’s Lifting-Line Theory** using horseshoe vortex elements to analyse 3D wing aerodynamics.

## Features
- Iterative circulation and downwash solver with numerical relaxation.
- Configurable wing geometry and airfoil polar data.
- Spanwise integration of lift and drag coefficients.
- Single angle-of-attack analysis, angle-of-attack sweeps, and discretisation convergence studies.
- Automated export of plots and numerical results.

## Repository Structure

```text
├── main.m
├── hshoe.m
├── N65-2-415.dat
├── output examples/ # Directory containing generated plots (.png) and data tables (.txt)
└── README.md
```

## Usage

1. Open the project folder in MATLAB.
2. Ensure the required helper functions and toolboxes are available.
3. Run:

   ```matlab
   main
   ```

4. Follow the prompts to configure the wing geometry, flow conditions, and analysis mode.

## Dependencies
- MATLAB
- `hshoe.m` — horseshoe vortex induced-velocity function.
- Aerospace Toolbox (`atmosisa`), if required by the selected analysis mode.
