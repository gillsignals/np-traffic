# np-traffic
MATLAB model of intracellular trafficking of nanoparticle-mRNA complexes. Written for Cellular Engineering final project.

*To Run:* Download all files from repo. Open `driver.m`. Set the flag value to the desired scenario and run the program. All other files will be automatically called inside of the driver.

- `colorblind_colors.m` - Defines a color palette for plots that is colorblind-accessible.

- `driver.m` - *Main file.* Runs a variety of scenarios using switch-case architecture, including basic simulations, local univariate sensitivity analysis, and global univariate sensitivity analysis.

- `equations.m` - ODEs describing the temporal evolution of the molecular species.

- `figures.m` - Create figure output for each scenario using switch-case architecture.

- `main_ode.m` - Runs the ODE solver using the equations file and parameters specified in driver.

- `parameters.m` - Define parameters and their initial values p0.
