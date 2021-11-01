function [T, Y, cmax, tmax] = main_ode(tspan, y0, sp, p)

% MAIN_ODE [T,Y] = main_ode(tspan, y0, p)
% Solves the system of ODEs specified in equations.m based on parameter values in driver.m.
% Inputs: time span, initial conditions, and parameters.
% Outputs: time vector and concentration matrix from the ODE solver.
%
% INPUT:
%   tspan = vector of beginning and end of time span over which to solve ODEs 
%   y0 = vector of initial conditions for every molecular species in the model
%   p = vector of reaction parameters - in this case, reaction rate constants
%
% OUTPUT:
%   T = vector of time values returned by ODE solver. Each time value matches to a row in the
%   concentration matrix Y.
%   Y = matrix of concentration values for each molecule species returned by ODE solver. Each row i is
%   a time point corresponding to the ith element in T. Each column j is a molecular species corresponding
%   to the jth species in sp. 
%   cmax = maximum concentration of each species


% Set ODE solver options
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-8, 'NonNegative', 1:length(fields(sp)), 'InitialStep', 1e-2);

% (option explanations from Christy Pickering's matlab-modeling repo)
% AbsTol = absolute error tolerance; if solution is smaller than this, solver will not try to get
% correct digits for this value
% RelTol = relative error tolerance; error tolerance relative to magnitude of each species. The max
% acceptable error is max(AbsTol, RelTol*abs(y(i))).
% NonNegative = species numbers required to be positive; all species here because all species
% represent physical numbers of chemical species, which must be non-negative.
% InitialStep = upper bound of initial time step size tried by the ODE solver

% Run the basic ODE solver on eqns in equations.m
[T,Y]       = ode45(@equations, tspan, y0, options, sp, p);

% extract maximum value for each species and its index
[cmax.c_x cmax_ind.c_x]     = max(Y(:, sp.c_x));
[cmax.c_e cmax_ind.c_e]     = max(Y(:, sp.c_e));
[cmax.c_c cmax_ind.c_c]     = max(Y(:, sp.c_c));
[cmax.m_c cmax_ind.m_c]     = max(Y(:, sp.m_c));
[cmax.p_c cmax_ind.p_c]     = max(Y(:, sp.p_c));

% find time when maximum value is reached for each species
tmax.c_x = T(cmax_ind.c_x);
tmax.c_e = T(cmax_ind.c_e);
tmax.c_c = T(cmax_ind.c_c);
tmax.m_c = T(cmax_ind.m_c);
tmax.p_c = T(cmax_ind.p_c);

return;