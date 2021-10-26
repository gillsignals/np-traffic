function [dydt] = equations(t, y, sp, p)
% equations for basic nanoparticle-mRNA trafficking and expression model
% t = time range over which to solve ODEs
% y = current values of molecular species
% sp = species index values
% p = parameter values (reaction rate constants)

% Equations
dydt(sp.c_x)	= -p.k_int * y(sp.c_x) - p.k_deg_x * y(sp.c_x);                           % d/dt nm_x
dydt(sp.c_e)	= p.k_int * y(sp.c_x) - p.k_escape * y(sp.c_e) - p.k_lys * y(sp.c_e);    % d/dt nm_e
dydt(sp.c_c)	= p.k_escape * y(sp.c_e) - p.k_deg_np * y(sp.c_c);                        % d/dt nm_c
dydt(sp.m_c)	= p.k_deg_np * y(sp.c_c) - p.k_deg_m * y(sp.m_c);                          % d/dt m_c
dydt(sp.p_c)	= p.k_expr * y(sp.m_c) - p.k_deg_p * y(sp.p_c);                            % d/dt p_c

% convert equations to column vector
dydt = dydt';

return;