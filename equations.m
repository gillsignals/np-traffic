function [dydt] = equations(t, y, p)
% equations for basic nanoparticle-mRNA trafficking and expression model
% t = time range over which to solve ODEs
% y = current values of molecular species
% p = parameter values (reaction rate constants)

% Variables
nm_x	= y(1);     % extracellular nanoparticle-mRNA complex concentration
nm_e	= y(2);     % endosomal nanoparticle-mRNA complex concentration
nm_c	= y(3);     % cytosolic nanoparticle-mRNA complex concentration
m_c     = y(4);     % cytosolic free mRNA concentration
p_c     = y(5);     % cytosolic protein concentration


% Parameters
k_deg_x     = p(1);     % extracellular degradation rate constant for nm_x             
k_int       = p(2);     % internalization rate constant for nm_x
k_lys       = p(3);     % lysosomal degradation rate constant for nm_e
k_escape    = p(4);     % endosomal escape rate constant for nm_e
k_deg_np	= p(5);     % unpackaging/degradation rate constant for nm_c (release mRNA, degrade np)
k_deg_m 	= p(6);     % degradation rate constant for m_c
k_expr      = p(7);     % translation/expression rate constant for p_c
k_deg_p     = p(8);     % degradation rate constant for p_c


% Equations
dydt(1)     = -k_int * nm_x - k_deg_x * nm_x;                   % d/dt nm_x
dydt(2)     = k_int * nm_x - k_escape * nm_e - k_lys * nm_e;    % d/dt nm_e
dydt(3)     = k_escape * nm_e - k_deg_np * nm_c;                % d/dt nm_c
dydt(4)     = k_deg_np * nm_c - k_deg_m * m_c;                  % d/dt m_c
dydt(5)     = -k_expr * m_c - k_deg_p * p_c;                    % d/dt p_c

return;