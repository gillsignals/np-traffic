%% Parameters - initial values

% extracellular degradation rate constant for nm_x 
k_deg_x     = 0;        % Assume no NP degradation in interstitium   

% internalization rate constant for nm_x
k_int       = 6e-3;     % min^-1. From ref (1), k_bind_uptake = 5-7e-3 min^-1

% lysosomal degradation rate constant for nm_e
k_lys       = 1.5e-2;	% min^-1. From ref (1), k_deg_vesicle = 1-2e-2 min^-1

% endosomal escape rate constant for nm_e
k_escape    = 1e-3;     % min^-1. From ref (1), k_escape highly variable: 1e-2 - 1e-5 min^-1

% unpackaging/degradation rate constant for nm_c (release mRNA, degrade np)
k_deg_np	= 1e-2;     % min^-1. From ref (1), k_unpack = 1e-2 - 5e-1 min^-1

% degradation rate constant for m_c
t_m         = 9 * 60;           % min. From ref (3), median mRNA half-life = 9h
k_deg_m 	= log(2) / t_m;     % convert half-life to degradation rate constant (~1e-3)

% translation/expression rate constant for p_c
k_expr      = 1e-2;     % min^-1 . From ref (2), k_protein = translation rate of GFP mRNA from plasmid

% degradation rate constant for p_c
t_p         = 46 * 60;          % min. From ref (3), median protein half-life = 46h
k_deg_p     = log(2) / t_p;     % convert half-life to degradation rate constant (~2.5e-4)


%% Citations for initial parameter values

% (1) Varga et al., 2005, Table 1 - lsq fit parameter values for gene delivery with adenovirus,
% PEI nanoparticles with different sizes/mods, lipofectamine, and naked DNA

% (2) Varga et al., 2001, Table 1 - Fit to FACS data

% (3) Schwanhausser et al., 2011 - converted from median protein/mRNA half-life values in 3T3 cells
