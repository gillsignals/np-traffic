%% Parameters - initial values

% extracellular degradation rate constant for nm_x 
pbae_t_nps  = [4.4 5.1 1.2 3.6 6.1 5.3 5.5 6.5 1.6 4.6 3.9];    % From ref (1), half-lives (1/h) of several PBAE NPs
t_np        = median(pbae_t_nps) * 60;  % Find median half-life and convert to minutes
p.k_deg_x     = log(2) / t_np;        % Convert half-life to degradation rate constant (2.5e-3)  

% internalization rate constant for nm_x
p.k_int       = 6e-3;     % min^-1. From ref (2), k_bind_uptake = 5-7e-3 min^-1

% lysosomal degradation rate constant for nm_e
p.k_lys       = 1.5e-2;	% min^-1. From ref (2), k_deg_vesicle = 1-2e-2 min^-1

% endosomal escape rate constant for nm_e
p.k_escape    = 1e-3;     % min^-1. From ref (2), k_escape highly variable: 1e-2 - 1e-5 min^-1

% unpackaging/degradation rate constant for nm_c (release mRNA, degrade np)
p.k_deg_np	= 1e-2;     % min^-1. From ref (2), k_unpack = 1e-2 - 5e-1 min^-1

% degradation rate constant for m_c
t_m         = 9 * 60;           % min. From ref (3), median mRNA half-life = 9h
p.k_deg_m 	= log(2) / t_m;     % convert half-life to degradation rate constant (~1e-3)
%p.k_deg_m     = .062 / 60;        % From ref (4b), gamma = 0.062/h

% translation/expression rate constant for p_c
%p.k_expr      = 1e-2;     % min^-1 . From ref (4a), k_protein = translation rate of GFP mRNA from plasmid
p.k_expr      = 170 / 60; % min^-1. From ref (4b), k_TL = 170/h...really depends how much mRNA gets in...

% degradation rate constant for p_c
t_p         = 46 * 60;          % min. From ref (3), median protein half-life = 46h
p.k_deg_p     = log(2) / t_p;     % min^-1. convert half-life to degradation rate constant (~2.5e-4 min^-1)
%o,k_deg_p     = .056 / 60;        % From ref (4b), gamma = 0.056/h


%% Citations for initial parameter values

% (1) Sunshine et al., 2012, Table 1 - half-life of 11 different PBAE nanoparticles in PBS at 37C

% (2) Varga et al., 2005, Table 1 - lsq fit parameter values for gene delivery with adenovirus,
% PEI nanoparticles with different sizes/mods, lipofectamine, and naked DNA

% (3) Schwanhausser et al., 2011 - converted from median protein/mRNA half-life values in 3T3 cells

% (4a) Varga et al., 2001, Table 1 - Fit to FACS data

% (4b) Leonhardt et al., 2014, Supplement