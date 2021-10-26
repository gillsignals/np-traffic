clear all;
close all;

%% Import colorblind-friendly colors for plotting
colorblind_colors;

%% Assign molecular species names to index values

sp.c_x     = 1;     % extracellular nanoparticle-mRNA complex (pg/cell)
sp.c_e     = 2;     % endosomal nanoparticle-mRNA complex concentration (pg/cell)
sp.c_c     = 3;     % cytosolic nanoparticle-mRNA complex concentration (pg/cell)
sp.m_c      = 4;     % cytosolic free mRNA number (#/cell)
sp.p_c      = 5;     % cytosolic protein number (#/cell)

%% Initial conditions

% initialize all species to 0
y0 = zeros(length(fields(sp)), 1);

y0(sp.c_x)     = 620;  % initialize extracellular nanoparticle amount to 620 pg/cell

%% Parameters - specify base values

% extracellular degradation rate constant for nm_x 
pbae_t_nps  = [4.4 5.1 1.2 3.6 6.1 5.3 5.5 6.5 1.6 4.6 3.9];    % From ref (1), half-lives (1/h) of several PBAE NPs
t_np        = median(pbae_t_nps) * 60;  % Find median half-life and convert to minutes
p.k_deg_x     = log(2) / t_np;        % Convert half-life to degradation rate constant (2.5e-3)  

% internalization rate constant for c_x
p.k_int       = 6e-3;     % min^-1. From ref (1), k_bind_uptake = 5-7e-3 min^-1

% lysosomal degradation rate constant for c_e
p.k_lys       = 1.5e-2;	% min^-1. From ref (1), k_deg_vesicle = 1-2e-2 min^-1

% endosomal escape rate constant for c_e
p.k_escape    = 1e-3;     % min^-1. From ref (1), k_escape highly variable: 1e-2 - 1e-5 min^-1

% unpackaging/degradation rate constant for c_c (release mRNA, degrade np)
p.k_deg_np	= 1e-2;     % min^-1. From ref (1), k_unpack = 1e-2 - 5e-1 min^-1

% degradation rate constant for m_c
t_m         = 9 * 60;           % min. From ref (2), median mRNA half-life = 9h
p.k_deg_m 	= log(2) / t_m;     % convert half-life to degradation rate constant (~1e-3)
%p.k_deg_m     = .062 / 60;        % From ref (3b), gamma = 0.062/h

% translation/expression rate constant for p_c
%p.k_expr      = 1e-2;     % min^-1 . From ref (3a), k_protein = translation rate of GFP mRNA from plasmid
p.k_expr      = 170 / 60; % min^-1. From ref (3b), k_TL = 170/h...really depends how much mRNA gets in...

% degradation rate constant for p_c
t_p         = 46 * 60;          % min. From ref (2), median protein half-life = 46h
p.k_deg_p     = log(2) / t_p;     % min^-1. convert half-life to degradation rate constant (~2.5e-4 min^-1)
%p.k_deg_p     = .056 / 60;        % From ref (3b), gamma = 0.056/h

%% SIMULATION 1

% simulate for 72h
[T, Y] = main_ode([0 72*60], y0, sp, p);

% calculate maximum concentrations of each species
cmax_c_x = max(Y(:, sp.c_x));
cmax_c_e = max(Y(:, sp.c_e));
cmax_c_c = max(Y(:, sp.c_c));
cmax_m_c = max(Y(:, sp.m_c));
cmax_p_c = max(Y(:, sp.p_c));    % max protein conc is a logical thing to look at for sensitivity analysis

% create figure for simulation 1
figure;

% plot C_x for sim 1
ax1 = subplot(3,2,1);
plot(ax1, T / 60, Y(:, sp.c_x), 'LineWidth', 2, 'Color', colors(cb.black))
title("Extracellular NP-mRNA complex (C_x)")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("C_x (pg/cell)")

ax2 = subplot(3,2,3);
plot(ax2, T / 60, Y(:, sp.c_e), 'LineWidth', 2, 'Color', colors(cb.red))
title("Endosomal NP-mRNA complex (C_e)")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("C_e (pg/cell)")

ax3 = subplot(3,2,5);
plot(ax3, T / 60, Y(:, sp.c_c), 'LineWidth', 2, 'Color', colors(cb.skyblue))
title("Cytosolic NP-mRNA complex (C_c)")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("C_c (pg/cell)")

ax4 = subplot(3,2,2);
plot(ax4, T / 60, Y(:, sp.m_c), 'LineWidth', 2, 'Color', colors(cb.orange))
title("Cytosolic free mRNA (M_c)")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("M_c (#/cell)")

ax5 = subplot(3,2,4);
plot(ax5, T / 60, Y(:, sp.p_c), 'LineWidth', 2, 'Color', colors(cb.green))
title("Cytosolic protein (P_c)")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("P_c (#/cell)")

ax6 = subplot(3,2,6);
plot(ax6, T / 60, Y(:, sp.c_x) / cmax_c_x, 'LineWidth', 2, 'Color', colors(cb.black))
hold on;
plot(ax6, T / 60, Y(:, sp.c_e) / cmax_c_e, 'LineWidth', 2, 'Color', colors(cb.red))
plot(ax6, T / 60, Y(:, sp.c_c) / cmax_c_c, 'LineWidth', 2, 'Color', colors(cb.skyblue))
plot(ax6, T / 60, Y(:, sp.m_c) / cmax_m_c, 'LineWidth', 2, 'Color', colors(cb.orange))
plot(ax6, T / 60, Y(:, sp.p_c) / cmax_p_c, 'LineWidth', 2, 'Color', colors(cb.green))
hold off;
title("Combined scaled time course")
xlabel("Time (h)")
xlim([0 72]);
xticks(linspace(0,72,7));
ylabel("Fraction of maximum amount")
legend(["C_x", "C_e", "C_c", "M_c", "P_c"], "Location", "eastoutside")




%% Citations for initial parameter values

% (1) Sunshine et al., 2012, Table 1 - half-life of 11 different PBAE nanoparticles in PBS at 37C

% (2) Varga et al., 2005, Table 1 - lsq fit parameter values for gene delivery with adenovirus,
% PEI nanoparticles with different sizes/mods, lipofectamine, and naked DNA

% (3) Schwanhausser et al., 2011 - converted from median protein/mRNA half-life values in 3T3 cells

% (4a) Varga et al., 2001, Table 1 - Fit to FACS data

% (4b) Leonhardt et al., 2014, Supplement
