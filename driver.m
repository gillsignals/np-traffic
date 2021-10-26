clear all;
close all;

%% Specify simulation type
flag = 3;
% 1 = simple simulation with basic parameter values for 6h
% 2 = simple simulation with basic parameter values for 72h
% 3 = figure out sensitivity analysis

%% Import colorblind-friendly colors for plotting
colorblind_colors;

%% Assign molecular species names to index values

sp.c_x	= 1;     % extracellular nanoparticle-mRNA complex (pg/cell)
sp.c_e	= 2;     % endosomal nanoparticle-mRNA complex concentration (pg/cell)
sp.c_c  = 3;     % cytosolic nanoparticle-mRNA complex concentration (pg/cell)
sp.m_c	= 4;     % cytosolic free mRNA number (#/cell)
sp.p_c	= 5;     % cytosolic protein number (#/cell)

%% Initial conditions

% initialize all species to 0
y0  = zeros(length(fields(sp)), 1);

y0(sp.c_x)  = 620;  % initialize extracellular nanoparticle amount to 620 pg/cell

%% Parameters - specify base values

% import parameters file to define p0
parameters;


switch flag
    
    %% SIMULATION 1 - initial parameter values, 6h
    case 1
        
        % simulate for 6h
        sim_length  = 6*60;  % minutes
        [T, Y]      = main_ode([0 sim_length], y0, sp, p0);  % run ODEs

        % calculate maximum concentrations of each species
        cmax_c_x    = max(Y(:, sp.c_x));
        cmax_c_e    = max(Y(:, sp.c_e));
        cmax_c_c    = max(Y(:, sp.c_c));
        cmax_m_c    = max(Y(:, sp.m_c));
        cmax_p_c    = max(Y(:, sp.p_c));    % max protein conc is a logical thing to look at for sensitivity analysis
        
        figures;

    %% SIMULATION 2 - initial parameter values, 72h
    case 2

        % simulate for 72h
        sim_length  = 72*60;  % minutes
        [T, Y]      = main_ode([0 sim_length], y0, sp, p); % run ODES

        % calculate maximum concentrations of each species
        cmax_c_x    = max(Y(:, sp.c_x));
        cmax_c_e    = max(Y(:, sp.c_e));
        cmax_c_c    = max(Y(:, sp.c_c));
        cmax_m_c    = max(Y(:, sp.m_c));
        cmax_p_c    = max(Y(:, sp.p_c));    % max protein conc is a logical thing to look at for sensitivity analysis
        
        figures;
        
    %% SIMULATION 3 - figure out local sensitivity analysis
    case 3
end