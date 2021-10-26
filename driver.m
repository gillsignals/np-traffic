clear all;
close all;

%% Specify simulation type
flag = 1;
% 1 = simple simulation with basic parameter values for 6h
% 2 = simple simulation with basic parameter values for 72h
% 3 = figure out sensitivity analysis

%% Import colorblind-friendly colors for plotting
colorblind_colors;

%% Assign molecular species names to index values

sp.c_x	= 1;     % c_x = extracellular nanoparticle-mRNA complex (pg/cell)
sp.c_e	= 2;     % c_e = endosomal nanoparticle-mRNA complex concentration (pg/cell)
sp.c_c  = 3;     % c_c = cytosolic nanoparticle-mRNA complex concentration (pg/cell)
sp.m_c	= 4;     % m_c = cytosolic free mRNA number (#/cell)
sp.p_c	= 5;     % p_c = cytosolic protein number (#/cell)

%% Initial conditions

% initialize all species to 0
y0  = zeros(length(fields(sp)), 1);

y0(sp.c_x)  = 620;  % initialize extracellular nanoparticle amount to 620 pg/cell

% import parameters file to define p0
parameters;


%% Select and run simulation

switch flag
    
    %% SIMULATION 1 - initial parameter values, 6h
    case 1
        
        % simulate for 6h
        sim_length  = 6*60;  % minutes
        [T, Y, cmax]      = main_ode([0 sim_length], y0, sp, p0);  % run ODEs
        
        figures;

    %% SIMULATION 2 - initial parameter values, 72h
    case 2

        % simulate for 72h
        sim_length  = 72*60;  % minutes
        [T, Y, cmax]      = main_ode([0 sim_length], y0, sp, p0); % run ODES
        
        figures;
        
    %% SIMULATION 3 - figure out local sensitivity analysis
    case 3
        
        
end