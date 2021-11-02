clear all;
close all;

%% Specify simulation type

flag = 4;
% 1 = simple simulation with basic parameter values for 6h
% 2 = simple simulation with basic parameter values for 72h
% 3 = local sensitivity analysis of max protein amount to 10% change
% 4 = local sensitivity analysis of max protein amount to various % change

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
        
        % announce program start
        tic
        fprintf('Simulating scenario %d...\n', flag);
        
        % simulate for 6h
        sim_mins  = 6*60;  % minutes simulated
        [T, Y, cmax, tmax]      = main_ode([0 sim_mins], y0, sp, p0);  % run ODEs
        
        % create figures
        fprintf('Creating figures...\n');
        figures;
        
        % announce program end
        fprintf('Scenario %d complete!\n', flag);
        toc

    %% SIMULATION 2 - initial parameter values, 72h
    case 2

        % announce program start
        tic
        fprintf('Simulating scenario %d...\n', flag);
        
        % simulate for 72h
        sim_days    = 3;	% days simulated 
        sim_mins    = sim_days * 24 *60;	% minutes simulated
        [T, Y, cmax, tmax]      = main_ode([0 sim_mins], y0, sp, p0);   % run ODES
        
        % create figures
        fprintf('Creating figures...\n');
        figures;
        
        % announce program end
        fprintf('Scenario %d complete!\n', flag);
        toc
        
    %% SIMULATION 3 - local sensitivity analysis of max protein amount to 10% change
    case 3
        
        % announce program start
        tic
        fprintf('Simulating scenario %d...\n', flag);
        
        % set variables
        delta = 0.1;  % test sensitivity to a 10% increase in parameter value
        sim_days    = 3;	% days simulated 
        sim_mins    = sim_days * 24 *60;	% minutes simulated
        p0_names = fieldnames(p0);  % extract parameter names
        sens_rel_norm   = zeros(length(p0_names), 1);
        
        % simulate base case
        [T0, Y0, cmax0, tmax0] = main_ode([0 sim_mins], y0, sp, p0);
        p_c0    = cmax0.p_c;    % max protein amount
        t_p_c0 = tmax0.p_c;     % time of max protein amount
        
        % perform sensitivity analysis for all parameters
        for i = 1:length(p0_names)
            
            fprintf('Analyzing parameter %s...\n', p0_names{i});
            
            % simulate test case
            p = p0;     % set initial parameter values
            p.(p0_names{i}) = p0.(p0_names{i}) * (1 + delta);	% increase ith parameter value by fraction delta
            [T1, Y1, cmax1, tmax1] = main_ode([0 sim_mins], y0, sp, p);   % run ODEs
            p_c = cmax1.p_c;    % new max protein amount
            t_p_c = tmax1.p_c;

            % calculate and store sensitivity
            sens_rel_not_norm   = (p_c - p_c0)/(p_c0);  % relative sensitivity (not normalized to parameter value)
            param_change	= (p.(p0_names{i}) - p0.(p0_names{i})) / p.(p0_names{i});   % magnitude of parameter change
            sens_rel_norm(i)   = sens_rel_not_norm / param_change;     % normalized relative sensitivity

        end
        
        % create figures
        fprintf('Creating figures...\n');
        figures;
        
        % announce program end
        fprintf('Scenario %d complete!\n', flag);
        toc
        
   %% SIMULATION 4 - local sensitivity analysis of max protein amount to various % change
    case 4
        
        % announce program start
        tic
        fprintf('Simulating scenario %d...\n', flag);
        
        % set variables
        test_deltas     = [.01 .05 .1 .25 .5 1];    % define deltas (% change in parameter) to test
        sim_days    = 3;	% days simulated 
        sim_mins    = sim_days * 24 *60;	% minutes simulated
        p0_names    = fieldnames(p0);               % extract parameter names
        sens_rel_norm   = zeros(length(p0_names), 1);
        combined_sens   = zeros(length(test_deltas), length(p0_names));
        
        % simulate base case
        [T0, Y0, cmax0, tmax0]     = main_ode([0 sim_mins], y0, sp, p0);
        p_c0    = cmax0.p_c;    % max protein amount in base case
        t_p_c0  = tmax0.p_c;    % time when max protein amount is reached
        
        % perform local univariate sensitivity analysis for each delta
        for j = 1:length(test_deltas)
           delta = test_deltas(j);
           fprintf('Exploring %d percent increase...\n', delta*100);
           
           % perform sensitivity analysis for all parameters
            for i = 1:length(p0_names)

                % simulate test case
                p = p0;     % set initial parameter values
                p.(p0_names{i}) = p0.(p0_names{i}) * (1 + delta);	% increase ith parameter value by fraction delta
                [T1, Y1, cmax1, tmax1] = main_ode([0 sim_mins], y0, sp, p);   % run ODEs
                p_c = cmax1.p_c;    % new max protein amount
                t_p_c = tmax1.p_c;  % new time when max protein amount reached

                % calculate sensitivity
                sens_rel_not_norm   = (p_c - p_c0)/(p_c0);  % relative sensitivity (not normalized to parameter value)
                param_change	= (p.(p0_names{i}) - p0.(p0_names{i})) / p.(p0_names{i});   % magnitude of parameter change
                sens_rel_norm(i)   = sens_rel_not_norm / param_change;     % normalized relative sensitivity

            end

            % specify parameter labels in initial order
            names_array = {'k_{deg_x}', 'k_{int}', 'k_{lys}', 'k_{escape}', ...
                'k_{deg_{np}}', 'k_{deg_m}', 'k_{expr}', 'k_{deg_p}'};
            
            % add sensitivity for this delta to the combined matrix
            combined_sens(j, :) = sens_rel_norm;
           
        end
        
        % create figures
        fprintf('Creating figures...\n');
        figures;
        
        % announce program end
        fprintf('Scenario %d complete!\n', flag);
        toc
        
    %% SIMULATION 5 - global univariate sensitivity analysis of max protein amount
    case 5
        
        % announce program start
        tic
        fprintf('Simulating scenario %d...\n', flag);
        
        % define variables
        exps = -3:3;                        % range of exponents to scan for each parameter
        outvars = ["Cmax_p", "Tmax_p"];     % output features to analyze
        sim_days = 15;                      % days simulated (15 req'd to reach Tmax in all conditions)
        sim_mins  = sim_days * 24 * 60;     % minutes simulated
        p0_names = fieldnames(p0);          % extract parameter names
        n_params = length(p0_names);        % number of parameters to test
        n_exps = length(exps);              % number of exponents per parameter
        n_outs = length(outvars);           % number of output features to analyze    
       
        % initialize empty output matrix m
        % format m(j,x,y) = m(parameter, exponent, output)
        % m = matrix of relative sensitivities, not normalized to param change size
        for y = 1:n_outs
            m(:, :, y) = zeros(n_params, n_exps);
        end
        
        % initialize empty output matrix m_norm with same dimensions
        % m_norm = matrix of relative sensitivities, normalized to size of parameter change
        m_norm = m;
        
        % run sensitivity analysis for every parameter across a log scale range
        for j = 1:1:n_params
            
            fprintf('Running simulation %d of %d', j, n_params)
        
            p = p0;     % set initial parameter values
            
            for x = 1:n_exps
                p.(p0_names{j}) = p0.(p0_names{j}) * 10^exps(x);    % scale param j by 10^exps(x)
                [T1, Y1, cmax, tmax] = main_ode([0 sim_mins], y0, sp, p);   % run ODEs
                outs = [cmax.p_c, tmax.p_c];    % return Cmax and Tmax of cytosolic protein (p_c)
                
                % store sensitivity of each output variable in m
                for y = 1:n_outs
                    m(j, x, y) = outs(y);
                end
                
                fprintf('.')
            end
           
            fprintf('\n')
            
        end
        
        %% TO DO: calculate relative sensitivity, but do it inside the loop
        % sensitivity relative to base case
        % normalized relative sensitivity
        
        % CODE COPYPASTA FROM 4
        % calculate sensitivity
        %        sens_rel_not_norm   = (p_c - p_c0)/(p_c0);  % relative sensitivity (not normalized to parameter value)
        %        param_change	= (p.(p0_names{i}) - p0.(p0_names{i})) / p.(p0_names{i});   % magnitude of parameter change
        %        sens_rel_norm(i)   = sens_rel_not_norm / param_change;     % normalized relative sensitivity
        
        % CODE SANDBOX FOR CALCULATING SENSITIVITIES
        % normalize output matrix to base case
        %sens_not_norm = m;
        %for y = 1:n_outs
        %    sens_not_norm(:, :, y) = ( m(:, :, y) - outs0(y) ) ./ outs0(y);
        %end
        
        % create figures
        fprintf('Creating figures...\n');
        figures;
        
        % announce program end
        fprintf('Scenario %d complete!\n', flag);
        toc
        
        % show warning if any sim hasn't reached Tmax
        if ismember(sim_mins, m)
            warning('Time is too short to reach Tmax - increase simulation length.')
        end
end