switch flag     % specified in driver
    % 1 = simple simulation with basic parameter values for 6h
    % 2 = simple simulation with basic parameter values for 72h
    % 3 = local sensitivity analysis of max protein conc to 10% change
    
    %% CASE 1 = simple simulation with basic parameter values for 6h
    case 1

        % create figure for simulation 1
        figure;

        % plot C_x for sim 1
        ax1 = subplot(3,2,1);
        plot(ax1, T / 60, Y(:, sp.c_x), 'LineWidth', 2, 'Color', colors(cb.black))
        title("Extracellular NP-mRNA complex (C_x)")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("C_x (pg/cell)")

        ax2 = subplot(3,2,3);
        plot(ax2, T / 60, Y(:, sp.c_e), 'LineWidth', 2, 'Color', colors(cb.red))
        title("Endosomal NP-mRNA complex (C_e)")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("C_e (pg/cell)")

        ax3 = subplot(3,2,5);
        plot(ax3, T / 60, Y(:, sp.c_c), 'LineWidth', 2, 'Color', colors(cb.skyblue))
        title("Cytosolic NP-mRNA complex (C_c)")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("C_c (pg/cell)")

        ax4 = subplot(3,2,2);
        plot(ax4, T / 60, Y(:, sp.m_c), 'LineWidth', 2, 'Color', colors(cb.orange))
        title("Cytosolic free mRNA (M_c)")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("M_c (#/cell)")

        ax5 = subplot(3,2,4);
        plot(ax5, T / 60, Y(:, sp.p_c), 'LineWidth', 2, 'Color', colors(cb.green))
        title("Cytosolic protein (P_c)")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("P_c (#/cell)")

        ax6 = subplot(3,2,6);
        plot(ax6, T / 60, Y(:, sp.c_x) / cmax.c_x, 'LineWidth', 2, 'Color', colors(cb.black))
        hold on;
        plot(ax6, T / 60, Y(:, sp.c_e) / cmax.c_e, 'LineWidth', 2, 'Color', colors(cb.red))
        plot(ax6, T / 60, Y(:, sp.c_c) / cmax.c_c, 'LineWidth', 2, 'Color', colors(cb.skyblue))
        plot(ax6, T / 60, Y(:, sp.m_c) / cmax.m_c, 'LineWidth', 2, 'Color', colors(cb.orange))
        plot(ax6, T / 60, Y(:, sp.p_c) / cmax.p_c, 'LineWidth', 2, 'Color', colors(cb.green))
        hold off;
        title("Combined scaled time course")
        xlabel("Time (h)")
        xlim([0 6]);
        xticks(linspace(0,6,7));
        ylabel("Fraction of maximum amount")
        legend(["C_x", "C_e", "C_c", "M_c", "P_c"], "Location", "eastoutside")
        
    %% CASE 2 = simple simulation with basic parameter values for 72h
    case 2

        % create figure for simulation 2
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
        plot(ax6, T / 60, Y(:, sp.c_x) / cmax.c_x, 'LineWidth', 2, 'Color', colors(cb.black))
        hold on;
        plot(ax6, T / 60, Y(:, sp.c_e) / cmax.c_e, 'LineWidth', 2, 'Color', colors(cb.red))
        plot(ax6, T / 60, Y(:, sp.c_c) / cmax.c_c, 'LineWidth', 2, 'Color', colors(cb.skyblue))
        plot(ax6, T / 60, Y(:, sp.m_c) / cmax.m_c, 'LineWidth', 2, 'Color', colors(cb.orange))
        plot(ax6, T / 60, Y(:, sp.p_c) / cmax.p_c, 'LineWidth', 2, 'Color', colors(cb.green))
        hold off;
        title("Combined scaled time course")
        xlabel("Time (h)")
        xlim([0 72]);
        xticks(linspace(0,72,7));
        ylabel("Fraction of maximum amount")
        legend(["C_x", "C_e", "C_c", "M_c", "P_c"], "Location", "eastoutside")
        
	%% CASE 3 = local sensitivity analysis of max protein conc to 10% change
    
    case 3

        % specify parameter labels in initial order
        names_array = {'k_{deg_x}', 'k_{int}', 'k_{lys}', 'k_{escape}', ...
            'k_{deg_{np}}', 'k_{deg_m}', 'k_{expr}', 'k_{deg_p}'};
        
        % find indices to sort sensitivity vector in descending order
        [~, ind] = sort(sens_rel_norm, 'ascend');        
             
        % sort parameter names
        names_array_sorted = names_array(ind);
               
        % horizontal, sorted version of bar graph
        figure;
        barh(sens_rel_norm(ind));
        ylabel("Parameter")
        xlabel("Sensitivity") % % change in max(P_c) per % change in parameter
        title("Sensitivity of max(P_c) to 10% increase in parameter")
        yticklabels(names_array_sorted)
end
