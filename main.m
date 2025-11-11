clear; clc; close all; 

%{
1. jiust dont even include over 35 AWG
2. make good good plots
3. add core permeability w/ some saturation addedd
4. 
%}


% --- Define Parameters --- % 
params.dens_cu = 8960; 
params.res_cu = 1.68e-8; 
params.dens_core = 8120; 
params.r_core = 0.002; 
params.L_core = 0.09; 
params.vol_core = pi*params.r_core^2*params.L_core; 
params.r_max = 0.0125; 
params.vol_max = pi*params.r_max^2*params.L_core; 
params.mass_core = params.vol_core * params.dens_core; 
params.mass_max = 0.5; 
params.I_battery_max = 0.5;
params.V_in = 4; 

% --- Get Magnetic Wire Data --- % 
T = readtable('awg.csv', 'VariableNamingRule', 'preserve'); 
AWG = T{:,'AWG'}; 
d_con_mm = T{:,'Bare Diameter (mm)'}; 
r_km = T{:, 'Resistance / KM'}; 
max_wire_current = T{:, 'Max Current (A)'}; 
d_ins_mm = T{:, 'Single Build Diameter (mm)'};

% --- Determine M_core --- % 
mu_r = 3500; % cold-rolled 1008 low carbon steel (guess)
x = params.L_core / params.r_core; 
Nd = (4*(log(x) - 1)) / ((x)^2 - 4*log(x));
gain = 1 + (mu_r - 1) / (1 + (mu_r - 1)*Nd);
for i = 1:height(AWG)
    d_ins = d_ins_mm(i) / 1000; 
    d_con = d_con_mm(i) / 1000; 
    awg = AWG(i); 
    % - max wire current is 80 % of what max wire can handle - % 
    I_wire_max = 0.8*max_wire_current(i); 

    N = floor(params.L_core / d_ins); 
    D = wrap(params, d_ins, d_con, N, awg, I_wire_max, gain); 
    D_results{i} = D; 
end 

% --- Write all designs to big table --- % 

AllData = []; 

for i = 1:length(D_results)
    D = D_results{i}; 
     headers = {'M_core','M_air','m_total','R_new','L_new', ...
                       'A_layer','I_new','Layer','r_outer','mass_new', ...
                       'AWG','d_ins','fail_code'};
     D = array2table(D, 'VariableNames',headers); 
     AllData = [AllData; D]; 
     D_results{i} = D; 
end 

writetable(AllData, 'D_results.csv'); 



% --- plot time --- % 

figure; hold on; grid on; 
cmap = turbo(length(D_results)); 
for i = 1:length(D_results)
    D = D_results{i};
    if ~isempty(D)
        % --- Plot magnetic moment vs outer radius --- %
        plot(D.r_outer, D.M_air, 'LineWidth', 2, 'Color', cmap(i,:));

        % --- Label the line with AWG --- %
        awg = D.AWG(1);  % from table column

        x_label = D.r_outer(end);
        y_label = D.M_air(end);

        text(x_label, y_label, sprintf('AWG %d', awg), ...
             'FontWeight','bold', 'FontSize',10, ...
             'VerticalAlignment','middle', 'HorizontalAlignment','left', ...
             'Color', cmap(i,:));
    end
end

xlabel('Outer Radius [m]');
ylabel('Magnetic Moment [A·m²]');



% figure; hold on; grid on;
% 
% % Optional: colormap for visual distinction
% cmap = turbo(length(D_results));
% 
% for i = 1:length(D_results)
%     D = D_results{i};
%     if ~isempty(D)
%         % Plot line
%         plot(D.r_outer)
% 
%         % Get AWG number (column 10)
%         awg = D(1,10);
% 
%         % Pick label location (end of the line)
%         x_label = D(end,9);
%         y_label = D(end,2);
% 
%         % Add text label near the end
%         text(x_label, y_label, sprintf('AWG %d', awg), ...
%              'FontWeight','bold', 'FontSize',10, ...
%              'VerticalAlignment','middle', 'HorizontalAlignment','left', ...
%              'Color', cmap(i,:));
%     end
% end
% 
% xlabel('Outer Radius [m]');
% ylabel('Magnetic Moment [A·m^2]');
% title('Magnetic Moment vs Outer Radius by Wire Gauge');

















% figure; hold on; grid on; 
% for i = 1:length(D_results)
%     D = D_results{i}; 
%     if ~isempty(D)
%         plot(D(:,10), D(:,2), 'LineWidth', 2); 
%     end 
% end 
% xlabel('AWG'); 
% ylabel('Moment'); 
% 
% figure; hold on; grid on; 
% legends = {}; 
% cmap = turbo(length(D_results)); 
% for i = 1:length(D_results)
%     D = D_results{i}; 
%     if ~isempty(D)
%         plot(D(:,9), D(:,2), 'LineWidth',2, 'Color', cmap(i,:)); 
%         legends{end+1} = sprintf('AWG %d', D(1,10)); 
%     end 
% end
% xlabel('Outer Radius'); 
% ylabel('Magnetic Moment');
% legend(legends); 
% 
% 
% 
% figure; hold on; grid on;
% 
% % Optional: colormap for visual distinction
% cmap = turbo(length(D_results));
% 
% for i = 1:length(D_results)
%     D = D_results{i};
%     if ~isempty(D)
%         % Plot line
%         plot(D(:,9), D(:,2), 'LineWidth', 2, 'Color', cmap(i,:));
% 
%         % Get AWG number (column 10)
%         awg = D(1,10);
% 
%         % Pick label location (end of the line)
%         x_label = D(end,9);
%         y_label = D(end,2);
% 
%         % Add text label near the end
%         text(x_label, y_label, sprintf('AWG %d', awg), ...
%              'FontWeight','bold', 'FontSize',10, ...
%              'VerticalAlignment','middle', 'HorizontalAlignment','left', ...
%              'Color', cmap(i,:));
%     end
% end
% 
% xlabel('Outer Radius [m]');
% ylabel('Magnetic Moment [A·m^2]');
% title('Magnetic Moment vs Outer Radius by Wire Gauge');
% 
% 
% % fix this plot 
% figure; hold on; grid on; 
% legends = {}; 
% cmap = turbo(length(D_results)); 
% for i = 1:length(D_results)
%     D = D_results{i}; 
%     if ~isempty(D)
%         plot(D(:,8), D(:,7), 'LineWidth',2, 'Color', cmap(i,:)); 
%         legends{end+1} = sprintf('AWG %d', D(1,10)); 
%     end 
% end
% xlabel('Number of Wraps'); 
% ylabel('Current');
% legend(legends); 