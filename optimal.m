clear; clc; close all; 

%{ 
“The Schleswig-Holstein question is so complicated, only three men in 
Europe have ever understood it. One was Prince Albert, who is dead. The 
second was a German professor who became mad. I am the third and I have 
forgotten all about it.” 
    - British Prime Minister Lord Palmerston
%}


%{
TODO
1. Add ferromagnetic alloy - need to iterate over rod sizes
    a. Iron-cobalt
    b. nickel-iron
    c. permalloy (78% nickel, 22% iron)
    d. permendur (50% cobalt, 50% iron)
    e. CK30
        Needs to have a linear relationship between dipole moment and
        current and small dipole moment produced when torquers are off. 
2. Make cooler plots
3. AWG Wires
    b. resistance per km 
%}


AWG = (-3:47)';  % numeric mapping

d_mm = [ ...
11.6840, 10.40384, 9.26592, 8.25246, 7.34822, 6.54304, 5.82676, 5.18922, ...
4.62126, 4.11582, 3.66522, 3.26390, 2.90576, 2.58826, 2.30378, 2.05232, ...
1.82880, 1.62814, 1.45034, 1.29032, 1.15062, 1.02362, 0.91186, 0.81280, ...
0.72390, 0.64516, 0.57404, 0.51054, 0.45466, 0.40386, 0.36068, 0.32004, ...
0.28702, 0.254, 0.22606, 0.2032, 0.18034, 0.16002, 0.14224, 0.12700, ...
0.11430, 0.10160, 0.08890, 0.07874, 0.07112, 0.06350, 0.05588, 0.05080, ... 
0.04470, 0.03988, 0.03556];

r_km = [ ...
0.160720, 0.202704, 0.255512, 0.322424, 0.406392, 0.512664, 0.646160, ... 
0.815080, 1.027624, 1.295928, 1.634096, 2.060496, 2.598088, 3.276392, ...
4.132800, 5.208640, 6.569840, 8.282000, 10.44352, 13.17248, 16.60992, ...
20.94280, 26.40728, 33.29200, 41.98400, 52.93920, 66.78080, 84.19760, ...
106.1736, 133.8568, 168.8216, 212.8720, 268.4024, 338.4960, 426.7280, ...
538.2480, 678.6320, 855.7520, 1079.120, 1360, 1715, 2163, 2728, 3442, ...
4341, 5443, 7031, 8507, 10984, 13802, 17359]; 

AWG_Table = table(AWG, d_mm', r_km', ...
    'VariableNames', {'AWG','d_mm','r_km'});

min_l = 0.01; 
max_l = 0.1; 
inc_l = 0.01; 

params.res_cu = 1.68e-8;   % ohm·m @ 20°C
params.rho_cu = 8960;      % kg/m^3ch
params.r_core = 0.01;      % m
params.m_core = 0.1;       % kg
params.mu_r = 2000;        % something chat gave me
params.m_max  = 0.5;       % kg (max coil + core)
params.r_max  = 0.05;      % m (max outer radius)
params.P_max  = 2;         % W
params.V_bus  = 3;         % V

 

l_values = min_l:inc_l:max_l; % range of lengths
% d_values = min_d:inc_d:max_d; % range of diameters

coil_template = struct( ...
    'd_wire', 0, ...
    'l_core', 0, ...
    'n_wraps', 0, ...
    'r_outer', params.r_core, ...
    'm_total', params.m_core, ...
    'res_per_m', 0, ...
    'res_total', 0, ...
    'current', 0, ...
    'M_dipole', 0, ...
    'M_9', 0);


coil = repmat(coil_template, length(l_values), height(AWG_Table)); % grid of coil templates
M = zeros(length(l_values), height(AWG_Table));
M_9 = zeros(length(l_values), height(AWG_Table)); 
P = zeros(length(l_values), height(AWG_Table)); 
m = zeros(length(l_values), height(AWG_Table)); 
r = zeros(length(l_values), height(AWG_Table)); 

for i = 1:length(l_values)
    l_core = l_values(i); 
    for j = 1:height(AWG_Table)
        d_wire = AWG_Table.d_mm(j)/1000; 
        res_per_m = AWG_Table.r_km(j)/1000; 
        % initialize coil 
        c = coil_template; 
        c.l_core = l_core; 
        c.d_wire = d_wire; 
        c.res_per_m = res_per_m; 
        % wrap until breaks constraints
        while true
            [c, wrapping] = wrapNext(c, params); 
            if ~wrapping
                break; 
            end 
        end
        % calculate magnetic dipole 
        N_total = floor(l_core / d_wire) * c.n_wraps; 
        r_mean = (params.r_core + c.r_outer) / 2; 
        A_loop = pi*r_mean^2; 
        c.M_dipole = N_total * c.current * A_loop; 

        % equation 9
        r_core = params.r_core; 
        Nd = demag(l_core, r_core); 
        gain = 1+(params.mu_r - 1)/(1 + (params.mu_r)*Nd); 
        A_core = pi * r_core^2; 
        c.M_9 = N_total*c.current*A_core*gain; 

        % store results
        coil(i,j) = c; 
        M(i,j) = c.M_dipole; % magnetorquer dipole 
        M_9(i,j) = c.M_9;
        P(i,j) = c.current*params.V_bus; 
        m(i,j) = c.m_total; 
        r(i,j) = params.r_core + c.r_outer; 

    end 
end 

% plotting

[DD, LL] = meshgrid(AWG, l_values);
M_plot = M_9;
M_plot(M_plot <= 0) = NaN; 
M_air = M; 
M_air(M_air <= 0) = NaN; 

figure; 
mesh(DD, LL, M_plot);
xlabel('Wire Diameter [m]');
ylabel('Core Length [m]');
zlabel('Magnetic Dipole Moment [A·m^2]');
title('Magnetorquer Design Space');
shading interp; colorbar; grid on;

% figure;
% P_plot = P; 
% scatter(P_plot, M_plot); 
% xlabel('Power [W]'); 
% ylabel('Magnetic Dipole Moment [A*m^2]'); 
% grid on; 
% 
% figure; 
% scatter(DD, M_plot); 
% xlabel('Diameter [m]'); 
% ylabel('Magnetic Dipole Moment [A*m^2]'); 
% grid on; 
% 
% figure; 
% scatter(m, M_plot); 
% xlabel('mass [kg]'); 
% ylabel('Magnetic Dipole Moment [A*m^2]'); 
% grid on; 
% 
% figure; 
% scatter(DD, r, 'o');
% set(gca, 'XScale', 'log');
% xlabel('Diameter [m]');
% ylabel('Radius [m]');
% grid on;






% helper functions

function [coil, wrapping] = wrapNext(coil, params)
    d_wire = coil.d_wire;
    n_turns = floor(coil.l_core / d_wire);

    % geometry for new layer
    r_inner = coil.r_outer; 
    r_outer = r_inner + d_wire; 
    r_mid = 0.5*(r_inner+r_outer); 

    % properties of new layer
    L_layer = 2*pi*r_mid*n_turns; 
    A_wire = pi * (d_wire^2)/4;
    R_layer = coil.res_per_m * L_layer / A_wire; % need to switch to specific AWG resistance
    m_layer = params.rho_cu * L_layer * A_wire; 

    % new totals
    R_new = coil.res_total + R_layer; 
    m_new = coil.m_total + m_layer; 
    P_new = (params.V_bus^2)/R_new; 

    % check constraints
    if (r_outer > params.r_max) || (m_new > params.m_max) || (P_new > params.P_max)
        wrapping = false;
        return
    end

    % wrap new layer
    coil.n_wraps = coil.n_wraps + 1; 
    coil.res_total = R_new; 
    coil.m_total = m_new; 
    coil.r_outer = r_outer; 
    coil.current = params.V_bus / coil.res_total; 
    wrapping = true; 
end

function Nd = demag(l_core, r_core)
    x = l_core/r_core; 
    Nd = 4*(log(x) - 1) / (x^2 - 4*log(x));
    Nd = min(max(Nd, 0), 1);
end