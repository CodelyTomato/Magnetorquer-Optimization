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
    
%}


min_l = 0.01; 
max_l = 0.1; 
min_d = 7.874e-5; 
max_d = 0.011684; 
inc_l = 0.01; 
inc_d = 0.0001; 

params.res_cu = 1.68e-8;   % ohm·m @ 20°C
params.rho_cu = 8960;      % kg/m^3
params.r_core = 0.01;      % m
params.m_core = 0.1;       % kg
params.mu_r = 2000;        % something chat gave me
params.m_max  = 0.5;       % kg (max coil + core)
params.r_max  = 0.05;      % m (max outer radius)
params.P_max  = 2;         % W
params.V_bus  = 3;         % V

 

l_values = min_l:inc_l:max_l; % range of lengths
d_values = min_d:inc_d:max_d; % range of diameters

coil_template = struct( ...
    'd_wire', 0, ...
    'l_core', 0, ...
    'n_wraps', 0, ...
    'r_outer', params.r_core, ...
    'm_total', params.m_core, ...
    'res_total', 0, ...
    'current', 0, ...
    'M_dipole', 0, ...
    'M_9', 0);


coil = repmat(coil_template, length(l_values), length(d_values)); % grid of coil templates
M = zeros(length(l_values), length(d_values));
M_9 = zeros(length(l_values), length(d_values)); 

for i = 1:length(l_values)
    l_core = l_values(i); 
    for j = 1:length(d_values)
        d_wire = d_values(j); 
        % initialize coil 
        c = coil_template; 
        c.l_core = l_core; 
        c.d_wire = d_wire; 
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
    end 
end 

% plotting

[DD, LL] = meshgrid(d_values, l_values);

M_plot = M_9;
M_plot(M_plot <= 0) = NaN; 
figure; 
surf(DD, LL, M_plot);
xlabel('Wire Diameter [m]');
ylabel('Core Length [m]');
zlabel('Magnetic Dipole Moment [A·m^2]');
title('Magnetorquer Design Space');
shading interp; colorbar; grid on;

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
    R_layer = params.res_cu * L_layer / A_wire; 
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
    x = max(l_core / max(r_core, eps), 1+1e-6);  % avoid l<=r singularity
        Nd = 4*(log(x) - 1) / (x^2 - 4*log(x));
        Nd = min(max(Nd, 0), 1);
end