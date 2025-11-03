clear; clc; close all; 
%{
Magnetorquer Sizing Optimization Program
Written by Aidan Moriarty
Assisted by Elias Dahl (roommate) and John Babineau (MOPS lead)
%}

% -------------------------------------------------------------------------

%{ 
“The Schleswig-Holstein question is so complicated, only three men in 
Europe have ever understood it. One was Prince Albert, who is dead. The 
second was a German professor who became mad. I am the third and I have 
forgotten all about it.” 
    - British Prime Minister Lord Palmerston
%}

% -------------------------------------------------------------------------

%{
TODO
    11. Machinable parameters
    12. add heat to power failure criteria 
    13. add in more choices for materials
    14. change outputs to optimize design rather than select largest
    magnetic moment. 
%}

% -------------------------------------------------------------------------

% Information taken from table of AWG standard data. 

AWG = (-3:47)';  
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

%--------------------------------------------------------------------------

% Get user inputs
disp('Please specify a set of constraints.'); 
disp('For default constraint, enter 0'); 
prompt = "What is the maximum core radius? [m] ";
txt = input(prompt); 
if txt == 0
    max_r = 0.01; 
else
    max_r = txt; 
end 
prompt = "What is the length? [m] "; 
txt = input(prompt); 
if txt == 0
    l_core = 0.2; 
else
    l_core = txt; 
end 
prompt = "What is the maximum volume? [m^3] "; 
txt = input(prompt); 
if txt == 0
    params.v_max = pi*(max_r*2)^2*l_core; 
else
    params.v_max = txt; 
end 
prompt = "What is the maximum mass? [kg] ";
txt = input(prompt); 
if txt == 0
    params.m_max = 1; 
else
    params.m_max = txt; 
end 
prompt = "What is the supplied bus voltage? [V] "; 
txt = input(prompt); 
if txt == 0
    params.V_bus = 3; 
else
    params.V_bus = txt; 
end 
prompt = "What is the maximum current? [A] "; 
txt = input(prompt); 
if txt == 0
    params.I_max = 0.25; % A
else
    params.I_max = txt; 
end 

%--------------------------------------------------------------------------
% Parameters

params.res_cu = 1.68e-8;   % ohm*m @ 20°C
params.rho_cu = 8960;      % kg/m^3

min_r = 0.005; 
inc_r = 0.0005; 
r_values = min_r:inc_r:max_r; % range of radii

coil_template = struct( ...
    'd_wire', 0, ...
    'l_core', 0, ...
    'n_wraps', 0, ...
    'r_outer', 0, ...
    'res_per_m', 0, ...
    'm_total', 0, ...
    'm_coil', 0, ...
    'res_total', 0, ...
    'current', 0, ...
    'M_dipole', 0, ...
    'M_9', 0, ...
    'v', 0, ...
    'fail_code', 0, ...
    'power', 0, ...
    'awg', 0 ...
);

core_template = struct( ...
    'name', '', ...
    'rho', [], ...
    'mass', [], ...
    'mu_r', [] ...
);

alloys = ["Vacoflux 50", "something"]; 

coil = repmat(coil_template, length(r_values), height(AWG_Table)); % grid of coil templates

sz = [1 12]; 
varTypes = ["double", "double", "string", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double"]; 
varNames = ["Magnetic Moment [A*m^2]", "Air Cored Magnetic Moment [A*m^2]", ...
    "Alloy", "Wire Gauge", "Core Length [m]", "Core Radius [m]", ...
    "Outer Radius [m]", "Volume [m^3]", "Number of Wraps", "Current [A]", ...
    "Power [W]", "Total Mass [kg]"];
results = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 

n = 1; 
good = 0; 
vol_fail = 0; 
mass_fail = 0; 
curr_fail = 0; 


%{
Loops through constraints and saves results to results. 
%}
for k = 1:numel(alloys)
    name = alloys(k); 
    properties = getProperties(name); 
    core = core_template; 
    core.name = name; 
    core.mu_r = properties.mu_r; 
    core.rho = properties.rho; 
    for i = 1:length(r_values) % switch to being maximum length and iterate over the radius ( 8 cm ) 
        r_core = r_values(i); 
        core.mass = core.rho*pi*r_core^2*l_core;
        for j = 1:height(AWG_Table)
            d_wire = (AWG_Table.d_mm(j))/1000; 
            res_per_m = (AWG_Table.r_km(j))/1000; 
            % initialize coil 
            c = coil_template; 
            c.l_core = l_core; 
            c.d_wire = d_wire; 
            c.res_per_m = res_per_m; 
            c.r_outer = r_core; 
            c.awg = AWG_Table.AWG(j); 
            % wrap until breaks constraints
            while true
                [c, wrapping] = wrapNext(core, c, params); 
                if ~wrapping
                    break; 
                end 
            end
            % calculate magnetic dipole 
            N_total = floor(l_core / d_wire) * c.n_wraps; 
            r_mean = (r_core + c.r_outer) / 2; 
            A_loop = pi*r_mean^2; 
            c.M_dipole = N_total * c.current * A_loop; 
    
            % equation 9
            Nd = demag(l_core, r_core); 
            gain = 1+(core.mu_r - 1)/(1 + (core.mu_r)*Nd); 
            A_core = pi * r_core^2; 
            c.M_9 = N_total*c.current*A_core*gain; 
            c.power = params.V_bus * c.current; 
            if c.M_9 > 0
                results(n,:) = {c.M_9, c.M_dipole, core.name, c.awg, l_core, r_core, c.r_outer, c.v, c.n_wraps, c.current, c.power, c.m_total}; 
                n = n+1; 
                good = good + 1; 
            end
            if c.fail_code == 1
                vol_fail = vol_fail + 1; 
            elseif c.fail_code ==2
                mass_fail = mass_fail + 1; 
            elseif c.fail_code == 3
                curr_fail = curr_fail + 1; 
            end 
        end 
    end 
end 

% basic plots to be optimized in future release

figure; 
scatter(results, "Power [W]", "Magnetic Moment [A*m^2]")

figure; 
scatter(results, "Total Mass [kg]", "Magnetic Moment [A*m^2]") 

figure; 
scatter(results,"Volume [m^3]", "Magnetic Moment [A*m^2]")

figure; 
x = ["Usable Designs" "Volume Fail" "Mass Fail" "Current Fail"]; 
goobers = [good vol_fail mass_fail curr_fail]; 
bar(x, goobers); 

%{
Checks if another layer of wire can be added and if so, wraps it. 
%}
function [coil, wrapping] = wrapNext(core, coil, params)
    d_wire = coil.d_wire;
    n_turns = floor(coil.l_core / d_wire);

    % geometry for new layer
    r_inner = coil.r_outer; 
    r_outer = r_inner + d_wire; 
    r_mid = 0.5*(r_inner+r_outer); 

    % properties of new layer
    L_layer = 2*pi*r_mid*n_turns; 
    A_wire = pi * (d_wire^2)/4;
    R_layer = coil.res_per_m * L_layer; 
    m_layer = params.rho_cu * L_layer * A_wire; 

    % new totals
    R_new = coil.res_total + R_layer; 
    I_new = params.V_bus / R_new;
    m_coil_new = coil.m_coil + m_layer;  
    m_total_new = core.mass + m_coil_new; 
    % P_new = params.V_bus * I_new; 
    v_new = pi*(r_outer^2)*coil.l_core; 

    % check constraints
    if v_new > params.v_max
        coil.fail_code = 1; % volume fail
        wrapping = false; 
        return
    elseif m_total_new > params.m_max
        coil.fail_code = 2; % mass fail 
        wrapping = false; 
        return 
    elseif I_new > params.I_max
        coil.fail_code = 3; % power fail 
        wrapping = false; 
        return; 
    end 

    % wrap new layer
    coil.v = v_new; 
    coil.n_wraps = coil.n_wraps + 1; 
    coil.res_total = R_new; 
    coil.m_total = m_total_new; 
    coil.m_coil = m_coil_new; 
    coil.r_outer = r_outer; 
    coil.current = params.V_bus / coil.res_total; 
    wrapping = true; 
end

% Demagnetizing calculation
function Nd = demag(l_core, r_core)
    x = l_core/r_core; 
    Nd = 4*(log(x) - 1) / (x^2 - 4*log(x));
end

% Properties of different alloy core choices. 
function props = getProperties(name)
    switch(name)
        case 'Vacoflux 50'
            props.mu_r = 5000; 
            props.rho = 8120; % kg/m^3
        case 'something' % TODO
            props.mu_r = 0; 
            props.rho = 0;
        otherwise 
            error('Unknown alloy: %s', name)
    end 
end 