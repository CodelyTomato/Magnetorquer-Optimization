clear; clc; close all; 
%{
Magnetorquer Sizing Optimization Program
Written by Aidan Moriarty
Assisted by Elias Dahl and John Babineau 
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
    1. Add alloys 
        - magnetic permeability will change with grade, temperature, and
        thickness
    2. add comments
    3. add volume constraint
    4. add inputs
    5. thermals 
    6. power output 
    7. fail rationale 
    8. saturation cap 
    9. computation of inductance 
    10. wire specific currents? 
%}

% -------------------------------------------------------------------------

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

params.res_cu = 1.68e-8;   % ohm*m @ 20°C
params.rho_cu = 8960;      % kg/m^3ch
params.m_max  = 1;         % kg (max coil + core)
params.r_max  = 0.05;      % m (max outer radius)
params.I_max  = 0.25;         % W
params.V_bus  = 3;         % V
params.v_max = 1e-5;       % m^3

l_values = min_l:inc_l:max_l; % range of lengths
r_values = linspace(0.001, 0.01, 10); % range of radius 

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
    'fail_code', 0 ...
);

core_template = struct( ...
    'name', '', ...
    'rho', [], ...
    'mass', [], ...
    'mu_r', [] ...
);

alloys = ["Vacoflux 50", "something"]; 

coil = repmat(coil_template, length(l_values), height(AWG_Table)); % grid of coil templates

M = zeros(length(l_values), height(AWG_Table), numel(alloys));
M_9 = zeros(length(l_values), height(AWG_Table), numel(alloys)); 
P = zeros(length(l_values), height(AWG_Table), numel(alloys)); 
m = zeros(length(l_values), height(AWG_Table), numel(alloys)); 
r = zeros(length(l_values), height(AWG_Table), numel(alloys)); 
d = zeros(length(l_values), height(AWG_Table), numel(alloys)); 

results = table([], [], [], [], [], [], 'VariableNames', {'Alloy','m_total','M_9', 'M_dipole', 'r_core', 'l_core'});

n = 1; 

for k = 1:numel(alloys)
    name = alloys(k); 
    properties = getProperties(name); 
    core = core_template; 
    core.name = name; 
    core.mu_r = properties.mu_r; 
    core.rho = properties.rho; 
    for l = 1:length(r_values)
        r_core = r_values(l); 

        for i = 1:length(l_values)
            l_core = l_values(i); 
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
        
                % store results
                if c.M_9 > 0
                    row = table(string(name), c.m_total, c.M_9, c.M_dipole, r_core, l_core, ...
                        'VariableNames', {'Alloy', 'm_total', 'M_9', 'M_dipole', 'r_core', 'l_core'});
                    results = [results; row];
                end 
            end 
        end 
    end 
end 

% alloyNames = unique(results.Alloy); 
% colors = lines(numel(alloyNames)); 
% 
% figure; hold on; grid on;
% for a = 1:numel(alloyNames)
%     mask = results.Alloy == alloyNames(a);
%     scatter(results.m_total(mask), results.M_9(mask), ...
%         40, colors(a,:), 'filled', 'DisplayName', alloyNames(a));
% end
% xlabel('Mass [kg]');
% ylabel('Magnetic Dipole [A·m^2]');
% legend show;
% title('Magnetic Dipole vs Mass by Alloy');


alloyNames = unique(results.Alloy);
colors = lines(numel(alloyNames));

figure; hold on; grid on; box on;
for a = 1:numel(alloyNames)
    mask = results.Alloy == alloyNames(a);
    scatter3( results.r_core(mask)*1e3, ...   % convert to mm
              results.l_core(mask)*1e3, ...
              results.M_9(mask), ...
              40, colors(a,:), 'filled', ...
              'DisplayName', alloyNames(a) );
end
xlabel('Core Radius [mm]');
ylabel('Core Length [mm]');
zlabel('Magnetic Dipole Moment, M_9 [A·m^2]');
title('Magnetorquer Design Space: Core Geometry vs Dipole Moment');
legend('Location','best');
view(45,30);




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
    m_coil_new = coil.m_coil + m_layer;  
    m_total_new = core.mass + m_coil_new; 
    I_new = params.V_bus/R_new; 
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

%{
fucntion does thing
what parameters are
what it retuns
%}
function Nd = demag(l_core, r_core)
    x = l_core/r_core; 
    Nd = 4*(log(x) - 1) / (x^2 - 4*log(x));
    % Nd = max(min(Nd, 0.9999), 1e-4); - not sure if this is necessary? 
end

function props = getProperties(name)
%GETPROPERTIES Returns properties of alloys
    switch(name)
        case 'Vacoflux 50'
            props.mu_r = 5000; 
            props.rho = 8120; % kg/m^3
        case 'something'
            props.mu_r = 1; 
            props.rho = 1; 
        otherwise 
            error('Unknown alloy: %s', name)
    end 
end 