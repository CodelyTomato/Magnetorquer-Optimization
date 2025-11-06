clear; clc; close all; 
%{
Magnetorquer Sizing Optimization Program
Written by Aidan Moriarty
Assisted by Elias Dahl (roommate) and John Babineau (MOPS lead)
%}

% -------------------------------------------------------------------------

%{ 
"The Schleswig-Holstein question is so complicated, only three men in 
Europe have ever understood it. One was Prince Albert, who is dead. The 
second was a German professor who became mad. I am the third and I have 
forgotten all about it." 
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

% Get AWG information 
T = readtable('awg.csv'); 
AWG = T{:,'AWG'}; 
d_conductor_mm = T{:,'BareDiameter_mm_'}; 
r_km = T{:, 'Resistance_KM'}; 
max_wire_current = T{:, 'MaxCurrent_A_'}; 
d_jacket_mm = T{:, 'SingleBuildDiameter_mm_'}; 

% Get user inputs
disp('Please specify a set of constraints.'); 
disp('For default constraint, enter 0'); 

prompt = "What is the maximum core radius? [m] ";
txt = input(prompt); 
if txt == 0
    max_r = 0.005; 
else
    max_r = txt; 
end 

prompt = "What is the maximum outer radius? [m]" ;
txt = input(prompt); 
if txt == 0
    params.r_outer_max = 0.05; % guess
else
    params.r_outer_max = txt;
end 

prompt = "What is the length? [m] "; 
txt = input(prompt); 
if txt == 0
    l_core = 0.1; 
else
    l_core = txt; 
end 

prompt = "What is the maximum mass? [kg] ";
txt = input(prompt); 
if txt == 0
    params.m_max = 0.300; 
else
    params.m_max = txt; 
end 

prompt = "What is the supplied voltage? [V] "; 
txt = input(prompt); 
if txt == 0
    params.V_bus = 5; % V
    % disp('Using 5 V maximum supplied power constraint.')
else
    params.V_bus = txt; 
end 

%--------------------------------------------------------------------------
% Parameters

params.res_cu = 1.68e-8;   % ohm*m @ 20°C
params.rho_cu = 8960;      % kg/m^3
% 
% min_P = 0.1; % W
% inc_P = 0.05; % W - adjust? 
% P_values = min_P:inc_P:params.P_max; 

min_r = 0.0001; 
inc_r = 0.0001; 
r_values = min_r:inc_r:max_r; % range of radii

coil_template = struct( ...
    'd_conductor', 0, ...
    'd_jacket', 0, ...
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
    'I_max', 0, ...
    'sum', 0, ...
    'awg', 0 ...
);

core_template = struct( ...
    'name', '', ...
    'rho', [], ...
    'mass', [], ...
    'mu_r', [], ...
    'radius', 0 ...
);

alloys = ["4140", "something"]; 

coil = repmat(coil_template, length(r_values), height(AWG)); % grid of coil templates

sz = [1 12]; 
varTypes = ["double", "double", "double", "string", "double", "double", "double", ...
    "double", "double", "double", "double", "double"]; 
varNames = ["Design Number", "Magnetic Moment [A*m^2]", "Air Cored Magnetic Moment [A*m^2]", ...
    "Alloy", "Wire Gauge", "Core Length [m]", "Core Radius [m]", ...
    "Outer Radius [m]", "Volume [m^3]", "Number of Wraps", "Current [A]", ...
    "Total Mass [kg]"];
results = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 

n = 1; 
good = 0; 
vol_fail = 0; 
mass_fail = 0; 
curr_fail = 0; 

% -------------------------------------------------------------------------

%{
Loops through design constraints provided and saves results. Results are 
saved if they provide > 0 magnetic moment. 
%}
for k = 1:numel(alloys)

    % Initialize core parameters
    name = alloys(k); 
    properties = getProperties(name); 
    core = core_template; 
    core.name = name; 
    core.mu_r = properties.mu_r; 
    core.rho = properties.rho; 

    % for P = 1:length(P_values)
    % 
    %     power = P_values(P); 

        for i = 1:length(r_values)

            r_core = r_values(i); 
            core.mass = core.rho*pi*r_core^2*l_core;
            core.radius = r_core; 

            for j = 1:height(AWG)

                % Initialize coil parameters 
                c = coil_template; 
                % c.power = power; 
                c.l_core = l_core; 
                c.r_outer = r_core; 
                c.d_conductor = d_conductor_mm(j)/1000; 
                c.d_jacket = d_jacket_mm(j)/1000; 
                c.res_per_m = r_km(j)/1000; 
                c.awg = AWG(j); 
                c.I_max = max_wire_current(j); 

                while true
                    [c, wrapping] = wrapNext(core, c, params); 
                    if ~wrapping
                        break; 
                    end 
                end

                    

                I_final = params.V_bus / c.res_total; 
                if I_final > c.I_max
                    curr_fail = curr_fail + 1; 
                    continue
                end 
                c.current = I_final; 

                % Determination of air core's magnetic moment contribution.
                A_loop = pi*c.sum; 
                c.M_dipole = c.current * A_loop; 

                % Determination of core's magnetic moment contribution. 
                Nd = demag(c.l_core, r_core); 
                gain = 1 + (core.mu_r - 1) / (1 + (core.mu_r - 1)*Nd);
                c.M_9 = c.M_dipole * gain;

                fprintf('AWG %2d | n_wraps=%3d | R=%.2f Ω | I=%.3f A | sum(r^2)=%.3e | M_air=%.3f | gain=%.1f | M=%.3f\n', ...
    c.awg, c.n_wraps, c.res_total, c.current, c.sum, c.M_dipole, gain, c.M_9);

                if c.M_9 > 0 
                    results(n,:) = {n, c.M_9, c.M_dipole, core.name, c.awg, l_core, r_core, c.r_outer, c.v, c.n_wraps, c.current, c.m_total}; 
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
    % end 
   
end 

writetable(results, 'results.csv')

results.("Moment Sensitivity [A*m^2 / A]") = results.("Magnetic Moment [A*m^2]") ./ results.("Current [A]"); 
figure; 
scatter(results, "Current [A]", "Moment Sensitivity [A*m^2 / A]", "filled");

figure; 
scatter(results, "Core Radius [m]", "Magnetic Moment [A*m^2]", "filled"); 
text(results.("Core Radius [m]"), results.("Magnetic Moment [A*m^2]"), ...
     string(results.("Design Number")), 'VerticalAlignment','top', ...
     'HorizontalAlignment','left')


% basic plots to be optimized in future release

figure; 
scatter(results, "Wire Gauge", "Magnetic Moment [A*m^2]", "filled")
text(results.("Wire Gauge"), results.("Magnetic Moment [A*m^2]"), ...
     string(results.("Design Number")), 'VerticalAlignment','top', ...
     'HorizontalAlignment','left')

figure; 
scatter(results, "Total Mass [kg]", "Magnetic Moment [A*m^2]", 'filled') 
text(results.("Total Mass [kg]"), results.("Magnetic Moment [A*m^2]"), ...
    string(results.("Design Number")), 'VerticalAlignment','top', ...
    'HorizontalAlignment','left')


figure; 
scatter(results,"Current [A]", "Magnetic Moment [A*m^2]", 'filled')
text(results.("Current [A]"), results.("Magnetic Moment [A*m^2]"), ...
    string(results.("Design Number")), 'VerticalAlignment','top', ...
    'HorizontalAlignment','left')

figure; 
x = ["Usable Designs" "Volume Fail" "Mass Fail" "Current Fail"]; 
goobers = [good vol_fail mass_fail curr_fail]; 
bar(x, goobers); 

figure; 
scatter(results, "Wire Gauge", "Current [A]", "filled")

%{
Checks if another layer of wire can be added and if so, wraps it. 
%}
function [coil, wrapping] = wrapNext(core, coil, params)

    n_turns = floor(coil.l_core / coil.d_jacket);

    % geometry for new layer
    r_inner = coil.r_outer;
    r_outer = r_inner + coil.d_jacket;
    r_mid   = r_inner + 0.5*coil.d_jacket;

    % properties of new layer (CONDUCTOR dia for mass/R)
    L_layer = 2*pi*r_mid*n_turns;
    A_wire  = pi*(coil.d_conductor/2)^2;
    R_layer = coil.res_per_m * L_layer;
    m_layer = params.rho_cu * L_layer * A_wire;

    % new totals
    R_new        = coil.res_total + R_layer;
    m_coil_new   = coil.m_coil + m_layer;
    m_total_new  = core.mass + m_coil_new;
    v_new        = pi*(r_outer^2)*coil.l_core;

    % geometry / mass constraints only
    if r_outer > params.r_outer_max
        coil.fail_code = 1; wrapping=false; return
    elseif m_total_new > params.m_max
        coil.fail_code = 2; wrapping=false; return
    end


    % accept layer
    coil.n_wraps   = coil.n_wraps + 1;
    coil.sum       = coil.sum + n_turns * (r_mid^2);   % SUM of r^2 over ALL turns
    coil.r_outer   = r_outer;

    coil.v         = v_new;
    coil.res_total = R_new;
    coil.m_total   = m_total_new;
    coil.m_coil    = m_coil_new;

    wrapping = true;
end


% Demagnetizing calculation
function Nd = demag(l_core, r_core)
    x = l_core/r_core; 
    Nd = (4*(log(x) - 1)) / ((x^2 - 4*log(x)));
end

% Properties of different alloy core choices. 
function props = getProperties(name)
    switch(name)
        case '4140'
            props.mu_r = 400; 
            props.rho = 7850; % kg/m^3
            % props.minmach = x; % minimum machinable radius 
        case 'something' % TODO
            props.mu_r = 0; 
            props.rho = 0;
        otherwise 
            error('Unknown alloy: %s', name)
    end 
end 