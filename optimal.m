
min_l = 0.01; 
max_l = 0.1; 
min_d = 0.0005;
max_d = 0.005;
inc_l = 0.01; 
inc_d = 0.0005; 

params.m_core = 0.1; 
params.r_core = 0.01; 
params.rho_cu = 8960; 
params.m_max = 0.3; 
params.r_max = 0.05;   
params.I = 0.25;       % current in amps

 

l_values = min_l:inc_l:max_l; % range of lengths
d_values = min_d:inc_d:max_d; % range of diameters

coil_template = struct( ...
    'd_wire', 0, ...
    'l_core', 0, ...
    'n_wraps', 0, ...
    'r_outer', params.r_core, ...
    'm_total', params.m_core, ...
    'M_dipole', 0 ...
);

coil = repmat(coil_template, length(l_values), length(d_values));
M = zeros(length(l_values), length(d_values));

for i = 1:length(l_values)
    l_core = l_values(i); 
    for j = 1:length(d_values)
        d_wire = d_values(j); 
        
        % initialize coil 
        coil(i,j).d_wire = d_wire; 
        coil(i,j).l_core = l_core; 
        coil(i,j).n_wraps = 0; 
        coil(i,j).r_outer = params.r_core; 
        coil(i,j).m_total = params.m_core; 
        % wrap until breaks constraints
        while canWrap(coil(i,j), params)
            coil(i,j) = wrap(coil(i,j), params); 
        end
        % calculate magnetic dipole for dipole plot
        n_turns = floor(l_core / d_wire);
        N_total = n_turns * coil(i,j).n_wraps;
        r_avg = (coil(i,j).r_outer + params.r_core) / 2;
        A = pi * r_avg^2;
        coil(i,j).M_dipole = N_total * params.I * A;
        M(i,j) = coil(i,j).M_dipole;
    end 
end 

function can = canWrap(coil, params)
    d_wire = coil.d_wire;
    n_layers = coil.n_wraps + 1; 
    n_turns = floor(coil.l_core / d_wire);
    r_outer_new = params.r_core + n_layers * d_wire;
    r_avg = (r_outer_new + params.r_core) / 2;
    l_wire = 2 * pi * r_avg * n_turns * n_layers;
    A_wire = pi * (d_wire^2)/4;
    m_wire_total = params.rho_cu * l_wire * A_wire;

    if (r_outer_new > params.r_max)
        can = false;
    elseif (m_wire_total + params.m_core > params.m_max)
        can = false;
    else
        can = true;
    end
end

function coil = wrap(coil, params)
    coil.n_wraps = coil.n_wraps + 1;
    d_wire = coil.d_wire;
    n_turns = floor(coil.l_core / d_wire);
    r_outer_new = params.r_core + coil.n_wraps * d_wire;
    r_avg = (r_outer_new + params.r_core) / 2;
    l_wire = 2 * pi * r_avg * n_turns * coil.n_wraps;
    A_wire = pi * (d_wire^2)/4;
    m_wire_total = params.rho_cu * l_wire * A_wire;
    coil.r_outer = r_outer_new;
    coil.m_total = m_wire_total + params.m_core;
end 

[DD, LL] = meshgrid(d_values, l_values);
figure;
surf(DD, LL, M);
xlabel('Wire Diameter [m]');
ylabel('Core Length [m]');
zlabel('Magnetic Dipole Moment [A·m^2]');
title('Magnetorquer Design Space');
shading interp;
colorbar;
grid on;

       

       




