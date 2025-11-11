

params.dens_cu = 8960; % kg/m^3
params.res_cu = 1.68e-8; % ohm*m @ 20 C




dens_core = 8120; % 4140 Steel Alloy 
r_core = 0.004; % m 
L_core = 0.9; % m 
vol_core = pi*r_core^2*L_core; 
vol_max = pi*0.01^2*L_core;  
mass_core = vol_core * dens_core; 
I_wire_max = 0.226; % A
I_battery_max = 2; % A
d_ins = 0.348/1000; % mm
d_con = 0.32/1000; % mm 
N = floor(L_core / d_ins); 
mass_max = 1; % kg 
V_in = 12; % V 

D = wrap(r_core, mass_core, vol_core, params, N, mass_max, ...
    vol_max, I_wire_max, I_battery_max, d_ins, d_con, V_in); 

T = array2table(D);

% 2. Define your desired column titles (must be 7 titles for 7 columns)
columnTitles = {'M_core', 'M_air', 'm', 'Resistance', 'Length', 'Area', 'Current'};

% 3. Assign the titles to the table's VariableNames property
T.Properties.VariableNames = columnTitles;

% Display the new table with titles
disp(T);




x = D(:, 7); 
y = D(:, 2); 

plot(x, y, 'o-'); 
xlabel('current'); 
ylabel('magnetic moment');


function D = wrap(r_core, mass_core, vol_core, params, N, mass_max, ...
    vol_max, I_wire_max, I_battery_max, d_ins, d_con, V_in) 

    CSA_wire = pi*d_con^2*0.25; % cross-sectional area of just conductor

    % --- initialize first layer --- 
    r_outer = r_core; 
    R_total = 0; 
    L_total = 0; 
    mass_total = mass_core; 
    vol_total = vol_core; 
    M_core = 0;
    M_air = 0; 
    D = []; 
   
    while true
        
        r_inner = r_outer; 
        r_outer = r_inner + d_ins; 
        r_ctr = r_inner + 0.5*d_ins; % centerline radius of wire
        
        A_layer = pi*(r_ctr^2); % this area is for moment calculation 
        L_layer = 2*pi*r_ctr*N; 
        R_layer = params.res_cu*L_layer / CSA_wire; 
        mass_layer = params.dens_cu*CSA_wire*L_layer; 
        vol_layer = 0.25*pi*d_ins^2*L_layer; 

        % -- new totals -- 
        R_new = R_total + R_layer;  
        L_new = L_total + L_layer; 
        I_new = V_in / R_new; 
        mass_new = mass_total + mass_layer; 
        vol_new = vol_total + vol_layer; 

        % --- check new totals against constraints --- 
        if mass_new > mass_max || vol_new > vol_max
            break; 
        else
            % --- Add new layer ---
            if I_new > I_wire_max
                I_new = min(I_battery_max, I_wire_max);
            end 
            m = N*A_layer; 
            M_air = I_new * m; 
            new_row = [M_core, M_air, m, R_new, L_new, A_layer, I_new]; 
            D = [D; new_row]; 

            % --- Update for next loop ---
            R_total = R_new; 
            L_total = L_new; 
            mass_total = mass_new; 
            vol_total = vol_new; 
        end 
    end 
end 