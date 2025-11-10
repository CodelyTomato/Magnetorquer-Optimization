

params.dens_wire = 8960; % kg/m^3


function D = wrap()

    CSA_wire = pi*d_con^2*0.25; % cross-sectional area of just conductor

    % --- initialize first layer --- 
    r_outer = r_core; 
    R_total = 0; 
    L_total = 0; 
    mass_total = mass_core; 
    vol_total = vol_core; 
    M_core = NaN;
    M_air = NaN; 
    D = []; 
   
    while true
        
        r_inner = r_outer; 
        r_outer = r_inner + d_ins; 
        r_ctr = r_inner + 0.5*d_ins; % centerline radius of wire
        
        A_layer = pi*(r_ctr^2); % this area is for moment calculation 
        L_layer = 2*pi*r_ctr*N; 
        R_layer = res*L_layer / CSA_wire; 
        mass_layer = params.dens_wire*CSA_wire*L_layer; 
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