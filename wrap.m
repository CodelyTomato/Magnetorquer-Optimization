function [D,F] = wrap(p, d_ins, d_con, N, awg, I_wire_max, gain)
    % --- initialize first layer --- % 
    k = 1; 
    r_outer = p.r_core; 
    R_total = 0; 
    L_total = 0; 
    CSA_wire = (pi/4)*d_con^2; 
    M_core = 0; 
    fail_code = 0; 
    D = []; % make huge to preinitialize ?

    while true 
        % --- new geometry --- % 
        r_inner = r_outer; 
        r_outer = r_inner + d_ins; 
        r_ctr = r_inner + 0.5*d_ins; 

        % --- new layer --- % 
        A_layer = pi*r_ctr^2; 
        L_layer = 2*pi*r_ctr*N; 
        R_layer = p.res_cu*L_layer / CSA_wire; 
        % mass_layer = p.dens_cu*CSA_wire*L_layer; 

        % --- new totals --- % 
        R_new = R_total + R_layer; 
        L_new = L_total + L_layer; 
        I_new = p.V_in / R_new; 
        mass_new = p.mass_core + (p.dens_cu * pi * (r_outer^2 - p.r_core^2)*p.L_core); 
        vol_new = pi*r_outer^2*p.L_core; 
        
        % --- check against constraints --- % 
        if mass_new > p.mass_max 
            fail_code = 1; 
            D(end,13) = fail_code; 
            break; 
        elseif vol_new > p.vol_max
            fail_code = 2;
            D(end,13) = fail_code; 
            break; 
        else
            if I_new > I_wire_max || I_new > p.I_battery_max
                I_new = min(p.I_battery_max, I_wire_max); % don't need to recaluclate all the time
            end
                m_i = N*A_layer; 
                if isempty(D)
                    m_total = m_i; 
                else
                    m_prev = D(end, 3); 
                    m_total = m_prev + m_i; 
                end 
                M_air = m_total * I_new;
                M_core = gain * M_air; 
                new_row = [M_core, M_air, m_total, R_new, L_new, A_layer, I_new, k, r_outer, mass_new, awg, d_ins, fail_code]; 
                D = [D; new_row]; 

                k = k+1; 
                R_total = R_new; 
                L_total = L_new; 
        end 
    end 


end

