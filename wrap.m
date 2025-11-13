%{
----- COIL WRAPPER -----
function wrap() takes in params and wire information and returns a
table of designs 'D', which contains rows of possible designs with that wire. 
%}
function D = wrap(p, d_ins, d_con, N, awg, I_wire_max, gain)
    % --- initialize first layer --- % 
    k = 1; % counts number of wraps in the longitudinal direction. 
    r_outer = p.r_core; % m 
    R_total = 0; % omh*m
    L_total = 0; % m
    CSA_wire = (pi/4)*d_con^2; % cross-sectional area, m^2 
    fail_code = 0; % keeps tracks of what caused a design to fail 
    D = []; 

    while true 
        % --- new torquer geometry --- % 
        r_inner = r_outer; % m 
        r_outer = r_inner + d_ins; % m 
        r_ctr = r_inner + 0.5*d_ins; % wire centerline radius 

        % --- new wire layer --- % 
        A_layer = pi*r_ctr^2; % m^2
        L_layer = 2*pi*r_ctr*N; % m 
        R_layer = p.res_cu*L_layer / CSA_wire; % ohm*m 

        % --- new totals --- % 
        R_new = R_total + R_layer; % ohm*m
        L_new = L_total + L_layer; % m 
        I_new = p.V_in / R_new; % A 
        mass_new = p.mass_core + (p.dens_cu * pi * (r_outer^2 - p.r_core^2)*p.L_core); % kg 
        vol_new = pi*r_outer^2*p.L_core; % m^3
        
        % --- check this against constraints --- % 
        % if the new mass is above the max allowable sets fail code to 1
        % and breaks. If the new volume is above the max allowable sets 
        % fail code to 2 and breaks.  
        if mass_new > p.mass_max 
            fail_code = 1; 
            D(end,13) = fail_code; 
            break; 
        elseif vol_new > p.vol_max
            fail_code = 2;
            D(end,13) = fail_code; 
            break; 
        else
            % - check current draw - % 
            % If the wire's new current draw is greater than the maximum
            % the battery or wire can support, set the new current draw to
            % be the smaller of the two. 
            if I_new > I_wire_max || I_new > p.I_battery_max
                % change to be function parameter
                I_new = min(p.I_battery_max, I_wire_max); 
            end
                m_i = N*A_layer; % little m 
                if isempty(D) % first wrap  
                    m_total = m_i;  
                else
                    % N*A_layer scales as more wraps are added by summing 
                    % m_i values for each design D(i,:);  
                    m_prev = D(end, 3); 
                    m_total = m_prev + m_i; 
                end 
                M_air = m_total * I_new; % A*m^2
                M_core = gain * M_air; % A*m^2
                % Append new design to 'D' 
                new_row = [M_core, M_air, m_total, R_new, L_new, ...
                    A_layer, I_new, k, r_outer, mass_new, awg, ...
                    d_ins, fail_code]; 
                D = [D; new_row]; 

                % - update totals - % 
                k = k+1; 
                R_total = R_new; 
                L_total = L_new; 
        end 
    end 
end

