function idx_set = myIdentifyParam(param_id)

param_names = {'c_as', 'c_vs', 'c_ap', 'c_vp', 'c_l', 'c_r', 'R_l', 'R_r', ...
    'kappa', 'alpha_l', 'alpha_r', 'beta_l', 'beta_r', 'gamma_l', 'gamma_r',...
    'M_O2', 'M_CO2', 'rho_O2', 'rho_CO2', 'V_tot', 'R_p_rest', 'R_p_exer', 'A_pesk_rest', ...
    'A_pesk_exer', 'P_IO2', 'P_ICO2', 'V_AO2', 'V_ACO2', 'V_TO2', 'V_TCO2', 'K_CO2', 'k_CO2',...
    'K_a1', 'K_a2', 'q_as', 'q_aco2', 'w1', 'w2', 'alpha', 'beta'};

idx_set = zeros(length(param_id),1);
for i=1:length(param_id)
    for j=1:length(param_names)
        if strcmp(param_id(i),param_names(j))==1
            idx_set(i) = j;
            break;
        end
    end
end

