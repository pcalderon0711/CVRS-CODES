function dy = myModelWithControl(t,Y,u,paramsToEstimate,paramsNotToEstimate,lookupToEstimate,lookupNotToEstimate,ts)

% ========= PARAMETERS ========== %
for i=1:37
    if ismember(i, lookupToEstimate)
        idx = find(lookupToEstimate==i,1);
        if i==1
            c_as = paramsToEstimate(idx);
        elseif i==2
            c_vs = paramsToEstimate(idx);
        elseif i==3
            c_ap = paramsToEstimate(idx);
        elseif i==4
            c_vp = paramsToEstimate(idx);
        elseif i==5
            c_l = paramsToEstimate(idx);
        elseif i==6
            c_r = paramsToEstimate(idx);
        elseif i==7
            R_l = paramsToEstimate(idx);
        elseif i==8
            R_r = paramsToEstimate(idx);
        elseif i==9
            kappa = paramsToEstimate(idx);
        elseif i==10
            alpha_l = paramsToEstimate(idx);
        elseif i==11
            alpha_r = paramsToEstimate(idx);
        elseif i==12
            beta_l = paramsToEstimate(idx);
        elseif i==13
            beta_r = paramsToEstimate(idx);
        elseif i==14
            gamma_l = paramsToEstimate(idx);
        elseif i==15
            gamma_r = paramsToEstimate(idx);
        elseif i==16
            M_O2 = paramsToEstimate(idx);
        elseif i==17
            M_CO2 = paramsToEstimate(idx);
        elseif i==18
            rho_O2 = paramsToEstimate(idx);
        elseif i==19
            rho_CO2 = paramsToEstimate(idx);
%         elseif i==20
%             q_as = paramsToEstimate(idx);
%         elseif i==21
%             V_tot = paramsToEstimate(idx);
        elseif i==22
            R_p_rest = paramsToEstimate(idx);
        elseif i==23
            A_pesk_rest = paramsToEstimate(idx);
        elseif i==24
            P_IO2 = paramsToEstimate(idx);
        elseif i==25
            P_ICO2 = paramsToEstimate(idx);
        elseif i==26
            V_AO2 = paramsToEstimate(idx);
        elseif i==27
            V_ACO2 = paramsToEstimate(idx);
        elseif i==28
            V_TO2 = paramsToEstimate(idx);
        elseif i==29
            V_TCO2 = paramsToEstimate(idx);
        elseif i==30
            K_CO2 = paramsToEstimate(idx);
        elseif i==31
            k_CO2 = paramsToEstimate(idx);
        elseif i==32
            K_a1 = paramsToEstimate(idx);
        elseif i==33
            K_a2 = paramsToEstimate(idx);
        elseif i==34
            W_rest = paramsToEstimate(idx);
        elseif i==35
            A_pesk_exer = paramsToEstimate(idx);
        elseif i==36
            R_p_exer = paramsToEstimate(idx);
        elseif i==37
            W_exer = paramsToEstimate(idx);
        end
    else
        idx = find(lookupNotToEstimate==i,1);
        if i==1
            c_as = paramsNotToEstimate(idx);
        elseif i==2
            c_vs = paramsNotToEstimate(idx);
        elseif i==3
            c_ap = paramsNotToEstimate(idx);
        elseif i==4
            c_vp = paramsNotToEstimate(idx);
        elseif i==5
            c_l = paramsNotToEstimate(idx);
        elseif i==6
            c_r = paramsNotToEstimate(idx);
        elseif i==7
            R_l = paramsNotToEstimate(idx);
        elseif i==8
            R_r = paramsNotToEstimate(idx);
        elseif i==9
            kappa = paramsNotToEstimate(idx);
        elseif i==10
            alpha_l = paramsNotToEstimate(idx);
        elseif i==11
            alpha_r = paramsNotToEstimate(idx);
        elseif i==12
            beta_l = paramsNotToEstimate(idx);
        elseif i==13
            beta_r = paramsNotToEstimate(idx);
        elseif i==14
            gamma_l = paramsNotToEstimate(idx);
        elseif i==15
            gamma_r = paramsNotToEstimate(idx);
        elseif i==16
            M_O2 = paramsNotToEstimate(idx);
        elseif i==17
            M_CO2 = paramsNotToEstimate(idx);
        elseif i==18
            rho_O2 = paramsNotToEstimate(idx);
        elseif i==19
            rho_CO2 = paramsNotToEstimate(idx);
%         elseif i==20
%             q_as = paramsNotToEstimate(idx);
%         elseif i==21
%             V_tot = paramsNotToEstimate(idx);
        elseif i==22
            R_p_rest = paramsNotToEstimate(idx);
        elseif i==23
            A_pesk_rest = paramsNotToEstimate(idx);
        elseif i==24
            P_IO2 = paramsNotToEstimate(idx);
        elseif i==25
            P_ICO2 = paramsNotToEstimate(idx);
        elseif i==26
            V_AO2 = paramsNotToEstimate(idx);
        elseif i==27
            V_ACO2 = paramsNotToEstimate(idx);
        elseif i==28
            V_TO2 = paramsNotToEstimate(idx);
        elseif i==29
            V_TCO2 = paramsNotToEstimate(idx);
        elseif i==30
            K_CO2 = paramsNotToEstimate(idx);
        elseif i==31
            k_CO2 = paramsNotToEstimate(idx);
        elseif i==32
            K_a1 = paramsNotToEstimate(idx);
        elseif i==33
            K_a2 = paramsNotToEstimate(idx);
        elseif i==34
            W_rest = paramsNotToEstimate(idx);
        elseif i==35
            A_pesk_exer = paramsNotToEstimate(idx);
        elseif i==36
            R_p_exer = paramsNotToEstimate(idx);
        elseif i==37
            W_exer = paramsNotToEstimate(idx);
        end
    end
end

% ===== STATE VARIABLES ===== %
P_as = Y(1);
P_vs = Y(2);
P_ap = Y(3);
P_vp = Y(4);
S_l = Y(5);
sigma_l = Y(6);
S_r = Y(7);
sigma_r = Y(8);
H = Y(9);
P_aCO2 = Y(10);
P_aO2 = Y(11);
C_vCO2 = Y(12);
C_vO2 = Y(13);
dotV_A = Y(14);

%%%% CONTROLS %%%%
u1 = u(1);
u2 = u(2);

% ========== COMPUTED QUANTITIES ========= %
W = @(v)myTrans(v,ts,W_rest+W_exer-W_exer,W_exer+W_rest-W_rest);
A_pesk = @(v)myTrans(v,ts,A_pesk_rest+A_pesk_exer-A_pesk_exer,A_pesk_exer+A_pesk_rest-A_pesk_rest);
R_p = @(v)myTrans(v,ts,R_p_rest+R_p_exer-R_p_exer,R_p_exer+R_p_rest-R_p_rest);

R_s = A_pesk(t) * C_vO2;
t_d = (H^(-0.5))*((H^(-0.5)) - kappa);
k_l = exp(-((c_l*R_l)^-1)*t_d);
a_l = 1 - k_l;
k_r = exp(-((c_r*R_r)^-1)*t_d);
a_r = 1 - k_r;

F_s = (R_s^-1) * (P_as - P_vs);
F_p = (R_p(t)^-1) * (P_ap - P_vp);
Q_l = H*((c_l*a_l*P_vp*S_l)*(a_l*P_as + k_l*S_l)^-1);
Q_r = H*((c_r*a_r*P_vs*S_r)*(a_r*P_ap + k_r*S_r)^-1);

% MR_O2 = M_O2 + rho_O2*W;
% MR_CO2 = M_CO2 + rho_CO2*W;

MR_O2 = M_O2 + rho_O2*W(t);
MR_CO2 = M_CO2 + rho_CO2*W(t);

% ========== DERIVATIVES ========== %
dy1 = (c_as^-1) * (Q_l - F_s);
dy2 = (c_vs^-1) * (F_s - Q_r);
dy3 = (c_ap^-1) * (Q_r - F_p);
dy4 = (c_vp^-1) * (F_p - Q_l);
dy5 = sigma_l;
dy6 = -gamma_l * sigma_l - alpha_l * S_l + beta_l * H;
dy7 = sigma_r;
dy8 = -gamma_r * sigma_r - alpha_r * S_r+ beta_r * H;
dy9 = 0;
dy10 = (V_ACO2^-1) * (863*F_p*(C_vCO2-K_CO2*P_aCO2 - k_CO2) + dotV_A * (P_ICO2 - P_aCO2));
dy11 = (V_AO2^-1) * (863*F_p*(C_vO2-K_a1*(1-exp(-K_a2*P_aO2))^2) + dotV_A * (P_IO2 - P_aO2));
dy12 = (V_TCO2^-1) * (MR_CO2 + F_s*(K_CO2*P_aCO2 + k_CO2 - C_vCO2));
dy13 = (V_TO2^-1) * (-MR_O2 + F_s*(K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2));
dy14 = 0;

dy = [dy1;dy2;dy3;dy4;dy5;dy6;dy7;dy8;dy9;dy10;dy11;dy12;dy13;dy14]; 
