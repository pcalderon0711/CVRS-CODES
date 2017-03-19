%% This ODE represents the HIV model in Section 4.2
function dy=ODEmodel(t,Y,X,run_num,tspan,ts)

%% PARAMETERS %%
Parameter_settings_EFAST;

c_as = X(run_num,1);
c_vs = X(run_num,2);
c_ap = X(run_num,3);
c_vp = X(run_num,4);
c_l = X(run_num,5);
c_r = X(run_num,6);
R_l = X(run_num,7);
R_r = X(run_num,8);
kappa = X(run_num,9);
alpha_l = X(run_num,10);
alpha_r = X(run_num,11);
beta_l = X(run_num,12);
beta_r = X(run_num,13);
gamma_l = X(run_num,14);
gamma_r = X(run_num,15);
M_O2 = X(run_num,16);
M_CO2 = X(run_num,17);
rho_O2 = X(run_num,18);
rho_CO2 = X(run_num,19);
R_p_rest = X(run_num,20);
A_pesk_rest = X(run_num,21);
P_IO2 = X(run_num,22);
P_ICO2 = X(run_num,23);
V_AO2 = X(run_num,24);
V_ACO2 = X(run_num,25);
V_TO2 = X(run_num,26);
V_TCO2 = X(run_num,27);
K_CO2 = X(run_num,28);
k_CO2 = X(run_num,29);
K_a1 = X(run_num,30);
K_a2 = X(run_num,31);
A_pesk_exer = X(run_num,32);
R_p_exer = X(run_num,33);
H_rest = X(run_num,34);
H_exer = X(run_num,35);
dot_VA_rest = X(run_num,36);
dot_VA_exer = X(run_num,37);

% ===== STATE VARIABLES ===== %
P_as = Y(1);
P_vs = Y(2);
P_ap = Y(3);
P_vp = Y(4);
S_l = Y(5);
sigma_l = Y(6);
S_r = Y(7);
sigma_r = Y(8);
P_aCO2 = Y(9);
P_aO2 = Y(10);
C_vCO2 = Y(11);
C_vO2 = Y(12);

% ========== COMPUTED QUANTITIES ========= %
% W = @(v)myTrans(v,ts,W_rest,W_exer);
% A_pesk = @(v)myTrans(v,ts,A_pesk_rest,A_pesk_exer);
% R_p = @(v)myTrans(v,ts,R_p_rest,R_p_exer);
% H = @(v)myTrans(v,ts,H_rest,H_exer);
% dot_VA = @(v)myTrans(v,ts,dot_VA_rest,dot_VA_exer);

W = @(v)myTrans(v,ts,W_rest,mySinwave(v-ts,W_rest,W_exer));
A_pesk = @(v)myTrans(t,ts,A_pesk_rest,mySinwave(v-ts,A_pesk_rest,A_pesk_exer)); %A_pesk function
R_p = @(v)myTrans(t,ts,R_p_rest,mySinwave(v-ts,R_p_rest,R_p_exer));
H = @(v)myTrans(v,ts,H_rest,mySinwave(v-ts,H_rest,H_exer));
dot_VA = @(v)myTrans(v,ts,dot_VA_rest,mySinwave(v-ts,dot_VA_rest,dot_VA_exer));

R_s = A_pesk(t) * C_vO2;
t_d = (H(t)^(-0.5))*((H(t)^(-0.5)) - kappa);
k_l = exp(-((c_l*R_l)^-1)*t_d);
a_l = 1 - k_l;
k_r = exp(-((c_r*R_r)^-1)*t_d);
a_r = 1 - k_r;

F_s = (R_s^-1) * (P_as - P_vs);
F_p = (R_p(t)^-1) * (P_ap - P_vp);
Q_l = H(t)*((c_l*a_l*P_vp*S_l)*(a_l*P_as + k_l*S_l)^-1);
Q_r = H(t)*((c_r*a_r*P_vs*S_r)*(a_r*P_ap + k_r*S_r)^-1);

% MR_O2 = M_O2 + rho_O2*W;
% MR_CO2 = M_CO2 + rho_CO2*W;

MR_O2 = M_O2 + rho_O2*W(t);
MR_CO2 = M_CO2 + rho_CO2*W(t);

% ========== DERIVATIVES ========== %
dy(1) = (c_as^-1) * (Q_l - F_s);
dy(2) = (c_vs^-1) * (F_s - Q_r);
dy(3) = (c_ap^-1) * (Q_r - F_p);
dy(4) = (c_vp^-1) * (F_p - Q_l);
dy(5) = sigma_l;
dy(6) = -gamma_l * sigma_l - alpha_l * S_l + beta_l * H(t);
dy(7) = sigma_r;
dy(8) = -gamma_r * sigma_r - alpha_r * S_r+ beta_r * H(t);
dy(9) = (V_ACO2^-1) * (863*F_p*(C_vCO2-K_CO2*P_aCO2 - k_CO2) + dot_VA(t) * (P_ICO2 - P_aCO2));
dy(10) = (V_AO2^-1) * (863*F_p*(C_vO2-K_a1*(1-exp(-K_a2*P_aO2))^2) + dot_VA(t) * (P_IO2 - P_aO2));
dy(11) = (V_TCO2^-1) * (MR_CO2 + F_s*(K_CO2*P_aCO2 + k_CO2 - C_vCO2));
dy(12) = (V_TO2^-1) * (-MR_O2 + F_s*(K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2));