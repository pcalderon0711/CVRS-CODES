function dy=myModel(t,Y,LHSmatrix,x,ts)
%% PARAMETERS %%
Parameter_settings_LHS;

c_as = LHSmatrix(x,1);
c_vs = LHSmatrix(x,2);
c_ap = LHSmatrix(x,3);
c_vp = LHSmatrix(x,4);
c_l = LHSmatrix(x,5);
c_r = LHSmatrix(x,6);
R_l = LHSmatrix(x,7);
R_r = LHSmatrix(x,8);
kappa = LHSmatrix(x,9);
alpha_l = LHSmatrix(x,10);
alpha_r = LHSmatrix(x,11);
beta_l = LHSmatrix(x,12);
beta_r = LHSmatrix(x,13);
gamma_l = LHSmatrix(x,14);
gamma_r = LHSmatrix(x,15);
M_O2 = LHSmatrix(x,16);
M_CO2 = LHSmatrix(x,17);
rho_O2 = LHSmatrix(x,18);
rho_CO2 = LHSmatrix(x,19);
R_p_rest = LHSmatrix(x,20);
A_pesk_rest = LHSmatrix(x,21);
P_IO2 = LHSmatrix(x,22);
P_ICO2 = LHSmatrix(x,23);
V_AO2 = LHSmatrix(x,24);
V_ACO2 = LHSmatrix(x,25);
V_TO2 = LHSmatrix(x,26);
V_TCO2 = LHSmatrix(x,27);
K_CO2 = LHSmatrix(x,28);
k_CO2 = LHSmatrix(x,29);
K_a1 = LHSmatrix(x,30);
K_a2 = LHSmatrix(x,31);
A_pesk_exer = LHSmatrix(x,32);
R_p_exer = LHSmatrix(x,33);
H_rest = LHSmatrix(x,34);
H_exer = LHSmatrix(x,35);
dot_VA_rest = LHSmatrix(x,36);
dot_VA_exer = LHSmatrix(x,37);

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
W = @(v)myTrans(v,ts,W_rest,W_exer);
A_pesk = @(v)myTrans(v,ts,A_pesk_rest,A_pesk_exer);
R_p = @(v)myTrans(v,ts,R_p_rest,R_p_exer);
H = @(v)myTrans(v,ts,H_rest,H_exer);
dot_VA = @(v)myTrans(v,ts,dot_VA_rest,dot_VA_exer);

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