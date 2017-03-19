function dy = myModelWithControlAndTimeVaryingWorkload(t,Y,x,LHSmatrix,t1,W,A, R_p)

%% myModelWithControlAndTimeVaryingWorkload
% Calculates the rhs of dot(t,Y) = f(t,Y) at time t and state Y
% Inputs:
%        t      : time
%        Y      : vector of state variable values
%        params : vector of parameters
%        W      : workload function on t
%        A      : A_pesk function on t
%        R_p    : R_p function on t
% Output:
%        dy     : value of f at time t and state Y
%
%% ========= PARAMETERS ========== %
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
V_tot = 5.0582;

%% ===== STATE VARIABLES ===== %
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

%% ========== COMPUTED QUANTITIES ========= %
R_s = A(t) * C_vO2;
t_d = 1/H - kappa/sqrt(abs(H));
k_l = exp(-t_d/(c_l*R_l));
a_l = 1 - k_l;
k_r = exp(-t_d/(c_r*R_r));
a_r = 1 - k_r;

F_s = (P_as - P_vs)/R_s;
F_p = (P_ap - P_vp)/R_p(t);
Q_l = c_l*a_l*P_vp*S_l*H/(a_l*P_as + k_l*S_l);
Q_r = c_r*a_r*P_vs*S_r*H/(a_r*P_ap + k_r*S_r);

MR_O2 = M_O2 + rho_O2*W(t);
MR_CO2 = M_CO2 + rho_CO2*W(t);

dy = zeros(14,1);

%% ========== DERIVATIVES ========== %
dy(1) = (Q_l - F_s)/c_as;
dy(2) = (F_s - Q_r)/c_vs;
dy(3) = (Q_r - F_p)/c_ap;
dy(4) = (F_p - Q_l)/c_vp;
dy(5) = sigma_l;
dy(6) = -gamma_l * sigma_l - alpha_l * S_l + beta_l * H;
dy(7) = sigma_r;
dy(8) = -gamma_r * sigma_r - alpha_r * S_r+ beta_r * H;
dy(9) = 0;
dy(10) = (863*F_p*(C_vCO2-K_CO2*P_aCO2 - k_CO2) + dotV_A * (P_ICO2 - P_aCO2))/V_ACO2;
dy(11) = (863*F_p*(C_vO2-K_a1*(1-exp(-K_a2*P_aO2))^2) + dotV_A * (P_IO2 - P_aO2))/V_AO2;
dy(12) = (MR_CO2 + F_s*(K_CO2*P_aCO2 + k_CO2 - C_vCO2))/V_TCO2;
dy(13) = (-MR_O2 + F_s*(K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2))/V_TO2;
dy(14) = 0;
