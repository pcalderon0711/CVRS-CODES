function jacobian = myDerivativesTimeVarying(t,Y,x,LHSmatrix,t1,A,R_p)

%% myDerivativesTimeVarying
% Calculates the jacobian of f at time t and state Y
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
%% initialize the columns of the jacobian
df_Pas = zeros(14,1);
df_Pvs = zeros(14,1);
df_Pap = zeros(14,1);
df_Pvp = zeros(14,1);
df_Sl = zeros(14,1);
df_sigmal = zeros(14,1);
df_Sr = zeros(14,1);
df_sigmar = zeros(14,1);
df_H = zeros(14,1);
df_PaCO2 = zeros(14,1);
df_PaO2 = zeros(14,1);
df_CvCO2 = zeros(14,1);
df_CvO2 = zeros(14,1);
df_dotVA = zeros(14,1);
    
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
S_r = Y(7);
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

N_l = a_l*P_as+k_l*S_l;
N_r = a_r*P_ap+k_r*S_r;

%% ===== P_as ===== %

dFs_Pas = 1/R_s;
dQl_Pas = -c_l*a_l^2*P_vp*S_l*H/N_l^2;

df_Pas(1) = (dQl_Pas - dFs_Pas)/c_as;
df_Pas(2) = dFs_Pas/c_vs;
df_Pas(3) = 0;
df_Pas(4) = -dQl_Pas/c_vp;
df_Pas(5) = 0;
df_Pas(6) = 0;
df_Pas(7) = 0;
df_Pas(8) = 0;
df_Pas(9) = 0;
df_Pas(10) = 0;
df_Pas(11) = 0;
df_Pas(12) = (K_CO2*P_aCO2 + k_CO2 - C_vCO2)*dFs_Pas/V_TCO2;
df_Pas(13) = (K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2)*dFs_Pas/V_TO2;
df_Pas(14) = 0;

%% ===== P_vs ===== %
dFs_Pvs = -1/R_s;
dQr_Pvs = c_r*a_r*S_r*H/N_r;

df_Pvs(1) = -dFs_Pvs/c_as;
df_Pvs(2) = (dFs_Pvs - dQr_Pvs)/c_vs;
df_Pvs(3) = dQr_Pvs/c_ap;
df_Pvs(4) = 0;
df_Pvs(5) = 0;
df_Pvs(6) = 0;
df_Pvs(7) = 0;
df_Pvs(8) = 0;
df_Pvs(9) = 0;
df_Pvs(10) = 0;
df_Pvs(11) = 0;
df_Pvs(12) = (K_CO2*P_aCO2 + k_CO2 - C_vCO2)*dFs_Pvs/V_TCO2;
df_Pvs(13) = (K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2)*dFs_Pvs/V_TO2;
df_Pvs(14) = 0;

%% ===== P_ap ===== %
dFp_Pap = 1/R_p(t);
dQr_Pap = -c_r*a_r^2*P_vs*S_r*H/N_r^2;

df_Pap(1) = 0;
df_Pap(2) = -dQr_Pap/c_vs;
df_Pap(3) = (dQr_Pap - dFp_Pap)/c_ap;
df_Pap(4) = dFp_Pap/c_vp;
df_Pap(5) = 0;
df_Pap(6) = 0;
df_Pap(7) = 0;
df_Pap(8) = 0;
df_Pap(9) = 0;
df_Pap(10) = 863*(C_vCO2 - K_CO2*P_aCO2 - k_CO2)*dFp_Pap/V_ACO2;
df_Pap(11) = 863*(C_vO2 - K_a1*(1-exp(-K_a2*P_aO2))^2)*dFp_Pap/V_AO2;
df_Pap(12) = 0;
df_Pap(13) = 0;
df_Pap(14) = 0;

%% ===== P_vp ===== %
dFp_Pvp = -1/R_p(t);
dQl_Pvp = c_l*a_l*S_l*H/N_l;

df_Pvp(1) = dQl_Pvp/c_as;
df_Pvp(2) = 0;
df_Pvp(3) = -dFp_Pvp/c_ap;
df_Pvp(4) = (dFp_Pvp - dQl_Pvp)/c_vp;
df_Pvp(5) = 0;
df_Pvp(6) = 0;
df_Pvp(7) = 0;
df_Pvp(8) = 0;
df_Pvp(9) = 0;
df_Pvp(10) = 863*(C_vCO2 - K_CO2*P_aCO2 - k_CO2)*dFp_Pvp/V_ACO2;
df_Pvp(11) = 863*(C_vO2 - K_a1*(1-exp(-K_a2*P_aO2))^2)*dFp_Pvp/V_AO2;
df_Pvp(12) = 0;
df_Pvp(13) = 0;
df_Pvp(14) = 0;

%% ===== S_l ===== %
dQl_Sl = c_l*a_l^2*P_as*P_vp*H/N_l^2;

df_Sl(1) = dQl_Sl/c_as;
df_Sl(2) = 0;
df_Sl(3) = 0;
df_Sl(4) = -dQl_Sl/c_vp;
df_Sl(5) = 0;
df_Sl(6) = -alpha_l;
df_Sl(7) = 0;
df_Sl(8) = 0;
df_Sl(9) = 0;
df_Sl(10) = 0;
df_Sl(11) = 0;
df_Sl(12) = 0;
df_Sl(13) = 0;
df_Sl(14) = 0;

%% ===== sigma_l ===== %
df_sigmal(1) = 0;
df_sigmal(2) = 0;
df_sigmal(3) = 0;
df_sigmal(4) = 0;
df_sigmal(5) = 1;
df_sigmal(6) = -gamma_l;
df_sigmal(7) = 0;
df_sigmal(8) = 0;
df_sigmal(9) = 0;
df_sigmal(10) = 0;
df_sigmal(11) = 0;
df_sigmal(12) = 0;
df_sigmal(13) = 0;
df_sigmal(14) = 0;

%% ===== S_r ===== %
dQr_Sr = c_r*a_r^2*P_ap*P_vs*H/N_r^2;

df_Sr(1) = 0;
df_Sr(2) = -dQr_Sr/c_vs;
df_Sr(3) = dQr_Sr/c_ap;
df_Sr(4) = 0;
df_Sr(5) = 0;
df_Sr(6) = 0;
df_Sr(7) = 0;
df_Sr(8) = -alpha_r;
df_Sr(9) = 0;
df_Sr(10) = 0;
df_Sr(11) = 0;
df_Sr(12) = 0;
df_Sr(13) = 0;
df_Sr(14) = 0;

%% ===== sigma_r ===== %
df_sigmar(1) = 0;
df_sigmar(2) = 0;
df_sigmar(3) = 0;
df_sigmar(4) = 0;
df_sigmar(5) = 0;
df_sigmar(6) = 0;
df_sigmar(7) = 1;
df_sigmar(8) = -gamma_r;
df_sigmar(9) = 0;
df_sigmar(10) = 0;
df_sigmar(11) = 0;
df_sigmar(12) = 0;
df_sigmar(13) = 0;
df_sigmar(14) = 0;

%% ===== H ===== %
t_prime = -(1-.5*kappa*sqrt(abs(H)))/H^2;
kl_prime = -k_l*t_prime/(c_l*R_l);
kr_prime = -k_r*t_prime/(c_r*R_r);
dQl_H = c_l*P_vp*S_l*(a_l*N_l-kl_prime*S_l*H)/N_l^2;
dQr_H = c_r*P_vs*S_r*(a_r*N_r-kr_prime*S_r*H)/N_r^2;

df_H(1) = dQl_H/c_as;
df_H(2) = -dQr_H/c_vs;
df_H(3) = dQr_H/c_ap;
df_H(4) = -dQl_H/c_vp;
df_H(5) = 0;
df_H(6) = beta_l;
df_H(7) = 0;
df_H(8) = beta_r;
df_H(9) = 0;
df_H(10) = 0;
df_H(11) = 0;
df_H(12) = 0;
df_H(13) = 0;
df_H(14) = 0;

%% ===== P_aCO2 ===== %
df_PaCO2(1) = 0;
df_PaCO2(2) = 0;
df_PaCO2(3) = 0;
df_PaCO2(4) = 0;
df_PaCO2(5) = 0;
df_PaCO2(6) = 0;
df_PaCO2(7) = 0;
df_PaCO2(8) = 0;
df_PaCO2(9) = 0;
df_PaCO2(10) = -(863*K_CO2*F_p + dotV_A)/V_ACO2;
df_PaCO2(11) = 0;
df_PaCO2(12) = K_CO2*F_s/V_TCO2;
df_PaCO2(13) = 0;
df_PaCO2(14) = 0;

%% ===== P_aO2 ===== %
df_PaO2(1) = 0;
df_PaO2(2) = 0;
df_PaO2(3) = 0;
df_PaO2(4) = 0;
df_PaO2(5) = 0;
df_PaO2(6) = 0;
df_PaO2(7) = 0;
df_PaO2(8) = 0;
df_PaO2(9) = 0;
df_PaO2(10) = 0;
df_PaO2(11) = -(1726*K_a1*K_a2*F_p*(1-exp(-K_a2*P_aO2))*exp(-K_a2*P_aO2) + dotV_A)/V_AO2;
df_PaO2(12) = 0;
df_PaO2(13) = 2*K_a1*K_a2*F_s*(1-exp(-K_a2*P_aO2))*exp(-K_a2*P_aO2)/V_TO2;
df_PaO2(14) = 0;

%% ===== C_vCO2 ===== %
df_CvCO2(1) = 0;
df_CvCO2(2) = 0;
df_CvCO2(3) = 0;
df_CvCO2(4) = 0;
df_CvCO2(5) = 0;
df_CvCO2(6) = 0;
df_CvCO2(7) = 0;
df_CvCO2(8) = 0;
df_CvCO2(9) = 0;
df_CvCO2(10) = 863*F_p/V_ACO2;
df_CvCO2(11) = 0;
df_CvCO2(12) = -F_s/V_TCO2;
df_CvCO2(13) = 0;
df_CvCO2(14) = 0;

%% ===== C_vO2 ===== %
df_CvO2(1) = A(t)*(P_as - P_vs)/(c_as*R_s^2);
df_CvO2(2) = -A(t)*(P_as - P_vs)/(c_vs*R_s^2);
df_CvO2(3) = 0;
df_CvO2(4) = 0;
df_CvO2(5) = 0;
df_CvO2(6) = 0;
df_CvO2(7) = 0;
df_CvO2(8) = 0;
df_CvO2(9) = 0;
df_CvO2(10) = 0;
df_CvO2(11) = 863*F_p/V_AO2;
% df_CvO2(12) = 0;
% df_CvO2(13) = -F_s/V_TO2;

caco2=K_CO2*P_aCO2 + k_CO2;
cao2=K_a1*(1-exp(-K_a2*P_aO2))^2;
dfdcvo2=-A(t)*(P_as - P_vs)/R_s^2;

df_CvO2(12) = ((caco2 - C_vCO2) * dfdcvo2)/V_TCO2;
df_CvO2(13) = (-F_s + (cao2 - C_vO2) * dfdcvo2)/V_TO2;
df_CvO2(14) = 0;

%% ===== dotV_A ===== %
df_dotVA(1) = 0;
df_dotVA(2) = 0;
df_dotVA(3) = 0;
df_dotVA(4) = 0;
df_dotVA(5) = 0;
df_dotVA(6) = 0;
df_dotVA(7) = 0;
df_dotVA(8) = 0;
df_dotVA(9) = 0;
df_dotVA(10) = (P_ICO2 - P_aCO2)/V_ACO2;
df_dotVA(11) = (P_IO2 - P_aO2)/V_AO2;
df_dotVA(12) = 0;
df_dotVA(13) = 0;
df_dotVA(14) = 0;

%%
jacobian = [df_Pas, df_Pvs, df_Pap, df_Pvp, df_Sl, ...
        df_sigmal, df_Sr, df_sigmar, df_H, ...
        df_PaCO2, df_PaO2, df_CvCO2, df_CvO2, ...
        df_dotVA];
end

