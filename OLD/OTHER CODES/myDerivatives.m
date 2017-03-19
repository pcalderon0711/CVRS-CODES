function jacobian = myDerivatives(Y,params)

% this is independent of W

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
    
% ========= PARAMETERS ========== %
c_as = params(1);
c_vs = params(2);
c_ap = params(3);
c_vp = params(4);
c_l = params(5);
c_r = params(6);
R_l = params(7);
R_r = params(8);
kappa = params(9);
alpha_l = params(10);
alpha_r = params(11);
beta_l = params(12);
beta_r = params(13);
gamma_l = params(14);
gamma_r = params(15);
M_O2 = params(16);
M_CO2 = params(17);
rho_O2 = params(18);
rho_CO2 = params(19);
q_as = params(20);
V_tot = params(21);
R_p = params(22);
A_pesk = params(23);
P_IO2 = params(24);
P_ICO2 = params(25);
V_AO2 = params(26);
V_ACO2 = params(27);
V_TO2 = params(28);
V_TCO2 = params(29);
K_CO2 = params(30);
k_CO2 = params(31);
K_a1 = params(32);
K_a2 = params(33);
W = params(34);

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

% ========== COMPUTED QUANTITIES ========= %
R_s = A_pesk * C_vO2;
t_d = (1/H^0.5)*((1/H^0.5) - kappa);
k_l = exp(-(1/(c_l*R_l))*t_d);
a_l = 1 - k_l;
k_r = exp(-(1/(c_r*R_r))*t_d);
a_r = 1 - k_r;

F_s = (1/R_s) * (P_as - P_vs);
F_p = (1/R_p) * (P_ap - P_vp);
Q_l = H*((c_l*a_l*P_vp*S_l)/(a_l*P_as + k_l*S_l));
Q_r = H*((c_r*a_r*P_vs*S_r)/(a_r*P_ap + k_r*S_r));

MR_O2 = M_O2 + rho_O2*W;
MR_CO2 = M_CO2 + rho_CO2*W;

% ===== P_as ===== %

dFs_Pas = 1/R_s;
dQl_Pas = -(c_l*(a_l^2)*P_vp*S_l*H)/(a_l*P_as+k_l*S_l)^2;

df_Pas(1) = (1/c_as)*(dQl_Pas - dFs_Pas);
df_Pas(2) = (1/c_vs)*dFs_Pas;
df_Pas(3) = 0;
df_Pas(4) = -(1/c_vp)*dQl_Pas;
df_Pas(5) = 0;
df_Pas(6) = 0;
df_Pas(7) = 0;
df_Pas(8) = 0;
df_Pas(9) = 0;
df_Pas(10) = 0;
df_Pas(11) = 0;
df_Pas(12) = (1/V_TCO2)*(K_CO2*P_aCO2 + k_CO2 - C_vCO2)*dFs_Pas;
df_Pas(13) = (1/V_TO2)*(K_a1*((1-exp(-K_a2*P_aO2))^2) - C_vO2)*dFs_Pas;
df_Pas(14) = 0;

% ===== P_vs ===== %
dFs_Pvs = -1/R_s;
dQr_Pvs = (c_r*a_r*S_r*H)/(a_r*P_ap + k_r*S_r);

df_Pvs(1) = -(1/c_as)*dFs_Pvs;
df_Pvs(2) = (1/c_vs)*(dFs_Pvs - dQr_Pvs);
df_Pvs(3) = (1/c_ap)*dQr_Pvs;
df_Pvs(4) = 0;
df_Pvs(5) = 0;
df_Pvs(6) = 0;
df_Pvs(7) = 0;
df_Pvs(8) = 0;
df_Pvs(9) = 0;
df_Pvs(10) = 0;
df_Pvs(11) = 0;
df_Pvs(12) = (1/V_TCO2)*(K_CO2*P_aCO2 + k_CO2 - C_vCO2)*dFs_Pvs;
df_Pvs(13) = (1/V_TO2)*(K_a1*((1-exp(-K_a2*P_aO2))^2) - C_vO2)*dFs_Pvs;
df_Pvs(14) = 0;

% ===== P_ap ===== %
dFp_Pap = 1/R_p;
dQr_Pap = -(c_r*(a_r^2)*P_vs*S_r*H)/(a_r*P_ap + k_r*S_r)^2;

df_Pap(1) = 0;
df_Pap(2) = -(1/c_vs)*dQr_Pap;
df_Pap(3) = (1/c_ap)*(dQr_Pap - dFp_Pap);
df_Pap(4) = (1/c_vp)*dFp_Pap;
df_Pap(5) = 0;
df_Pap(6) = 0;
df_Pap(7) = 0;
df_Pap(8) = 0;
df_Pap(9) = 0;
df_Pap(10) = (863/V_ACO2)*(C_vCO2 - K_CO2*P_aCO2 - k_CO2)*dFp_Pap;
df_Pap(11) = (863/V_AO2)*(C_vO2 - K_a1*(1-exp(-K_a2*P_aO2))^2)*dFp_Pap;
df_Pap(12) = 0;
df_Pap(13) = 0;
df_Pap(14) = 0;

% ===== P_vp ===== %
dFp_Pvp = -1/R_p;
dQl_Pvp = (c_l*a_l*S_l*H)/(a_l*P_as + k_l*S_l);

df_Pvp(1) = (1/c_as)*dQl_Pvp;
df_Pvp(2) = 0;
df_Pvp(3) = -(1/c_ap)*dFp_Pvp;
df_Pvp(4) = (1/c_vp)*(dFp_Pvp - dQl_Pvp);
df_Pvp(5) = 0;
df_Pvp(6) = 0;
df_Pvp(7) = 0;
df_Pvp(8) = 0;
df_Pvp(9) = 0;
df_Pvp(10) = (863/V_ACO2)*(C_vCO2 - K_CO2*P_aCO2 - k_CO2)*dFp_Pvp;
df_Pvp(11) = (863/V_AO2)*(C_vO2 - K_a1*(1-exp(-K_a2*P_aO2))^2)*dFp_Pvp;
df_Pvp(12) = 0;
df_Pvp(13) = 0;
df_Pvp(14) = 0;

% ===== S_l ===== %
dQl_Sl = (c_l*(a_l^2)*P_as*P_vp*H)/(a_l*P_as + k_l*S_l)^2;

df_Sl(1) = (1/c_as)*dQl_Sl;
df_Sl(2) = 0;
df_Sl(3) = 0;
df_Sl(4) = -(1/c_vp)*dQl_Sl;
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

% ===== sigma_l ===== %
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

% ===== S_r ===== %
dQr_Sr = (c_r*(a_r^2)*P_ap*P_vs*H)/(a_r*P_ap + k_r*S_r)^2;

df_Sr(1) = 0;
df_Sr(2) = -(1/c_vs)*dQr_Sr;
df_Sr(3) = (1/c_ap)*dQr_Sr;
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

% ===== sigma_r ===== %
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

% ===== H ===== %
t_prime = -(1/H^2)*(1-(kappa/2)*H^0.5);
kl_prime = -(1/(c_l*R_l))*k_l*t_prime;
al_prime = -kl_prime;
kr_prime = -(1/(c_r*R_r))*k_r*t_prime;
ar_prime = -kr_prime;
dQl_H = c_l*P_vp*S_l*((a_l*(a_l*P_as + k_l*S_l))-kl_prime*S_l*H)/((a_l*P_as + k_l*S_l)^2);
dQr_H = c_r*P_vs*S_r*((a_r*(a_r*P_ap + k_r*S_r))-kr_prime*S_r*H)/((a_r*P_ap + k_r*S_r)^2);

df_H(1) = (1/c_as)*dQl_H;
df_H(2) = -(1/c_vs)*dQr_H;
df_H(3) = (1/c_ap)*dQr_H;
df_H(4) = -(1/c_vp)*dQl_H;
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

% ===== P_aCO2 ===== %
df_PaCO2(1) = 0;
df_PaCO2(2) = 0;
df_PaCO2(3) = 0;
df_PaCO2(4) = 0;
df_PaCO2(5) = 0;
df_PaCO2(6) = 0;
df_PaCO2(7) = 0;
df_PaCO2(8) = 0;
df_PaCO2(9) = 0;
df_PaCO2(10) = (-863*K_CO2*F_p/V_ACO2) - (dotV_A/V_ACO2);
df_PaCO2(11) = 0;
df_PaCO2(12) = (K_CO2*F_s)/V_TCO2;
df_PaCO2(13) = 0;
df_PaCO2(14) = 0;

% ===== P_aO2 ===== %
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
df_PaO2(11) = -(1726*K_a1*K_a2*F_p/V_AO2)*(1-exp(-K_a2*P_aO2))*exp(-K_a2*P_aO2) - dotV_A/V_AO2;
df_PaO2(12) = 0;
df_PaO2(13) = (2*K_a1*K_a2*F_s/V_TO2)*(1-exp(-K_a2*P_aO2))*exp(-K_a2*P_aO2);
df_PaO2(14) = 0;

% ===== C_vCO2 ===== %
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
df_CvCO2(12) = -(F_s/V_TCO2);
df_CvCO2(13) = 0;
df_CvCO2(14) = 0;

% ===== C_vO2 ===== %
df_CvO2(1) = (A_pesk/(c_as*R_s^2))*(P_as - P_vs);
df_CvO2(2) = -(A_pesk/(c_vs*R_s^2))*(P_as - P_vs);
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
df_CvO2(12) = ((K_CO2*P_aCO2 + k_CO2 - C_vCO2) * (-(A_pesk/R_s^2)*(P_as - P_vs)))/V_TCO2;
df_CvO2(13) = (-F_s + ((K_a1*(1-exp(-K_a2*P_aO2))^2) - C_vO2) * (-(A_pesk/R_s^2)*(P_as - P_vs)))/V_TO2;
df_CvO2(14) = 0;

% ===== dotV_A ===== %
df_dotVA(1) = 0;
df_dotVA(2) = 0;
df_dotVA(3) = 0;
df_dotVA(4) = 0;
df_dotVA(5) = 0;
df_dotVA(6) = 0;
df_dotVA(7) = 0;
df_dotVA(8) = 0;
df_dotVA(9) = 0;
df_dotVA(10) = (1/V_ACO2)*(P_ICO2 - P_aCO2);
df_dotVA(11) = (1/V_AO2)*(P_IO2 - P_aO2);
df_dotVA(12) = 0;
df_dotVA(13) = 0;
df_dotVA(14) = 0;

%%%%%%%%%%%%

jacobian = [df_Pas, df_Pvs, df_Pap, df_Pvp, df_Sl, ...
        df_sigmal, df_Sr, df_sigmar, df_H, ...
        df_PaCO2, df_PaO2, df_CvCO2, df_CvO2, ...
        df_dotVA];
end

