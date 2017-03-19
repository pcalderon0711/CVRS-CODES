%% PARAMETER INITIALIZATION
% set up max and mix matrices
c_as = 0.01016;
c_vs = 0.6500;
c_ap = 0.03608;
c_vp = 0.1408;
c_l = 0.02305;
c_r = 0.04413;
R_l = 0.2671;
R_r = 0.04150;
kappa = 0.05164;
alpha_l = 30.5587;
alpha_r = 28.6785;
beta_l = 25.0652;
beta_r = 1.4132;
gamma_l = -1.6744;
gamma_r = -1.8607;
M_O2 = 0.35;
M_CO2 = 0.28;
rho_O2 = 0.011;
rho_CO2 = 0.009;
q_as = 163.047;
V_tot = 5.0582;
R_p_rest = 1.5446;
A_pesk_rest = 177.682;
P_IO2 = 150;
P_ICO2 = 0;
V_AO2 = 2.5;
V_ACO2 = 3.2;
V_TO2 = 6.0;
V_TCO2 = 15.0;
K_CO2 = 0.0057;
k_CO2 = 0.224;
K_a1 = 0.2;
K_a2 = 0.05;
W_rest = 0;
A_pesk_exer = 270;
R_p_exer = 0.3;
W_exer = 75;
H_rest = 78.5;
H_exer = 107;
dot_VA_rest = 6.04;
dot_VA_exer = 20.6;
dummy = 1;

pmean = [c_as, c_vs, c_ap, c_vp, c_l, c_r, R_l, R_r, kappa, alpha_l, alpha_r,...
    beta_l, beta_r, gamma_l, gamma_r, M_O2, M_CO2, rho_O2, rho_CO2, R_p_rest, A_pesk_rest,...
    P_IO2, P_ICO2, V_AO2, V_ACO2, V_TO2, V_TCO2, K_CO2, k_CO2, K_a1, K_a2, A_pesk_exer,...
    R_p_exer, H_rest, H_exer, dot_VA_rest,dot_VA_exer, dummy]';
pmin = 0.8*pmean;
pmax = 1.2*pmean;

% pmin=[1e-2, 
% 0.4, 
% 0.02, 
% 0.1, 
% 0.01, 
% 0.03, 
% 0.1,
% 0,
% 0.03,
% 25,
% 25,
% 20,
% 1,
% -2,
% -2,
% 0.1,
% 0.1,
% 0.005,
% 0.005,
% 1.2,
% 160,
% 140,
% 0,
% 2,
% 2,
% 4,
% 13,
% 0.001,
% 0.2,
% 0.1,
% 0,
% 250,
% 0.1,
% 1]; % dummy
% 
% pmax=[2e-2, 
% 0.6,
% 0.05,
% 0.2,
% 0.03, 
% 0.05,
% 0.05,
% 0.1, % muV
% 0.07,
% 35,
% 35,
% 30,
% 2,
% -1,
% -1,
% 1,
% 1,
% 0.015,
% 0.015,
% 1.8,
% 190,
% 160,
% 1,
% 3,
% 4,
% 8,
% 17,
% 0.01,
% 0.3,
% 0.3,
% 0.1,
% 300,
% 0.5,
% 10]; % dummy
% 
% pmean = [0.01016,
% 0.6500,
% 0.03608,
% 0.1408,
% 0.02305,
% 0.04413,
% 0.2671,
% 0.04150,
% 0.05164,
% 30.5587,
% 28.6785,
% 25.0652,
% 1.4132,
% -1.6744,
% -1.8607,
% 0.35,
% 0.28,
% 0.011,
% 0.009,
% 1.5446,
% 177.682,
% 150,
% 0,
% 2.5,
% 3.2,
% 6.0,
% 15.0,
% 0.0057,
% 0.224,
% 0.2,
% 0.05,
% 270,
% 0.3,
% 1];
% 
pstd = [2e-3,
0.05,
0.005,
0.05,
0.005,
0.005,
0.05,
0.01,
0.005,
5,
5,
5,
0.5,
0.5,
0.5,
0.05,
0.05,
0.0005,
0.0005,
0.1,
10,
10,
0.05,
0.5,
0.5,
1,
5,
0.05,
0.05,
0.05,
0.05,
10,
0.05,
0.1];

efast_var={'c_{as}', 'c_{vs}', 'c_{ap}', 'c_{vp}', ...
    'c_l','c_r', 'R_l','R_r', '\kappa','\alpha_l', ...
    '\alpha_r', '\beta_l', '\beta_r', '\gamma_l', ...
    '\gamma_r', 'M_{O2}', 'M_{CO2}', '\rho_{O2}', ...
    '\rho_{CO2}', 'R_{p,rest}', 'A_{pesk,rest}', ...
    'P_{IO2}', 'P_{ICO2}', 'V_{AO2}', 'V_{ACO2}',...
    'V_{TO2}', 'V_{TCO2}', 'K_{CO2}', 'k_{CO2}', ...
    'K_{a1}', 'K_{a2}', 'A_{pesk,exer}', 'R_{p,exer}',...
    'H_{rest}', 'H_{exer}', 'dot_{VA,rest}', 'dot_{VA,exer}','dummy'};% 

%% TIME SPAN OF THE SIMULATION
t_end=15; % length of the simulations
tspan=(0:0.1:t_end);   % time points where the output is calculated
time_points=[40 60 120]; % time points of interest for the US analysis

% Variables Labels
y_var_label={'P_{as}','P_{vs}', 'P_{ap}', 'P_{vp}', 'S_l', '\sigma_l', 'S_r', ...
    '\sigma_r', 'P_{aCO2}', 'P_{aO2}', 'C_{vCO2}', 'C_{vO2}'};
