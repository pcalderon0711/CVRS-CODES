% PARAMETER BASELINE VALUES
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

% Parameter Labels 
PRCC_var={'c_{as}', 'c_{vs}', 'c_{ap}', 'c_{vp}', ...
    'c_l','c_r', 'R_l','R_r', '\kappa','\alpha_l', ...
    '\alpha_r', '\beta_l', '\beta_r', '\gamma_l', ...
    '\gamma_r', 'M_{O_2}', 'M_{CO_2}', '\rho_{O_2}', ...
    '\rho_{CO_2}', 'R_{p,rest}', 'A_{pesk,rest}', ...
    'P_{IO_2}', 'P_{ICO_2}', 'V_{AO_2}', 'V_{ACO_2}',...
    'V_{TO_2}', 'V_{TCO_2}', 'K_{CO_2}', 'k_{CO_2}', ...
    'K_{a1}', 'K_{a2}', 'A_{pesk,exer}', 'R_{p,exer}', ...
    'H_{rest}', 'H_{exer}', 'dot_{VA,rest}', 'dot_{VA,exer}'};% 

%% TIME SPAN OF THE SIMULATION
t_end=15; % length of the simulations
tspan=(0:0.1:t_end);   % time points where the output is calculated
time_points=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]; % time points of interest for the US analysis time points 4,6,12

% Variables Labels
y_var_label={'P_{as}','P_{vs}', 'P_{ap}', 'P_{vp}', 'S_l', '\sigma_l', 'S_r', ...
    '\sigma_r', 'P_{aCO2}', 'P_{aO2}', 'C_{vCO2}', 'C_{vO2}'};
