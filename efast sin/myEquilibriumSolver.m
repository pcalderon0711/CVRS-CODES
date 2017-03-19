function equilibrium = myEquilibriumSolver(LHSmatrix,x,P_aCO2,mode)

%% myEquilibriumSolver
% Solve for equilibrium values (x: dot(x) = 0) for given H and PaCO2.
% Inputs:
%       params       : vector of parameters
%       H            : value of heart rate
%       P_aCO2       : value of PaCO2
% Output:
%       equilibrium  : vector of equilibrium values

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
H = LHSmatrix(x,34);
V_tot = 5.0582;
if mode == 0 % rest
    A_pesk = LHSmatrix(x,21);
    R_p = LHSmatrix(x,20);
    W = 0;
else
    A_pesk = LHSmatrix(x,32);
    R_p = LHSmatrix(x,33);
    W = 75;
end
%% Solve for the equilibrium values (refer to Habib)

% Metabolic rates
MR_O2 = M_O2 + rho_O2*W;
MR_CO2 = M_CO2 + rho_CO2*W;

% contractilities
sigma_l = 0;
sigma_r = 0;
S_l = beta_l*H/alpha_l;
S_r = beta_r*H/alpha_r;

dotV_A = 863*MR_CO2/(P_aCO2-P_ICO2);
P_aO2 = P_IO2 - MR_O2*(P_aCO2-P_ICO2)/MR_CO2;
R_s = @(F)(A_pesk*(K_a1*((1-exp(-K_a2*P_aO2))^2) - MR_O2/F));

t_d = 1/H - kappa/(H^0.5);
k_l = exp(-t_d/(c_l*R_l));
k_r = exp(-t_d/(c_r*R_r));
a_l = 1 - k_l;
a_r = 1 - k_r;

mu_r = beta_r*c_r*a_r*H^2/alpha_r;
mu_l = beta_l*c_l*a_l*H^2/alpha_l;
lambda_r = beta_r*k_r*H/alpha_r;
lambda_l = beta_l*k_l*H/alpha_l;

g = @(F)(V_tot*(a_l*a_r*F^2 - mu_l*mu_r) + ...
    c_as*(mu_l*(lambda_r+mu_r*R_s(F)) + a_r*(lambda_l+mu_l*R_p)*F)*F + ...
    c_vs*(mu_l*(lambda_r+a_r*R_p*F)+a_r*(lambda_l+a_l*R_s(F)*F)*F)*F + ...
	c_ap*(mu_r*(lambda_l+mu_l*R_p)+a_l*(lambda_r+mu_r*R_s(F))*F)*F + ...
    c_vp*(mu_r*(lambda_l+a_l*R_s(F)*F)+a_l*(lambda_r+a_r*R_p*F)*F)*F);

F0=H^2*sqrt(c_l*c_r*beta_l*beta_r/(alpha_l*alpha_r));
% solve for the root of g near F0
F = fzero(g,F0);

C_vCO2 = K_CO2*P_aCO2 + k_CO2 + MR_CO2/F;
C_vO2 = K_a1*(1-exp(-K_a2*P_aO2))^2 - MR_O2/F;

D = a_l*a_r*F^2 - mu_l*mu_r;
P_vs = -(mu_l*(lambda_r+a_r*R_p*F)+a_r*(lambda_l+a_l*R_s(F)*F)*F)*F/D;
P_vp = -(mu_r*(lambda_l+a_l*R_s(F)*F)+a_l*(lambda_r+a_r*R_p*F)*F)*F/D;
P_as = R_s(F)*F + P_vs;
P_ap = R_p*F + P_vp;

equilibrium = [P_as; P_vs; P_ap; P_vp; S_l; sigma_l; S_r; sigma_r; P_aCO2; P_aO2; C_vCO2; C_vO2];
end

