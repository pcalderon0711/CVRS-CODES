%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear all;
close all;

%% Sample size N
runs=1000;

%% LHS MATRIX  %%
Parameter_settings_LHS;

c_as_LHS=LHS_Call(.8*c_as, c_as, 1.2*c_as, 2e-3,runs,'unif'); % baseline = 10
c_vs_LHS=LHS_Call(.8*c_vs, c_vs, 1.2*c_vs, 0.05 ,runs,'unif'); % baseline = 2e-2
c_ap_LHS=LHS_Call(.8*c_ap, c_ap, 1.2*c_ap, 0.005, runs,'unif'); % baseline = 3e-2
c_vp_LHS=LHS_Call(.8*c_vp, c_vp, 1.2*c_vp,0.01 ,runs,'unif'); % baseline = 2.4e-5
c_l_LHS=LHS_Call(.8*c_l, c_l, 1.2*c_l, 0.005, runs,'unif'); % baseline = 3e-3
c_r_LHS= LHS_Call(.8*c_r , c_r , 1.2*c_r , 0.005 , runs , 'unif');  % baseline = 0.24
R_l_LHS=LHS_Call(.8*R_l,R_l, 1.2*R_l,0.3 ,runs,'unif'); % baseline = 1200
R_r_LHS=LHS_Call(.8*R_r,R_r, 1.2*R_r,0.1,runs,'unif'); % baseline = 2.4
kappa_LHS=LHS_Call(.8*kappa,kappa,1.2*kappa, 0.005 ,runs,'unif'); % dummy parameter
alpha_l_LHS=LHS_Call(.8*alpha_l, alpha_l, 1.2*alpha_l, 5 ,runs,'unif'); % baseline = 10
alpha_r_LHS=LHS_Call(.8*alpha_r, alpha_r, 1.2*alpha_r, 5 ,runs,'unif'); % baseline = 2e-2
beta_l_LHS=LHS_Call(.8*beta_l, beta_l, 1.2*beta_l, 5, runs,'unif'); % baseline = 3e-2
beta_r_LHS=LHS_Call(.8*beta_r, beta_r, 1.2*beta_r, 0.5 ,runs,'unif'); % baseline = 2.4e-5
gamma_l_LHS=LHS_Call(.8*gamma_l, gamma_l, 1.2*gamma_l, 0.5, runs,'unif'); % baseline = 3e-3
gamma_r_LHS= LHS_Call(.8*gamma_r, gamma_r , 1.2*gamma_r , 0.5 , runs , 'unif');  % baseline = 0.24
M_O2_LHS=LHS_Call(.8*M_O2, M_O2,1.2*M_O2, 0.05 ,runs,'unif'); % baseline = 1200
M_CO2_LHS=LHS_Call(.8*M_O2, M_CO2,1.2*M_CO2, 0.05 ,runs,'unif'); % baseline = 2.4
rho_O2_LHS=LHS_Call(.8*rho_O2,rho_O2,1.2*rho_O2, 0.0005 ,runs,'unif'); % dummy parameter
rho_CO2_LHS=LHS_Call(.8*rho_CO2,rho_CO2, 1.2*rho_CO2, 0.0005 ,runs,'unif'); % baseline = 10
R_p_rest_LHS=LHS_Call(.8*R_p_rest,R_p_rest,1.2*R_p_rest, 0.1 ,runs,'unif'); % baseline = 2.4e-5
A_pesk_rest_LHS=LHS_Call(.8*A_pesk_rest, A_pesk_rest, 1.2*A_pesk_rest, 10, runs,'unif'); % baseline = 3e-3
P_IO2_LHS= LHS_Call(.8*P_IO2 , P_IO2 , 1.2*P_IO2 , 10 , runs , 'unif');  % baseline = 0.24
P_ICO2_LHS=LHS_Call(.8*P_ICO2 ,P_ICO2,1.2*P_ICO2, 1 ,runs,'unif'); % baseline = 1200
V_AO2_LHS=LHS_Call(.8*V_AO2,V_AO2,1.2*V_AO2, 0.5 ,runs,'unif'); % baseline = 2.4
V_ACO2_LHS=LHS_Call(.8*V_ACO2,V_ACO2,1.2*V_ACO2, 0.5 ,runs,'unif'); % dummy parameter
V_TO2_LHS=LHS_Call(.8*V_TO2, V_TO2, 1.2*V_TO2, 1 ,runs,'unif'); % baseline = 10
V_TCO2_LHS=LHS_Call(.8*V_TCO2, V_TCO2, 1.2*V_TCO2, 5 ,runs,'unif'); % baseline = 2e-2
K_CO2_LHS=LHS_Call(.8*K_CO2, K_CO2, 1.2*K_CO2, 0.05, runs,'unif'); % baseline = 3e-2
k_CO2_LHS=LHS_Call(.8*k_CO2,k_CO2,1.2*k_CO2, 0.05 ,runs,'unif'); % baseline = 2.4e-5
K_a1_LHS=LHS_Call(.8*K_a1, K_a1,1.2*K_a1, 0.05, runs,'unif'); % baseline = 3e-3
K_a2_LHS= LHS_Call(.8*K_a2 , K_a2, 1.2*K_a2 , 0.05 , runs , 'unif');  % baseline = 0.24
A_pesk_exer_LHS=LHS_Call(.8*A_pesk_exer,A_pesk_exer,1.2*A_pesk_exer, 10 ,runs,'unif'); % baseline = 2.4
R_p_exer_LHS=LHS_Call(.8*R_p_exer,R_p_exer,1.2*R_p_exer, 0.05 ,runs,'unif'); % dummy parameter
H_rest_LHS=LHS_Call(.8*H_rest,H_rest,1.2*H_rest, 5 ,runs,'unif');
H_exer_LHS=LHS_Call(.8*H_exer,H_exer,1.2*H_exer, 5 ,runs,'unif');
dot_VA_rest_LHS=LHS_Call(.8*dot_VA_rest,dot_VA_rest,1.2*dot_VA_rest, 1 ,runs,'unif');
dot_VA_exer_LHS=LHS_Call(.8*dot_VA_exer,dot_VA_exer,1.2*dot_VA_exer, 1 ,runs,'unif');

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix=[c_as_LHS c_vs_LHS c_ap_LHS c_vp_LHS c_l_LHS c_r_LHS ...
              R_l_LHS R_r_LHS kappa_LHS  alpha_l_LHS alpha_r_LHS ...
              beta_l_LHS beta_r_LHS gamma_l_LHS gamma_r_LHS M_O2_LHS M_CO2_LHS ...
              rho_O2_LHS rho_CO2_LHS R_p_rest_LHS A_pesk_rest_LHS P_IO2_LHS ...
              P_ICO2_LHS V_AO2_LHS V_ACO2_LHS V_TO2_LHS V_TCO2_LHS K_CO2_LHS ...
              k_CO2_LHS K_a1_LHS K_a2_LHS A_pesk_exer_LHS R_p_exer_LHS, ...
              H_rest_LHS, H_exer_LHS, dot_VA_rest_LHS, dot_VA_exer_LHS];

for i=1:runs
i
ff=myEquilibriumSolver(LHSmatrix,i,40,0);
ff(1)
end

ts = 5;
opts = optimset('Diagnostics','off', 'Display','off');
h = 0.1;

for r=1:runs %Run solution x times choosing different values
    f=@ODE_LHS;
    r
    y0 = myEquilibriumSolver(LHSmatrix,r,40,0);
    A = myBWEulerSolver(@(t,y) f(t,y,LHSmatrix,r,ts), h, tspan,y0);

    cd 'output'
    x = plot(tspan,A(:,1));
    x.LineWidth = 4;
    xlabel('t');
    ylabel('P_{as}');
    set(gca, 'FontSize', 15)
    savefig(sprintf('output_run_%d',r));
    print(sprintf('output_run_%d',r),'-dpng');
    clf    
    cd ..
    
%     [t,y]=ode15s(@(t,y)f(t,y,LHSmatrix,x,runs),tspan,y0,[]); 
%     A=[t y]; % [time y]

    %% Save the outputs at ALL time points [tspan]
    %T_lhs(:,x)=Anew(:,1);
    %CD4_lhs(:,x)=Anew(:,2);
    %T1_lhs(:,x)=Anew(:,3);
    %T2_lhs(:,x)=Anew(:,4);
    %V_lhs(:,x)=Anew(:,5);
    
    %% Save only the outputs at the time points of interest [time_points]:
    %% MORE EFFICIENT
    Pas_lhs(:,r)=A(time_points+1,1);
end
%% Save the workspace
save Model_LHS.mat;
% CALCULATE PRCC
[prcc sign sign_label]=PRCC(LHSmatrix,Pas_lhs,time_points,PRCC_var,0.05);
for i=1:3
    PRCC_PLOT(LHSmatrix,Pas_lhs,i,PRCC_var,0.05);
end
close all